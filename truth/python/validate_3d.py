#!/usr/bin/env python3
"""validate_3d.py -- E0.6: imaging and points against the true depos (pr_ml doc 19 sec 7).

Reads, per event, the products of the native truth route:
  * blobtruth archives written by BlobDepoLabelTensor (pr_ml/cfg/wct-blobtruth.jsonnet or an in-chain sink): the
    blobtruth (M,6) and blobtruth_sum (B,7) tensors with the completeness denominator in the metadata
    (q_depos_window = |q| of the depos inside the anode's sensitive volume and the slices' time span);
  * optionally the cluster archives the blobtruth was made from (clusters-*.tar.gz): the solved charge per blob
    (row-aligned with blobtruth_sum: both are in graph vertex order) for the ghost fraction by charge;
  * optionally a point-tree dump (TensorFileSink of a clustering job run with TruthPointLabeler): the "truth" local
    point cloud of every live blob node, the tree parentage (cluster membership) and the "3d" positions;
  * optionally the truth file truth.npz (class of a trackid from its pdg when no cloud is given).

Writes <out>/summary.json, <out>/per_event.json, <out>/summary.md and PNGs.  Metrics (sec 7.1, 7.2, 7.4):
  per blob: completeness = sum q_true(blobs) / q_depos_window; ghost fraction by count (q_true <= 0) and by solved
            charge; purity = distribution of max_k f_bk per class; mixed fraction (max f < shared_frac); per-slice
            occupancy of mixed blobs;
  per point: labelled / valid fractions; d_traj bands (tau_on = 2 pitches, tau_off = 5 pitches) per class; nearest-
            neighbour tid coherence (kNN on the 3d positions); points per cm of true trajectory where the table gives lengths;
  matching: the intersection table shared[cluster, tid] by POINT COUNT (charge per point is not in the tree) ->
            SPINE purity/efficiency/IoU, Pandora's cut, pimm PQ (instances < 10 points excluded, stated), ARI with LE and
            shared points masked.
usage:
  validate_3d.py --out <dir> --det <name> --pitch-mm <p> [--shared-frac 0.5] \
      --event <tag> --blobtruth <archive>[,<archive>...] [--clusters <archive>[,...]] [--pctree <dump>] [--truth truth.npz] \
      [--event ... (repeat the group per event)]
"""
import argparse, io, json, os, re, sys, tarfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import truth_table as T  # noqa: E402

CLS_NAMES = {0: "track", 1: "shower", 2: "michel", 3: "delta", 4: "le", -1: "unknown"}
TRACK_PDGS = {13, 211, 321, 2212, 3222, 3112}


def load_tensors(fn):
    """{datapath: (metadata, array)} of a TensorFileSink archive."""
    md, arr = {}, {}
    with tarfile.open(fn) as t:
        for m in t.getmembers():
            if not m.isfile():
                continue
            d = t.extractfile(m).read()
            if m.name.endswith("metadata.json"):
                md[m.name] = json.loads(d)
            elif m.name.endswith("array.npy"):
                arr[m.name] = np.load(io.BytesIO(d), allow_pickle=False)
    out = {}
    for k, m in md.items():
        if "datapath" not in m:
            continue
        a = arr.get(k.replace("metadata.json", "array.npy"))
        out[m["datapath"]] = (m, a)
    return out


def by_datatype(tens, dt):
    return [(m, a) for dp, (m, a) in sorted(tens.items()) if m.get("datatype") == dt]


def cluster_blob_values(fn):
    """solved charge per blob, graph vertex order (the *_bnodes rows of a ClusterFileSink numpy archive, tar or npz)."""
    vals = []
    if not fn or not os.path.exists(fn):
        return np.zeros(0)
    if fn.endswith(".npz"):
        z = np.load(fn, allow_pickle=False)
        for n in sorted(k for k in z.files if k.endswith("_bnodes")):
            r = z[n]; vals.append(r.reshape(-1, r.shape[-1])[:, 2].astype(np.float64))
        return np.concatenate(vals) if vals else np.zeros(0)
    with tarfile.open(fn) as t:
        names = sorted(m.name for m in t.getmembers() if m.name.endswith("_bnodes.npy"))
        for n in names:
            r = np.load(io.BytesIO(t.extractfile(n).read()), allow_pickle=False)
            vals.append(r.reshape(-1, r.shape[-1])[:, 2].astype(np.float64))
    return np.concatenate(vals) if vals else np.zeros(0)


def pdg_class(pdg):
    a = abs(int(pdg))
    if a in TRACK_PDGS or a > 1000000000:
        return 0
    if a in (11, 22):
        return 1
    return -1


class PCTree:
    """A live point-tree dump: per-node sizes of the 3d / truth clouds, the cluster of every node, the arrays."""

    def __init__(self, tens, datapath=None):
        tops = [(dp, m, a) for dp, (m, a) in tens.items() if m.get("datatype") == "pctree" and dp.endswith("/live")]
        if datapath:
            tops = [x for x in tops if x[0] == datapath]
        if not tops:
            raise RuntimeError("no live pctree in the dump")
        dp, m, parentage = tops[0]
        self.dp = dp
        self.parentage = np.asarray(parentage).astype(int)
        lp = tens[m["lpcmaps"]][0]["arrays"]
        self.sizes = {k: np.asarray(tens[v][1]).astype(int) for k, v in lp.items()}
        pcs = tens[m["pointclouds"]][0]["items"]
        self.arrays = {}
        for name, path in pcs.items():
            amd = tens[path][0]
            self.arrays[name] = {k: np.asarray(tens[v][1]) for k, v in amd["arrays"].items()}
        # node -> cluster node (the child of the root on its path)
        n = len(self.parentage)
        self.cluster_of = np.full(n, -1, dtype=int)
        for i in range(1, n):
            p = self.parentage[i]
            self.cluster_of[i] = i if p == 0 else self.cluster_of[p]

    def point_cluster(self, pcname="3d"):
        """cluster node index of every point of the concatenated cloud, and the node index."""
        sizes = self.sizes[pcname]
        node = np.repeat(np.arange(len(sizes)), sizes)
        return self.cluster_of[node], node


def per_blob(bt_files, cl_files, shared_frac, truth_tab):
    rows_sum, rows_tab, q_window, ndepos, solved = [], [], 0.0, 0, []
    for i, fn in enumerate(bt_files):
        tens = load_tensors(fn)
        for m, a in by_datatype(tens, "blobtruth_sum"):
            rows_sum.append(np.asarray(a, dtype=np.float64).reshape(-1, len(m["columns"])))
            q_window += float(m.get("q_depos_window", m.get("q_depos", 0.0)))
            ndepos += int(m.get("ndepos_window", m.get("ndepos", 0)))
            cols_sum = m["columns"]
        for m, a in by_datatype(tens, "blobtruth"):
            rows_tab.append(np.asarray(a, dtype=np.float64).reshape(-1, len(m["columns"])))
            cols_tab = m["columns"]
        if cl_files and i < len(cl_files) and cl_files[i]:
            solved.append(cluster_blob_values(cl_files[i]))
    if not rows_sum:
        # an event whose imaging made no blob at all (R2: a PDHD K+ that decayed before the TPC): no table, no metrics
        return {"nblobs": 0, "nfilled": 0, "nghost": 0, "note": "no blobtruth_sum table in the archives"}, np.zeros(0)
    S = np.concatenate(rows_sum)
    Tb = np.concatenate(rows_tab) if rows_tab else np.zeros((0, 6))
    c = {k: cols_sum.index(k) for k in cols_sum}
    q_true, q_top, tid_top, ntr = S[:, c["q_true"]], S[:, c["q_top"]], S[:, c["tid_top"]].astype(int), S[:, c["ntracks"]].astype(int)
    filled = q_true > 0
    fmax = np.where(filled, q_top / np.maximum(q_true, 1e-12), 0.0)
    mixed = filled & (fmax < shared_frac)
    res = {
        "nblobs": int(len(S)), "nfilled": int(filled.sum()), "nghost": int((~filled).sum()),
        "ghost_frac_count": float((~filled).mean()) if len(S) else 0.0,
        "ndepos_window": ndepos, "q_depos_window": q_window, "q_true_in_blobs": float(q_true.sum()),
        "completeness": float(q_true.sum() / q_window) if q_window > 0 else None,
        "nmulti": int((ntr >= 2).sum()), "nmixed": int(mixed.sum()),
        "mixed_frac_of_filled": float(mixed.sum() / max(filled.sum(), 1)),
        "fmax_hist": np.histogram(fmax[filled], bins=np.linspace(0, 1, 11))[0].tolist(),
        "fmax_median": float(np.median(fmax[filled])) if filled.any() else None,
    }
    if solved:
        v = np.concatenate(solved)
        if len(v) == len(S):
            res["solved_total"] = float(v.sum())
            res["ghost_frac_charge"] = float(v[~filled].sum() / v.sum()) if v.sum() > 0 else None
            res["true_q_frac_median"] = float(np.median(q_true[filled] / np.maximum(v[filled], 1e-12))) if filled.any() else None
        else:
            res["solved_mismatch"] = [int(len(v)), int(len(S))]
    # purity per class of the top track (needs the table for pdg)
    if truth_tab is not None and filled.any():
        pf = truth_tab["truth_pf"]
        tid2pdg = dict(zip(pf[:, T.pf_col("trackid")].astype(int), pf[:, T.pf_col("pdg")].astype(int)))
        cls = np.array([pdg_class(tid2pdg.get(int(t), 0)) for t in tid_top])
        res["fmax_median_by_class"] = {CLS_NAMES[k]: float(np.median(fmax[filled & (cls == k)])) for k in (0, 1, -1) if (filled & (cls == k)).any()}
        res["mixed_frac_by_class"] = {CLS_NAMES[k]: float(mixed[filled & (cls == k)].mean()) for k in (0, 1, -1) if (filled & (cls == k)).any()}
    # mixed blobs per slice occupancy
    if "slice_ident" in c and filled.any():
        sl = S[:, c["slice_ident"]].astype(int)
        occ = {}
        for s in np.unique(sl):
            occ[int(s)] = int((sl == s).sum())
        occ_arr = np.array([occ[int(s)] for s in sl])
        bins = [1, 2, 4, 8, 16, 1e9]
        res["mixed_frac_by_slice_occupancy"] = {}
        for lo, hi in zip(bins[:-1], bins[1:]):
            m = filled & (occ_arr >= lo) & (occ_arr < hi)
            if m.any():
                res["mixed_frac_by_slice_occupancy"]["%d-%s" % (lo, "inf" if hi > 1e8 else int(hi - 1))] = float(mixed[m].mean())
    # per-track visible charge (for the table's nvis_q comparison)
    if len(Tb):
        ct = {k: cols_tab.index(k) for k in cols_tab}
        tids, q = Tb[:, ct["trackid"]].astype(int), Tb[:, ct["q_bk"]]
        res["ntracks_with_charge"] = int(len(np.unique(tids)))
        res["q_bk_total"] = float(q.sum())
    return res, fmax[filled]


def per_point(pctree_file, pitch_mm, truth_tab, shared_frac):
    tens = load_tensors(pctree_file)
    try:
        tree = PCTree(tens)
    except RuntimeError as err:   # an empty dump (an event without blobs)
        return {"error": str(err)}, None
    if "truth" not in tree.arrays:
        return {"error": "no truth cloud in the dump"}, None
    tr, p3 = tree.arrays["truth"], tree.arrays["3d"]
    tid, valid, cls = tr["tid"].astype(int), tr["instance_valid"].astype(int), tr["cls"].astype(int)
    d = tr["d_traj"].astype(float)
    n = len(tid)
    lab = tid >= 0
    res = {"npoints": int(n), "labelled_frac": float(lab.mean()), "ghost_frac": float((tid == T.TID_GHOST).mean()),
           "shared_frac": float((tid == T.TID_SHARED).mean()), "valid_frac": float(valid.mean()),
           "ntracks": int(len(np.unique(tid[lab]))), "cls_counts": {CLS_NAMES.get(int(k), str(k)): int(v) for k, v in zip(*np.unique(cls, return_counts=True))}}
    tau_on, tau_off = 2 * pitch_mm, 5 * pitch_mm
    res["d_traj"] = {}
    for k in (0, 1, 2, 3, 4):
        m = lab & (cls == k)
        if not m.any():
            continue
        res["d_traj"][CLS_NAMES[k]] = {"n": int(m.sum()), "median_mm": float(np.median(d[m])), "p90_mm": float(np.percentile(d[m], 90)),
                                       "frac_within_tau_on": float((d[m] <= tau_on).mean()), "frac_within_tau_off": float((d[m] <= tau_off).mean())}
    # nearest-neighbour tid coherence
    try:
        from scipy.spatial import cKDTree
        xyz = np.stack([p3["x"], p3["y"], p3["z"]], axis=1).astype(float)
        if n > 1:
            tree_k = cKDTree(xyz)
            _, idx = tree_k.query(xyz, k=2)
            nn = idx[:, 1]
            m = lab & (tid[nn] >= 0)
            res["nn_tid_coherence"] = float((tid[m] == tid[nn][m]).mean()) if m.any() else None
    except ImportError:
        res["nn_tid_coherence"] = "scipy missing"
    # points per cm of true trajectory
    if truth_tab is not None:
        pf = truth_tab["truth_pf"]
        c = T.pf_col
        L = {}
        for r in pf:
            t = int(r[c("trackid")])
            L[t] = float(np.hypot(np.hypot(r[c("end_x")] - r[c("start_x")], r[c("end_y")] - r[c("start_y")]), r[c("end_z")] - r[c("start_z")]))
        u, cnt = np.unique(tid[lab], return_counts=True)
        dens = [cnt[i] / L[int(t)] for i, t in enumerate(u) if L.get(int(t), 0) > 5.0]
        res["points_per_cm_median"] = float(np.median(dens)) if dens else None
    # matching: cluster <-> tid by point counts
    cl, node = tree.point_cluster("3d")
    inter = {}
    for ci, ti in zip(cl[lab], tid[lab]):
        inter[(int(ci), int(ti))] = inter.get((int(ci), int(ti)), 0) + 1
    reco_tot, true_tot = {}, {}
    for (ci, ti), w in inter.items():
        reco_tot[ci] = reco_tot.get(ci, 0) + w
        true_tot[ti] = true_tot.get(ti, 0) + w
    for ci in np.unique(cl):
        reco_tot.setdefault(int(ci), 0)
    pairs = []
    for (ci, ti), w in inter.items():
        pur, eff = w / reco_tot[ci], w / true_tot[ti]
        pairs.append((ci, ti, w, pur, eff, w / (reco_tot[ci] + true_tot[ti] - w)))
    pairs.sort(key=lambda p: (-p[5], p[0], p[1]))
    used_r, used_t, matched = set(), set(), []
    min_pts = 10
    for ci, ti, w, pur, eff, iou in pairs:
        if iou < 0.5 or ci in used_r or ti in used_t:
            continue
        used_r.add(ci); used_t.add(ti); matched.append(iou)
    recos = [ci for ci, v in reco_tot.items() if v >= min_pts]
    trues = [ti for ti, v in true_tot.items() if v >= min_pts]
    tp = len(matched); fp = len(recos) - len([ci for ci in used_r if reco_tot[ci] >= min_pts]); fn_ = len(trues) - len([ti for ti in used_t if true_tot[ti] >= min_pts])
    sq = float(np.mean(matched)) if matched else 0.0
    rq = tp / (tp + 0.5 * max(fp, 0) + 0.5 * max(fn_, 0)) if tp + fp + fn_ > 0 else 0.0
    best = {}
    for ci, ti, w, pur, eff, iou in pairs:
        if ci not in best or w > best[ci][1]:
            best[ci] = (ti, w, pur, eff)
    res["matching_counts"] = {
        "nclusters": int(len(reco_tot)), "ntracks": int(len(true_tot)), "min_points": min_pts,
        "cluster_best_purity_median": float(np.median([b[2] for b in best.values()])) if best else None,
        "cluster_best_efficiency_median": float(np.median([b[3] for b in best.values()])) if best else None,
        "pandora_good_clusters": int(sum(1 for b in best.values() if b[1] >= 5 and b[2] >= 0.5 and b[3] >= 0.1)),
        "pq": sq * rq, "sq": sq, "rq": rq, "tp": tp, "fp": int(max(fp, 0)), "fn": int(max(fn_, 0)),
    }
    # ARI over valid points: cluster label vs tid
    try:
        mask = lab & (valid == 1)
        res["ari_valid"] = ari(cl[mask], tid[mask])
    except Exception as e:  # noqa: BLE001
        res["ari_valid"] = str(e)
    return res, (d[lab], cls[lab])


def ari(a, b):
    from collections import Counter
    n = len(a)
    if n < 2:
        return 1.0
    nij = Counter(zip(a.tolist(), b.tolist())); ai = Counter(a.tolist()); bj = Counter(b.tolist())
    c2 = lambda x: x * (x - 1) / 2.0  # noqa: E731
    sij = sum(c2(v) for v in nij.values()); sa = sum(c2(v) for v in ai.values()); sb = sum(c2(v) for v in bj.values())
    cn = c2(n); exp = sa * sb / cn; mx = 0.5 * (sa + sb)
    return 1.0 if mx == exp else float((sij - exp) / (mx - exp))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True); ap.add_argument("--det", required=True); ap.add_argument("--pitch-mm", type=float, required=True)
    ap.add_argument("--shared-frac", type=float, default=0.5)
    ap.add_argument("--event", action="append", required=True, help="event tag; the following --blobtruth/--clusters/--pctree/--truth belong to it")
    ap.add_argument("--blobtruth", action="append", default=[], help="comma list of blobtruth archives of the event")
    ap.add_argument("--clusters", action="append", default=[], help="comma list of the cluster archives the blobtruth came from ('' = none)")
    ap.add_argument("--pctree", action="append", default=[], help="pctree dump with the truth cloud ('' = none)")
    ap.add_argument("--truth", action="append", default=[], help="truth.npz ('' = none)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    per_event, fmax_all, dtraj_all = {}, [], []
    for i, ev in enumerate(a.event):
        bt = a.blobtruth[i].split(",") if i < len(a.blobtruth) and a.blobtruth[i] else []
        cl = a.clusters[i].split(",") if i < len(a.clusters) and a.clusters[i] else []
        pc = a.pctree[i] if i < len(a.pctree) else ""
        tf = a.truth[i] if i < len(a.truth) else ""
        tab = T.read_truth_npz(tf) if tf else None
        r = {"event": ev}
        if bt:
            r["blob"], fm = per_blob(bt, cl, a.shared_frac, tab)
            fmax_all.append(fm)
        if pc:
            r["point"], dd = per_point(pc, a.pitch_mm, tab, a.shared_frac)
            if dd is not None:
                dtraj_all.append(dd)
        per_event[ev] = r
        print(ev, json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if not isinstance(vv, (dict, list))}) for k, v in r.items() if k != "event"})[:600])
    # summary
    def agg(key, sub):
        vals = [e[key][sub] for e in per_event.values() if key in e and isinstance(e[key].get(sub), (int, float))]
        return {"mean": float(np.mean(vals)), "min": float(np.min(vals)), "max": float(np.max(vals)), "n": len(vals)} if vals else None
    summary = {"detector": a.det, "pitch_mm": a.pitch_mm, "shared_frac": a.shared_frac, "nevents": len(per_event),
               "completeness": agg("blob", "completeness"), "ghost_frac_count": agg("blob", "ghost_frac_count"),
               "ghost_frac_charge": agg("blob", "ghost_frac_charge"), "mixed_frac_of_filled": agg("blob", "mixed_frac_of_filled"),
               "fmax_median": agg("blob", "fmax_median"), "labelled_frac": agg("point", "labelled_frac"),
               "valid_frac": agg("point", "valid_frac"), "nn_tid_coherence": agg("point", "nn_tid_coherence")}
    json.dump(summary, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    json.dump(per_event, open(os.path.join(a.out, "per_event.json"), "w"), indent=1)
    with open(os.path.join(a.out, "summary.md"), "w") as f:
        f.write("| metric | mean | min | max | n |\n|---|---|---|---|---|\n")
        for k, v in summary.items():
            if isinstance(v, dict):
                f.write("| %s | %.4f | %.4f | %.4f | %d |\n" % (k, v["mean"], v["min"], v["max"], v["n"]))
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        if fmax_all:
            fm = np.concatenate(fmax_all)
            plt.figure(figsize=(5, 3.5)); plt.hist(fm, bins=20, range=(0, 1)); plt.xlabel("max_k f_bk per filled blob"); plt.ylabel("blobs")
            plt.title("%s label purity (n=%d)" % (a.det, len(fm))); plt.tight_layout(); plt.savefig(os.path.join(a.out, "fmax.png"), dpi=120); plt.close()
        if dtraj_all:
            d = np.concatenate([x[0] for x in dtraj_all]); c = np.concatenate([x[1] for x in dtraj_all])
            plt.figure(figsize=(5, 3.5))
            for k in (0, 1, 4):
                if (c == k).any():
                    plt.hist(d[c == k], bins=25, range=(0, 5 * a.pitch_mm), histtype="step", label=CLS_NAMES[k])
            for x, lbl in ((2 * a.pitch_mm, "tau_on"), (5 * a.pitch_mm, "tau_off")):
                plt.axvline(x, ls="--", lw=0.8, color="gray")
            plt.xlabel("d_traj (mm)"); plt.ylabel("points"); plt.legend(); plt.title("%s point-to-trajectory distance" % a.det)
            plt.tight_layout(); plt.savefig(os.path.join(a.out, "d_traj.png"), dpi=120); plt.close()
    except Exception as e:  # noqa: BLE001
        print("no plots:", e)
    print("wrote", a.out, json.dumps({k: v for k, v in summary.items() if not isinstance(v, dict)}))


if __name__ == "__main__":
    main()
