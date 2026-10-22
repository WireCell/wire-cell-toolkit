#!/usr/bin/env python3
"""Per-blob parity of the native truth labels against larwirecell's wclsTensorSetLabeler (pr_ml doc 19 sec 7.3).

The LArSoft labeler writes the dominant G4 track id of every live blob into the blob's "scalar" PC and a Bee points
set "truth_trackid" (raw x, y, z in cm, cluster_id = the blob's track id, -1 = no matching deposit).  The native
TruthPointLabeler writes the "truth" PC of every live blob node.  When both run on the same imaging, the sampled
points coincide one for one, so the join is exact: every Bee point is matched to the native point at the same
position (distance ~0), and through it to the native blob node.

Inputs: --pctree = the clustering dump of the native labeler run with refine_by_nearest_depo=false (every point
carries its blob's argmax tid, -2 when the largest fraction is below shared_frac, -1 ghost); optional --pctree-point
= the default run (per-point tid of the nearest contributing depo), used only to list which tracks the native
association sees in a blob; --bee = truth_trackid.json; --truth = the event's truth.npz (MC tree for the relation of
two disagreeing ids).

Per blob: L = LArSoft id (one value per blob, checked), N = native id.  Categories: both labelled and equal / both
labelled and different / native only / LArSoft only / neither.  Agreement = equal / both labelled, by blob count,
by points and by Bee charge.  A disagreement is classed by the MC relation of N and L (N ancestor of L, descendant,
same parent, other) and whether L is among the tracks the native per-point run sees in the blob.
"""

import argparse, json, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_3d import load_tensors, PCTree  # noqa: E402
import truth_table as T  # noqa: E402

CM = 10.0  # Bee cm -> WCT mm


def live_points(fn):
    """x, y, z (mm), tid, pid_frac, global blob-node index of every live point of the dump."""
    tens = load_tensors(fn)
    tops = sorted(dp for dp, (m, _) in tens.items() if m.get("datatype") == "pctree" and dp.endswith("/live"))
    out = {k: [] for k in ("x", "y", "z", "tid", "pid_frac", "node")}
    off = 0
    for dp in tops:
        tree = PCTree(tens, dp)
        if "3d" not in tree.arrays or "truth" not in tree.arrays:
            continue
        p3, tr = tree.arrays["3d"], tree.arrays["truth"]
        _, node = tree.point_cluster("3d")
        for k in ("x", "y", "z"):
            out[k].append(np.asarray(p3[k], dtype=float))
        out["tid"].append(np.asarray(tr["tid"]).astype(np.int64))
        out["pid_frac"].append(np.asarray(tr["pid_frac"], dtype=float))
        out["node"].append(node + off)
        off += len(tree.parentage)
    if not out["x"]:
        raise RuntimeError("no live tree with 3d + truth clouds in " + fn)
    return {k: np.concatenate(v) for k, v in out.items()}


def mc_parents(truth_npz):
    if not truth_npz:
        return {}
    pf = np.load(truth_npz)["truth_pf"]
    c = T.pf_col
    return {int(r[c("trackid")]): int(r[c("mother_trackid")]) for r in pf}


def ancestors(t, par, nmax=200):
    out = []
    while t in par and len(out) < nmax:
        t = par[t]
        out.append(t)
    return out


def relation(n, l, par):
    """MC relation of the native id n to the LArSoft id l (mother_trackid chain)."""
    if not par:
        return "no table"
    if l in ancestors(n, par):
        return "native is a descendant of LArSoft's"
    if n in ancestors(l, par):
        return "native is an ancestor of LArSoft's"
    if par.get(n) is not None and par.get(n) == par.get(l):
        return "siblings"
    return "unrelated"


def mc_tree_nodes(mc_json):
    """{id: parent id in the tree (0 = root)} of a LArSoft Bee "mc" particle-flow tree."""
    out = {}

    def walk(x, parent):
        out[int(x["id"])] = parent
        for ch in x.get("children", []):
            walk(ch, int(x["id"]))
    for x in json.load(open(mc_json)):
        walk(x, 0)
    return out


def native_pf_tree(truth_npz, ke_min_gev, nu_only):
    """The same tree from the native truth table: particles with KE > ke_min, nested under the nearest kept ancestor
    by mother_trackid (the LArSoft labeler's pf_ke_min / Mother() rule); nu_only keeps the rows whose truth_nu row
    is a beam neutrino (nu_idx >= 0; the cosmic generator's row has nu_idx -1)."""
    tab = np.load(truth_npz)
    pf, nu = tab["truth_pf"], tab["truth_nu"]
    c = T.pf_col
    is_nu = {i for i, r in enumerate(nu) if r[T.nu_col("nu_idx")] >= 0}
    mother = {int(r[c("trackid")]): int(r[c("mother_trackid")]) for r in pf}
    kept = {int(r[c("trackid")]) for r in pf if r[c("KE")] > ke_min_gev and (not nu_only or int(r[c("nu_row")]) in is_nu)}
    out = {}
    for t in kept:
        m = mother.get(t, 0)
        while m in mother and m not in kept:
            m = mother[m]
        out[t] = m if m in kept else 0
    return out


def compare_pf(mc_json, truth_npz, ke_min_gev, nu_only):
    L = mc_tree_nodes(mc_json)
    N = native_pf_tree(truth_npz, ke_min_gev, nu_only)
    common = set(L) & set(N)
    return {"nu_only": nu_only, "n_larsoft": len(L), "n_native": len(N), "common": len(common),
            "only_larsoft": sorted(set(L) - set(N))[:20], "only_native": sorted(set(N) - set(L))[:20],
            "parent_agree": int(sum(1 for t in common if L[t] == N[t])),
            "parent_differ": [(t, L[t], N[t]) for t in sorted(common) if L[t] != N[t]][:20]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--event", required=True)
    ap.add_argument("--pctree", required=True, help="native dump, refine_by_nearest_depo=false")
    ap.add_argument("--pctree-point", default="", help="native dump, labeler defaults (optional)")
    ap.add_argument("--bee", required=True, help="LArSoft truth_trackid.json")
    ap.add_argument("--truth", default="", help="truth.npz of the event")
    ap.add_argument("--tol-mm", type=float, default=0.05, help="max Bee-to-native point distance of a match")
    ap.add_argument("--out", required=True)
    ap.add_argument("--label", default="", help="name of the LArSoft variant (for the tables)")
    ap.add_argument("--mc", default="", help="LArSoft Bee mc.json (particle-flow tree) to compare with the native table")
    ap.add_argument("--mc-nu-only", action="store_true", help="the mc tree was made with pf_nu_only")
    ap.add_argument("--pf-ke-min-mev", type=float, default=10.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    from scipy.spatial import cKDTree

    P = live_points(a.pctree)
    B = json.load(open(a.bee))
    bx = np.asarray(B["x"], float) * CM
    by = np.asarray(B["y"], float) * CM
    bz = np.asarray(B["z"], float) * CM
    bl = np.asarray(B["cluster_id"]).astype(np.int64)
    bq = np.asarray(B.get("q", np.ones(len(bx))), float)
    res = {"event": a.event, "label": a.label, "inputs": vars(a), "n_native_points": int(len(P["x"])),
           "n_bee_points": int(len(bx))}

    # Neighbouring blobs can share a sample point at the same position (a few per mille of the points); such a point
    # cannot be attributed to one blob by position, so it is left out (the label is per blob, the blob keeps its
    # other points).  Both dumps list the same points in a different order.
    kd = cKDTree(np.stack([P["x"], P["y"], P["z"]], 1))
    d2, j2 = kd.query(np.stack([bx, by, bz], 1), k=2)
    d, j = d2[:, 0], j2[:, 0]
    ambiguous = d2[:, 1] <= a.tol_mm
    ok = (d <= a.tol_mm) & ~ambiguous
    nodes_all = np.unique(P["node"])
    res["match"] = {"matched_frac": float((d <= a.tol_mm).mean()), "d_median_mm": float(np.median(d)), "d_max_mm": float(d.max()),
                    "ambiguous_points": int(ambiguous.sum()), "one_to_one": bool(len(np.unique(j[ok])) == ok.sum()),
                    "native_blob_nodes": int(len(nodes_all)),
                    "blobs_without_unambiguous_point": int(len(nodes_all) - len(np.unique(P["node"][j[ok]])))}

    # per native blob node: LArSoft id (must be unique per blob), native id, points, Bee charge
    node = P["node"][j[ok]]
    L = bl[ok]
    N = P["tid"][j[ok]]
    q = bq[ok]
    order = np.argsort(node, kind="stable")
    node, L, N, q = node[order], L[order], N[order], q[order]
    starts = np.r_[0, np.flatnonzero(np.diff(node)) + 1]
    ends = np.r_[starts[1:], len(node)]
    blobs = []
    nonunique = 0
    for s, e in zip(starts, ends):
        lu = np.unique(L[s:e])
        nu = np.unique(N[s:e])
        if len(lu) > 1 or len(nu) > 1:
            nonunique += 1
        lmode = lu[np.bincount(np.searchsorted(lu, L[s:e])).argmax()]
        nmode = nu[np.bincount(np.searchsorted(nu, N[s:e])).argmax()]
        blobs.append((int(node[s]), int(lmode), int(nmode), int(e - s), float(q[s:e].sum())))
    res["nblobs"] = len(blobs)
    res["blobs_with_mixed_ids"] = nonunique

    # categories.  larg4 keys the deposits of untracked EM-shower daughters on -(saved ancestor id); the LArSoft labeler
    # keeps that negative id as the blob's label (its own check counts it as unlabelled), the native labeler folds it to
    # the ancestor (fold_negative_ids).  "folded" compares |L| for L < -1 (the convention-normalised parity), "raw"
    # counts L < 0 as unlabelled (LArSoft's own labelled fraction).
    def categorize(fold):
        cats = {"agree": [], "disagree": [], "native_only": [], "larsoft_only": [], "neither": []}
        for b in blobs:
            nd, l, n, npt, qq = b
            if fold and l < -1:
                l = -l
            nl, ll = n >= 0, l >= 0
            if nl and ll:
                cats["agree" if n == l else "disagree"].append((nd, l, n, npt, qq))
            elif nl:
                cats["native_only"].append(b)
            elif ll:
                cats["larsoft_only"].append(b)
            else:
                cats["neither"].append(b)
        return cats

    def tot(lst, i):
        return float(sum(x[i] for x in lst)) if i else float(len(lst))

    def summarize(cats):
        summ = {}
        for w, i in (("blobs", 0), ("points", 3), ("charge", 4)):
            both = tot(cats["agree"], i) + tot(cats["disagree"], i)
            allb = sum(tot(v, i) for v in cats.values())
            summ[w] = {k: tot(v, i) for k, v in cats.items()}
            summ[w]["agreement"] = tot(cats["agree"], i) / both if both else None
            summ[w]["native_labelled_frac"] = (tot(cats["agree"], i) + tot(cats["disagree"], i) + tot(cats["native_only"], i)) / allb
            summ[w]["larsoft_labelled_frac"] = (tot(cats["agree"], i) + tot(cats["disagree"], i) + tot(cats["larsoft_only"], i)) / allb
        return summ

    res["summary_raw"] = summarize(categorize(False))
    cats = categorize(True)
    res["summary"] = summ = summarize(cats)
    negb = [b for b in blobs if b[1] < -1]
    res["larsoft_negative_id_blobs"] = {"n": len(negb), "native_equals_abs": int(sum(1 for b in negb if b[2] == -b[1])),
                                        "native_unlabelled": int(sum(1 for b in negb if b[2] < 0))}
    # native unlabelled split: ghost (-1) vs shared (-2), for the LArSoft-only blobs
    res["larsoft_only_native_code"] = {str(k): int(v) for k, v in zip(*np.unique([b[2] for b in cats["larsoft_only"]], return_counts=True))} if cats["larsoft_only"] else {}

    # disagreement diagnosis
    par = mc_parents(a.truth)
    seen = {}
    if a.pctree_point:
        Q = live_points(a.pctree_point)
        if len(Q["x"]) == len(P["x"]) and np.allclose(Q["x"], P["x"]):
            idx = np.argsort(Q["node"], kind="stable")
            qn, qt = Q["node"][idx], Q["tid"][idx]
            st = np.r_[0, np.flatnonzero(np.diff(qn)) + 1]
            en = np.r_[st[1:], len(qn)]
            seen = {int(qn[s]): set(int(t) for t in np.unique(qt[s:e])) for s, e in zip(st, en)}
        else:
            res["point_run_warning"] = "point-mode dump does not have the same points; contributor check skipped"
    pidf = {}
    for nd, pf in zip(P["node"], P["pid_frac"]):
        pidf.setdefault(int(nd), float(pf))
    rel = {}
    dis_rows = []
    for nd, l, n, npt, qq in cats["disagree"]:
        r = relation(n, l, par)
        inn = (l in seen.get(nd, set())) if seen else None
        key = r + ("" if inn is None else (", L among native contributors" if inn else ", L not seen natively"))
        rel[key] = rel.get(key, 0) + 1
        dis_rows.append({"node": nd, "larsoft": l, "native": n, "points": npt, "charge": qq, "native_pid_frac": pidf.get(nd),
                         "relation": r, "larsoft_in_native_contributors": inn})
    res["disagree_relations"] = rel
    fr = np.array([x["native_pid_frac"] for x in dis_rows if x["native_pid_frac"] is not None])
    res["disagree_native_pid_frac"] = {"median": float(np.median(fr)), "frac_below_0p7": float((fr < 0.7).mean())} if len(fr) else {}
    agree_fr = np.array([pidf.get(b[0], np.nan) for b in cats["agree"]])
    res["agree_native_pid_frac_median"] = float(np.nanmedian(agree_fr)) if len(agree_fr) else None
    if a.mc and a.truth:
        res["pf_tree"] = compare_pf(a.mc, a.truth, a.pf_ke_min_mev / 1000.0, a.mc_nu_only)
        print("pf tree:", res["pf_tree"])
    json.dump(res, open(os.path.join(a.out, a.event + ".json"), "w"), indent=1)
    json.dump(dis_rows, open(os.path.join(a.out, a.event + "_disagree.json"), "w"), indent=1)
    sb = summ["blobs"]
    sr = res["summary_raw"]["blobs"]
    print("%s %s: match %.4f (d_max %.3g mm, %d ambiguous points, 1:1 %s, mixed-id blobs %d, %d of %d blob nodes without a point) blobs %d | folded: both %d agree %.4f (pts %.4f, q %.4f) native-only %d larsoft-only %d neither %d labelled native %.3f larsoft %.3f | raw: both %d agree %.4f larsoft labelled %.3f | neg-id blobs %s | disagree %s"
          % (a.event, a.label, res["match"]["matched_frac"], res["match"]["d_max_mm"], res["match"]["ambiguous_points"],
             res["match"]["one_to_one"], nonunique, res["match"]["blobs_without_unambiguous_point"], res["match"]["native_blob_nodes"], len(blobs), sb["agree"] + sb["disagree"],
             sb["agreement"] or 0, summ["points"]["agreement"] or 0, summ["charge"]["agreement"] or 0,
             sb["native_only"], sb["larsoft_only"], sb["neither"], sb["native_labelled_frac"], sb["larsoft_labelled_frac"],
             sr["agree"] + sr["disagree"], sr["agreement"] or 0, sr["larsoft_labelled_frac"], res["larsoft_negative_id_blobs"], rel))


if __name__ == "__main__":
    main()
