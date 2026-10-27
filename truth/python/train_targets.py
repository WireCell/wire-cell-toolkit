#!/usr/bin/env python3
"""Training targets from the native per-point truth labels (wcp-porting-img pr_ml doc 25 sec 7.4; owner decisions
2026-10-24 / 2026-10-25).  Every training reads its per-point targets through targets(), so the rules below are applied
once, the same way, with the parameters recorded (PARAMS, version).

targets() is a pure function on per-point and per-particle arrays; read_event() builds those arrays from a clustering
point-tree dump labelled by TruthPointLabeler (pctree-evt.tar.gz), the event's truth table (truth.npz) and the wires
file.  The rules, in this order (the distance cut must see the delta / Michel nodes before the fold removes them):

  1. frame (post-Q/L, the light-less stand-in = the true t0): every cluster moves along the drift by -dir * v * t0.
     t0 = the start time of the cluster's dominant primary (the Geant4 ancestor without a mother, weighted by the
     points' charge); dir = the drift direction of the point's face (+1 when the face's U plane lies at larger x than
     its W plane, i.e. the drift volume is on the +x side); v = the label-frame drift speed.  A cluster without a
     labelled point keeps its position (frame_valid 0).
  2. em-daughter distance cut (Michel and delta nodes alike, R = radius_mm): a point of a Michel or delta node whose
     particle (orig_tid) is not the node's root particle and that lies farther than R from the root's own points
     (in the corrected frame) leaves the node.  Its branch -- the root's child it descends from -- becomes a shower
     instance of its own when the branch's visible charge reaches compton_q (the PFMerger rule), else LE (instance
     not valid).  The branch's visible charge is the truth nvis_q of its Geant4 subtree; a particle not in the truth
     table (larg4 dropped it, SBND / FD-HD) uses the image charge of its detached points instead.  A node whose root
     has no point of its own is left as it is (counted in root_missing).
  3. delta fold: a delta node's points take the track class and the instance of the parent node (the node of the
     root's first ancestor with points in another node); without one, the track class and their own instance.
  4. d_traj mask: track-labelled points with d_traj >= dtraj_cut_mm (the labeler's 25 mm cap) are left out of the class
     and instance losses; ghosts are never masked by it (they carry the cap by construction).
  5. orphan_le (default off): a labelled point without a PF node (tid >= 0, pf_node -1: a nucleus with no visible ancestor
     in the truth table -- its mother not stored, or an all-invisible chain such as K0 -> K0L -> alpha) takes class LE with
     no instance.  Off, its labeler class -1 passes through and reads as GHOST (pr_ml doc 26 sec 7.3, owner 2026-10-26).

Per point out: x / y / z (mm, corrected frame), cls (-1 ghost, -2 no label; else 0 track, 1 shower, 2 Michel, 3 delta,
4 LE after the rules), inst (-1 none), mask_cls, mask_inst, frame_valid, detached.

usage: train_targets.py --pctree <pctree-evt.tar.gz> --truth <truth.npz> --wires <wires.json.bz2> --speed-mmus <v>
                        [--out targets.npz] [--json summary.json]
       train_targets.py --selftest
"""
import argparse, bz2, json, os, sys
import numpy as np
from scipy.spatial import cKDTree

PARAMS = dict(version=1, frame="true_t0_dominant_primary", radius_mm=100.0, em_cut_classes=[2, 3], delta_fold=True,
              dtraj_cut_mm=25.0, compton_q=2.0e5, orphan_le=False)
TRACK, SHOWER, MICHEL, DELTA, LE = 0, 1, 2, 3, 4
GHOST, NOLABEL = -1, -2


def _ancestors(t, mother, nmax=100000):
    out, seen = [], set()
    t = mother.get(t, 0)
    while t > 0 and t not in seen and len(out) < nmax:
        out.append(t); seen.add(t); t = mother.get(t, 0)
    return out


def _primary(t, mother):
    a = _ancestors(t, mother)
    return a[-1] if a else t


def targets(P, truth, speed_mmus, params=None):
    """P: per-point arrays X (N,3) mm, cl (cluster), dir (+-1), q (image charge), tid, orig_tid, pf_node, cls,
    instance_valid, d_traj (mm).  truth: dicts mother {tid: mother tid}, start_t_ns {tid: ns}, nvis_q {tid: e}.
    Returns (out arrays, summary)."""
    p = dict(PARAMS, **(params or {}))
    X = np.asarray(P["X"], float).copy()
    cl, dirx, q = np.asarray(P["cl"]), np.asarray(P["dir"], float), np.asarray(P["q"], float)
    tid, otid, pfn = np.asarray(P["tid"]).astype(np.int64), np.asarray(P["orig_tid"]).astype(np.int64), np.asarray(P["pf_node"]).astype(np.int64)
    cls, iv, dtr = np.asarray(P["cls"]).astype(np.int64).copy(), np.asarray(P["instance_valid"]).astype(bool), np.asarray(P["d_traj"], float)
    mother, st, nvis = truth["mother"], truth["start_t_ns"], truth["nvis_q"]
    n = len(tid)
    lab = tid >= 0
    S = dict(npoints=int(n), nlabelled=int(lab.sum()), params=p)

    # 1. frame: per cluster, the dominant primary's start time
    frame_valid = np.zeros(n, bool)
    t0_cl, dom_share = {}, []
    prim = np.array([_primary(int(t), mother) if t >= 0 else -1 for t in tid])
    for c in np.unique(cl):
        s = (cl == c) & lab
        if not s.any():
            continue
        w, wt = {}, {}
        for pr, qq in zip(prim[s], np.maximum(q[s], 1e-9)):
            w[int(pr)] = w.get(int(pr), 0.0) + float(qq)
            if int(pr) in st:   # primaries of one interaction share its time: the mixing measure groups by time (10 ns)
                b = round(st[int(pr)] / 10.0)
                wt[b] = wt.get(b, 0.0) + float(qq)
        best = max(sorted(w), key=lambda k: w[k])
        if best not in st:
            continue
        t0_cl[int(c)] = st[best] / 1000.0   # us
        dom_share.append(max(wt.values()) / sum(wt.values()) if wt else 1.0)
        sel = cl == c
        X[sel, 0] -= dirx[sel] * speed_mmus * t0_cl[int(c)]
        frame_valid[sel] = True
    S["frame"] = dict(clusters=int(len(np.unique(cl))), shifted=len(t0_cl),
                      dominant_share_p10=float(np.percentile(dom_share, 10)) if dom_share else None,
                      clusters_dominant_below_0p9=int(sum(d < 0.9 for d in dom_share)))

    # node bookkeeping (pr_truth parent logic): majority node per track key, node root
    inst = np.where(lab & (pfn >= 0), pfn, -1).astype(np.int64)
    key2node, best_n = {}, {}
    pairs, cnt = (np.unique(np.column_stack([tid[lab & (pfn >= 0)], pfn[lab & (pfn >= 0)]]), axis=0, return_counts=True)
                  if (lab & (pfn >= 0)).any() else (np.zeros((0, 2), int), np.zeros(0, int)))
    for (t, k), c in zip(pairs.tolist(), cnt.tolist()):
        if c > best_n.get(t, 0):
            key2node[t], best_n[t] = k, c
    keys_of = {}
    for t, k in key2node.items():
        keys_of.setdefault(k, []).append(t)
    root_of = {}
    for k, ks in keys_of.items():
        kset = set(ks)
        roots = [t for t in ks if not any(a in kset for a in _ancestors(t, mother))] or ks
        if len(roots) == 1:
            root_of[k] = roots[0]
        else:
            # several top particles: the node's root is their lowest common ancestor (every particle of a PF node
            # descends from its root), possibly a particle without points of its own
            chains = [[t] + _ancestors(t, mother) for t in roots]
            common = set(chains[0]).intersection(*map(set, chains[1:]))
            lca = next((a for a in chains[0] if a in common), None)
            root_of[k] = lca if lca is not None else max(roots, key=lambda t: (best_n[t], -t))
    node_cls = {}
    for k in keys_of:
        s = lab & (pfn == k)
        node_cls[k] = int(np.bincount(np.clip(cls[s], 0, 4)).argmax())

    # 2. em-daughter distance cut on Michel / delta nodes
    detached = np.zeros(n, bool)
    children = {}
    for t, m in mother.items():
        children.setdefault(m, []).append(t)

    def subtree_q(b):
        tot, stack, seen = 0.0, [b], set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x); tot += nvis.get(x, 0.0); stack.extend(children.get(x, ()))
        return tot

    next_inst = int(inst.max()) + 1 if n and inst.max() >= 0 else 0
    root_missing = cut_nodes = 0
    em_nodes = []
    for k in sorted(keys_of):
        if node_cls[k] not in p["em_cut_classes"]:
            continue
        s = np.flatnonzero(lab & (pfn == k))
        root = root_of[k]
        own = s[otid[s] == root]
        dau = s[otid[s] != root]
        if not len(dau):
            continue
        if not len(own):
            root_missing += 1
            continue
        d = cKDTree(X[own]).query(X[dau])[0]
        far = dau[d > p["radius_mm"]]
        em_nodes.append(dict(node=int(k), cls=node_cls[k], root=int(root), points=int(len(s)), detached=int(len(far)),
                             max_dist_mm=round(float(d.max()), 1)))
        if not len(far):
            continue
        cut_nodes += 1
        branch_of = {}
        for o in np.unique(otid[far]).tolist():
            chain = [o] + _ancestors(o, mother)
            b = o
            for a, nxt in zip(chain, chain[1:] + [0]):
                if nxt == root:
                    b = a; break
            branch_of[o] = b
        br = np.array([branch_of[int(o)] for o in otid[far]])
        for b in np.unique(br).tolist():
            sel = far[br == b]
            qb = subtree_q(b) if b in mother else float(q[sel].sum())
            if qb >= p["compton_q"]:
                cls[sel] = SHOWER; inst[sel] = next_inst; iv[sel] = True; next_inst += 1
            else:
                cls[sel] = LE; inst[sel] = -1; iv[sel] = False
            detached[sel] = True
    S["em_cut"] = dict(nodes_cut=cut_nodes, root_missing=root_missing, points_detached=int(detached.sum()),
                       detached_to_shower=int((detached & (cls == SHOWER)).sum()), detached_to_le=int((detached & (cls == LE)).sum()))
    for nm, kc in (("michel", MICHEL), ("delta", DELTA)):
        pts = int(sum(e["points"] for e in em_nodes if e["cls"] == kc)); dt = int(sum(e["detached"] for e in em_nodes if e["cls"] == kc))
        tot = int(sum((lab & (pfn == k)).sum() for k in keys_of if node_cls[k] == kc))
        S["em_cut"][nm] = dict(node_points=tot, detached=dt, detached_share=round(dt / max(tot, 1), 4))
    S["em_nodes"] = em_nodes

    # 3. delta fold
    folded = 0
    if p["delta_fold"]:
        for k in sorted(keys_of):
            if node_cls[k] != DELTA:
                continue
            s = lab & (pfn == k) & ~detached
            parent = None
            for a in _ancestors(root_of[k], mother):
                n2 = key2node.get(a)
                if n2 is not None and n2 != k:
                    parent = n2; break
            cls[s] = TRACK
            if parent is not None:
                inst[s] = parent
            folded += int(s.sum())
    S["delta_fold"] = dict(points=folded)

    # 5. orphan nuclei (no PF node) -> LE, no instance
    orphan = lab & (pfn < 0)
    if p["orphan_le"]:
        cls[orphan] = LE; inst[orphan] = -1; iv[orphan] = False
    S["orphan"] = dict(points=int(orphan.sum()), to_le=bool(p["orphan_le"]))

    # 4. masks
    out_cls = np.where(tid == -1, GHOST, np.where(tid < -1, NOLABEL, cls))
    far_traj = lab & (dtr >= p["dtraj_cut_mm"])
    mask_cls = (out_cls != NOLABEL) & ~far_traj
    mask_inst = lab & iv & ~far_traj & (inst >= 0)
    S["masks"] = dict(dtraj_masked=int(far_traj.sum()), dtraj_masked_share_of_labelled=float(far_traj.sum() / max(lab.sum(), 1)),
                      class_points=int(mask_cls.sum()), instance_points=int(mask_inst.sum()))
    out = dict(x=X[:, 0], y=X[:, 1], z=X[:, 2], cls=out_cls.astype(np.int16), inst=inst.astype(np.int64),
               mask_cls=mask_cls, mask_inst=mask_inst, frame_valid=frame_valid, detached=detached)
    return out, S


# ---------------------------------------------------------------- readers
def face_dirs(wires_file):
    """{(anode ident, face index): +-1}: +1 when the face's U plane lies at larger x than its W plane."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from xres_check import unwrap
    S = json.load(bz2.open(wires_file))["Store"]
    out = {}
    for a in S["anodes"]:
        a = unwrap(a)
        for fidx, fi in enumerate(a["faces"]):
            f = unwrap(S["faces"][fi])
            xs = []
            for pi in f["planes"]:
                pl = unwrap(S["planes"][pi]); w = unwrap(S["wires"][pl["wires"][0]])
                xs.append(unwrap(S["points"][w["tail"]])["x"])
            out[(a["ident"], fidx)] = 1 if xs[0] > xs[-1] else -1
    return out


def read_event(pctree, truth_npz, wires):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from steiner_truth import load_tree, per_node
    import truth_table as T
    tree = load_tree(pctree)
    P3, TR, SC = per_node(tree, "3d"), per_node(tree, "truth"), per_node(tree, "scalar")
    co = tree["cluster_of"]
    fd = face_dirs(wires)
    cols = {k: [] for k in ("x", "y", "z", "cl", "dir", "q", "tid", "orig_tid", "pf_node", "cls", "instance_valid", "d_traj")}
    for nd in sorted(P3):
        p3 = P3[nd]; m = len(p3["x"]); tr = TR.get(nd); sc = SC.get(nd)
        g = (lambda k, dflt: np.asarray(tr[k]) if tr is not None else np.full(m, dflt))
        wpid = int(sc["wpid"][0]) if sc is not None else 0
        apa, face = wpid >> 4, (wpid >> 3) & 1
        qb = float(sc["charge"][0]) / max(int(sc["npoints"][0]), 1) if sc is not None else 0.0
        for k, v in (("x", p3["x"]), ("y", p3["y"]), ("z", p3["z"]), ("cl", np.full(m, int(co[nd]))), ("dir", np.full(m, fd.get((apa, face), 1))),
                     ("q", np.full(m, qb)), ("tid", g("tid", -3)), ("orig_tid", g("orig_tid", -3)), ("pf_node", g("pf_node", -1)),
                     ("cls", g("cls", -1)), ("instance_valid", g("instance_valid", 0)), ("d_traj", g("d_traj", 0.0))):
            cols[k].append(np.asarray(v))
    P = {k: np.concatenate(v) if v else np.zeros(0) for k, v in cols.items()}
    P["X"] = np.column_stack([P.pop("x"), P.pop("y"), P.pop("z")]).astype(float)
    z = np.load(truth_npz)
    pf = z["truth_pf"]; c = T.pf_col
    tids = pf[:, c("trackid")].astype(np.int64)
    truth = dict(mother=dict(zip(tids.tolist(), pf[:, c("mother_trackid")].astype(np.int64).tolist())),
                 start_t_ns=dict(zip(tids.tolist(), pf[:, c("start_t")].astype(float).tolist())),
                 nvis_q=dict(zip(tids.tolist(), pf[:, c("nvis_q")].astype(float).tolist())))
    return P, truth


# ---------------------------------------------------------------- self test
def selftest():
    """Synthetic event: a muon (1) at t0 = 100 us stopping at the origin; its Michel (2) with a near photon (3) and a far
    photon (4) in one em_fold node; a delta (5) of the muon with a far, large photon (6); a ghost point; drift +x on
    cluster 0 and -x on cluster 1 (the far photon of the Michel converted across the cathode)."""
    v, t0 = 1.6, 100.0
    shift = v * t0
    mother = {1: 0, 2: 1, 3: 2, 4: 2, 5: 1, 6: 5, 7: 0}
    start = {1: t0 * 1000, 2: t0 * 1000 + 2000, 3: t0 * 1000 + 2000, 4: t0 * 1000 + 2000, 5: t0 * 1000, 6: t0 * 1000, 7: t0 * 1000}
    nvis = {1: 1e6, 2: 3e5, 3: 1e4, 4: 1e4, 5: 5e4, 6: 4e5, 7: 3e3}
    pts = []   # (true x, y, z, cluster, dir, tid, orig_tid, pf_node, cls, iv, d_traj)
    for i in range(20):
        pts.append((-200.0 + 10 * i, 0, 0, 0, +1, 1, 1, 0, TRACK, 1, 3.0))       # muon, node 0
    for i in range(5):
        pts.append((5.0 * i, 5, 0, 0, +1, 2, 2, 1, MICHEL, 1, 2.0))              # Michel own, node 1
    pts.append((20.0, 20, 0, 0, +1, 3, 3, 1, MICHEL, 1, 2.0))                     # near photon (< R)
    for i in range(3):
        pts.append((-300.0 - i, 400, 0, 1, -1, 4, 4, 1, MICHEL, 1, 2.0))         # far photon (> R), other drift side
    for i in range(4):
        pts.append((-100.0, 10 + 3 * i, 0, 0, +1, 5, 5, 2, DELTA, 1, 2.0))       # delta, node 2
    for i in range(3):
        pts.append((-100.0, 400 + i, 0, 0, +1, 6, 6, 2, DELTA, 1, 2.0))          # delta's far photon, 4e5 e -> shower
    pts.append((-150.0, 0, 30, 0, +1, 1, 1, 0, TRACK, 1, 30.0))                   # muon point far off its trajectory
    pts.append((0.0, -50, 0, 0, +1, -1, -1, -1, -1, 0, 25.0))                    # ghost (d_traj at the cap)
    for i in range(2):
        pts.append((-120.0, 30 + i, 0, 0, +1, 7, 7, -1, -1, 1, 5.0))             # a nucleus (7) without a PF node
    A = np.array(pts, dtype=float)
    Xtrue = A[:, :3]
    X = Xtrue.copy(); X[:, 0] += A[:, 4] * shift                                 # what the raw frame shows
    P = dict(X=X, cl=A[:, 3].astype(int), dir=A[:, 4], q=np.full(len(A), 1e4), tid=A[:, 5].astype(int), orig_tid=A[:, 6].astype(int),
             pf_node=A[:, 7].astype(int), cls=A[:, 8].astype(int), instance_valid=A[:, 9].astype(int), d_traj=A[:, 10])
    out, S = targets(P, dict(mother=mother, start_t_ns=start, nvis_q=nvis), v)
    Xo = np.column_stack([out["x"], out["y"], out["z"]])
    ok = []
    ok.append(("frame: both drift sides back to the true positions", np.allclose(Xo, Xtrue)))
    ok.append(("frame: every cluster shifted", out["frame_valid"].all()))
    mich_own = np.flatnonzero((P["tid"] == 2))
    ok.append(("Michel own points stay Michel", (out["cls"][mich_own] == MICHEL).all() and not out["detached"][mich_own].any()))
    near = np.flatnonzero(P["tid"] == 3)
    ok.append(("near Michel photon stays in the Michel", (out["cls"][near] == MICHEL).all() and (out["inst"][near] == 1).all()))
    far = np.flatnonzero(P["tid"] == 4)
    ok.append(("far Michel photon (1e4 e) detached to LE, no instance", (out["cls"][far] == LE).all() and (~out["mask_inst"][far]).all()))
    dl = np.flatnonzero(P["tid"] == 5)
    ok.append(("delta folded into the muon: track class, muon instance", (out["cls"][dl] == TRACK).all() and (out["inst"][dl] == 0).all()))
    dfar = np.flatnonzero(P["tid"] == 6)
    ok.append(("far delta photon (4e5 e) detached to its own shower instance",
               (out["cls"][dfar] == SHOWER).all() and len(set(out["inst"][dfar])) == 1 and out["inst"][dfar][0] not in (0, 1, 2)))
    offt = np.flatnonzero(P["d_traj"] == 30.0)
    ok.append(("track point beyond the d_traj cut masked", (~out["mask_cls"][offt]).all() and (~out["mask_inst"][offt]).all()))
    gh = np.flatnonzero(P["tid"] == -1)
    ok.append(("ghost stays ghost and in the class loss", (out["cls"][gh] == GHOST).all() and out["mask_cls"][gh].all()))
    orp = np.flatnonzero(P["tid"] == 7)
    ok.append(("orphan nucleus, orphan_le off: labeler class -1 passes through (reads as GHOST), counted",
               (out["cls"][orp] == GHOST).all() and S["orphan"]["points"] == 2 and not S["orphan"]["to_le"]))
    outo, So = targets(P, dict(mother=mother, start_t_ns=start, nvis_q=nvis), v, dict(orphan_le=True))
    ok.append(("orphan nucleus, orphan_le on: LE, in the class loss, no instance",
               (outo["cls"][orp] == LE).all() and outo["mask_cls"][orp].all() and (~outo["mask_inst"][orp]).all() and (outo["inst"][orp] == -1).all()))
    nonorp = np.flatnonzero(P["tid"] != 7)
    ok.append(("orphan_le on changes nothing else", all(np.array_equal(out[k][nonorp], outo[k][nonorp]) for k in out)))
    # root without own points: the Michel node is left as it is
    P2 = {k: (v.copy() if hasattr(v, "copy") else v) for k, v in P.items()}
    keep = P2["tid"] != 2
    P2 = {k: v[keep] for k, v in P2.items()}
    out2, S2 = targets(P2, dict(mother=mother, start_t_ns=start, nvis_q=nvis), v)
    far2 = np.flatnonzero(P2["tid"] == 4)
    ok.append(("root without own points: node unchanged, counted", (out2["cls"][far2] == MICHEL).all() and S2["em_cut"]["root_missing"] == 1))
    # SBND-style folded ids: the photons carry the Michel's tid, their own orig_tid, and are not in the truth table
    P3 = {k: (v.copy() if hasattr(v, "copy") else v) for k, v in P.items()}
    fol = np.isin(P3["orig_tid"], [3, 4])
    P3["tid"] = np.where(fol, 2, P3["tid"])
    m3 = {k: v for k, v in mother.items() if k not in (3, 4)}
    out3, _ = targets(P3, dict(mother=m3, start_t_ns=start, nvis_q=nvis), v)
    far3 = np.flatnonzero(P3["orig_tid"] == 4)
    near3 = np.flatnonzero(P3["orig_tid"] == 3)
    ok.append(("folded ids (orig_tid path): far photon detached, near kept", (out3["cls"][far3] == LE).all() and (out3["cls"][near3] == MICHEL).all()))
    for name, good in ok:
        print(("PASS " if good else "FAIL ") + name)
    print(json.dumps({k: v for k, v in S.items() if k != "params"}))
    return all(g for _, g in ok)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--pctree"); ap.add_argument("--truth"); ap.add_argument("--wires"); ap.add_argument("--speed-mmus", type=float)
    ap.add_argument("--radius-mm", type=float, default=PARAMS["radius_mm"]); ap.add_argument("--dtraj-cut-mm", type=float, default=PARAMS["dtraj_cut_mm"])
    ap.add_argument("--orphan-le", action="store_true", help="rule 5: nucleus points without a PF node -> LE (default off)")
    ap.add_argument("--out"); ap.add_argument("--json")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    P, truth = read_event(a.pctree, a.truth, a.wires)
    out, S = targets(P, truth, a.speed_mmus, dict(radius_mm=a.radius_mm, dtraj_cut_mm=a.dtraj_cut_mm, orphan_le=a.orphan_le))
    if a.out:
        np.savez_compressed(a.out, **out, params=json.dumps(S["params"]))
    txt = json.dumps(S, indent=1)
    if a.json:
        open(a.json, "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
