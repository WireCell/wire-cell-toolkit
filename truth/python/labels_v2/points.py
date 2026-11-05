"""labels_v2.points -- from deposits to the PR points (pr_ml doc 54 sec 0, choice 3).

A PR point of the pctree carries the G4 ids of the deposit the C++ labeler chose for it (tid, orig_tid: TruthPointLabeler
refine = nearest eligible deposit) but not that deposit's index.  v2 links each point to the nearest deposit of the same
particle in the reco frame: first by orig_tid among the deposits whose own G4 id is that particle, else by tid among the
deposits folded to that id, else (flagged) the nearest deposit of any particle.  The deposit's labels are then copied to
the point.  Ghost points (tid == -1) keep class ghost and store the distance to the nearest labelled point; points with
tid < -1 (shared / no truth) are "other"; points farther than far_mm from their own trajectory (lab d_traj) keep their
label with the far bit set (50b 2.7, "Other": masked, out of the loss, label kept).  A point whose own instance holds < 2,400 e
is re-owned to the nearest deposit of a visible instance (Q16 b, Xin 2026-11-05), with other_reason bit 128 set.
"""
import numpy as np
from scipy.spatial import cKDTree

from .rules import OBJ, CLS7

FAR_MM = 25.0
# other_reason bits
OR_GHOST, OR_SHARED, OR_NOTRUTH, OR_FAR, OR_NODEPO, OR_FALLBACK_TID, OR_FALLBACK_ANY, OR_INVISIBLE_OWNER = 1, 2, 4, 8, 16, 32, 64, 128


def link(ev, R):
    """-> dict of per-point arrays"""
    pts = ev.pts
    X = pts["X"]; npt = len(X)
    tid, otid, d_traj = pts["tid"], pts["orig_tid"], pts["d_traj"]
    dep = np.full(npt, -1, np.int64); mdist = np.full(npt, np.inf)
    reason = np.zeros(npt, np.int16)
    lab = tid >= 0
    Xr = ev.Xreco; okdep = np.isfinite(Xr[:, 0])
    # deposits grouped by own G4 id (otid) and by stored id (tid)
    def group_match(keys_dep, keys_pt, sel_pt, bit):
        """for the selected points, nearest deposit sharing the key"""
        if not sel_pt.any():
            return
        kd = keys_dep.copy(); kd[~okdep] = -1
        od = np.argsort(kd, kind="stable"); kds = kd[od]
        kp = keys_pt[sel_pt]; ip = np.flatnonzero(sel_pt)
        op = np.argsort(kp, kind="stable"); kps = kp[op]; ips = ip[op]
        bp = np.r_[0, np.flatnonzero(np.diff(kps)) + 1, len(kps)]
        for a, b in zip(bp[:-1], bp[1:]):
            k = kps[a]
            lo, hi = np.searchsorted(kds, k), np.searchsorted(kds, k, side="right")
            if hi <= lo:
                continue
            dd = od[lo:hi]; pp = ips[a:b]
            d, j = cKDTree(Xr[dd]).query(X[pp])
            dep[pp] = dd[j]; mdist[pp] = d; reason[pp] |= bit
    group_match(ev.otid, otid, lab & (dep < 0), 0)
    group_match(ev.tid, tid, lab & (dep < 0), OR_FALLBACK_TID)
    rest = lab & (dep < 0)
    if rest.any() and okdep.any():
        d, j = cKDTree(Xr[okdep]).query(X[rest]); dd = np.flatnonzero(okdep)[j]
        dep[rest] = dd; mdist[rest] = d; reason[rest] |= OR_FALLBACK_ANY
    reason[lab & (dep < 0)] |= OR_NODEPO
    # copy the deposit labels
    has = dep >= 0
    out = dict(inst=np.full(npt, -1, np.int64), cls_obj=np.full(npt, OBJ["other"], np.int8), cls7=np.full(npt, CLS7["other"], np.int8),
               piece_id=np.full(npt, -1, np.int64), piece_type=np.zeros(npt, np.int8), d_own_seg=np.full(npt, -1.0, np.float32),
               overlap=np.zeros(npt, np.int16), stretch=np.zeros(npt, np.int8), dep=dep, match_dist=np.where(np.isfinite(mdist), mdist, -1.0).astype(np.float32),
               other_reason=reason, ghost_dist=np.full(npt, -1.0, np.float32), lab_d_traj=d_traj.astype(np.float32), lab_tid=tid, lab_orig_tid=otid,
               lab_cls=pts["cls"], lab_pf_node=pts["pf_node"], cl=pts["cl"])
    for k in ("inst", "cls_obj", "cls7", "piece_id", "piece_type", "d_own_seg", "overlap", "stretch"):
        src = getattr(R, k)
        out[k][has] = src[dep[has]]
    # ghost, shared, no truth, far
    g = tid == -1
    out["cls_obj"][g] = OBJ["ghost"]; out["cls7"][g] = CLS7["ghost"]; reason[g] |= OR_GHOST
    sh = tid == -2; nt = tid < -2
    reason[sh] |= OR_SHARED; reason[nt] |= OR_NOTRUTH
    far = lab & (d_traj >= FAR_MM)
    reason[far] |= OR_FAR
    if g.any() and (~g & lab).any():
        out["ghost_dist"][g] = cKDTree(X[~g & lab]).query(X[g])[0]
    # the point's own instance holds less than one hit of charge (a nuclear recoil, a sub-hit photon): v1 gave the point to it by the
    # nearest-deposit refinement; flagged, label kept (doc 54 sec 4, to rule)
    from .rules import FL_VISIBLE
    fl = np.asarray(R.table["inst_flags"]); io_ = out["inst"]
    inv = (io_ >= 0) & ((fl[np.maximum(io_, 0)] & FL_VISIBLE) == 0)
    reason[inv] |= OR_INVISIBLE_OWNER
    # Q16 b (Xin 2026-11-05): such a point is re-owned to the nearest deposit of a visible instance (>= 2,400 e); the flag stays
    vis_inst = (fl & FL_VISIBLE) > 0
    vd = okdep & (R.inst >= 0) & vis_inst[np.maximum(R.inst, 0)]
    if inv.any() and vd.any():
        d, j = cKDTree(Xr[vd]).query(X[inv]); dd = np.flatnonzero(vd)[j]
        dep[inv] = dd; out["match_dist"][inv] = d.astype(np.float32)
        for k in ("inst", "cls_obj", "cls7", "piece_id", "piece_type", "d_own_seg", "overlap", "stretch"):
            out[k][inv] = getattr(R, k)[dd]
    out["dep"] = dep
    out["other_reason"] = reason
    return out
