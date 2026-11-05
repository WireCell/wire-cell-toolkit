"""labels_v2.keypoints -- G5, typed key points per instance (50b sec 2.7 G5, G7; 50i; Q15).

  KP_TRACK_END   1  start / end of a Track, when physical: the particle's first / last own deposit lies within 10 mm of its
                    G4 start / end (it did not enter or leave the active volume there; a stop in a dead region reads as
                    physical, 50i caveat), not inside a shower-hidden stretch, and the track is not buried
  KP_STUB_START  2  the start of a Stub (physical), not swallowed, not buried; only instances with >= THR of charge (visible)
                    give key points or join a vertex
  KP_SHOWER      3  the start of a shower: the first-in-time deposit of its earliest segment
  KP_MICHEL      4  the first own deposit of a resolved Michel (hidden start: none)
  KP_DELTA       5  the attachment point of a resolved delta (its first own deposit; hidden start: none)
  KP_NU_VERTEX   6  the neutrino vertex: median of the primaries' first deposits per interaction
  KP_VERTEX      7  a secondary vertex: >= 2 visible instances (tracks, stubs, showers, deltas, Michels) start within 10 mm
Positions are the reco-frame positions of deposits (x from the arrival time on the face), the frame of the PR points.
Per point: kp_type = the type of the nearest key point within radius_mm (10), kp_dist its distance.
"""
import numpy as np
from scipy.spatial import cKDTree

from .objects import KIND
from .rules import OBJ, FL_BURIED, FL_SWALLOWED, FL_HIDDEN_START, FL_HIDDEN_END, FL_HIDDEN_ATTACH, FL_RESOLVED, FL_VISIBLE

KP = {"track_end": 1, "stub_start": 2, "shower_start": 3, "michel_start": 4, "delta_attach": 5, "nu_vertex": 6, "vertex": 7}
PHYS_MM = 10.0
RADIUS_MM = 10.0
VERTEX_MM = 10.0


def build(ev, R, table):
    T = R.T; nr = T.nrow
    Xr = ev.Xreco
    o_r, st_r = R.o_r, R.st_r
    kp_x, kp_kind, kp_inst = [], [], []
    starts = []                                              # (position, instance) of every visible instance start, for vertices
    inst_of_root = {int(r): i for i, r in enumerate(R.inst_root.tolist())}

    def own_sorted(r):
        mem = o_r[st_r[r]:st_r[r + 1]]
        own = mem[(ev.otid[mem] == T.tt[r]) & ev.placed[mem] & np.isfinite(Xr[mem, 0])]
        if not len(own):
            own = mem[np.isfinite(Xr[mem, 0])]
        if not len(own):
            return own
        return own[np.argsort(ev.T0[own], kind="stable")]

    def add(p, kind, inst):
        kp_x.append(np.asarray(p, np.float64)); kp_kind.append(kind); kp_inst.append(inst)

    NI = len(table["inst_q"])
    for i in range(NI):
        r = int(table["inst_root_row"][i])
        if r < 0:
            continue
        cls = int(table["inst_cls_obj"][i]); fl = int(table["inst_flags"][i]); k = T.kind[r]
        if not (fl & FL_VISIBLE):                        # below one hit of charge: no key point, no vertex membership
            continue
        if k == KIND["TRACK"]:
            own = own_sorted(r)
            if not len(own):
                continue
            p0, p1 = Xr[own[0]], Xr[own[-1]]
            ps = np.linalg.norm(ev.P3[own[0]] - T.S0[r]) <= PHYS_MM
            pe = np.linalg.norm(ev.P3[own[-1]] - T.S1[r]) <= PHYS_MM
            if fl & (FL_BURIED | FL_SWALLOWED):
                continue
            if cls == OBJ["track"]:
                if ps and not (fl & FL_HIDDEN_START):
                    add(p0, KP["track_end"], i); starts.append((p0, i))
                if pe and not (fl & FL_HIDDEN_END):
                    add(p1, KP["track_end"], i)
            elif cls == OBJ["stub"]:
                if ps:
                    add(p0, KP["stub_start"], i); starts.append((p0, i))
        elif k in (KIND["DELTA"], KIND["MICHEL"]):
            if not (fl & FL_RESOLVED):
                continue
            own = own_sorted(r)
            if not len(own):
                continue
            p0 = Xr[own[0]]
            if cls == OBJ["shower"]:                        # a resolved delta with >= 2 segments: a shower start
                e = R.eidx[r]; j = R.seg_start_dep[e] if e >= 0 else -1
                pj = Xr[j] if j >= 0 and np.isfinite(Xr[j, 0]) else p0
                if not (fl & FL_HIDDEN_ATTACH):
                    add(pj, KP["shower_start"], i)
                starts.append((pj, i))
            else:
                if not (fl & FL_HIDDEN_ATTACH):
                    add(p0, KP["michel_start"] if k == KIND["MICHEL"] else KP["delta_attach"], i)
                starts.append((p0, i))
        elif k == KIND["EM"]:
            if cls != OBJ["shower"]:
                continue
            e = R.eidx[r]; j = R.seg_start_dep[e] if e >= 0 else -1
            if j < 0 or not np.isfinite(Xr[j, 0]):
                continue
            add(Xr[j], KP["shower_start"], i); starts.append((Xr[j], i))
    # neutrino vertices: median of the primaries' first deposits per interaction (all truth_nu rows; as steiner_truth)
    nus = sorted(set(int(v) for v in T.nu_row if v >= 0))
    for nu in nus:
        prim = np.flatnonzero((T.nu_row == nu) & (T.mot == 0))
        S = []
        for r in prim.tolist():
            own = own_sorted(r) if st_r[r + 1] > st_r[r] else np.zeros(0, np.int64)
            if len(own):
                S.append(Xr[own[0]])
        if S:
            add(np.median(np.array(S), axis=0), KP["nu_vertex"], -1)
    # secondary vertices: >= 2 instance starts within VERTEX_MM (single linkage), away from a nu vertex
    if len(starts) >= 2:
        SP = np.array([s[0] for s in starts]); SI = np.array([s[1] for s in starts])
        tree = cKDTree(SP); pairs = tree.query_pairs(VERTEX_MM, output_type="ndarray")
        if len(pairs):
            from scipy.sparse import coo_matrix
            from scipy.sparse.csgraph import connected_components
            A = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(SP), len(SP)))
            _, lab = connected_components(A, directed=False)
            nuv = np.array([x for x, k in zip(kp_x, kp_kind) if k == KP["nu_vertex"]]) if kp_x else np.zeros((0, 3))
            for c in np.unique(lab):
                m = np.flatnonzero(lab == c)
                if len(np.unique(SI[m])) < 2:
                    continue
                p = SP[m].mean(0)
                if len(nuv) and np.linalg.norm(nuv - p, axis=1).min() <= VERTEX_MM:
                    continue
                add(p, KP["vertex"], -1)
    KX = np.array(kp_x) if kp_x else np.zeros((0, 3))
    KK = np.array(kp_kind, np.int8) if kp_x else np.zeros(0, np.int8)
    KI = np.array(kp_inst, np.int64) if kp_x else np.zeros(0, np.int64)
    # per point
    X = ev.pts["X"]; npt = len(X)
    kp_type = np.zeros(npt, np.int8); kp_dist = np.full(npt, -1.0, np.float32); kp_index = np.full(npt, -1, np.int64)
    if len(KX) and npt:
        d, j = cKDTree(KX).query(X)
        near = d <= RADIUS_MM
        kp_type[near] = KK[j[near]]; kp_dist[:] = d; kp_index[near] = j[near]
    return dict(kp_x=KX[:, 0] if len(KX) else np.zeros(0), kp_y=KX[:, 1] if len(KX) else np.zeros(0), kp_z=KX[:, 2] if len(KX) else np.zeros(0),
                kp_kind=KK, kp_inst=KI), dict(kp_type=kp_type, kp_dist=kp_dist, kp_index=kp_index)
