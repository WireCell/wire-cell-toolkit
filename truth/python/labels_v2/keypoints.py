"""labels_v2.keypoints -- G5, typed key points per instance (50b sec 2.7 G5, G7; 50i; Q15).

  KP_TRACK_END   1  start / end of a Track, when physical: the particle's first / last own deposit lies within 10 mm of its
                    G4 start / end (it did not enter or leave the active volume there; a stop in a dead region reads as
                    physical, 50i caveat), not inside a shower-hidden stretch, and the track is not buried
  KP_STUB_START  2  the start of a Stub (physical), not swallowed, not buried; only instances with >= THR of charge (visible)
                    give key points or join a vertex
  KP_SHOWER      3  the start of a shower: the first-in-time deposit of its earliest segment
  KP_MICHEL      4  the first own deposit of a resolved Michel (hidden start: none)
  KP_DELTA       5  the attachment point of a resolved delta (its first own deposit; hidden start: none)
  KP_NU_VERTEX   6  the neutrino vertex: median of the primaries' first deposits per interaction row with a neutrino pdg (cosmic
                    rain and particle-gun rows carry pdg 0: none), where >= 2 primaries start within 10 mm of it (or one primary)
  Every key point's deposit must be seen (inside the active volume and the readout window, in a live pixel): a boundary
  point is not a target (G7).  A delta that is a shower gets its attachment point and, when > 10 mm away, its shower start.
  KP_VERTEX      7  a secondary vertex: >= 2 visible instances (tracks, stubs, showers, deltas, Michels) start within 10 mm
Positions are the reco-frame positions of deposits (x from the arrival time on the face), the frame of the PR points.
Per point: kp_type = the type of the nearest key point within radius_mm (10), kp_dist its distance.
"""
import numpy as np
from scipy.spatial import cKDTree

from .objects import KIND
from .rules import OBJ, FL_BURIED, FL_SWALLOWED, FL_HIDDEN_START, FL_HIDDEN_END, FL_HIDDEN_ATTACH, FL_RESOLVED, FL_VISIBLE, FL_UNRESOLVED

KP = {"track_end": 1, "stub_start": 2, "shower_start": 3, "michel_start": 4, "delta_attach": 5, "nu_vertex": 6, "vertex": 7}
PHYS_MM = 10.0
RADIUS_MM = 10.0
VERTEX_MM = 10.0
KINK_W_MM = 30.0
KINK_DEG = 15.0


def kinks(ev, R, table, w_mm=KINK_W_MM, deg=KINK_DEG):
    """Hard scatters of ONE particle (pr_ml doc 71 sec 6: a segment boundary allowed there, not required).  Only called with the
    option kink_points; its members are written beside the key-point table, which is untouched (no new kp_kind, pt_kp_* as before).

    Per visible, resolved instance of class track: its own deposits with a true position, ordered in time, equal-time runs ordered
    along the track (steiner_truth.tie_order; the depo file stores float32 times), zero-length steps dropped, are its trajectory.  At
    every vertex the angle between the chord from the point w_mm behind to it and from it to the point w_mm ahead; a kink = a local
    maximum of that angle within +- w_mm at or above deg, with both chords at their full length (at least w_mm of trajectory on
    either side), placed at the nearest own deposit within RADIUS_MM along the trajectory that is seen (G7) and has a reco
    position (the kink's own deposit when it is seen); none: no kink.
    -> kink_x / kink_y / kink_z (reco frame, mm), kink_inst, kink_angle (deg), kink_to_end (mm of trajectory after the kink)."""
    from .rules import tie_order
    T = R.T; Xr = ev.Xreco; o_r, st_r = R.o_r, R.st_r; seen = R.P.seen.any(1); out = []
    for i in range(len(table["inst_q"])):
        r = int(table["inst_root_row"][i])
        if r < 0 or T.kind[r] != KIND["TRACK"] or int(table["inst_cls_obj"][i]) != OBJ["track"]:
            continue
        fl = int(table["inst_flags"][i])
        if not (fl & FL_VISIBLE) or fl & (FL_BURIED | FL_SWALLOWED | FL_UNRESOLVED):
            continue
        mem = o_r[st_r[r]:st_r[r + 1]]; own = mem[(ev.otid[mem] == T.tt[r]) & np.isfinite(ev.P3[mem]).all(1)]
        if len(own) < 3:
            continue
        own = own[np.argsort(ev.T0[own], kind="stable")]; own = own[tie_order(ev.T0[own], ev.P3[own])]; V = ev.P3[own]
        keep = np.r_[True, np.linalg.norm(V[1:] - V[:-1], axis=1) > 0]; own, V = own[keep], V[keep]; n = len(V)
        if n < 3:
            continue
        s = np.r_[0.0, np.cumsum(np.linalg.norm(V[1:] - V[:-1], axis=1))]; m = np.arange(1, n - 1)
        if s[-1] < 2 * w_mm:
            continue
        j = np.clip(np.searchsorted(s, s[m] - w_mm, side="right") - 1, 0, None); k = np.clip(np.searchsorted(s, s[m] + w_mm, side="left"), None, n - 1)
        b, a = V[m] - V[j], V[k] - V[m]; nb, na = np.linalg.norm(b, axis=1), np.linalg.norm(a, axis=1)
        ok = (s[m] >= w_mm) & (s[-1] - s[m] >= w_mm) & (nb > 0) & (na > 0)
        ang = np.zeros(n); ang[m[ok]] = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", a, b)[ok] / (na * nb)[ok], -1, 1)))
        for c in np.flatnonzero(ang >= deg).tolist():
            w = np.arange(np.searchsorted(s, s[c] - w_mm, side="left"), np.searchsorted(s, s[c] + w_mm, side="right")); w = w[w != c]
            if np.any(ang[w] > ang[c]) or np.any((ang[w] == ang[c]) & (w < c)):
                continue
            near = np.arange(np.searchsorted(s, s[c] - RADIUS_MM, side="left"), np.searchsorted(s, s[c] + RADIUS_MM, side="right"))
            near = near[seen[own[near]] & np.isfinite(Xr[own[near]]).all(1)]      # the kink's own deposit may sit in an unseen pixel of a seen track:
            if len(near):                                                        # take the nearest seen own deposit within RADIUS_MM along the trajectory
                d = int(own[near[np.argmin(np.abs(s[near] - s[c]))]])
                out.append((Xr[d, 0], Xr[d, 1], Xr[d, 2], i, ang[c], s[-1] - s[c]))
    A = np.array(out, np.float64).reshape(-1, 6)
    return dict(kink_x=A[:, 0].astype(np.float32), kink_y=A[:, 1].astype(np.float32), kink_z=A[:, 2].astype(np.float32),
                kink_inst=A[:, 3].astype(np.int32), kink_angle=A[:, 4].astype(np.float32), kink_to_end=A[:, 5].astype(np.float32))


def build(ev, R, table):
    T = R.T; nr = T.nrow
    Xr = ev.Xreco
    o_r, st_r = R.o_r, R.st_r
    kp_x, kp_kind, kp_inst = [], [], []
    starts = []                                              # (position, instance) of every visible instance start, for vertices
    inst_of_root = {int(r): i for i, r in enumerate(R.inst_root.tolist())}

    seen = R.P.seen.any(1)

    def own_sorted(r):
        mem = o_r[st_r[r]:st_r[r + 1]]
        own = mem[(ev.otid[mem] == T.tt[r]) & ev.placed[mem] & np.isfinite(Xr[mem, 0])]
        if not len(own):
            own = mem[np.isfinite(Xr[mem, 0])]
        if not len(own):
            return own
        return own[np.argsort(ev.T0[own], kind="stable")]

    def vis(j):
        """a key point is a target only where its deposit is seen (inside the active volume and the readout, in a live pixel): G7"""
        return j >= 0 and bool(seen[j])

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
            if fl & (FL_BURIED | FL_SWALLOWED | FL_UNRESOLVED):      # never resolved: no key points (Q17)
                continue
            if cls == OBJ["track"]:
                if ps and not (fl & FL_HIDDEN_START) and vis(own[0]):
                    add(p0, KP["track_end"], i); starts.append((p0, i))
                if pe and not (fl & FL_HIDDEN_END) and vis(own[-1]):
                    add(p1, KP["track_end"], i)
            elif cls == OBJ["stub"]:
                if ps and vis(own[0]):
                    add(p0, KP["stub_start"], i); starts.append((p0, i))
        elif k in (KIND["DELTA"], KIND["MICHEL"]):
            if not (fl & FL_RESOLVED):
                continue
            own = own_sorted(r)
            if not len(own):
                continue
            p0 = Xr[own[0]]
            if not (fl & FL_HIDDEN_ATTACH) and vis(own[0]):          # the attachment point (delta) / start (Michel), also for a delta that is a shower
                add(p0, KP["michel_start"] if k == KIND["MICHEL"] else KP["delta_attach"], i)
            starts.append((p0, i))
            if cls == OBJ["shower"]:                        # a resolved delta with >= 2 segments: also its shower start when distinct
                e = R.eidx[r]; j = R.seg_start_dep[e] if e >= 0 else -1
                if j >= 0 and np.isfinite(Xr[j, 0]) and vis(j) and np.linalg.norm(Xr[j] - p0) > RADIUS_MM and not (fl & FL_HIDDEN_ATTACH):
                    add(Xr[j], KP["shower_start"], i)
        elif k == KIND["EM"]:
            if cls != OBJ["shower"]:
                continue
            e = R.eidx[r]; j = R.seg_start_dep[e] if e >= 0 else -1
            if j < 0 or not np.isfinite(Xr[j, 0]) or not vis(j):     # a shower starting outside the active volume / readout: a boundary point, no target (G7)
                continue
            add(Xr[j], KP["shower_start"], i); starts.append((Xr[j], i))
    # neutrino vertices: median of the primaries' first deposits per interaction (all truth_nu rows; as steiner_truth)
    # only interaction rows with a neutrino pdg (cosmic-rain and particle-gun rows carry pdg 0), and only where >= 2 primaries
    # start together (within VERTEX_MM of the median): a scattered "rain" group gives no vertex
    import truth_table as TT
    nu_tab = ev.truth.get("truth_nu", np.zeros((0, len(TT.NU_COLUMNS))))
    nu_pdg = {int(row[TT.nu_col("nu_idx")]) if TT.NU_COLUMNS[0] == "nu_idx" else k_: int(abs(row[TT.nu_col("pdg")])) for k_, row in enumerate(np.atleast_2d(nu_tab))} if len(nu_tab) else {}
    nus = sorted(set(int(v) for v in T.nu_row if v >= 0))
    for nu in nus:
        pdg_nu = nu_pdg.get(nu, nu_pdg.get(int(nu), 0)) if nu_pdg else 0
        if len(nu_tab) and 0 <= nu < len(np.atleast_2d(nu_tab)):
            pdg_nu = int(abs(np.atleast_2d(nu_tab)[nu, TT.nu_col("pdg")]))
        if pdg_nu not in (12, 14, 16):
            continue
        prim = np.flatnonzero((T.nu_row == nu) & (T.mot == 0))
        S = []
        for r in prim.tolist():
            own = own_sorted(r) if st_r[r + 1] > st_r[r] else np.zeros(0, np.int64)
            if len(own) and seen[own[0]]:
                S.append(Xr[own[0]])
        if len(S) >= 1:
            c = np.median(np.array(S), axis=0)
            near = sum(1 for p in S if np.linalg.norm(p - c) <= VERTEX_MM)
            if near >= 2 or len(S) == 1:
                add(c, KP["nu_vertex"], -1)
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
