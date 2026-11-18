"""labels_v2.c7 -- the E-1 target C7 (pr_ml docs 85 / 85a / 86; Xin's rulings of 2026-11-16 and 2026-11-17): label what is
visible, with key-point PLACES and visible PIECES.  Called only with the option c7 (which implies c6), after the points are
linked to their deposits; the rules and thresholds are those of the view C7v5 (pr_ml doc 86, 86_results/scripts/d86_c7.py),
with the labeller's own data where the view had to approximate.

Per counted point (class above ghost, an instance, not far, a matched deposit):
  S1  a stub whose free run (inst_run_mm) is >= 5 mm stays a stub: own piece, a place at its start
  S2  any other stub: its covered points (overlap bit 1) take the class and the piece of the HOST, flag 1.  An object is folded
      whole into ONE host piece (the piece most of its covered points name, when at least half of them name one), or not at all
  S3  ... its free points become EM blip, flag 8, in a piece outside the graph
  U1  a track flagged unresolved (and not buried): every point takes the host's class and piece, flag 2
  H1  a hidden END of a track (stretch 2; runs with < 4 steps of seen path between them are one) longer than 8 steps: the points
      in that stretch take class EM segment and the shower's piece, flag 4; the track ends there at a shower-boundary place.  A
      traversed stretch (stretch 1) stays track, flag 32 (Xin 2026-11-17: a through-going track is one piece)
  H3  a buried track: every point takes class EM segment and the shower's piece, flag 4
  D2  (Xin 2026-11-17, doc 86a D2; replaces U1, H1 and H3) none of U1, H1, H3 is applied: an unresolved track, a hidden start or end
      and a buried track stay track in their own pieces (flags 128 unresolved kept, 64 hidden start / end or buried kept)
  P1  the rest of a track: its points fall into RUNS, the single-linkage groups at G = 30 mm in the image (a gap of more than
      30 mm with no point of the instance ends a piece, as in the view; a frame change of an out-of-time track at the cathode is
      such a gap); in a run the points are ordered by the path length over the particle's own deposits in time order; every run
      is one visible piece
  M1  a track piece with a seen length under 2 steps is no piece: covered points to the host (flag 2), free points EM blip (flag 8)
  E1  a delta, Michel or shower instance with fewer than 5 free points, or whose free points span less than 2 steps, keeps its
      class but is outside the graph and gives no place
  Z1  flag 16: a point of a track piece within 15 mm of another track piece's centre line (the instance is not scored there)
The HOST of a covered point is the instance the overlap test (Projection.clean_share) found in the ring of its deposit's main
pixel, per view; of a hidden end / a buried track the shower instance of the shower-ring test.  The named instance is taken only
when it has an eligible point within 6 mm of the point in the image (a shower: 30 mm): the test is made per wire view, and an
instance that shares a pixel in one view can be far away in space.  Otherwise: the nearest eligible point within 6 mm (a shower:
10 mm, then 30 mm); none: no host (a stub's or a short piece's point becomes EM blip, an unresolved track's point stays track in
a piece outside the graph).

Candidates of key points, merged into places:
  K1  both ends of every visible track piece.  The seen-deposit end is the reco position of its first / last matched own
      deposits (median over the deposits within 5 mm of path of the extreme one); for an end tagged physical the candidate is
      the reco position of the particle's first / last own deposit, seen or not (the true end of its charge, as today's track-end
      key point).  Candidates are merged, and a place that several candidates share sits, at these positions.  The TIP of an end is where the piece's points stop: the outward
      direction runs from the piece's own deposits 10 to 30 mm of path inside to the seen-deposit end; of the piece's points
      matched to the last 5 mm of path (a point matched to another particle of the instance counts with its nearest own
      deposit) and within G = 30 mm of the seen-deposit end, those beyond it are ranked by their distance from it, and the tip
      is the median of the points within 5 mm of the farthest; with no point beyond it, the median of the points within 5 mm
      (along the direction) of the most outward one.  A place whose only member is one track-piece end (a free end) sits at
      the tip; pc_e0 / pc_e1 are the tips.  (Image points reach past the last deposit, and an absorbed electron or a short
      daughter can prolong the visible track; the blind reads of doc 86, 2026-11-17.  At a vertex the tips of the tracks
      spread over the vertex blob, so a shared place stays at the deposits.)  pc_len is the path between the first and the
      last seen own deposit.  Tag 1 physical (within 10 mm of path of the particle's true start / end),
      4 shower boundary (H1), 2 edge (the particle's next own deposits beyond the end are outside the active volume, or there are
      none), 8 break (the rest)
  K3 / K4 / K5  today's key points of a kept stub (start), a shower (start), a resolved delta (attach) and Michel (start), when the
      object is in the graph
  merge: single linkage at 10 mm; and ONE PLACE PER TRUE VERTEX: two candidates whose true vertex positions (the particle's G4
  start for a start-side end and for a stub / Michel / delta start, its G4 end for an end-side end, the true position of a
  shower's first segment deposit) lie within 5 mm, and whose reco positions lie within 30 mm, are one place.  A piece end is start
  (end) side when the path from the particle's true start (end) to it is at most 30 mm.
The implied graph: a track piece is a chain of segments through the places at its ends and those within 6 mm of its centre line;
a stub, shower, Michel or delta hangs on its place(s); a shower is also joined to the boundary places of the hidden ends it took.

Kinds of a candidate: 1 track-piece end, 2 stub start, 3 shower start, 4 Michel start, 5 delta attach.
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from steiner_truth import tie_order

from .rules import OBJ, CLS7, STEP, FL_BURIED, FL_UNRESOLVED, OV_NOTFREE, STR_TRAVERSED, STR_HIDDEN_END
from .points import OR_FAR

FL_FSTUB, FL_FTRK, FL_HID, FL_LE, FL_ZONE, FL_THRU, FL_HIDKEEP, FL_UNRKEEP = 1, 2, 4, 8, 16, 32, 64, 128
TG_PHYS, TG_EDGE, TG_SHB, TG_BREAK, TG_VIS = 1, 2, 4, 8, 16
PRM = dict(L_stub=5.0, L_hid_steps=8.0, hid_join_steps=4.0, G=30.0, link=12.0, R_merge=10.0, R_zone=15.0, host_r=6.0, sh_r=10.0, end_w=5.0, L_trk_steps=2.0,
           free_n=5, free_ext_steps=2.0, vertex_true=5.0, vertex_reco=30.0, phys=10.0, pass_r=6.0, bw=10.0, min_pts=3, end_in=(10.0, 30.0),
           tip_nmin=3, tip_axis=(5.0, 25.0))
NEW_MEMBERS = ("pt_c7_piece", "pt_c7_flag", "pt_c7_place_dist", "pc_type", "pc_inst", "pc_link", "pc_graph", "pc_len", "pc_e0", "pc_e1", "pc_place0", "pc_place1",
               "pl_xyz", "pl_kinds", "pl_tag", "pl_member_place", "pl_member_piece", "seg_a", "seg_b", "seg_piece")


def _resample(C, step=2.0):
    if len(C) < 2:
        return C.copy()
    arc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))]
    if arc[-1] <= 0:
        return C[:1].copy()
    u = np.linspace(0.0, arc[-1], max(int(np.ceil(arc[-1] / step)) + 1, 2))
    return np.column_stack([np.interp(u, arc, C[:, k]) for k in range(3)])


def runs_at(Xp, G, link):
    """the single-linkage groups of the points at G (mm): the groups of a graph with links up to `link`, joined nearest group
    first from the group of the first point; a join longer than G starts a new run (as 86_results/scripts/d86_c7.py along).
    -> run id per point"""
    n = len(Xp); run = np.zeros(n, np.int64)
    if n < 2:
        return run
    T = cKDTree(Xp); M = T.sparse_distance_matrix(T, min(link, G), output_type="coo_matrix")
    nc, lab = connected_components(coo_matrix((np.ones(len(M.row)), (M.row, M.col)), shape=(n, n)), directed=False)
    if nc == 1:
        return run
    run[:] = -1; done = lab == lab[0]; run[done] = 0; nrun = 1
    while not done.all():
        rest = np.flatnonzero(~done); dn = np.flatnonzero(done); d_, j = cKDTree(Xp[dn]).query(Xp[rest]); k = int(np.argmin(d_)); g = lab == lab[rest[k]]
        if d_[k] > G:
            run[g] = nrun; nrun += 1
        else:
            run[g] = run[dn[j[k]]]
        done |= g
    return run


def build(ev, R, pt, kps, det):
    """-> (dict of the new members, log dict); pt["cls7"] is changed in place"""
    V = PRM; step = STEP[det]; T = R.T; nr = T.nrow; tab = R.table
    X = ev.pts["X"]; npt = len(X); Xr = ev.Xreco; P3 = ev.P3
    cls = pt["cls7"]; inst = pt["inst"]; dep = pt["dep"]; ov = pt["overlap"].astype(np.int64); st = pt["stretch"].astype(np.int64)
    co = np.asarray(tab["inst_cls_obj"]).astype(np.int64); fl = np.asarray(tab["inst_flags"]).astype(np.int64); ni = len(co)
    run_mm = np.asarray(tab["inst_run_mm"], np.float64); g4 = np.asarray(tab["inst_g4_id"]).astype(np.int64); root = np.asarray(R.inst_root).astype(np.int64)
    free = (ov & OV_NOTFREE) == 0
    cnt = (cls >= 1) & (inst >= 0) & (dep >= 0) & ((pt["other_reason"].astype(np.int64) & OR_FAR) == 0)
    piece = np.full(npt, -1, np.int64); flag = np.zeros(npt, np.int64); rule = np.zeros(npt, np.int64)
    # the points of every instance
    order = np.argsort(np.where(cnt, inst, ni), kind="stable"); ii = np.where(cnt, inst, ni)[order]
    ui, s0_ = np.unique(ii, return_index=True); e0_ = np.r_[s0_[1:], len(ii)]
    rows_of = {int(u): order[a:b] for u, a, b in zip(ui.tolist(), s0_.tolist(), e0_.tolist()) if u < ni}
    insts = sorted(rows_of)
    stub_keep = [i for i in insts if co[i] == OBJ["stub"] and run_mm[i] >= V["L_stub"]]
    stub_fold = [i for i in insts if co[i] == OBJ["stub"] and run_mm[i] < V["L_stub"]]
    buried_kept = [i for i in insts if co[i] == OBJ["track"] and (fl[i] & FL_BURIED)]
    unres_kept = [i for i in insts if co[i] == OBJ["track"] and (fl[i] & FL_UNRESOLVED) and not (fl[i] & FL_BURIED)]
    buried, unres = [], []                                # D2 (Xin 2026-11-17, doc 86a): no H3 and no U1 -- a buried or an unresolved track is a track like any other
    tracks = [i for i in insts if co[i] == OBJ["track"]]
    others = [i for i in insts if co[i] not in (OBJ["track"], OBJ["stub"])]
    Pk = dict(kind=[], inst=[], n=[], length=[], e0=[], e1=[], graph=[], tag0=[], tag1=[], vtx0=[], vtx1=[], cl=[])
    nan3 = np.full(3, np.nan)

    def new_piece(kind, i, rows, length=0.0, e0=nan3, e1=nan3, graph=True, tag0=0, tag1=0, vtx0=nan3, vtx1=nan3, cl=None):
        k = len(Pk["kind"])
        for key, v in (("kind", kind), ("inst", i), ("n", len(rows)), ("length", length), ("e0", e0), ("e1", e1), ("graph", graph), ("tag0", tag0), ("tag1", tag1),
                       ("vtx0", vtx0), ("vtx1", vtx1), ("cl", cl)):
            Pk[key].append(v)
        piece[rows] = k
        return k

    # ---- tracks: the path over the own deposits, the hidden ends (H1), the visible pieces (P1), the short pieces (M1)
    o_r, st_r = R.o_r, R.st_r
    to_shower = []; short = []                            # (rows, instance, rule)
    TIP = {}                                              # per graph track piece: where its points stop at either end
    arc = np.full(ev.n, np.nan)                           # path length of every own deposit of a track instance
    for i in tracks:
        r = int(root[i])
        if r >= nr:
            new_piece(OBJ["track"], i, rows_of[i], graph=False); continue
        mem = o_r[st_r[r]:st_r[r + 1]]
        own = mem[(ev.otid[mem] == T.tt[r]) & ev.placed[mem] & np.isfinite(Xr[mem, 0])]
        P = rows_of[i]
        if len(own) < 2:
            new_piece(OBJ["track"], i, P, graph=False); continue
        own = own[np.argsort(ev.T0[own], kind="stable")]; own = own[tie_order(ev.T0[own], P3[own])]
        sseg = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(P3[own], axis=0), axis=1))]; arc[own] = sseg
        pre = float(np.linalg.norm(P3[own[0]] - T.S0[r])); post = float(np.linalg.norm(P3[own[-1]] - T.S1[r])); stot = float(sseg[-1])
        s = arc[dep[P]]; m_ = ~np.isfinite(s)
        if m_.any():                                      # a point matched to a deposit of the instance that is not the particle's own
            s[m_] = sseg[cKDTree(Xr[own]).query(X[P[m_]])[1]]
        rn = runs_at(X[P], V["G"], V["link"])            # the runs: no point of the instance within G in the image ends a piece
        first = np.full(int(rn.max()) + 1, np.inf); np.minimum.at(first, rn, s)
        rn = np.argsort(np.argsort(first, kind="stable"), kind="stable")[rn]          # numbered by their first path coordinate
        o = np.lexsort((s, rn)); P, s, runid = P[o], s[o], rn[o]                       # run by run, by the path in a run
        taken = np.zeros(len(P), bool); cuts = []
        hi = np.zeros(0, np.int64)                        # D2 (Xin 2026-11-17): no H1 -- a hidden start or end stays in its track's piece
        if len(hi):
            b = np.r_[0, np.flatnonzero((np.diff(s[hi]) > V["hid_join_steps"] * step) | (np.diff(runid[hi]) != 0)) + 1, len(hi)]
            for a_, b_ in zip(b[:-1], b[1:]):
                lo, up = s[hi[a_]], s[hi[b_ - 1]]
                if up - lo > V["L_hid_steps"] * step:
                    taken |= (s >= lo) & (s <= up) & (runid == runid[hi[a_]]); cuts.append((lo, up, int(runid[hi[a_]])))
        if taken.any():
            to_shower.append((P[taken], i, 4))
        keep = np.flatnonzero(~taken)
        if not len(keep):
            continue
        sk = s[keep]; rk = runid[keep]
        nb = (np.diff(rk) != 0) | np.array([any(sk[j] < lo and sk[j + 1] > up and rk[j] == rc_ for lo, up, rc_ in cuts) for j in range(len(sk) - 1)], bool) if len(sk) > 1 else np.zeros(0, bool)
        b = np.r_[0, np.flatnonzero(nb) + 1, len(sk)]
        ownset = np.zeros(ev.n, bool); ownset[own] = True
        for a_, b_ in zip(b[:-1], b[1:]):
            rr = P[keep[a_:b_]]; ss = sk[a_:b_]
            du = np.unique(dep[rr]); du = du[ownset[du]]
            if len(du):
                sa = arc[du]; q = np.argsort(sa, kind="stable"); du, sa = du[q], sa[q]
            else:                                         # no own deposit matched: the points' coordinate
                du = np.zeros(0, np.int64); sa = np.array([ss[0], ss[-1]])
            lo_, up_ = float(sa[0]), float(sa[-1]); length = up_ - lo_
            e0 = np.median(Xr[du[sa <= lo_ + V["end_w"]]], axis=0) if len(du) else X[rr[0]].copy()
            e1 = np.median(Xr[du[sa >= up_ - V["end_w"]]], axis=0) if len(du) else X[rr[-1]].copy()
            if len(rr) < V["min_pts"] or length < V["L_trk_steps"] * step:
                short.append((rr, i)); continue
            if len(du):                                   # K1: the place of an end is where the piece's points stop
                def tip(d, inner, m):
                    """D1 (Xin 2026-11-17, doc 86a D1 (a), rule C of doc 86's tip diagnostic, detail C0): where the piece's dots stop
                    along the end's line, put on the centre line of its dots.  The dots: the piece's points within G of the end; the stop:
                    the median coordinate along the line of the last end_w of them (at least tip_nmin dots, so a lone dot past a gap is
                    averaged with the band's end); the centre line: the median offset across the line of the piece's points tip_axis
                    inside the end."""
                    u = d - inner; n = np.linalg.norm(u)
                    if n < 1e-6 or not m.any():
                        return d
                    u = u / n; Q = X[rr]; Q = Q[np.linalg.norm(Q - d, axis=1) <= V["G"]]
                    if not len(Q):
                        return d
                    t = (Q - d) @ u
                    def stop(tk):
                        sl = tk >= tk.max() - V["end_w"]
                        if sl.sum() < V["tip_nmin"]:
                            sl = np.zeros(len(tk), bool); sl[np.argsort(tk)[-V["tip_nmin"]:]] = True
                        return float(np.median(tk[sl]))
                    st_ = stop(t)                         # C0 (Xin 2026-11-17): no isolated-dot drop
                    b_ = (t >= -V["tip_axis"][1]) & (t <= -V["tip_axis"][0])
                    off = np.median((Q[b_] - d) - np.outer(t[b_], u), axis=0) if b_.sum() >= 3 else np.zeros(3)
                    return d + off + st_ * u
                i0 = du[(sa >= lo_ + V["end_in"][0]) & (sa <= lo_ + V["end_in"][1])]; i1 = du[(sa <= up_ - V["end_in"][0]) & (sa >= up_ - V["end_in"][1])]
                t0 = tip(e0, np.median(Xr[i0], axis=0) if len(i0) else e1, ss <= lo_ + V["end_w"])
                t1 = tip(e1, np.median(Xr[i1], axis=0) if len(i1) else e0, ss >= up_ - V["end_w"])
            else:
                t0, t1 = e0, e1
            h0 = any(0 <= ss[0] - up <= V["hid_join_steps"] * step and rk[a_] == rc_ for lo, up, rc_ in cuts); h1 = any(0 <= lo - ss[-1] <= V["hid_join_steps"] * step and rk[a_] == rc_ for lo, up, rc_ in cuts)
            def tag(side, hid):
                if hid:
                    return TG_SHB
                d_ = (pre + lo_) if side == 0 else (post + stot - up_)
                if d_ <= V["phys"]:
                    return TG_PHYS
                bey = own[(sseg < lo_) & (sseg >= lo_ - V["G"])] if side == 0 else own[(sseg > up_) & (sseg <= up_ + V["G"])]
                return TG_EDGE if (not len(bey) or not ev.inside[bey].any()) else TG_BREAK
            tg0, tg1 = tag(0, h0), tag(1, h1)
            if tg0 == TG_PHYS:                            # the particle's own start is seen: the candidate sits at its first own deposit, the true end of its charge
                e0 = Xr[own[0]].astype(np.float64)
            if tg1 == TG_PHYS:
                e1 = Xr[own[-1]].astype(np.float64)
            v0 = T.S0[r].astype(np.float64) if pre + lo_ <= V["vertex_reco"] else nan3
            v1 = T.S1[r].astype(np.float64) if post + stot - up_ <= V["vertex_reco"] else nan3
            if len(du) >= 2:                              # the centre line: the matched own deposits in bins of bw along the path
                bi = np.floor((sa - lo_) / V["bw"]).astype(np.int64); ub, inv = np.unique(bi, return_inverse=True)
                C = np.zeros((len(ub), 3)); np.add.at(C, inv, Xr[du]); C /= np.bincount(inv)[:, None]
                C = np.vstack([e0, C[1:-1], e1]) if len(ub) > 2 else np.vstack([e0, e1])
            else:
                C = np.vstack([e0, e1])
            TIP[new_piece(OBJ["track"], i, rr, length, e0, e1, True, tg0, tg1, v0, v1, C)] = (t0, t1)
    # ---- every other instance that keeps its own piece (P2, S1), with E1
    for i in sorted(others + stub_keep):
        rows = rows_of[i]; g = co[i] in (OBJ["stub"], OBJ["delta"], OBJ["michel"], OBJ["shower"])
        if g and co[i] != OBJ["stub"]:
            rf = rows[free[rows]]; g = len(rf) >= V["free_n"] and float(np.linalg.norm(X[rf].max(0) - X[rf].min(0))) >= V["free_ext_steps"] * step
        new_piece(int(co[i]), i, rows, graph=bool(g))
    # ---- hosts
    elig = piece >= 0; er = np.flatnonzero(elig); et = cKDTree(X[er]) if len(er) else None
    host = np.where(R.absorbed_into >= 0, R.absorbed_into, np.where(R.swallowed_into >= 0, R.swallowed_into, np.arange(nr)))
    def inst_of(x):                                       # the final instance of a present-instance id of the overlap test (a root row, or a pseudo id)
        x = int(x)
        if x < 0:
            return -1
        y = int(host[x]) if x < nr else x
        k = int(np.searchsorted(root, y))
        return k if k < len(root) and root[k] == y else -1
    itree = {}
    def nearest_of_inst(j, p):                            # the nearest eligible point of instance j to the point p
        if j not in itree:
            rj = rows_of.get(j, np.zeros(0, np.int64)); rj = rj[piece[rj] >= 0]
            itree[j] = (rj, cKDTree(X[rj]) if len(rj) else None)
        rj, t = itree[j]
        return int(rj[t.query(X[p])[1]]) if t is not None else -1
    cause1 = getattr(R, "cause1", None); n_exact = 0; n_near = 0; n_none = 0
    def host_point(p, causes, radius, sel=None, me=-1):
        """the eligible point whose class and piece the folded point p takes: by the overlap test, else the nearest within radius"""
        nonlocal n_exact, n_near, n_none
        if causes is not None and dep[p] >= 0:
            for c in causes[dep[p]].tolist():
                j = inst_of(c)
                if j >= 0 and j != me and (sel is None or co[j] == OBJ["shower"]):
                    q = nearest_of_inst(j, p)                 # named by the 2-D test: taken only when it is also near in the image
                    if q >= 0 and np.linalg.norm(X[q] - X[p]) <= (radius if sel is None else 3 * radius):
                        n_exact += 1; return q
        if et is not None:
            if sel is None:
                d_, j = et.query(X[p], distance_upper_bound=radius)
                if np.isfinite(d_):
                    n_near += 1; return int(er[j])
            else:
                for rad, m in ((radius, sel), (radius, None), (3 * radius, sel)):
                    cand = er[m[er]] if m is not None else er[np.isin(cls[er], (CLS7["segment"], CLS7["blip"], CLS7["michel"], CLS7["delta"]))]
                    if len(cand):
                        d_, j = cKDTree(X[cand]).query(X[p], distance_upper_bound=rad)
                        if np.isfinite(d_):
                            n_near += 1; return int(cand[j])
        n_none += 1
        return -1
    is_sh = np.zeros(npt, bool); is_sh[cnt] = co[inst[cnt]] == OBJ["shower"]
    ins = {}; n_h_nohost = 0; shb_host = []               # shb_host: (track instance, rows, shower piece)
    def move(p, q, fl_, rl_, to_cls=None):
        ins.setdefault(int(piece[q]), []).append(p); cls[p] = cls[q] if to_cls is None else to_cls; flag[p] |= fl_; rule[p] = rl_
    for rows, i, rl in to_shower + [(rows_of[i], i, 5) for i in buried]:           # H1, H3: to the shower
        tookp = {}
        for p in rows.tolist():
            q = host_point(p, R.causes, V["sh_r"], is_sh, i)
            if q >= 0:
                tookp[p] = q
        if len(tookp) * 2 < len(rows):                    # no shower found for most of it: the stretch stays track, an own piece outside the graph
            n_h_nohost += 1; new_piece(OBJ["track"], i, rows, graph=False); continue
        hp = np.bincount([piece[q] for q in tookp.values()]).argmax(); hq = np.flatnonzero(piece == hp)[0]
        for p in rows.tolist():
            move(p, tookp.get(p, hq), FL_HID, rl, CLS7["segment"])
        if rl == 4:
            shb_host.append((i, rows, int(hp)))
    def group_host(rows, me):
        """ONE host piece for a folded object: the piece that most of its rows name (host_point per row), when at least half of
        the rows name one; -> (piece or -1, per row the nearest point of that piece).  An object is folded whole or not at all."""
        if not len(rows):
            return -1, None
        q = np.array([host_point(p, cause1, V["host_r"], None, me) for p in rows.tolist()], np.int64); ok = q >= 0
        if not ok.any() or int(ok.sum()) * 2 < len(rows):
            return -1, None
        hp = int(np.bincount(piece[q[ok]]).argmax()); hr = np.flatnonzero(piece == hp)
        return hp, hr[cKDTree(X[hr]).query(X[rows])[1]]
    for i in unres:                                       # U1
        rows = rows_of[i]; hp, near = group_host(rows, i)
        if hp < 0:
            new_piece(OBJ["track"], i, rows, graph=False); continue
        for p, q in zip(rows.tolist(), near.tolist()):
            move(p, q, FL_FTRK, 3)
    for rows, i in short + [(rows_of[i], i) for i in stub_fold]:                  # M1, S2 / S3: the covered points to the one host, the rest EM blip
        stub = co[i] == OBJ["stub"]; nf = ~free[rows]; hp, near = group_host(rows[nf], i)
        if hp >= 0:
            for p, q in zip(rows[nf].tolist(), near.tolist()):
                move(p, q, FL_FSTUB if stub else FL_FTRK, 1 if stub else 6)
            rest = rows[~nf]
        else:
            rest = rows
        if len(rest):
            new_piece(OBJ["le"], i, rest, graph=False); cls[rest] = CLS7["blip"]; flag[rest] |= FL_LE; rule[rest] = 2 if stub else 7
    for hp in sorted(ins):
        rr = np.array(ins[hp], np.int64); piece[rr] = hp; Pk["n"][hp] += len(rr)
    flag[(piece >= 0) & (cls == CLS7["track"]) & (st == STR_TRAVERSED)] |= FL_THRU
    flag[(piece >= 0) & (cls == CLS7["track"]) & (st == STR_HIDDEN_END)] |= FL_HIDKEEP          # D2: a hidden start / end kept as track
    pki = np.array(Pk["inst"], np.int64)
    for i in buried_kept:
        r_ = rows_of[i]; flag[r_[(piece[r_] >= 0) & (cls[r_] == CLS7["track"])]] |= FL_HIDKEEP
    for i in unres_kept:                                  # D2: an unresolved track in its own piece
        r_ = rows_of[i]; own_ = piece[r_] >= 0; own_[own_] = pki[piece[r_][own_]] == i; flag[r_[own_]] |= FL_UNRKEEP
    # ---- candidates
    npc = len(Pk["kind"]); cx, ck, ct, cp, ce, cv, cex = [], [], [], [], [], [], []       # position, kind, tag, piece, end, true vertex, exact position
    for k in range(npc):
        if Pk["kind"][k] == OBJ["track"] and Pk["graph"][k]:
            for e, pos, tg, vt in ((0, Pk["e0"][k], Pk["tag0"][k], Pk["vtx0"][k]), (1, Pk["e1"][k], Pk["tag1"][k], Pk["vtx1"][k])):
                cx.append(pos); ck.append(1); ct.append(tg); cp.append(k); ce.append(e); cv.append(vt); cex.append(tg == TG_PHYS)
    pin = {}
    for k in range(npc):
        pin.setdefault(Pk["inst"][k], []).append(k)
    KX = np.c_[kps["kp_x"], kps["kp_y"], kps["kp_z"]].astype(np.float64).reshape(-1, 3); KK = np.asarray(kps["kp_kind"]).astype(np.int64); KI = np.asarray(kps["kp_inst"]).astype(np.int64)
    folded = set(stub_fold) | set(unres) | set(buried)
    for j in range(len(KX)):
        kd, i = int(KK[j]), int(KI[j])
        if kd not in (2, 3, 4, 5) or i < 0 or i in folded or i not in pin:
            continue
        ks = [k for k in pin[i] if Pk["kind"][k] == co[i]]
        if not ks or not Pk["graph"][ks[0]]:
            continue
        k = ks[0]; r = int(root[i])
        if kd == 3:
            rr = np.flatnonzero(piece == k)
            if not len(rr) or cKDTree(X[rr]).query(KX[j])[0] > 10.0:
                continue
            e_ = R.eidx[r] if r < nr else -1; d_ = R.seg_start_dep[e_] if e_ >= 0 else -1
            vt = P3[d_].astype(np.float64) if d_ >= 0 else nan3
        else:
            vt = T.S0[r].astype(np.float64) if r < nr else nan3
        cx.append(KX[j]); ck.append(kd); ct.append(0); cp.append(k); ce.append(0); cv.append(vt); cex.append(True)
    # D3 (Xin 2026-11-17, doc 86a D3 (a)): a piece in the graph that is not a track and has no candidate gets one at its first visible
    # own deposit: of its points of its own instance (no point moved in), told apart in a wire view, with a matched deposit, the deposit
    # with the smallest Geant4 time; the median of its points.  Tag 16, and 2 when the instance's earliest deposit is not inside.
    KD3 = {OBJ["stub"]: 2, OBJ["shower"]: 3, OBJ["michel"]: 4, OBJ["delta"]: 5}; has_c = set(cp); n_d3 = n_d3_edge = n_d3_none = 0; d3 = []
    ri = np.asarray(R.inst).astype(np.int64); dsrt = np.argsort(ri, kind="stable"); ok_ = ri[dsrt] >= 0          # the deposits of every instance, by index
    ui_, fi_, ni_ = np.unique(ri[dsrt][ok_], return_index=True, return_counts=True); dsrt = dsrt[ok_]; deps_of = {int(u): (int(a_), int(a_ + c_)) for u, a_, c_ in zip(ui_, fi_, ni_)}
    YZ = [(fg[0]["ymin"], fg[0]["ymax"], fg[0]["zmin"], fg[0]["zmax"]) for faces in ev.geo.values() for fg in faces.values()]
    def start_outside(i):                                 # the root particle starts outside the y-z extent of every wire-plane face (it came in)
        r = int(root[i])
        if r >= nr:
            return False
        y, z = float(T.S0[r][1]), float(T.S0[r][2])
        return not any(a_ <= y <= b_ and c_ <= z <= d_ for a_, b_, c_, d_ in YZ)
    moved = (flag & (FL_FSTUB | FL_FTRK | FL_HID)) > 0
    for k in range(npc):
        if not Pk["graph"][k] or Pk["kind"][k] not in KD3 or k in has_c:
            continue
        i = int(Pk["inst"][k]); rows = np.flatnonzero((piece == k) & (inst == i) & free & (dep >= 0) & ~moved)
        if not len(rows):
            n_d3_none += 1; continue
        dd = dep[rows]; a_, b_ = deps_of.get(i, (0, 0)); od = dsrt[a_:b_]                      # the instance's deposits in time, ties along the path (tie_order)
        od = od[np.argsort(ev.T0[od], kind="stable")]; od = od[tie_order(ev.T0[od], P3[od])] if len(od) > 1 else od
        hit = np.flatnonzero(np.isin(od, dd))
        d0 = int(od[hit[0]]) if len(hit) else int(dd[np.lexsort((dd, ev.T0[dd]))[0]]); pos = np.median(X[rows[dd == d0]], axis=0)
        e_ = int(od[0]) if len(od) else -1; tg = TG_VIS | (TG_EDGE if (e_ < 0 or not ev.inside[e_] or start_outside(i)) else 0)
        n_d3 += 1; n_d3_edge += int((tg & TG_EDGE) > 0); d3.append(len(cx))
        cx.append(pos); ck.append(KD3[Pk["kind"][k]]); ct.append(tg); cp.append(k); ce.append(0); cv.append(P3[d0].astype(np.float64)); cex.append(False)
    nc = len(cx); cx = np.array(cx, np.float64).reshape(-1, 3); ck = np.array(ck, np.int64); ct = np.array(ct, np.int64); cp = np.array(cp, np.int64); ce = np.array(ce, np.int64)
    cv = np.array(cv, np.float64).reshape(-1, 3); cex = np.array(cex, bool)
    # ---- places: single linkage at R_merge, one place per true vertex, a shower's start with the boundary it took
    link = []; shb_piece = {}
    if nc:
        hv = np.flatnonzero(np.isfinite(cv[:, 0]))
        if len(hv) > 1:
            for a_, b_ in cKDTree(cv[hv]).query_pairs(V["vertex_true"], output_type="ndarray").tolist():
                if np.linalg.norm(cx[hv[a_]] - cx[hv[b_]]) <= V["vertex_reco"]:
                    link.append((int(hv[a_]), int(hv[b_])))
        for j in np.flatnonzero((ck == 1) & (ct == TG_SHB)).tolist():
            best = (1e9, -1)
            for i_, rows, hp in shb_host:
                if i_ == Pk["inst"][cp[j]]:
                    d_ = float(np.linalg.norm(X[rows] - cx[j], axis=1).min())
                    if d_ < best[0]:
                        best = (d_, hp)
            if best[1] >= 0 and best[0] <= 40.0:
                shb_piece[j] = best[1]; js = np.flatnonzero((cp == best[1]) & (ck == 3))
                if len(js) and np.linalg.norm(cx[js[0]] - cx[j]) <= V["vertex_reco"]:
                    link.append((j, int(js[0])))
        pr = cKDTree(cx).query_pairs(V["R_merge"], output_type="ndarray").reshape(-1, 2)
        if link:
            pr = np.vstack([pr, np.array(link, np.int64).reshape(-1, 2)])
        npl, cpl = connected_components(coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(nc, nc)), directed=False)
    else:
        npl, cpl = 0, np.zeros(0, np.int64)
    PL = np.zeros((npl, 3)); pl_kinds = np.zeros(npl, np.int64); pl_tag = np.zeros(npl, np.int64)
    for u in range(npl):
        m = cpl == u; mm = m & cex; PL[u] = cx[mm].mean(0) if mm.any() else cx[m].mean(0)
        pl_kinds[u] = int(np.bitwise_or.reduce(1 << ck[m])); pl_tag[u] = int(np.bitwise_or.reduce(ct[m]))
        if m.sum() == 1 and ck[m][0] == 1 and int(cp[m][0]) in TIP:      # a free end of a track piece: the place is where its points stop
            PL[u] = TIP[int(cp[m][0])][int(ce[m][0])]
    # X2: a place shared by two or more pieces that lies more than host_r from every member piece's points moves onto them: to the median
    # of the member pieces' points within end_w of the member point nearest to it (doc 86 X2; 19 places on dev, at most 20 mm off)
    n_x2 = 0
    for u in range(npl):
        mem = np.unique(cp[cpl == u])
        if len(mem) < 2:
            continue
        rm = np.flatnonzero(np.isin(piece, mem))
        if not len(rm):
            continue
        dm = np.linalg.norm(X[rm] - PL[u], axis=1); jn = int(np.argmin(dm))
        if dm[jn] > V["host_r"]:
            PL[u] = np.median(X[rm[np.linalg.norm(X[rm] - X[rm[jn]], axis=1) <= V["end_w"]]], axis=0); n_x2 += 1
    d3only = np.zeros(npl, bool); d3par = {}
    if nc:
        isd3 = np.zeros(nc, bool); isd3[d3] = True; par = np.asarray(tab["inst_parent_inst"]).astype(np.int64)
        for u in range(npl):
            m = cpl == u
            if isd3[m].all():
                d3only[u] = True; d3par[u] = {int(par[Pk["inst"][int(c)]]) for c in cp[m]}
    n_d3_chain = n_d3_skip = 0
    # ---- the implied graph
    ea, eb, ep = [], [], []; nv = 0; p0 = np.full(npc, -1, np.int64); p1 = np.full(npc, -1, np.int64); CL = {}
    for k in range(npc):
        if not Pk["graph"][k]:
            continue
        mk = cp == k
        if Pk["kind"][k] == OBJ["track"]:
            j0 = np.flatnonzero(mk & (ce == 0)); j1 = np.flatnonzero(mk & (ce == 1)); u0, u1 = int(cpl[j0[0]]), int(cpl[j1[0]]); p0[k], p1[k] = u0, u1
            C = Pk["cl"][k].copy()
            if 2.0 < np.linalg.norm(PL[u0] - C[0]) <= 40.0:
                C = np.vstack([PL[u0], C])
            if 2.0 < np.linalg.norm(PL[u1] - C[-1]) <= 40.0:
                C = np.vstack([C, PL[u1]])
            C = _resample(C); CL[k] = C; ac = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(C, axis=0), axis=1))]
            d_, j = cKDTree(C).query(PL); near_u = [u for u in np.flatnonzero(d_ <= V["pass_r"]).tolist() if u not in (u0, u1)]
            for u in near_u:                              # a place of D3 candidates only joins the chain of its object's parent (D5: a crossing is not a vertex)
                if d3only[u]:
                    if Pk["inst"][k] in d3par[u]:
                        n_d3_chain += 1
                    else:
                        n_d3_skip += 1
            mid = sorted((float(ac[j[u]]), int(u)) for u in near_u if not d3only[u] or Pk["inst"][k] in d3par[u])
            chain = [u0] + [u for _, u in mid] + [u1]
            for a_, b_ in zip(chain[:-1], chain[1:]):
                ea.append(a_); eb.append(b_); ep.append(k)
        else:
            j = np.flatnonzero(mk); bnd = sorted({int(cpl[jj]) for jj, hp in shb_piece.items() if hp == k})
            if len(j):
                c0 = int(cpl[j[0]]); bnd = [u for u in bnd if u != c0]; p0[k] = c0
                if not bnd:
                    nv += 1; ea.append(c0); eb.append(-nv); ep.append(k)
                for u in bnd:
                    ea.append(c0); eb.append(u); ep.append(k)
            elif bnd:
                nv += 1; p0[k] = bnd[0]
                for u in bnd:
                    ea.append(-nv); eb.append(u); ep.append(k)
    # ---- the zone (Z1), the distance to the nearest place
    if CL:
        ks = sorted(CL); cs = np.vstack([CL[k] for k in ks]); cid = np.concatenate([np.full(len(CL[k]), k) for k in ks]); tk = np.flatnonzero((piece >= 0) & np.isin(piece, ks))
        if len(tk):
            kq = min(48, len(cs)); d_, j = cKDTree(cs).query(X[tk], k=kq, distance_upper_bound=V["R_zone"]); d_ = d_.reshape(len(tk), -1); j = j.reshape(len(tk), -1)
            flag[tk[(np.isfinite(d_) & (cid[np.minimum(j, len(cs) - 1)] != piece[tk][:, None])).any(1)]] |= FL_ZONE
    pd = cKDTree(PL).query(X)[0] if npl else np.full(npt, -1.0)
    mem = np.unique(np.c_[cpl, cp], axis=0) if nc else np.zeros((0, 2), np.int64)
    out = dict(pt_c7_piece=piece, pt_c7_flag=flag, pt_c7_place_dist=np.asarray(pd, np.float32),
               pc_type=np.array(Pk["kind"], np.int64), pc_inst=np.array(Pk["inst"], np.int64), pc_link=g4[np.array(Pk["inst"], np.int64)] if npc else np.zeros(0, np.int64),
               pc_graph=np.array(Pk["graph"], bool), pc_len=np.array(Pk["length"], np.float32), pc_e0=np.array([TIP[k][0] if k in TIP else Pk["e0"][k] for k in range(npc)], np.float32).reshape(-1, 3), pc_e1=np.array([TIP[k][1] if k in TIP else Pk["e1"][k] for k in range(npc)], np.float32).reshape(-1, 3),
               pc_place0=p0, pc_place1=p1, pl_xyz=PL.astype(np.float32), pl_kinds=pl_kinds, pl_tag=pl_tag, pl_member_place=mem[:, 0].astype(np.int64), pl_member_piece=mem[:, 1].astype(np.int64),
               seg_a=np.array(ea, np.int64), seg_b=np.array(eb, np.int64), seg_piece=np.array(ep, np.int64))
    rb = np.bincount(rule[cnt], minlength=8)
    log = dict(n_c7_pieces=int(npc), n_c7_graph_pieces=int(np.sum(Pk["graph"])), n_c7_track_pieces=int(sum(1 for k in range(npc) if Pk["kind"][k] == OBJ["track"] and Pk["graph"][k])),
               n_c7_places=int(npl), n_c7_segments=int(len(ea)), n_c7_stub_kept=int(len(stub_keep)), n_c7_stub_folded=int(len(stub_fold)), n_c7_unresolved=int(len(unres_kept)), n_c7_buried=int(len(buried_kept)),
               n_c7_pts_folded_stub=int(rb[1]), n_c7_pts_stub_le=int(rb[2]), n_c7_pts_folded_track=int(rb[3]), n_c7_pts_hidden_end=int(rb[4]), n_c7_pts_buried=int(rb[5]),
               n_c7_pts_short_host=int(rb[6]), n_c7_pts_short_le=int(rb[7]), n_c7_pts_through=int(((flag & FL_THRU) > 0).sum()), n_c7_pts_zone=int(((flag & FL_ZONE) > 0).sum()),
               n_c7_host_exact=int(n_exact), n_c7_host_near=int(n_near), n_c7_host_none=int(n_none), n_c7_hidden_no_shower=int(n_h_nohost),
               n_c7_d3=int(n_d3), n_c7_d3_edge=int(n_d3_edge), n_c7_d3_none=int(n_d3_none), n_c7_d3_chain=int(n_d3_chain), n_c7_d3_skip=int(n_d3_skip), n_c7_x2=int(n_x2))
    return out, log
