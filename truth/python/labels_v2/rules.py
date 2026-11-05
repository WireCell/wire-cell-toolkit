"""labels_v2.rules -- the ruled definitions applied to one event's deposits (50b secs 2.2, 2.3, 2.5, 2.6, 2.7, Q12 A, Q15).

Order (pr_ml doc 54 sec 1):
  pass 1   provisional instances = every root of objects.build (tracks, delta / Michel candidates, EM objects); 2.2 free
           flags for the candidates' deposits
  2.6      per candidate: own pieces (all its deposits) and free pieces (its free deposits), 2-of-3, extent over the deposits
           seen in >= 2 views; main free piece = the one holding most of the root electron's own charge; resolved when its
           extent >= L_tb.  Unresolved -> folded into the parent track (electron deposits), its photons root EM objects
  pass 2   final instances; 2.2 free flags for every deposit (the overlap flag)
  2.5      pieces of every EM instance (G3 objects, resolved deltas and Michels): segment / blip, nseg, distance to the
           nearest own segment; shower = >= 1 segment; a resolved delta is a shower when it has >= 2 own segments
  Q12 A    an EM object (not born on a track) with no free deposit and >= half its charge shower-ringed -> absorbed
  2.2/2.3  per track kind: Track / Stub by chord and longest resolved run; stretches against showers: traversed / hidden end
           / swallowed stub (absorbed, 10 mm vertex guard) / buried track (G8: 2 steps)
  Q15      an unresolved delta followed its parent (the fold); a resolved delta / Michel starting in a shower-hidden stretch
           of its parent is flagged hidden_start (no key point there)
Outputs per deposit: inst (compact), cls_obj, cls7, piece, piece_type, d_own_seg, overlap bits, stretch; and the instance table.
"""
import numpy as np
from scipy.spatial import cKDTree, ConvexHull
from scipy.spatial.distance import pdist

from . import objects
from .objects import KIND
from .pixels import Projection, pieces_2of3, THR
from steiner_truth import tie_order

STEP = {"sbnd": 6.0, "pdhd": 8.0, "pdvd": 8.0, "fdhd": 8.0, "fdvd": 8.0, "uboone": 6.0}
GAP_MM = 3.0
VERTEX_GUARD_MM = 10.0

# object classes (cls_obj) and the 7-class head (cls7)
OBJ = {"track": 0, "stub": 1, "delta": 2, "michel": 3, "shower": 4, "le": 5, "ghost": 6, "other": 7}
CLS7 = {"ghost": 0, "track": 1, "stub": 2, "segment": 3, "blip": 4, "michel": 5, "delta": 6, "other": -1}
# overlap bits (per deposit / point)
OV_NOTFREE, OV_SHRING, OV_FOLDED, OV_ABSORBED, OV_SWALLOWED, OV_WEAKFOLD = 1, 2, 4, 8, 16, 32
# stretch codes
STR_NONE, STR_TRAVERSED, STR_HIDDEN_END, STR_SWALLOWED, STR_BURIED = 0, 1, 2, 3, 4
# instance flags
FL_ORPHAN, FL_BURIED, FL_HIDDEN_START, FL_HIDDEN_END, FL_VERTEX_STUB, FL_ABSORBED, FL_SWALLOWED, FL_RESOLVED, FL_HIDDEN_ATTACH, FL_PRIMARY, FL_NODELTA_TRUTH, FL_VISIBLE = (
    1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
EDGE = {"none": 0, "delta": 1, "michel": 2, "daughter": 3, "absorbed": 4, "swallowed": 5}


def diameter(p):
    if len(p) < 2:
        return 0.0
    if len(p) > 60:
        try:
            p = p[ConvexHull(p).vertices]
        except Exception:
            pass
    if len(p) > 3000:
        p = p[np.random.default_rng(0).choice(len(p), 3000, replace=False)]
    return float(pdist(p).max())


def longest_run(r, s, P):
    """longest run of resolved deposits along the path: a run breaks only where an UNRESOLVED stretch spans > GAP_MM of path
    (d50d_real.longest_run, unchanged) -> (i0, i1, chord, nrun) or None"""
    idx = np.flatnonzero(r)
    if len(idx) == 0:
        return None
    runs = []; st = idx[0]; pv = idx[0]
    for j in idx[1:]:
        if j > pv + 1 and (s[j - 1] - s[pv + 1]) + 0.5 * ((s[pv + 1] - s[pv]) + (s[j] - s[j - 1])) > GAP_MM:
            runs.append((st, pv)); st = j
        pv = j
    runs.append((st, pv))
    a, b = max(runs, key=lambda ab: s[ab[1]] - s[ab[0]])
    return a, b, float(np.linalg.norm(P[b] - P[a])), len(runs)


def runs_of(mask):
    idx = np.flatnonzero(mask)
    if not len(idx):
        return []
    br = np.flatnonzero(np.diff(idx) > 1)
    return list(zip(np.r_[idx[0], idx[br + 1]].tolist(), np.r_[idx[br], idx[-1]].tolist()))


def csr(keys, n_keys):
    """rows grouped by key: order, start offsets"""
    o = np.argsort(keys, kind="stable")
    st = np.searchsorted(keys[o], np.arange(n_keys + 1))
    return o, st


def piece_stats(piece, dob, comp, q, P3, selfdep=None):
    """per piece (sorted by piece id): object, charge, self charge, extent over the deposits seen in >= 2 views; and the
    piece id list"""
    sel = np.flatnonzero(piece >= 0)
    if not len(sel):
        z = np.zeros(0)
        return dict(pid=np.zeros(0, np.int64), obj=np.zeros(0, np.int64), q=z, sq=z, ext=z, members=[])
    o = sel[np.argsort(piece[sel], kind="stable")]
    pid = piece[o]; b = np.r_[0, np.flatnonzero(np.diff(pid)) + 1, len(o)]
    obj = dob[o[b[:-1]]]; pq = np.add.reduceat(q[o], b[:-1])
    sq = np.add.reduceat(q[o] * selfdep[o], b[:-1]) if selfdep is not None else pq
    mv = (comp[o] >= 0).sum(1) >= 2
    ext = np.array([diameter(P3[o[b[i]:b[i + 1]][mv[b[i]:b[i + 1]]]]) for i in range(len(b) - 1)])
    members = [o[b[i]:b[i + 1]] for i in range(len(b) - 1)]
    return dict(pid=pid[b[:-1]], obj=obj, q=pq, sq=sq, ext=ext, members=members)


class Result:
    pass


def label_event(ev, det=None, params=None):
    det = det or ev.det
    step = STEP[det]; ltb = 2.0 * step
    prm = dict(step_mm=step, ltb_mm=ltb, thr_e=THR, eps=0.5, ring=3, gap_mm=GAP_MM, vertex_guard_mm=VERTEX_GUARD_MM,
               wait_ns=objects.WAIT_NS, michel_capture_ke_min_mev=objects.MICHEL_CAPTURE_KE_MIN, extent="2view", two_pass=True)
    if params:
        prm.update(params)
    n, q, P3 = ev.n, ev.q, ev.P3
    T = objects.build(ev, THR)
    P = Projection(ev)
    R = Result(); R.T, R.P, R.prm = T, P, prm
    log = R.log = {}
    nr = T.nrow

    # ------------------------------------------------------------------ pass 1 and 2.6: the candidates
    inst0 = T.inst_dep
    cand_rows = np.flatnonzero(((T.kind == KIND["DELTA"]) | (T.kind == KIND["MICHEL"])) & ~T.fold)
    cq = np.bincount(np.where(inst0 >= 0, inst0, 0), weights=q * (inst0 >= 0), minlength=nr + 1)
    cand_rows = cand_rows[cq[cand_rows] >= THR]
    K = len(cand_rows)
    cidx = np.full(nr + 1 + int(ev.tid.max() if n else 0) + 1, -1, np.int64)
    cidx[cand_rows] = np.arange(K)
    dob_c = np.where(inst0 >= 0, cidx[np.minimum(inst0, len(cidx) - 1)], -1)
    is_cand_dep = dob_c >= 0
    resolved = np.zeros(K, bool); main_ext = np.zeros(K); nseg_own = np.zeros(K, np.int64)
    if K:
        CS0, _ = P.clean_share(inst0, is_cand_dep)
        free0 = P.free(CS0) & is_cand_dep
        selfdep = is_cand_dep & (T.grow == inst0)
        comp_own = P.components(dob_c); p_own = pieces_2of3(dob_c, comp_own, q)
        so = piece_stats(p_own, dob_c, comp_own, q, P3, selfdep)
        np.add.at(nseg_own, so["obj"][so["ext"] >= ltb], 1)
        dfree = np.where(free0, dob_c, -1)
        comp_fr = P.components(dfree); p_fr = pieces_2of3(dfree, comp_fr, q)
        sf = piece_stats(p_fr, dfree, comp_fr, q, P3, selfdep)
        main_sq = np.zeros(K)
        for i in range(len(sf["obj"])):
            k = sf["obj"][i]
            if sf["sq"][i] > main_sq[k]:
                main_sq[k], main_ext[k] = sf["sq"][i], sf["ext"][i]
        resolved = main_ext >= ltb
        for k in np.flatnonzero(~resolved).tolist():
            T.fold_candidate(int(cand_rows[k]))
        T.root_of = T.roots()
        T.inst_dep = np.where(T.grow >= 0, T.root_of[np.maximum(T.grow, 0)], -1)
        lost = T.inst_dep < 0
        T.inst_dep[lost] = nr + ev.tid[lost]
    log.update(n_candidates=int(K), n_resolved=int(resolved.sum()), n_resolved_michel=int((resolved & (T.kind[cand_rows] == KIND["MICHEL"])).sum()) if K else 0)
    R.cand_rows, R.cand_resolved, R.cand_main_ext, R.cand_nseg_own = cand_rows, resolved, main_ext, nseg_own
    res_rows = cand_rows[resolved]
    res_shower = np.zeros(nr, bool); res_shower[res_rows[nseg_own[resolved] >= 2]] = True
    resolved_row = np.zeros(nr, bool); resolved_row[res_rows] = True

    # ------------------------------------------------------------------ pass 2: final instances, every deposit's overlap flag
    inst1 = T.inst_dep                                    # root rows (>= nr: pseudo rows by stored id)
    test_all = inst1 >= 0
    CS1, cause1 = P.clean_share(inst1, test_all)
    free1 = P.free(CS1) & test_all
    R.free, R.CS = free1, CS1

    # ------------------------------------------------------------------ 2.5: pieces of every EM instance
    is_row = inst1 < nr
    rk = np.where(is_row, T.kind[np.minimum(inst1, nr - 1)], KIND["OTHER"])
    em_root = is_row & ((rk == KIND["EM"]) | (resolved_row[np.minimum(inst1, nr - 1)] & is_row))
    em_rows = np.unique(inst1[em_root])
    M = len(em_rows); eidx = np.full(nr, -1, np.int64); eidx[em_rows] = np.arange(M)
    dob_e = np.where(em_root, eidx[np.minimum(inst1, nr - 1)], -1)
    piece_type = np.zeros(n, np.int8); piece_id = np.full(n, -1, np.int64); d_own_seg = np.full(n, -1.0, np.float32)
    nseg_e = np.zeros(M, np.int64); npiece_e = np.zeros(M, np.int64)
    seg_start_dep = np.full(M, -1, np.int64)              # the first-in-time deposit of the earliest segment
    if M:
        comp_e = P.components(dob_e); piece_e = pieces_2of3(dob_e, comp_e, q)
        se = piece_stats(piece_e, dob_e, comp_e, q, P3)
        piece_id = piece_e
        isseg = se["ext"] >= ltb
        np.add.at(nseg_e, se["obj"], 1 * 0 + 1); npiece_e = nseg_e.copy(); nseg_e[:] = 0
        np.add.at(nseg_e, se["obj"][isseg], 1)
        t0seg = np.full(M, np.inf)
        for i, mem in enumerate(se["members"]):
            piece_type[mem] = 1 if isseg[i] else 2
            if isseg[i]:
                k = se["obj"][i]; j = mem[np.argmin(ev.T0[mem])]
                if ev.T0[j] < t0seg[k]:
                    t0seg[k], seg_start_dep[k] = ev.T0[j], j
        # distance of every EM deposit to the nearest own-segment deposit (0 inside a segment)
        o_e, st_e = csr(np.where(dob_e >= 0, dob_e, M), M + 1)
        for k in np.flatnonzero(nseg_e > 0).tolist():
            mem = o_e[st_e[k]:st_e[k + 1]]
            sd = mem[piece_type[mem] == 1]
            d_own_seg[sd] = 0.0
            bl = mem[piece_type[mem] != 1]
            if len(bl) and len(sd):
                d_own_seg[bl] = cKDTree(P3[sd]).query(P3[bl])[0]
    R.em_rows, R.nseg_e, R.npiece_e = em_rows, nseg_e, npiece_e

    # ------------------------------------------------------------------ classes per root row (before absorption)
    cls_row = np.full(nr, OBJ["other"], np.int64)
    cls_row[T.kind == KIND["TRACK"]] = OBJ["stub"]                     # Track decided below
    for k, r in enumerate(em_rows.tolist()):
        if T.kind[r] == KIND["EM"]:
            cls_row[r] = OBJ["shower"] if nseg_e[k] >= 1 else OBJ["le"]
        elif T.kind[r] == KIND["MICHEL"]:
            cls_row[r] = OBJ["michel"]
        else:
            cls_row[r] = OBJ["shower"] if res_shower[r] else OBJ["delta"]
    is_shower_row = cls_row == OBJ["shower"]

    # shower-ringed shares: only shower instances can be present
    shower_dep = is_row & is_shower_row[np.minimum(inst1, nr - 1)]
    CSs, causes = P.clean_share(inst1, test_all, qmask=shower_dep.astype(np.float64))
    ringed = ((CSs <= 0.5 + 1e-9) & P.seen).any(1) & test_all
    R.ringed = ringed

    # ------------------------------------------------------------------ Q12 A: EM objects with no free deposit, shower-ringed
    absorbed_into = np.full(nr, -1, np.int64)
    o_r, st_r = csr(np.where(inst1 < nr, inst1, nr), nr + 1)
    for k, r in enumerate(em_rows.tolist()):
        if T.kind[r] != KIND["EM"]:
            continue
        mem = o_r[st_r[r]:st_r[r + 1]]
        if free1[mem].any():
            continue
        qr = q[mem]; rg = ringed[mem]
        if qr[rg].sum() < 0.5 * qr.sum() or not rg.any():
            continue
        cs = causes[mem][rg]; cs = cs[cs >= 0]
        cs = cs[(cs < nr)]
        if not len(cs):
            continue
        cs = cs[is_shower_row[cs]]
        if not len(cs):
            continue
        u, cnt = np.unique(cs, return_counts=True)
        absorbed_into[r] = u[np.argmax(cnt)]
    # follow chains (an absorbed object's host absorbed itself)
    for _ in range(5):
        m = (absorbed_into >= 0) & (absorbed_into[np.maximum(absorbed_into, 0)] >= 0)
        if not m.any():
            break
        absorbed_into[m] = absorbed_into[absorbed_into[m]]
    log["n_absorbed_em"] = int((absorbed_into >= 0).sum())
    q_root = np.bincount(np.where(inst1 < nr, inst1, nr), weights=q, minlength=nr + 1)
    log["n_absorbed_em_visible"] = int(((absorbed_into >= 0) & (q_root[:nr] >= THR)).sum())       # 50i's population: objects >= 2,400 e

    # ------------------------------------------------------------------ 2.2 / 2.3: track kinds
    trk_rows = np.unique(inst1[is_row & (rk == KIND["TRACK"])])
    stretch = np.zeros(n, np.int8)
    trk_info = {}
    sh_first = {}                                                     # shower row -> first-in-time deposit position (true frame)
    def shower_start(r):
        if r not in sh_first:
            mem = o_r[st_r[r]:st_r[r + 1]]
            sh_first[r] = P3[mem[np.argmin(ev.T0[mem])]] if len(mem) else None
        return sh_first[r]
    swallowed_into = np.full(nr, -1, np.int64)
    for r in trk_rows.tolist():
        mem = o_r[st_r[r]:st_r[r + 1]]
        own = mem[(ev.otid[mem] == T.tt[r]) & ev.placed[mem]]
        info = dict(track=False, chord=0.0, path=0.0, run=0.0, ratio=0.0, hfrac=0.0, flags=0, i_first=-1, i_last=-1, start_hidden=False, end_hidden=False)
        trk_info[r] = info
        if len(own) < 2:
            info["ratio"] = float(free1[own].sum() / max(len(own), 1))
            info["i_first"] = info["i_last"] = int(own[0]) if len(own) else -1
            continue
        own = own[np.lexsort((ev.T0[own],))]
        own = own[tie_order(ev.T0[own], P3[own])]
        Pp = P3[own]
        sseg = np.r_[0, np.cumsum(np.linalg.norm(np.diff(Pp, axis=0), axis=1))]
        chord = float(np.linalg.norm(Pp[-1] - Pp[0]))
        Rm = free1[own]; Hm = ~Rm & ringed[own]
        qq = q[own]; qt = max(qq.sum(), 1e-9)
        lr = longest_run(Rm, sseg, Pp)
        run = lr[2] if lr else 0.0
        track = chord >= 2 * step and run >= step
        info.update(track=bool(track), chord=chord, path=float(sseg[-1]), run=run, ratio=float(qq[Rm].sum() / qt), hfrac=float(qq[Hm].sum() / qt),
                    i_first=int(own[0]), i_last=int(own[-1]))
        good = [(a, b) for a, b in runs_of(Rm) if np.linalg.norm(Pp[b] - Pp[a]) >= step]
        if track:
            for a, b in runs_of(~Rm):
                hq = qq[a:b + 1][Hm[a:b + 1]].sum(); tq = qq[a:b + 1].sum()
                if hq < 0.5 * tq:
                    continue
                left = any(e_ < a for s_, e_ in good); right = any(s_ > b for s_, e_ in good)
                if not (left or right):
                    continue
                code = STR_TRAVERSED if (left and right) else STR_HIDDEN_END
                stretch[own[a:b + 1]] = code
                if code == STR_HIDDEN_END:
                    if a == 0:
                        info["start_hidden"] = True; info["flags"] |= FL_HIDDEN_START
                    if b == len(own) - 1:
                        info["end_hidden"] = True; info["flags"] |= FL_HIDDEN_END
        elif info["hfrac"] >= 0.5:
            if chord >= 2 * step:
                stretch[mem] = STR_BURIED; info["flags"] |= FL_BURIED
            else:
                cs = causes[own][Hm]; cs = cs[(cs >= 0) & (cs < nr)]
                cs = cs[is_shower_row[cs]] if len(cs) else cs
                if len(cs):
                    u, cnt = np.unique(cs, return_counts=True); host = int(u[np.argmax(cnt)])
                    hs = shower_start(host)
                    if hs is not None and np.linalg.norm(Pp[0] - hs) <= VERTEX_GUARD_MM:
                        info["flags"] |= FL_VERTEX_STUB
                    else:
                        swallowed_into[r] = host; stretch[mem] = STR_SWALLOWED; info["flags"] |= FL_SWALLOWED
        if track:
            cls_row[r] = OBJ["track"]
    log["n_track"] = int(sum(1 for i in trk_info.values() if i["track"])); log["n_stub"] = int(len(trk_info) - log["n_track"])
    log["n_swallowed"] = int((swallowed_into >= 0).sum()); log["n_buried"] = int(sum(1 for i in trk_info.values() if i["flags"] & FL_BURIED))
    R.trk_info = trk_info

    # ------------------------------------------------------------------ Q15: hidden start of a resolved delta / Michel
    hidden_attach = np.zeros(nr, bool)
    for r in res_rows.tolist():
        m = T.mrow[r]
        if m < 0 or m not in trk_info:
            continue
        pm = o_r[st_r[m]:st_r[m + 1]]
        pm = pm[ev.otid[pm] == T.tt[m]]
        if not len(pm):
            continue
        j = pm[np.argmin(np.linalg.norm(P3[pm] - T.S0[r], axis=1))]
        cj = causes[j][causes[j] >= 0]                     # the shower instances ringing the attachment deposit (per view)
        if stretch[j] in (STR_TRAVERSED, STR_HIDDEN_END, STR_BURIED) and (not len(cj) or (cj != r).any()):   # not shadowed by its own shower only
            hidden_attach[r] = True
    log["n_hidden_attach"] = int(hidden_attach.sum())

    # ------------------------------------------------------------------ final instance of every deposit, compact ids
    host = np.where(absorbed_into >= 0, absorbed_into, np.where(swallowed_into >= 0, swallowed_into, np.arange(nr)))
    inst_row = np.where(is_row, host[np.minimum(inst1, nr - 1)], inst1)       # root row after absorption, or pseudo id
    uroots, inst = np.unique(inst_row, return_inverse=True)
    inst = inst.ravel().astype(np.int64)
    NI = len(uroots)
    # per deposit class
    cls_obj = np.full(n, OBJ["other"], np.int8)
    rowc = np.where(is_row, cls_row[np.minimum(inst1, nr - 1)], OBJ["other"])
    cls_obj[:] = rowc
    abs_dep = is_row & (absorbed_into[np.minimum(inst1, nr - 1)] >= 0)
    swl_dep = is_row & (swallowed_into[np.minimum(inst1, nr - 1)] >= 0)
    cls_obj[abs_dep | swl_dep] = OBJ["shower"]
    cls7 = np.full(n, CLS7["other"], np.int8)
    cls7[cls_obj == OBJ["track"]] = CLS7["track"]; cls7[cls_obj == OBJ["stub"]] = CLS7["stub"]
    cls7[cls_obj == OBJ["delta"]] = CLS7["delta"]; cls7[cls_obj == OBJ["michel"]] = CLS7["michel"]
    sh = cls_obj == OBJ["shower"]
    cls7[sh] = np.where(piece_type[sh] == 1, CLS7["segment"], CLS7["blip"])
    cls7[cls_obj == OBJ["le"]] = CLS7["blip"]
    # overlap bits
    overlap = np.zeros(n, np.int16)
    overlap[~free1 & test_all] |= OV_NOTFREE
    overlap[ringed] |= OV_SHRING
    folded_dep = (T.grow >= 0) & T.fold[np.maximum(T.grow, 0)]
    overlap[folded_dep] |= OV_FOLDED
    overlap[abs_dep] |= OV_ABSORBED; overlap[swl_dep] |= OV_SWALLOWED
    R.inst, R.inst_root, R.cls_obj, R.cls7, R.piece_id, R.piece_type, R.d_own_seg, R.overlap, R.stretch = inst, uroots, cls_obj, cls7, piece_id, piece_type, d_own_seg, overlap, stretch
    R.absorbed_into, R.swallowed_into, R.hidden_attach, R.cls_row, R.seg_start_dep, R.eidx = absorbed_into, swallowed_into, hidden_attach, cls_row, seg_start_dep, eidx
    R.o_r, R.st_r = o_r, st_r
    # ------------------------------------------------------------------ the instance table
    R.table = instance_table(ev, T, R, NI, uroots, inst)
    log.update(n_instances=int(NI), n_em_objects=int(M), n_deposits=int(n), n_free=int(free1.sum()))
    return R


def instance_table(ev, T, R, NI, uroots, inst):
    nr = T.nrow; n = ev.n; q = ev.q
    cols = {k: np.zeros(NI, np.int64) for k in ("root_row", "g4_id", "kind", "cls_obj", "pdg", "parent_inst", "parent_g4", "edge", "interaction", "origin", "flags", "n_dep", "nseg", "npiece")}
    fcols = {k: np.zeros(NI, np.float64) for k in ("q", "e_dep_mev", "ke_mev", "vis_frac", "ratio", "chord_mm", "run_mm", "wait_ns", "main_ext_mm")}
    origins = []; oidx = {}
    root_inst = {int(r): i for i, r in enumerate(uroots.tolist())}
    qi = np.bincount(inst, weights=q, minlength=NI); ei = np.bincount(inst, weights=ev.energy, minlength=NI); ni = np.bincount(inst, minlength=NI)
    cand_pos = {int(r): k for k, r in enumerate(R.cand_rows.tolist())}
    for i, r in enumerate(uroots.tolist()):
        cols["n_dep"][i] = ni[i]; fcols["q"][i] = qi[i]; fcols["e_dep_mev"][i] = ei[i]
        if r >= nr:                                       # pseudo instance: a stored id with no table row
            cols["root_row"][i] = -1; cols["g4_id"][i] = r - nr; cols["kind"][i] = KIND["OTHER"]; cols["cls_obj"][i] = OBJ["other"]
            cols["parent_inst"][i] = -1; cols["parent_g4"][i] = -1; cols["interaction"][i] = -1; cols["origin"][i] = -1
            continue
        cols["root_row"][i] = r; cols["g4_id"][i] = T.tt[r]; cols["kind"][i] = T.kind[r] if not T.fold[r] else KIND["OTHER"]
        cols["cls_obj"][i] = R.cls_row[r]; cols["pdg"][i] = T.pdg[r]
        fcols["ke_mev"][i] = T.ke[r]; fcols["vis_frac"][i] = ei[i] / T.ke[r] if T.ke[r] > 0 else 0.0
        m = T.mrow[r]
        cols["parent_g4"][i] = T.tt[m] if m >= 0 else -1
        pr = T.root_of[m] if m >= 0 else -1
        host = R.absorbed_into[pr] if (pr >= 0 and R.absorbed_into[pr] >= 0) else (R.swallowed_into[pr] if (pr >= 0 and R.swallowed_into[pr] >= 0) else pr)
        cols["parent_inst"][i] = root_inst.get(int(host), -1) if pr >= 0 else -1
        k = T.kind[r]
        cols["edge"][i] = EDGE["delta"] if k == KIND["DELTA"] else EDGE["michel"] if k == KIND["MICHEL"] else (EDGE["daughter"] if m >= 0 else EDGE["none"])
        cols["interaction"][i] = T.interaction(r)
        og = T.origin(r)
        if og not in oidx:
            oidx[og] = len(origins); origins.append(og)
        cols["origin"][i] = oidx[og]
        fl = 0
        if qi[i] >= THR: fl |= FL_VISIBLE                 # >= one hit's worth of charge (2.5 / 2.7: an object below THR has no present pixel)
        if T.orphan[r]: fl |= FL_ORPHAN
        if T.primary[r]: fl |= FL_PRIMARY
        if r in R.trk_info:
            ti = R.trk_info[r]; fl |= ti["flags"]; fcols["ratio"][i] = ti["ratio"]; fcols["chord_mm"][i] = ti["chord"]; fcols["run_mm"][i] = ti["run"]
        if r in cand_pos:
            kk = cand_pos[r]; fcols["main_ext_mm"][i] = R.cand_main_ext[kk]; fcols["wait_ns"][i] = T.wait[r] if np.isfinite(T.wait[r]) else 0.0
            if R.cand_resolved[kk]: fl |= FL_RESOLVED
            if R.hidden_attach[r]: fl |= FL_HIDDEN_ATTACH
        if R.eidx[r] >= 0:
            cols["nseg"][i] = R.nseg_e[R.eidx[r]]; cols["npiece"][i] = R.npiece_e[R.eidx[r]]
        cols["flags"][i] = fl
    # the members absorbed into each instance (EM objects and swallowed stubs): count and list by root row
    absorbed_members = {}
    for src in np.flatnonzero(R.absorbed_into >= 0).tolist():
        absorbed_members.setdefault(root_inst.get(int(R.absorbed_into[src]), -1), []).append((src, EDGE["absorbed"]))
    for src in np.flatnonzero(R.swallowed_into >= 0).tolist():
        absorbed_members.setdefault(root_inst.get(int(R.swallowed_into[src]), -1), []).append((src, EDGE["swallowed"]))
    out = dict(**{"inst_" + k: v for k, v in cols.items()}, **{"inst_" + k: v for k, v in fcols.items()})
    out["origins"] = origins
    out["absorbed_members"] = absorbed_members
    return out
