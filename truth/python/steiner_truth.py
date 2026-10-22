#!/usr/bin/env python3
"""pr_ml doc 23 (WP0.2 "G" labels): truth labels on the Steiner graph and the skeleton metrics of the rule routing.

Input: the pattern-recognition job's point-tree dump (TensorFileSink, save_tensors) with the clouds that the
clus visitor SteinerTruthLabeler adds to every cluster carrying a Steiner graph:
    steiner_pc     the Steiner vertices (default-scope coordinates, flag_steiner_terminal)
    steiner_truth  per vertex: tid (majority key of the cluster points within 0.5 cm), nearest point, frame shift
    steiner_edges  src, dst, weight, length
plus the drifted depos and the truth table of the event.

The trajectories T_p (pr_ml doc 12 sec 4.1) are the polylines of the drifted depos of each Geant4 particle p (the
original id, sed_ids column 1), placed in the reconstructed frame of THEIR OWN face (x from the arrival time, as the
sampler places blobs) and ordered by the time of their undrifted prior.  The track key of a particle is its larg4 id
folded as the labeler folds it (an unsaved EM daughter's "-ancestor" -> the ancestor); labels are made per particle
(doc 12 sec 4.1: daughters keep their own id, so fragment labels exist) and the key is stored beside them, so the
merging can change offline.  Labels (doc 15 sec 2.2, doc 12 sec 4.3), tau_on / tau_off = 15 / 40 mm on every detector
(owner 2026-10-22, doc 23 sec 7.1; --tau-on-pitch / --tau-off-pitch give the per-pitch bands of doc 12, e.g. 2 / 5):

  node  kstar, dstar     nearest particle and distance (the distance field of doc 12 sec 4.1); key = its track key
        is_vertex        a true vertex within 1 cm: the start of a particle (length >= 1 cm) that another such
                         particle passes within 1 cm of
        is_nu_vertex     the neutrino vertex (median start of the primaries of an interaction) within 1 cm
  edge  cls4             0 on-trajectory: kstar(a) = kstar(b) = p, both within tau_on, every chord sample (3 mm)
                           within tau_on of T_p, and |s(b) - s(a)| <= 1.5 |ab| + 1 pitch (arc-length slack)
                         1 junction: different particles, both ends within tau_on, and the start of one of the two
                           lies within 1 cm of the other particle and within 1 cm of the chord
                         2 false bridge: different particles, both ends within tau_on, no such start
                         3 ghost: an end beyond tau_off
                        -1 don't care (the band cases)
        same_key         the two ends' particles carry the same track key (a parent and its folded daughter)
        flow             +1 / -1: sign of the trajectory tangent (time order) along src -> dst at the midpoint,
                         on-trajectory edges only; else 0

Metrics of the rule skeleton (doc 12 sec 4.4) over the clusters with a graph: coverage recall R_d and precision
P_d (d = 0.5, 1, 2 cm), edge-class counts and lengths, distance of the true particle ends to the nearest node and
terminal (the stored graph keeps cycles, so an end is rarely a leaf), junction recall
(true vertices with a node of degree >= 3 within d), terminal density, and the agreement of the C++ node key
(blob-level majority) with the key of kstar.  The true arc length counted by R_d is the part of the trajectories within
2 cm of a point of those clusters (the rest is out of scope or out of the image, not a routing miss).

Cross-check: kstar/dstar are recomputed at the cluster points and compared with the C++ per-point d_traj of
the "truth" cloud (same keys, same capping).

  steiner_truth.py --pctree <pr dump> --depos <sim-depos.tar.bz2> --truth <truth.npz> --wires <wires file>
                   --depo-toff-us T --reco-toff-us R --speed-mmus V --pitch-mm P --event NAME --out DIR
"""
import argparse, json, os, sys
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_3d import load_tensors  # noqa: E402
from xres_check import load_depos, load_faces, assign_depo_faces  # noqa: E402
import truth_table as T  # noqa: E402

US = 1000.0          # ns per us (WCT time unit is ns)
CM = 10.0            # mm per cm
NOTRUTH = -3


# ---------------------------------------------------------------------------- the dump
def load_tree(fn):
    """The live tree of a PR dump: node parentage, per-node sizes and concatenated arrays of every named cloud."""
    tens = load_tensors(fn)
    tops = sorted(dp for dp, (m, _) in tens.items() if m.get("datatype") == "pctree" and dp.endswith("/live"))
    if not tops:
        raise RuntimeError("no live pctree in " + fn)
    dp = tops[0]
    m, parentage = tens[dp]
    lp = tens[m["lpcmaps"]][0]["arrays"]
    sizes = {k: np.asarray(tens[v][1]).astype(int) for k, v in lp.items()}
    pcs = tens[m["pointclouds"]][0]["items"]
    arrays = {}
    for name, path in pcs.items():
        amd = tens[path][0]
        arrays[name] = {k: np.asarray(tens[v][1]) for k, v in amd["arrays"].items()}
    parentage = np.asarray(parentage).astype(int)
    n = len(parentage)
    cluster_of = np.full(n, -1, dtype=int)
    for i in range(1, n):
        p = parentage[i]
        cluster_of[i] = i if p == 0 else cluster_of[p]
    return dict(parentage=parentage, sizes=sizes, arrays=arrays, cluster_of=cluster_of)


def per_node(tree, name):
    """{node index: {key: rows}} of a named cloud."""
    if name not in tree["sizes"]:
        return {}
    s = tree["sizes"][name]
    off = np.r_[0, np.cumsum(s)]
    out = {}
    for i in np.flatnonzero(s > 0):
        out[int(i)] = {k: v[off[i]:off[i + 1]] for k, v in tree["arrays"][name].items()}
    return out


# ---------------------------------------------------------------------------- trajectories
class Polyline:
    def __init__(self, V):
        self.V = V
        self.A, self.B = V[:-1], V[1:]
        d = self.B - self.A
        self.L = np.linalg.norm(d, axis=1)
        self.T = np.where(self.L[:, None] > 0, d / np.maximum(self.L, 1e-12)[:, None], 0.0)
        self.s = np.r_[0.0, np.cumsum(self.L)]
        self.length = float(self.s[-1])
        self.kd = cKDTree(0.5 * (self.A + self.B)) if len(self.L) else None

    def nearest(self, Q, kc=8):
        """distance, arc length, unit tangent of the closest polyline point of every row of Q."""
        Q = np.atleast_2d(Q)
        if self.kd is None:
            d = np.linalg.norm(Q - self.V[0], axis=1)
            return d, np.zeros(len(Q)), np.zeros((len(Q), 3))
        k = min(kc, len(self.L))
        _, cand = self.kd.query(Q, k=k)
        cand = cand.reshape(len(Q), k)
        a, t, L = self.A[cand], self.T[cand], self.L[cand]
        u = np.clip(np.einsum("nkj,nkj->nk", Q[:, None, :] - a, t) / np.maximum(L, 1e-12), 0.0, 1.0)
        c = a + t * (u * L)[..., None]
        d = np.linalg.norm(Q[:, None, :] - c, axis=2)
        # ties: the earlier segment wins (as DepoPolyline)
        best = np.argmin(d + 1e-12 * cand, axis=1)
        r = np.arange(len(Q))
        seg = cand[r, best]
        return d[r, best], self.s[seg] + u[r, best] * L[r, best], self.T[seg]


class KeyTraj:
    """The trajectory T_k of one track key: one Polyline per original G4 id carrying the key.  A positive larg4 id
    is its own key with one unit; the unsaved EM daughters keyed "-ancestor" (doc 20 sec 6.1) are folded to the
    ancestor's key but keep their own unit (sed_ids column 1), so a shower is a set of daughter polylines rather
    than one time-ordered zigzag through all of them."""

    def __init__(self, key, units, origs, t0s):
        self.key = key
        self.units = units
        self.origs = origs
        V = [u.V for u in units]
        self.V = np.concatenate(V)
        self.vunit = np.concatenate([np.full(len(v), i) for i, v in enumerate(V)])
        self.kd = cKDTree(self.V)
        self.length = float(sum(u.length for u in units))
        main = [i for i, o in enumerate(origs) if o == key]
        self.main = main[0] if main else int(np.argmin(t0s))
        self.start = units[self.main].V[0]

    def nearest(self, Q, kc=8):
        """distance, arc length on its unit, tangent and unit index of the closest point of T_k."""
        Q = np.atleast_2d(Q)
        n = len(Q)
        d = np.full(n, np.inf); s = np.zeros(n); t = np.zeros((n, 3)); u = np.full(n, -1)
        if len(self.units) == 1:
            d, s, t = self.units[0].nearest(Q)
            return d, s, t, np.zeros(n, dtype=int)
        k = min(kc, len(self.V))
        _, ii = self.kd.query(Q, k=k)
        ii = ii.reshape(n, k)
        pairs = np.unique(np.column_stack([self.vunit[ii].ravel(), np.repeat(np.arange(n), k)]), axis=0)
        st = np.r_[0, np.flatnonzero(np.diff(pairs[:, 0])) + 1]
        en = np.r_[st[1:], len(pairs)]
        for s0, e0 in zip(st, en):
            ui, qs = int(pairs[s0, 0]), pairs[s0:e0, 1]
            dd, ss, tt = self.units[ui].nearest(Q[qs])
            better = dd < d[qs]
            q = qs[better]
            d[q], s[q], t[q], u[q] = dd[better], ss[better], tt[better], ui
        return d, s, t, u


def build_trajectories(depos, truth, wires, depo_toff_us, reco_toff_us, speed):
    D, I = load_depos(depos)
    n = len(D)
    act = np.flatnonzero(I[:, 2] == 0)
    # the prior of drifted row i is the generation-1 row whose child column is i
    prior_t = np.full(n, np.nan)
    pri = np.flatnonzero(I[:, 2] == 1)
    ch = I[pri, 3].astype(int)
    ok = (ch >= 0) & (ch < n)
    prior_t[ch[ok]] = D[pri[ok], 0]
    tkey = np.where(np.isfinite(prior_t[act]), prior_t[act], D[act, 0])
    ids = I[act, 0].astype(np.int64)
    if truth is not None:
        sed = np.load(truth)["sed_ids"]
        inrange = (ids >= 0) & (ids < len(sed))
        j = np.clip(ids, 0, len(sed) - 1)
        tid = np.where(inrange, sed[j, 0], ids)
        orig = np.where(inrange, sed[j, 1] if sed.shape[1] > 1 else sed[j, 0], ids)
    else:
        tid = ids; orig = ids
    nfold = int((tid < 0).sum())
    tid = np.abs(tid).astype(np.int64)   # larg4's -ancestor ids folded (doc 20 sec 6.1)
    orig = np.abs(orig).astype(np.int64)
    Da = D[act]
    faces = load_faces(wires)
    keys, dface = assign_depo_faces(Da, faces)
    x = np.full(len(Da), np.nan)
    for ki, k in enumerate(keys):
        sel = dface == ki
        if not sel.any():
            continue
        xo = faces[k]["xorig"]
        dirx = float(np.sign(np.median(Da[sel, 2] - xo))) or 1.0
        x[sel] = xo + dirx * (Da[sel, 0] / US + depo_toff_us + reco_toff_us) * speed
    P = np.column_stack([x, Da[:, 3], Da[:, 4]])
    good = np.isfinite(x)
    order = np.lexsort((tkey[good], orig[good], tid[good]))
    tg, og, Pg, tk = tid[good][order], orig[good][order], P[good][order], tkey[good][order]
    brk = np.flatnonzero((np.diff(tg) != 0) | (np.diff(og) != 0)) + 1
    st = np.r_[0, brk]
    en = np.r_[st[1:], len(tg)]
    units = {}
    for s0, e0 in zip(st, en):
        V = Pg[s0:e0]
        keep = np.r_[True, np.any(np.diff(V, axis=0) != 0, axis=1)]   # coincident consecutive vertices merged
        units.setdefault(int(tg[s0]), []).append((int(og[s0]), float(tk[s0]), Polyline(V[keep])))
    out = {k: KeyTraj(k, [u for _, _, u in us], [o for o, _, _ in us], [t for _, t, _ in us]) for k, us in units.items()}
    return out, {"ndepos": int(good.sum()), "ndepos_unplaced": int((~good).sum()), "ndepos_folded": nfold,
                 "nunits": int(sum(len(v.units) for v in out.values()))}


class Field:
    """The distance field over all trajectories: nearest track key and distance of query points."""

    def __init__(self, trajs, reach):
        self.trajs = trajs
        self.keys = np.array(sorted(trajs), dtype=np.int64)
        pts, lab = [], []
        for i, k in enumerate(self.keys):
            pts.append(trajs[k].V); lab.append(np.full(len(trajs[k].V), i))
        self.P = np.concatenate(pts); self.lab = np.concatenate(lab)
        self.kd = cKDTree(self.P)
        self.maxseg = max((float(u.L.max()) for t in trajs.values() for u in t.units if len(u.L)), default=0.0)
        self.reach = reach + 0.5 * self.maxseg

    def nearest(self, Q, kc=8):
        """candidate tracks = the tracks of the kc nearest depo vertices within reach (depo steps are sub-mm, so
        the nearest trajectory always owns one of them); exact polyline distance per candidate track."""
        Q = np.atleast_2d(Q)
        kstar = np.full(len(Q), -1, dtype=np.int64)
        dstar = np.full(len(Q), np.inf)
        k = min(kc, len(self.P))
        dd, ii = self.kd.query(Q, k=k, distance_upper_bound=self.reach)
        dd, ii = dd.reshape(len(Q), k), ii.reshape(len(Q), k)
        ok = np.isfinite(dd)
        qi = np.repeat(np.arange(len(Q)), k).reshape(len(Q), k)[ok]
        pairs = np.unique(np.column_stack([self.lab[ii[ok]], qi]), axis=0)
        if not len(pairs):
            return kstar, dstar
        st = np.r_[0, np.flatnonzero(np.diff(pairs[:, 0])) + 1]
        en = np.r_[st[1:], len(pairs)]
        for s0, e0 in zip(st, en):
            qs = pairs[s0:e0, 1]
            d = self.trajs[int(self.keys[pairs[s0, 0]])].nearest(Q[qs])[0]
            better = d < dstar[qs]
            dstar[qs[better]] = d[better]; kstar[qs[better]] = self.keys[pairs[s0, 0]]
        return kstar, dstar


def true_vertices(trajs, lmin, dv):
    """starts of tracks (length >= lmin) that another such track passes within dv of."""
    keys = [k for k in sorted(trajs) if trajs[k].length >= lmin]
    out = []
    for k in keys:
        c = trajs[k].start[None]
        for k2 in keys:
            if k2 != k and trajs[k2].nearest(c)[0][0] <= dv:
                out.append((k, trajs[k].start)); break
    return out


def nu_vertices(truth, trajs):
    if truth is None:
        return []
    z = np.load(truth)
    pf = z["truth_pf"]
    c = T.pf_col
    out = []
    for nu in sorted(set(int(v) for v in pf[:, c("nu_row")] if v >= 0)):
        prim = pf[(pf[:, c("nu_row")] == nu) & (pf[:, c("mother_trackid")] == 0)]
        S = [trajs[int(t)].start for t in prim[:, c("trackid")].astype(int) if int(t) in trajs]
        if S:
            out.append(np.median(np.array(S), axis=0))
    return out


def seg_dist(P, A, B):
    """distance of points P (n,3) to segments A-B (n,3) row by row."""
    d = B - A
    L2 = np.maximum((d * d).sum(1), 1e-12)
    u = np.clip(((P - A) * d).sum(1) / L2, 0, 1)
    return np.linalg.norm(P - (A + d * u[:, None]), axis=1)


# ---------------------------------------------------------------------------- labels
def label_cluster(sp, st, se, coords, trajs, field, prm, vtx, nuv):
    X = np.column_stack([np.asarray(sp[c], float) for c in coords])
    shift = np.column_stack([np.asarray(st[k], float) for k in ("dx_cor", "dy_cor", "dz_cor")])
    R = X - shift                                   # vertices in the raw (labeler) frame
    kstar, dstar = field.nearest(R)
    nv = len(R)
    isv = np.zeros(nv, dtype=int); isnu = np.zeros(nv, dtype=int)
    if vtx:
        dv, _ = cKDTree(np.array([p for _, p in vtx])).query(R)
        isv = (dv <= prm["dv"]).astype(int)
    if nuv:
        dn, _ = cKDTree(np.array(nuv)).query(R)
        isnu = (dn <= prm["dv"]).astype(int)
    nodes = dict(kstar=kstar, dstar=dstar, is_vertex=isv, is_nu_vertex=isnu)
    if se is None:
        return R, nodes, None
    a = np.asarray(se["src"]).astype(int); b = np.asarray(se["dst"]).astype(int)
    ne = len(a)
    cls4 = np.full(ne, -1, dtype=int); flow = np.zeros(ne, dtype=int)
    on, off = prm["tau_on"], prm["tau_off"]
    ka, kb, da, db = kstar[a], kstar[b], dstar[a], dstar[b]
    cls4[(da > off) | (db > off)] = 3
    both = (da < on) & (db < on)
    same = both & (ka == kb) & (ka >= 0)
    diff = both & (ka != kb) & (ka >= 0) & (kb >= 0)
    # same track: chord samples and arc length
    for k in np.unique(ka[same]):
        es = np.flatnonzero(same & (ka == k))
        tr = trajs[int(k)]
        A, B = R[a[es]], R[b[es]]
        L = np.linalg.norm(B - A, axis=1)
        ns = np.maximum(2, np.ceil(L / prm["step"]).astype(int) + 1)
        rep = np.repeat(np.arange(len(es)), ns)
        f = (np.arange(len(rep)) - np.repeat(np.r_[0, np.cumsum(ns)[:-1]], ns)) / (ns[rep] - 1)
        S = A[rep] + (B - A)[rep] * f[:, None]
        d, s, t, u = tr.nearest(S)
        worst = np.zeros(len(es)); np.maximum.at(worst, rep, d)
        first = np.r_[0, np.cumsum(ns)[:-1]]; last = np.cumsum(ns) - 1
        # arc length only along one G4 particle: ends on two daughters of a folded key skip the U-turn test
        one = u[first] == u[last]
        darc = np.where(one, np.abs(s[last] - s[first]), 0.0)
        ok = (worst < on) & (darc <= 1.5 * L + prm["pitch"])
        cls4[es[ok]] = 0
        # flow: tangent at the midpoint vs the edge direction, when both ends are on the same particle
        dm, sm, tm, um = tr.nearest(0.5 * (A + B))
        sg = np.sign(np.einsum("ij,ij->i", tm, B - A)).astype(int)
        flow[es[ok & one]] = sg[ok & one]
    # different tracks: junction or false bridge (a start of one track within dv of the other track and of the chord)
    if diff.any():
        ed = np.flatnonzero(diff)
        near = {}   # (kx, ky) -> distance of the start of kx to the trajectory of ky
        jn = np.zeros(len(ed), dtype=bool)
        for kx_arr, ky_arr in ((ka[ed], kb[ed]), (kb[ed], ka[ed])):
            for kx, ky in set(zip(kx_arr.tolist(), ky_arr.tolist())):
                if (kx, ky) not in near:
                    near[(kx, ky)] = trajs[ky].nearest(trajs[kx].start[None])[0][0]
            dpair = np.array([near[(x, y)] for x, y in zip(kx_arr.tolist(), ky_arr.tolist())])
            C = np.array([trajs[x].start for x in kx_arr.tolist()])
            dchord = seg_dist(C, R[a[ed]], R[b[ed]])
            jn |= (dpair <= prm["dv"]) & (dchord <= prm["dv"])
        cls4[ed] = np.where(jn, 1, 2)
    edges = dict(src=a, dst=b, cls4=cls4, flow=flow, length=np.linalg.norm(R[b] - R[a], axis=1))
    return R, nodes, edges


# ---------------------------------------------------------------------------- metrics
def sample_traj(tr, step):
    """points every <= step along the polyline and the arc length each one stands for."""
    if not len(tr.L):
        return tr.V[:1], np.zeros(1)
    n = np.maximum(1, np.ceil(tr.L / step).astype(int))
    rep = np.repeat(np.arange(len(tr.L)), n)
    f = (np.arange(len(rep)) - np.repeat(np.r_[0, np.cumsum(n)[:-1]], n) + 0.5) / n[rep]
    return tr.A[rep] + (tr.B - tr.A)[rep] * f[:, None], (tr.L / n)[rep]


def edge_samples(R, a, b, step):
    A, B = R[a], R[b]
    L = np.linalg.norm(B - A, axis=1)
    n = np.maximum(1, np.ceil(L / step).astype(int))
    rep = np.repeat(np.arange(len(a)), n)
    f = (np.arange(len(rep)) - np.repeat(np.r_[0, np.cumsum(n)[:-1]], n) + 0.5) / n[rep]
    return A[rep] + (B - A)[rep] * f[:, None], (L / n)[rep], rep


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pctree", required=True)
    ap.add_argument("--depos", required=True)
    ap.add_argument("--truth", default="")
    ap.add_argument("--wires", required=True)
    ap.add_argument("--depo-toff-us", type=float, required=True)
    ap.add_argument("--reco-toff-us", type=float, default=0.0)
    ap.add_argument("--speed-mmus", type=float, required=True)
    ap.add_argument("--pitch-mm", type=float, required=True)
    ap.add_argument("--tau-on-mm", type=float, default=15.0, help="tau_on in mm, one band for all detectors (doc 23 sec 7.1)")
    ap.add_argument("--tau-off-mm", type=float, default=40.0, help="tau_off in mm")
    ap.add_argument("--tau-on-pitch", type=float, default=None, help="tau_on in wire pitches instead (doc 12 sec 4.1: 2)")
    ap.add_argument("--tau-off-pitch", type=float, default=None, help="tau_off in wire pitches instead (doc 12: 5)")
    ap.add_argument("--d-cap-mm", type=float, default=25.0, help="the C++ labeler's d_traj cap (cross-check)")
    ap.add_argument("--event", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    tau_on = a.tau_on_mm if a.tau_on_pitch is None else a.tau_on_pitch * a.pitch_mm
    tau_off = a.tau_off_mm if a.tau_off_pitch is None else a.tau_off_pitch * a.pitch_mm
    prm = dict(pitch=a.pitch_mm, tau_on=tau_on, tau_off=tau_off, dv=1 * CM, step=3.0)
    truth = a.truth or None

    try:
        tree = load_tree(a.pctree)
    except RuntimeError as err:   # an event whose imaging made no blob (R2: PDHD beam2, doc 20)
        json.dump({"event": a.event, "error": str(err)}, open(os.path.join(a.out, a.event + ".json"), "w"), indent=1)
        print(a.event, "no tree:", err); return
    trajs, tinfo = build_trajectories(a.depos, truth, a.wires, a.depo_toff_us, a.reco_toff_us, a.speed_mmus)
    # edge classes and vertices at the G4-particle level (doc 12 sec 4.1: shower daughters keep their own id so
    # fragment-level labels exist); the merged track key of each particle is kept beside it
    ptrajs, pkey = {}, {}
    for k, tr in trajs.items():
        for o, u in zip(tr.origs, tr.units):
            ptrajs[o] = KeyTraj(o, [u], [o], [0.0]); pkey[o] = k
    field = Field(ptrajs, prm["tau_off"])
    vtx = true_vertices(ptrajs, 1 * CM, prm["dv"])
    nuv = nu_vertices(truth, trajs)
    SP, ST, SE = per_node(tree, "steiner_pc"), per_node(tree, "steiner_truth"), per_node(tree, "steiner_edges")
    res = {"event": a.event, "params": prm, "inputs": vars(a), **tinfo,
           "ntracks": len(trajs), "nparticles": len(ptrajs), "max_segment_mm": field.maxseg, "n_true_vertices": len(vtx), "n_nu_vertices": len(nuv)}
    if not ST:
        res["error"] = "no steiner_truth cloud in the dump"
        json.dump(res, open(os.path.join(a.out, a.event + ".json"), "w"), indent=1)
        print(a.event, res["error"]); return
    keys0 = list(SP[next(iter(SP))].keys())
    coords = [next(k for k in keys0 if k.startswith(ax) and k not in ("wpid",)) for ax in ("x", "y", "z")]
    res["coords"] = coords

    # cluster points (raw frame) of the clusters with a graph, and their C++ labels
    P3, TR = per_node(tree, "3d"), per_node(tree, "truth")
    clusters = sorted(ST)
    cp, ctid, cdt = [], [], []
    for node, pc in P3.items():
        if tree["cluster_of"][node] in ST and node in TR:
            cp.append(np.column_stack([pc["x"], pc["y"], pc["z"]]).astype(float))
            ctid.append(np.asarray(TR[node]["tid"]).astype(int)); cdt.append(np.asarray(TR[node]["d_traj"], float))
    CP = np.concatenate(cp) if cp else np.zeros((0, 3))
    CT = np.concatenate(ctid) if ctid else np.zeros(0, int)
    CD = np.concatenate(cdt) if cdt else np.zeros(0)
    res["npoints_in_scope"] = int(len(CP))

    # cross-check of the distance field against the C++ per-point labels
    if len(CP):
        lab = CT > 0
        dpy = np.full(len(CP), np.inf)
        for k in np.unique(CT[lab]):
            sel = np.flatnonzero(CT == k)
            if int(k) in trajs:
                dpy[sel] = trajs[int(k)].nearest(CP[sel])[0]
        cap = a.d_cap_mm
        multi = np.array([int(k) in trajs and len(trajs[int(k)].units) > 1 for k in CT], dtype=bool)
        res["cross_check_d_traj"] = {}
        # single-particle keys: the same polyline as the C++ labeler; folded keys: the C++ labeler orders all the
        # daughters' depos into one polyline (zigzag), so its d_traj is <= the per-particle distance
        for nm, sel0 in (("single", ~multi), ("folded", multi)):
            both = lab & np.isfinite(dpy) & sel0
            dd = np.abs(np.minimum(dpy[both], cap) - np.minimum(CD[both], cap))
            res["cross_check_d_traj"][nm] = {
                "npoints": int(both.sum()), "frac_within_0.1mm": float((dd <= 0.1).mean()) if both.any() else None,
                "frac_within_1mm": float((dd <= 1.0).mean()) if both.any() else None,
                "python_smaller_by_gt_1mm": int(((np.minimum(dpy[both], cap) + 1.0) < np.minimum(CD[both], cap)).sum()),
                "python_larger_by_gt_1mm": int(((np.minimum(dpy[both], cap) - 1.0) > np.minimum(CD[both], cap)).sum())}

    # labels per cluster
    allR, allE, nodes_out, edges_out = [], [], {}, {}
    agree = {"n_on": 0, "cpp_eq_kstar": 0}
    off = 0
    for cl in clusters:
        sp, st = SP[cl], ST[cl]
        se = SE.get(cl)
        R, nodes, edges = label_cluster(sp, st, se, coords, ptrajs, field, prm, vtx, nuv)
        nodes["key"] = np.array([pkey.get(int(p), -1) for p in nodes["kstar"]], dtype=np.int64)
        if edges is not None:
            edges["same_key"] = (nodes["key"][edges["src"]] == nodes["key"][edges["dst"]]).astype(int)
        cpp = np.asarray(st["tid"]).astype(int)
        on = nodes["dstar"] < prm["tau_on"]
        agree["n_on"] += int(on.sum()); agree["cpp_eq_kstar"] += int((cpp[on] == nodes["key"][on]).sum())
        nodes_out[cl] = dict(nodes, tid_cpp=cpp, terminal=np.asarray(sp["flag_steiner_terminal"]).astype(int), xyz_raw=R)
        allR.append(R)
        if edges is not None:
            edges_out[cl] = edges
            allE.append((R, edges))
        off += len(R)
    res["node_key_agreement"] = dict(agree, frac=agree["cpp_eq_kstar"] / max(agree["n_on"], 1))
    nn = sum(len(v["kstar"]) for v in nodes_out.values())
    res["nodes"] = {"n": nn, "n_terminal": int(sum(v["terminal"].sum() for v in nodes_out.values())),
                    "frac_within_tau_on": float(sum((v["dstar"] < prm["tau_on"]).sum() for v in nodes_out.values()) / max(nn, 1)),
                    "frac_beyond_tau_off": float(sum((v["dstar"] > prm["tau_off"]).sum() for v in nodes_out.values()) / max(nn, 1)),
                    "n_is_vertex": int(sum(v["is_vertex"].sum() for v in nodes_out.values())),
                    "n_is_nu_vertex": int(sum(v["is_nu_vertex"].sum() for v in nodes_out.values()))}
    names = {0: "on_trajectory", 1: "junction", 2: "false_bridge", 3: "ghost", -1: "dont_care"}
    ec = {v: {"n": 0, "length_cm": 0.0} for v in names.values()}
    for R, e in allE:
        for c, nm in names.items():
            m = e["cls4"] == c
            ec[nm]["n"] += int(m.sum()); ec[nm]["length_cm"] += float(e["length"][m].sum() / CM)
    res["edges"] = ec
    res["flow"] = {"plus": int(sum((e["flow"] > 0).sum() for _, e in allE)), "minus": int(sum((e["flow"] < 0).sum() for _, e in allE))}

    # skeleton metrics (doc 12 sec 4.4)
    if allE and len(CP):
        skel = np.concatenate([edge_samples(R, e["src"], e["dst"], prm["step"])[0] for R, e in allE])
        skw = np.concatenate([edge_samples(R, e["src"], e["dst"], prm["step"])[1] for R, e in allE])
        sk_kd = cKDTree(skel)
        cp_kd = cKDTree(CP)
        tot, cov = 0.0, {0.5: 0.0, 1.0: 0.0, 2.0: 0.0}
        per_track = []
        for k, tr in sorted(ptrajs.items()):
            sw = [sample_traj(un, prm["step"]) for un in tr.units]
            S, w = np.concatenate([x for x, _ in sw]), np.concatenate([y for _, y in sw])
            dcp, _ = cp_kd.query(S)
            ins = dcp <= 2 * CM
            if not ins.any():
                continue
            ds, _ = sk_kd.query(S[ins])
            L = float(w[ins].sum()); tot += L
            row = {"tid": k, "length_in_scope_cm": L / CM}
            for d in cov:
                c = float(w[ins][ds <= d * CM].sum()); cov[d] += c
                row["R_%g" % d] = c / L if L else None
            per_track.append(row)
        res["coverage_recall"] = {"R_%gcm" % d: (cov[d] / tot if tot else None) for d in cov}
        res["true_length_in_scope_cm"] = tot / CM
        kst, dst = field.nearest(skel)
        res["precision"] = {"P_%gcm" % d: float(skw[dst <= d * CM].sum() / skw.sum()) for d in (0.5, 1.0, 2.0)}
        # true track ends vs the nearest Steiner node / terminal (the stored graph keeps cycles, so an end is
        # rarely a leaf); true vertices vs nodes of degree >= 3
        allnodes, terms, junc = [], [], []
        for cl, v in nodes_out.items():
            allnodes.append(v["xyz_raw"]); terms.append(v["xyz_raw"][v["terminal"] > 0])
        for R, e in allE:
            deg = np.bincount(np.r_[e["src"], e["dst"]], minlength=len(R))
            junc.append(R[deg >= 3])
        allnodes, terms = np.concatenate(allnodes), np.concatenate(terms)
        junc = np.concatenate(junc) if junc else np.zeros((0, 3))
        ends = []
        for k, tr in trajs.items():
            for un in tr.units:
                if un.length < 3 * CM:
                    continue
                for p in (un.V[0], un.V[-1]):
                    if cp_kd.query(p)[0] <= 2 * CM:
                        ends.append(p)
        if ends:
            E = np.array(ends)
            dn = cKDTree(allnodes).query(E)[0]
            dt = cKDTree(terms).query(E)[0] if len(terms) else np.full(len(E), np.inf)
            res["endpoints"] = {"n": len(ends), "to_node_median_cm": float(np.median(dn) / CM), "to_node_le_1cm": float((dn <= CM).mean()),
                                "to_terminal_median_cm": float(np.median(dt) / CM), "to_terminal_le_1cm": float((dt <= CM).mean()),
                                "to_terminal_le_2cm": float((dt <= 2 * CM).mean())}
        vin = [p for _, p in vtx if cp_kd.query(p)[0] <= 2 * CM]
        if vin:
            dj = cKDTree(junc).query(np.array(vin))[0] if len(junc) else np.full(len(vin), np.inf)
            res["junction_recall"] = {"n_true_vertices": len(vin), "le_1cm": float((dj <= CM).mean()), "le_2cm": float((dj <= 2 * CM).mean())}
        nterm = res["nodes"]["n_terminal"]
        res["terminals_per_cm_true"] = nterm / (tot / CM) if tot else None
        json.dump(per_track, open(os.path.join(a.out, a.event + "_tracks.json"), "w"), indent=1)

    # labels file: one row per Steiner vertex / edge, cluster node index as the join key
    out = {}
    for cl, v in nodes_out.items():
        for k, x in v.items():
            out.setdefault("node_" + k, []).append(x)
        out.setdefault("node_cluster", []).append(np.full(len(v["kstar"]), cl))
        out.setdefault("node_row", []).append(np.arange(len(v["kstar"])))
    for cl, e in edges_out.items():
        for k, x in e.items():
            out.setdefault("edge_" + k, []).append(x)
        out.setdefault("edge_cluster", []).append(np.full(len(e["src"]), cl))
    np.savez_compressed(os.path.join(a.out, a.event + "_labels.npz"), **{k: np.concatenate(v) for k, v in out.items()})
    json.dump(res, open(os.path.join(a.out, a.event + ".json"), "w"), indent=1)
    cc = res.get("cross_check_d_traj", {}).get("single", {})
    print(a.event, "clusters", len(clusters), "nodes", nn, "on", "%.3f" % res["nodes"]["frac_within_tau_on"],
          "key agree %.4f" % res["node_key_agreement"]["frac"],
          "edges", {k: v["n"] for k, v in ec.items()},
          "R1cm", res.get("coverage_recall", {}).get("R_1cm"), "P1cm", res.get("precision", {}).get("P_1cm"),
          "xcheck<=0.1mm", cc.get("frac_within_0.1mm"))


if __name__ == "__main__":
    main()
