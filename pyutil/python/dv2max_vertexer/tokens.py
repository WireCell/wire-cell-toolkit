"""Event-level physics analysis and token helpers (truth-free).

    an = analyze(event)                       # EventAnalysis: trees, PR graph (oriented), cluster stats
    et = end_tokens(an)                       # one row per segment END: position, inward dir, END_FEATURES, p_start
    pos, kind = candidate_positions(an)       # PR cands, tag_nu, WC reco starts, PR nodes, segment ends,
                                              # cluster axis ends
    D = describe_positions(an, pos)           # (m, len(CAND_FEATURES)) local physics descriptors
    tt = trajectory_tokens(an)                # 3 cm oriented pieces of fit segments (dQ/dx mean/slope, rr)
    wt = wc_particle_tokens(an)               # Wire-Cell reco particles: PID, energy, start descriptors
    xyz, info = heuristic_vertex(an)          # hand-made truth-free vertex pick

`event` must be truth-free (event.strip_truth()); nothing here reads event.truth.
"""
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from . import fitgraph as fg
from . import local
from .geometry import boundary_distance

NO_TAG_CM = fg.NO_TAG_CM
DIRS = np.array([[a, b, c] for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) if (a, b, c) != (0, 0, 0)],
                np.float64)
DIRS /= np.linalg.norm(DIRS, axis=1, keepdims=True)
SHELL_CM = (2.0, 8.0)
DIR_MIN_Q = 3e4            # charge (~5 points) needed to call a direction bin occupied
END_R_CM = 3.0             # segment ends this close to a position count as attached prongs
END_MIN_LEN_CM = 1.0       # ignore 2-point stub segments of small clusters when counting prongs
POS_CHUNK = 128            # query positions per 10 cm-ball chunk (bounds memory)
KINDS = ("cand", "tag", "reco_start", "node", "end", "cluster_end")
CLUSTER_END_MIN_POINTS = 20   # clusters this big contribute their 2 axis extremes + upstream point

CAND_FEATURES = (
    "on_charge",     # log1p(distance to the nearest clustering_pr point)
    "n1", "n3", "n10",            # log1p(point count within 1/3/10 cm)
    "q1", "q3", "q10",            # log1p(charge within 1/3/10 cm / 1e4)
    "one_sided",     # |charge-weighted mean unit vector| of points in the 1.5-10 cm shell (1 = all on one side)
    "shell_dz",      # z component of that mean vector (> 0: charge lies downstream)
    "fwd_frac",      # fraction of the 10 cm ball charge downstream (z > z0)
    "n_dirs",        # occupied direction bins (26) in the 2-8 cm shell: a prong-count proxy
    "node_d",        # log1p(distance to the nearest PR-graph node)
    "node_deg",      # degree of that node if within 1 cm, else 0
    "n_ends3",       # segment ends within 3 cm (segments >= 1 cm)
    "n_out",         # sum over those ends of P(end is the segment start): prongs leaving
    "n_in",          # sum of 1 - P(start): prongs ending here (Bragg / stopping side)
    "n_out_bragg",   # ends with P(start) > 0.7 and L > 3 cm: Bragg-consistent outgoing prongs
    "len_out",       # log1p(sum of L * P(start) over those ends)
    "max_len",       # log1p(longest segment with an end within 3 cm)
    "bragg_in",      # max bragg_mean3 over ends within 3 cm (a Bragg peak arriving here)
    "is_wc",         # nearest point's cluster is a Wire-Cell reco-neutrino cluster
    "cl_size", "cl_q", "cl_rank",  # log1p(size), log1p(total q / 1e4), log1p(charge rank in event)
    "cl_faces",      # TPC faces touched by the cluster (0-6)
    "cl_frac_z",     # position along the cluster axis (0 = upstream end)
    "cl_dz",         # log1p(z - cluster z_min)
    "d_tag",         # log1p(distance to tag_nu); log1p(1000) if none
    "n_cand2",       # PR vertex candidates within 2 cm
    "d_wall",        # log1p(distance to the nearest crossable TPC boundary)
    "shw_frac3",     # shower-tagged fraction of the 3 cm charge
    "trk_frac3",     # track-tagged fraction of the 3 cm charge
    "wc_starts",     # Wire-Cell reco particle starts within 1 cm
    "wc_vtx",        # a Wire-Cell reco-neutrino vertex within 0.5 cm
)
CF = {n: k for k, n in enumerate(CAND_FEATURES)}

END_TOKEN_FIELDS = ("xyz", "dir", "seg", "end", "length", "straight", "log_med_dqdx", "p_start", "node",
                    "node_degree", "cluster", "is_wc", "feat")


@dataclass
class EventAnalysis:
    event: object
    X: np.ndarray              # (n,3) float64 clustering_pr points
    Q: np.ndarray              # (n,) float64
    tree: object               # cKDTree over X (None if no points)
    shower: np.ndarray         # (n,) int8 0/1/2 (none / track / shower tag)
    cl_point: dict             # per-point cluster stats (see local.cluster_stats)
    cl: dict                   # per-cluster stats
    cl_inv: np.ndarray         # (n,) cluster index
    graph: fg.FitGraph         # parsed, context-filled and oriented PR graph
    cand: np.ndarray           # (m,3) float64 PR candidates
    tag: np.ndarray            # (3,) float64 (NaN if none)
    wc_clusters: tuple
    wc_vtx: np.ndarray         # (k,3) Wire-Cell reco-nu vertices
    wc_starts: np.ndarray      # (p,3) Wire-Cell reco particle starts


def analyze(event):
    """Build the truth-free EventAnalysis for an Event (see module docstring)."""
    if getattr(event, "truth", None) is not None:
        raise ValueError("analyze() needs a truth-free event (use event.strip_truth())")
    X = np.asarray(event.points_xyz, np.float64).reshape(-1, 3)
    Q = np.maximum(np.asarray(event.points_q, np.float64), 0.0)
    tree = cKDTree(X) if len(X) else None
    shower = local.shower_tags(X, tree, event.shower_xyz, event.shower_q) if tree is not None else \
        np.zeros(0, np.int8)
    wc_clusters = tuple(int(w.cluster_id) for w in event.wc_reco)
    cl_point, cl, cl_inv = local.cluster_stats(X, Q, np.asarray(event.points_cluster_id), wc_clusters)
    g = fg.build_fit_graph(event.fit_xyz, event.fit_q, event.fit_cluster_id)
    fg.add_context(g, X, Q, tree, shower, event.cand_xyz, event.tag_nu)
    fg.orient(g)
    wv = [np.asarray(w.vtx, np.float64) for w in event.wc_reco]
    ws = [np.asarray(p.start, np.float64) for w in event.wc_reco for p in w.particles]
    return EventAnalysis(event, X, Q, tree, shower, cl_point, cl, cl_inv, g,
                         np.asarray(event.cand_xyz, np.float64).reshape(-1, 3),
                         np.asarray(event.tag_nu, np.float64).reshape(3), wc_clusters,
                         np.array(wv).reshape(-1, 3), np.array(ws).reshape(-1, 3))


def end_tokens(an):
    """One row per segment end (2 per segment).  dict of arrays; 'feat' is (2S, len(END_FEATURES)) with
    column names fitgraph.END_FEATURES; 'p_start' = P(this end is the start of its segment)."""
    g = an.graph
    S = g.n_seg
    ends = np.tile(np.array([0, 1]), S)
    seg = np.repeat(np.arange(S), 2)
    ps = np.stack([g.p_start, 1.0 - g.p_start], 1).reshape(-1) if S else np.zeros(0)
    node = g.end_node.reshape(-1)
    cl = g.seg_cluster[seg] if S else np.zeros(0, np.int64)
    return {"xyz": g.end_xyz.reshape(-1, 3), "dir": g.end_dir.reshape(-1, 3), "seg": seg, "end": ends,
            "length": g.length[seg], "straight": g.straight[seg],
            "log_med_dqdx": np.log(np.maximum(g.med_dqdx[seg], 1.0) / fg.MIP_DQDX), "p_start": ps, "node": node,
            "node_degree": g.node_degree[node] if S else np.zeros(0, np.int64), "cluster": cl,
            "is_wc": np.isin(cl, np.asarray(an.wc_clusters, np.int64)).astype(np.float64),
            "feat": g.end_feat.reshape(-1, len(fg.END_FEATURES))}


def cluster_end_positions(an):
    """Axis extremes and upstream (min z) point of every cluster with >= CLUSTER_END_MIN_POINTS points,
    de-duplicated: candidates on clusters that Wire-Cell PR did not process."""
    cl = an.cl
    if not cl:
        return np.zeros((0, 3))
    big = cl["size"] >= CLUSTER_END_MIN_POINTS
    idx = np.unique(np.concatenate([cl["end_lo"][big], cl["end_hi"][big], cl["z_lo"][big]]))
    return an.X[idx]


def candidate_positions(an, kinds=KINDS):
    """-> (xyz (m,3), kind (m,) int index into KINDS) for the requested token kinds."""
    g = an.graph
    parts = {"cand": an.cand, "tag": an.tag[None] if np.all(np.isfinite(an.tag)) else np.zeros((0, 3)),
             "reco_start": an.wc_starts, "node": g.node_xyz[g.node_degree >= 1] if g.n_node else np.zeros((0, 3)),
             "end": g.end_xyz.reshape(-1, 3), "cluster_end": cluster_end_positions(an)}
    xs, ks = [], []
    for k in kinds:
        p = np.asarray(parts[k], np.float64).reshape(-1, 3)
        xs.append(p)
        ks.append(np.full(len(p), KINDS.index(k), np.int64))
    return (np.concatenate(xs) if xs else np.zeros((0, 3))), (np.concatenate(ks) if ks else np.zeros(0, np.int64))


def _ball_descriptors(an, P, D):
    """Fill the charge-ball columns of D (rows for positions P) from points within 10 cm."""
    m = len(P)
    M = cKDTree(P).sparse_distance_matrix(an.tree, 10.0, output_type="ndarray")
    i, j, dist = M["i"], M["j"], M["v"]
    q = an.Q[j]
    for r, cn, qn in ((1.0, "n1", "q1"), (3.0, "n3", "q3"), (10.0, "n10", "q10")):
        k = dist <= r
        D[:, CF[cn]] = np.log1p(np.bincount(i[k], minlength=m))
        D[:, CF[qn]] = np.log1p(np.bincount(i[k], q[k], m) / local.Q_UNIT)
    dv = an.X[j] - P[i]
    q10 = np.maximum(np.bincount(i, q, m), 1e-9)
    D[:, CF["fwd_frac"]] = np.bincount(i, q * (dv[:, 2] > 0), m) / q10
    sh = dist > 1.5
    u = dv[sh] / dist[sh][:, None]
    qs = q[sh]
    qsum = np.maximum(np.bincount(i[sh], qs, m), 1e-9)
    mu = np.stack([np.bincount(i[sh], qs * u[:, a], m) for a in range(3)], 1) / qsum[:, None]
    D[:, CF["one_sided"]] = np.linalg.norm(mu, axis=1)
    D[:, CF["shell_dz"]] = mu[:, 2]
    sh2 = (dist >= SHELL_CM[0]) & (dist <= SHELL_CM[1])
    if np.any(sh2):
        ub = dv[sh2] / dist[sh2][:, None]
        b = np.argmax(ub @ DIRS.T, axis=1)
        occ = np.bincount(i[sh2] * len(DIRS) + b, q[sh2], m * len(DIRS)).reshape(m, len(DIRS))
        D[:, CF["n_dirs"]] = (occ >= DIR_MIN_Q).sum(1)
    k3 = dist <= 3.0
    q3 = np.maximum(np.bincount(i[k3], q[k3], m), 1e-9)
    tg = an.shower[j]
    D[:, CF["shw_frac3"]] = np.bincount(i[k3], q[k3] * (tg[k3] == 2), m) / q3
    D[:, CF["trk_frac3"]] = np.bincount(i[k3], q[k3] * (tg[k3] == 1), m) / q3
    return D


def describe_positions(an, P):
    """(m,3) query positions -> (m, len(CAND_FEATURES)) float64 local physics descriptors."""
    P = np.asarray(P, np.float64).reshape(-1, 3)
    m = len(P)
    D = np.zeros((m, len(CAND_FEATURES)))
    if m == 0:
        return D
    D[:, CF["d_tag"]] = np.log1p(np.linalg.norm(P - an.tag, axis=1)) if np.all(np.isfinite(an.tag)) else \
        np.log1p(NO_TAG_CM)
    D[:, CF["d_wall"]] = np.log1p(boundary_distance(P))
    if len(an.cand):
        D[:, CF["n_cand2"]] = cKDTree(an.cand).query_ball_point(P, 2.0, return_length=True)
    if len(an.wc_starts):
        D[:, CF["wc_starts"]] = cKDTree(an.wc_starts).query_ball_point(P, 1.0, return_length=True)
    if len(an.wc_vtx):
        D[:, CF["wc_vtx"]] = (cKDTree(an.wc_vtx).query(P)[0] < 0.5).astype(np.float64)
    ptree = cKDTree(P)
    if an.tree is not None:
        dmin, jn = an.tree.query(P)
        D[:, CF["on_charge"]] = np.log1p(dmin)
        c = an.cl_inv[jn]
        cl = an.cl
        D[:, CF["is_wc"]] = cl["is_wc"][c]
        D[:, CF["cl_size"]] = np.log1p(cl["size"][c])
        D[:, CF["cl_q"]] = np.log1p(cl["qtot"][c] / local.Q_UNIT)
        D[:, CF["cl_rank"]] = np.log1p(cl["rank"][c])
        D[:, CF["cl_faces"]] = cl["faces"][c]
        proj = np.einsum("ij,ij->i", P - cl["center"][c], cl["e1"][c])
        D[:, CF["cl_frac_z"]] = np.clip((proj - cl["pmin"][c]) / np.maximum(cl["ext1"][c], 1e-6), 0, 1)
        D[:, CF["cl_dz"]] = np.log1p(np.maximum(P[:, 2] - cl["zmin"][c], 0.0))
        for c0 in range(0, m, POS_CHUNK):
            sl = slice(c0, min(c0 + POS_CHUNK, m))
            D[sl] = _ball_descriptors(an, P[sl], D[sl])
    g = an.graph
    if g.n_node:
        dn, jn = cKDTree(g.node_xyz).query(P)
        D[:, CF["node_d"]] = np.log1p(dn)
        D[:, CF["node_deg"]] = np.where(dn < 1.0, g.node_degree[jn], 0)
    et = end_tokens(an) if g.n_seg else None
    if et is not None and np.any(et["length"] >= END_MIN_LEN_CM):
        keep = et["length"] >= END_MIN_LEN_CM
        et = {k: v[keep] for k, v in et.items()}
        M = ptree.sparse_distance_matrix(cKDTree(et["xyz"]), END_R_CM, output_type="ndarray")
        i, j = M["i"], M["j"]
        ps = et["p_start"][j]
        L = et["length"][j]
        D[:, CF["n_ends3"]] = np.bincount(i, minlength=m)
        D[:, CF["n_out"]] = np.bincount(i, ps, m)
        D[:, CF["n_in"]] = np.bincount(i, 1.0 - ps, m)
        D[:, CF["n_out_bragg"]] = np.bincount(i, ((ps > 0.7) & (L > 3.0)).astype(np.float64), m)
        D[:, CF["len_out"]] = np.log1p(np.bincount(i, ps * L, m))
        mx = np.zeros(m)
        np.maximum.at(mx, i, L)
        D[:, CF["max_len"]] = np.log1p(mx)
        bi = np.full(m, -np.inf)
        np.maximum.at(bi, i, et["feat"][j, fg.EF["bragg_mean3"]])
        D[:, CF["bragg_in"]] = np.where(np.isfinite(bi), bi, 0.0)
    return D


def trajectory_tokens(an, chunk_cm=3.0):
    """Oriented track_fit segments cut into `chunk_cm` pieces along the arc length.  dict of arrays, one row
    per chunk: xyz (mean position), dir (unit, likely start -> stop), seg, n (fit points), log_dqdx (mean log
    smoothed dQ/dx / MIP), dqdx_slope (d log dQ/dx per cm toward the stop end), rr_start / rr_stop (mean cm),
    length (segment), p (P(orientation right) >= 0.5), is_start_chunk / is_stop_chunk."""
    g = an.graph
    keys = ("xyz", "dir", "seg", "n", "log_dqdx", "dqdx_slope", "rr_start", "rr_stop", "length", "p",
            "is_start_chunk", "is_stop_chunk")
    idx = np.flatnonzero(g.seg >= 0)
    if len(idx) == 0:
        return {k: np.zeros((0, 3)) if k in ("xyz", "dir") else np.zeros(0) for k in keys}
    sp = fg.segment_point_features(g)
    sg = g.seg[idx]
    rs = sp["rr_start"][idx]
    nmax = int(np.floor(g.length.max() / chunk_cm)) + 2
    ck = sg * nmax + np.minimum(np.floor(rs / chunk_cm), np.floor(g.length[sg] / chunk_cm)).astype(np.int64)
    uk, inv = np.unique(ck, return_inverse=True)
    inv = inv.reshape(-1)
    m = len(uk)
    cnt = np.bincount(inv, minlength=m).astype(np.float64)
    xyz = np.stack([np.bincount(inv, g.xyz[idx, a], m) for a in range(3)], 1) / cnt[:, None]
    dsum = np.stack([np.bincount(inv, sp["dir"][idx, a], m) for a in range(3)], 1)
    y = np.log(np.maximum(g.dqdx_s[idx], 1.0) / fg.MIP_DQDX)
    ym = np.bincount(inv, y, m) / cnt
    xm = np.bincount(inv, rs, m) / cnt
    sxx = np.bincount(inv, (rs - xm[inv]) ** 2, m)
    sxy = np.bincount(inv, (rs - xm[inv]) * (y - ym[inv]), m)
    seg = uk // nmax
    L = g.length[seg]
    piece = uk % nmax
    last = np.floor(L / chunk_cm)
    return {"xyz": xyz, "dir": np.where(np.linalg.norm(dsum, axis=1, keepdims=True) > 0,
                                        dsum / np.maximum(np.linalg.norm(dsum, axis=1, keepdims=True), 1e-9), 0.0),
            "seg": seg, "n": cnt, "log_dqdx": ym, "dqdx_slope": np.where(sxx > 1e-6, sxy / np.maximum(sxx, 1e-6), 0.0),
            "rr_start": xm, "rr_stop": L - xm, "length": L,
            "p": np.bincount(inv, sp["p"][idx], m) / cnt, "is_start_chunk": (piece == 0).astype(np.float64),
            "is_stop_chunk": (piece == last).astype(np.float64)}


WC_PID_CLASSES = (11, 13, 22, 211, 2212)   # |pdg| one-hot for Wire-Cell reco particles; other -> last column


def wc_particle_tokens(an):
    """One row per Wire-Cell reco particle (reco rows: legitimate inputs).  dict: xyz (start), nu (reco nu
    index), pid (p, len(WC_PID_CLASSES)+1) one-hot of |pdg|, log_e (log1p energy/MeV), is_nu_vtx (start
    within 0.5 cm of its reco-nu vertex), desc (p, len(CAND_FEATURES)) = describe_positions at the start."""
    rows = [(w.nu, p, w.vtx) for w in an.event.wc_reco for p in w.particles]
    k = len(WC_PID_CLASSES) + 1
    if not rows:
        return {"xyz": np.zeros((0, 3)), "nu": np.zeros(0, np.int64), "pid": np.zeros((0, k)),
                "log_e": np.zeros(0), "is_nu_vtx": np.zeros(0), "desc": np.zeros((0, len(CAND_FEATURES)))}
    xyz = np.array([np.asarray(p.start, np.float64) for _, p, _ in rows]).reshape(-1, 3)
    pid = np.zeros((len(rows), k))
    for r, (_, p, _) in enumerate(rows):
        a = abs(int(p.pdg))
        pid[r, WC_PID_CLASSES.index(a) if a in WC_PID_CLASSES else k - 1] = 1.0
    vt = np.array([np.asarray(v, np.float64) for _, _, v in rows]).reshape(-1, 3)
    return {"xyz": xyz, "nu": np.array([n for n, _, _ in rows], np.int64), "pid": pid,
            "log_e": np.log1p(np.maximum([p.energy for _, p, _ in rows], 0.0)),
            "is_nu_vtx": (np.linalg.norm(xyz - vt, axis=1) < 0.5).astype(np.float64),
            "desc": describe_positions(an, xyz)}


# Hand-made vertex-end score over CAND_FEATURES (truth-free): prongs attached (longest attached segment,
# node degree), non-Bragg end (few stopping ends), upstream (charge downstream, early along the cluster
# axis), little shower charge, on charge, on the Wire-Cell neutrino cluster.  Physics-set start, then a
# coarse coordinate search on the dev cache (547 target events; in-sample 1 cm eff 0.565 nue / 0.577 BNB,
# from 0.484 / 0.515 at the start); out-of-sample numbers from scripts/validate_features.py.
HEURISTIC_WEIGHTS = {
    "max_len": 0.5, "node_deg": 0.3, "n_in": -1.0, "shw_frac3": -1.0, "shell_dz": 2.0, "fwd_frac": 1.0,
    "cl_frac_z": -1.0, "is_wc": 2.0, "q3": 0.3, "on_charge": -2.0,
}


def heuristic_score(D, weights=None):
    """Linear vertex-end score for descriptor rows D (m, len(CAND_FEATURES))."""
    w = HEURISTIC_WEIGHTS if weights is None else weights
    s = np.zeros(len(D))
    for k, v in w.items():
        s += v * D[:, CF[k]]
    return s


def heuristic_vertex(an, kinds=("cand", "tag", "node", "end"), weights=None):
    """Truth-free pick: the candidate position with the highest heuristic_score.
    Returns (xyz (3,) or NaN, info dict with the pool, kinds, descriptors and scores)."""
    P, K = candidate_positions(an, kinds)
    if len(P) == 0:
        return np.full(3, np.nan), {"pos": P, "kind": K, "D": np.zeros((0, len(CAND_FEATURES))), "score": np.zeros(0)}
    D = describe_positions(an, P)
    s = heuristic_score(D, weights)
    return P[int(np.argmax(s))], {"pos": P, "kind": K, "D": D, "score": s}
