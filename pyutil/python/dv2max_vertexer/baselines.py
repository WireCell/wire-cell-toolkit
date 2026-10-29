"""Hypothesis pool and its truth-free features: vendored from sbndvtx/eval/baselines.py (inference subset:
SOURCES / MERGE_CM constants, Pool, build_pool, cluster_end_points, basic_features, physics_features,
CacheFeatureAggregates, DEFAULT_FEATURE_FNS, pool_features).  Function bodies unchanged; dropped: the baselines,
oracles, pool tables and the OOF ranker (pandas, cache and r3 I/O).

Hypothesis pool per event (build_pool): tag_nu, PR candidates, Wire-Cell reco-nu vertices, Wire-Cell reco particle
starts, both ends of every PR-graph segment of the track_fit layer and ('clusend') the PCA-axis extremes and upstream
point of every clustering_pr cluster with >= 20 points.  Hypotheses closer than MERGE_CM are merged.
"""
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

SOURCES = ("tag", "cand", "wcnu", "wcp", "fitend", "clusend")   # merge priority order
SPEC_SOURCES = ("tag", "cand", "wcnu", "wcp", "fitend")         # the pool named in the eval spec
# Default = with cluster ends.  Dev cache (547 events): pool oracle 0.687 -> 0.762 at 1 cm; OOF ranker
# -0.3 to +3.5 pts depending on the feature set (not settled at dev size; the full-cache run decides).
DEFAULT_SOURCES = SOURCES
CLUSTER_END_MIN_POINTS = 20   # clusend: clusters this big give their 2 PCA-axis extremes + min-z point
MERGE_CM = 0.05
FEATURE_PREFIX = "f_"      # feature columns of the pool table (labels are 'label_*')
TRACK_PDG = {13, 211, 2212, 321}
SHOWER_PDG = {11, 22}
END_NEAR_CM = 3.0          # physics_features: nearest PR-graph segment end within this distance
DIRS = np.array([[a, b, c] for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1) if (a, b, c) != (0, 0, 0)],
                np.float64)
DIRS /= np.linalg.norm(DIRS, axis=1, keepdims=True)
AV_HALF = (200.0, 200.0)
AV_Z = (0.0, 500.0)


# ---------------------------------------------------------------------------- hypothesis pool
@dataclass
class Pool:
    xyz: np.ndarray          # (n,3) float64 merged hypotheses
    source: np.ndarray       # (n,) int index into SOURCES of the highest-priority merged source
    counts: np.ndarray       # (n, len(SOURCES)) raw hypotheses merged into each, per source
    cand_q: np.ndarray       # (n,) max PR-candidate q among merged (NaN if none)
    cand_index: np.ndarray   # (n,) first merged PR candidate index (-1 if none)
    wcp_track: np.ndarray    # (n,) merged Wire-Cell track-like particle starts (mu, pi, p, K)
    wcp_shower: np.ndarray   # (n,) merged Wire-Cell e/gamma starts
    wcp_emax: np.ndarray     # (n,) max Wire-Cell particle energy among merged (0 if none)
    raw_xyz: np.ndarray      # (m,3) all raw hypotheses before merging

    def __len__(self):
        return len(self.xyz)


def _fit_segment_ends(ev):
    """(2S,3) both ends of every PR-graph segment of the track_fit layer (fitgraph parse)."""
    if ev.n_fit == 0:
        return np.zeros((0, 3))
    from .fitgraph import build_fit_graph
    g = build_fit_graph(ev.fit_xyz, ev.fit_q, ev.fit_cluster_id)
    return np.asarray(g.end_xyz, np.float64).reshape(-1, 3)


def cluster_end_points(ev, min_points=CLUSTER_END_MIN_POINTS):
    """(m,3) PCA-axis extremes and the most upstream (min z) point of every clustering_pr cluster with
    >= min_points points: hypotheses on clusters Wire-Cell PR did not process."""
    P = np.asarray(ev.points_xyz, np.float64)
    if len(P) == 0:
        return np.zeros((0, 3))
    cid = np.asarray(ev.points_cluster_id)
    order = np.argsort(cid, kind="stable")
    uc, start, cnt = np.unique(cid[order], return_index=True, return_counts=True)
    out = []
    for a, c in zip(start, cnt):
        if c < min_points:
            continue
        idx = order[a:a + c]
        pts = P[idx]
        mu = pts.mean(0)
        ax = np.linalg.eigh(np.cov((pts - mu).T))[1][:, -1]
        s = (pts - mu) @ ax
        out += [pts[np.argmin(s)], pts[np.argmax(s)], pts[np.argmin(pts[:, 2])]]
    return np.unique(np.array(out).reshape(-1, 3), axis=0)


def build_pool(ev, merge_cm=MERGE_CM, sources=DEFAULT_SOURCES):
    """Hypothesis pool of a (truth-free) Event from the given sources (subset of SOURCES)."""
    bad = set(sources) - set(SOURCES)
    if bad:
        raise ValueError(f"unknown pool sources {sorted(bad)}")
    xs, src, q, cidx, pdg, en = [], [], [], [], [], []

    def add(v, s, qq=np.nan, ci=-1, pd_=0, e=0.0):
        v = np.asarray(v, np.float64).reshape(3)
        if np.all(np.isfinite(v)):
            xs.append(v), src.append(s), q.append(qq), cidx.append(ci), pdg.append(pd_), en.append(e)

    use = set(sources)
    if "tag" in use and ev.has_tag_nu:
        add(ev.tag_nu, 0)
    if "cand" in use:
        for i, (c, cq) in enumerate(zip(ev.cand_xyz, ev.cand_q)):
            add(c, 1, float(cq), i)
    if "wcnu" in use:
        for w in ev.wc_reco:
            add(w.vtx, 2)
    if "wcp" in use:
        for w in ev.wc_reco:
            for p in w.particles:
                add(p.start, 3, pd_=p.pdg, e=p.energy)
    if "fitend" in use:
        for v in _fit_segment_ends(ev):
            add(v, 4)
    if "clusend" in use:
        for v in cluster_end_points(ev):
            add(v, 5)
    m = len(xs)
    if m == 0:
        z = np.zeros(0)
        return Pool(np.zeros((0, 3)), np.zeros(0, np.int64), np.zeros((0, len(SOURCES)), np.int64), z,
                    np.zeros(0, np.int64), np.zeros(0, np.int64), np.zeros(0, np.int64), z, np.zeros((0, 3)))
    R = np.stack(xs)
    K = np.empty_like(R)
    assign = np.empty(m, np.int64)
    nk = 0
    for i in range(m):      # greedy merge in priority order
        if nk:
            d2 = ((K[:nk] - R[i]) ** 2).sum(1)
            j = int(np.argmin(d2))
            if d2[j] < merge_cm ** 2:
                assign[i] = j
                continue
        K[nk] = R[i]
        assign[i] = nk
        nk += 1
    src, q, cidx, pdg, en = map(np.asarray, (src, q, cidx, pdg, en))
    counts = np.zeros((nk, len(SOURCES)), np.int64)
    np.add.at(counts, (assign, src), 1)
    first_src = np.full(nk, len(SOURCES), np.int64)
    np.minimum.at(first_src, assign, src)
    cq = np.full(nk, -np.inf)
    np.maximum.at(cq, assign, np.where(np.isfinite(q), q, -np.inf))
    ci = np.full(nk, np.iinfo(np.int64).max)
    np.minimum.at(ci, assign, np.where(cidx >= 0, cidx, np.iinfo(np.int64).max))
    apdg = np.abs(pdg)
    trk = np.zeros(nk, np.int64)
    np.add.at(trk, assign, (src == 3) & np.isin(apdg, list(TRACK_PDG)))
    shw = np.zeros(nk, np.int64)
    np.add.at(shw, assign, (src == 3) & np.isin(apdg, list(SHOWER_PDG)))
    emax = np.zeros(nk)
    np.maximum.at(emax, assign, np.where(src == 3, en, 0.0))
    return Pool(K[:nk].copy(), first_src, counts, np.where(np.isfinite(cq), cq, np.nan),
                np.where(ci == np.iinfo(np.int64).max, -1, ci), trk, shw, emax, R)


# ---------------------------------------------------------------------------- features
def _pairs(tree_h, tree_p, r):
    M = tree_h.sparse_distance_matrix(tree_p, r, output_type="ndarray")
    return M["i"].astype(np.int64), M["j"].astype(np.int64), M["v"]


def _bincount(i, w, n):
    return np.bincount(i, weights=w, minlength=n)[:n] if len(i) else np.zeros(n)


BASIC_NAMES = (
    [f"n_src_{s}" for s in SOURCES]
    + ["cand_q", "cand_q_rank", "wcp_track", "wcp_shower", "wcp_emax", "d_tag", "dz_tag", "same_clus_tag",
       "x", "y", "z", "d_wall_x", "d_wall_y", "d_wall_z", "n_hyp_1", "n_hyp_3", "n_hyp_10", "n_pool", "n_cand_ev",
       "log_n_points_ev", "n_fit_ev", "n_clusters_ev", "wc_ok",
       "d_nn", "npt_1", "npt_3", "npt_10", "q_1", "q_3", "q_10", "n_dirs_occ", "u_mean", "u_mean_z",
       "uq_mean", "uq_mean_z", "q_down_frac_10",
       "clus_logn", "clus_logq", "clus_rank_q", "clus_ext", "clus_frac", "clus_dend", "clus_dz_min", "clus_zfrac",
       "clus_ymax", "clus_wall", "clus_is_wc0", "clus_is_wc", "d_wcclus", "dz_wcclus_min",
       "d_fit_nn", "shw_n_5", "shw_frac_5"])


def basic_features(ev, pool):
    """Default truth-free hypothesis features (names BASIC_NAMES); see module docstring."""
    n = len(pool)
    X = {k: np.zeros(n) for k in BASIC_NAMES}
    if n == 0:
        return np.zeros((0, len(BASIC_NAMES)), np.float32), list(BASIC_NAMES)
    H = pool.xyz
    for k, s in enumerate(SOURCES):
        X[f"n_src_{s}"] = pool.counts[:, k].astype(np.float64)
    X["cand_q"] = np.where(np.isfinite(pool.cand_q), pool.cand_q, -1.0)
    if ev.n_cand:
        order = np.argsort(np.argsort(-np.asarray(ev.cand_q, np.float64), kind="stable"), kind="stable")
        X["cand_q_rank"] = np.where(pool.cand_index >= 0, order[np.maximum(pool.cand_index, 0)] / ev.n_cand, -1.0)
    else:
        X["cand_q_rank"][:] = -1.0
    X["wcp_track"], X["wcp_shower"], X["wcp_emax"] = pool.wcp_track, pool.wcp_shower, pool.wcp_emax
    X["x"], X["y"], X["z"] = H[:, 0], H[:, 1], H[:, 2]
    X["d_wall_x"] = AV_HALF[0] - np.abs(H[:, 0])
    X["d_wall_y"] = AV_HALF[1] - np.abs(H[:, 1])
    X["d_wall_z"] = np.minimum(H[:, 2] - AV_Z[0], AV_Z[1] - H[:, 2])
    raw_tree = cKDTree(pool.raw_xyz)
    own = pool.counts.sum(1)
    for r in (1, 3, 10):
        X[f"n_hyp_{r}"] = raw_tree.query_ball_point(H, r, return_length=True) - own
    X["n_pool"][:] = n
    X["n_cand_ev"][:] = ev.n_cand
    X["log_n_points_ev"][:] = np.log1p(ev.n_points)
    X["n_fit_ev"][:] = ev.n_fit
    X["wc_ok"][:] = float(ev.wc_status == "ok")
    tag = ev.tag_nu.astype(np.float64) if ev.has_tag_nu else None
    if tag is not None:
        X["d_tag"], X["dz_tag"] = np.sqrt(((H - tag) ** 2).sum(1)), H[:, 2] - tag[2]
    else:
        X["d_tag"][:] = -1.0

    # ---- clustering_pr context
    P = ev.points_xyz.astype(np.float64)
    if len(P):
        Q = np.maximum(np.asarray(ev.points_q, np.float64), 0)
        C = np.asarray(ev.points_cluster_id)
        uc, inv, cnt = np.unique(C, return_inverse=True, return_counts=True)
        X["n_clusters_ev"][:] = len(uc)
        cq = np.bincount(inv, weights=Q)
        crank = np.argsort(np.argsort(-cq, kind="stable"), kind="stable")
        ptree = cKDTree(P)
        d_nn, j_nn = ptree.query(H)
        X["d_nn"] = d_nn
        hc = inv[j_nn]
        htree = cKDTree(H)
        i, j, d = _pairs(htree, ptree, 10.0)
        for r in (1, 3, 10):
            mk = d <= r
            X[f"npt_{r}"] = _bincount(i[mk], None, n)
            X[f"q_{r}"] = _bincount(i[mk], Q[j[mk]], n) / 1e4
        mk = (d > 2) & (d <= 8)
        if mk.any():
            u = (P[j[mk]] - H[i[mk]]) / d[mk, None]
            b = np.argmax(u @ DIRS.T, 1)
            occ = np.zeros((n, len(DIRS)))
            np.add.at(occ, (i[mk], b), 1)
            X["n_dirs_occ"] = (occ > 3).sum(1)
            c8 = np.maximum(_bincount(i[mk], None, n), 1)
            mu = np.stack([_bincount(i[mk], u[:, k], n) for k in range(3)], 1) / c8[:, None]
            X["u_mean"], X["u_mean_z"] = np.linalg.norm(mu, axis=1), mu[:, 2]
        mk = (d > 2) & (d <= 10)
        if mk.any():
            w = Q[j[mk]]
            u = (P[j[mk]] - H[i[mk]]) / d[mk, None]
            wq = np.maximum(_bincount(i[mk], w, n), 1e-9)
            mq = np.stack([_bincount(i[mk], w * u[:, k], n) for k in range(3)], 1) / wq[:, None]
            X["uq_mean"], X["uq_mean_z"] = np.linalg.norm(mq, axis=1), mq[:, 2]
        mk = d <= 10
        down = P[j[mk], 2] > H[i[mk], 2]
        X["q_down_frac_10"] = _bincount(i[mk], Q[j[mk]] * down, n) / np.maximum(_bincount(i[mk], Q[j[mk]], n), 1e-9)
        wc_ids = [w.cluster_id for w in ev.wc_reco]
        clus = {}
        for c in np.unique(hc):
            pts = P[inv == c]
            mu = pts.mean(0)
            ax = np.linalg.svd(pts - mu, full_matrices=False)[2][0] if len(pts) > 2 else np.array([0, 0, 1.0])
            s = (pts - mu) @ ax
            wall = np.minimum(AV_HALF[1] - np.abs(pts[:, 1]), np.minimum(pts[:, 2] - AV_Z[0], AV_Z[1] - pts[:, 2]))
            clus[c] = (mu, ax, s.min(), s.max(), pts[:, 2].min(), pts[:, 2].max(), pts[:, 1].max(), wall.min())
        for k in range(n):
            c = hc[k]
            mu, ax, smin, smax, zmin, zmax, ymax, wall = clus[c]
            sh = (H[k] - mu) @ ax
            ext = smax - smin
            X["clus_logn"][k] = np.log1p(cnt[c])
            X["clus_logq"][k] = np.log1p(cq[c] / 1e4)
            X["clus_rank_q"][k] = crank[c]
            X["clus_ext"][k] = ext
            X["clus_frac"][k] = (sh - smin) / max(ext, 1e-3)
            X["clus_dend"][k] = min(abs(sh - smin), abs(smax - sh))
            X["clus_dz_min"][k] = H[k, 2] - zmin
            X["clus_zfrac"][k] = (H[k, 2] - zmin) / max(zmax - zmin, 1e-3)
            X["clus_ymax"][k] = ymax
            X["clus_wall"][k] = wall
            X["clus_is_wc0"][k] = float(bool(wc_ids) and uc[c] == wc_ids[0])
            X["clus_is_wc"][k] = float(uc[c] in wc_ids)
        if tag is not None:
            X["same_clus_tag"] = (hc == inv[ptree.query(tag)[1]]).astype(np.float64)
        if wc_ids and np.any(C == wc_ids[0]):
            wpts = P[C == wc_ids[0]]
            X["d_wcclus"] = cKDTree(wpts).query(H)[0]
            X["dz_wcclus_min"] = H[:, 2] - wpts[:, 2].min()
        else:
            X["d_wcclus"][:] = -1.0
    else:
        X["d_nn"][:] = -1.0
        X["d_wcclus"][:] = -1.0

    # ---- track_fit proximity (Bragg / orientation cues come from physics_features)
    X["d_fit_nn"][:] = cKDTree(ev.fit_xyz.astype(np.float64)).query(H)[0] if ev.n_fit else -1.0

    # ---- shower_track tags (q = 15000 shower-tagged, 0 track-tagged)
    if ev.n_shower:
        stree = cKDTree(ev.shower_xyz.astype(np.float64))
        i, j, d = _pairs(cKDTree(H), stree, 5.0)
        X["shw_n_5"] = _bincount(i, None, n)
        X["shw_frac_5"] = _bincount(i, (np.asarray(ev.shower_q)[j] > 0).astype(float), n) / np.maximum(X["shw_n_5"], 1)
    return np.stack([np.broadcast_to(X[k], (n,)) for k in BASIC_NAMES], 1).astype(np.float32), list(BASIC_NAMES)


def physics_features(ev, pool):
    """sbndvtx.features.tokens.describe_positions at each hypothesis ('ph_*') plus the END_FEATURES,
    P(start), length and inward direction of the nearest PR-graph segment end within END_NEAR_CM
    ('pe_*'; NaN if none)."""
    from . import fitgraph, tokens
    an = tokens.analyze(ev)
    n = len(pool)
    D = tokens.describe_positions(an, pool.xyz)
    et = tokens.end_tokens(an)
    pe_names = ["d", "p_start", "log_len", "straight", "dir_z", "log_med_dqdx"] + list(fitgraph.END_FEATURES)
    E = np.full((n, len(pe_names)), np.nan)
    if n and len(et["xyz"]):
        d, j = cKDTree(et["xyz"]).query(pool.xyz)
        ok = d <= END_NEAR_CM
        j = j[ok]
        E[ok] = np.column_stack([d[ok], et["p_start"][j], np.log1p(et["length"][j]), et["straight"][j],
                                 et["dir"][j, 2], et["log_med_dqdx"][j], et["feat"][j]])
    names = [f"ph_{c}" for c in tokens.CAND_FEATURES] + [f"pe_{c}" for c in pe_names]
    return np.concatenate([D, E], 1).astype(np.float32), names


class CacheFeatureAggregates:
    """Cached sbndvtx.features around each hypothesis: 'point'/'fit'/'shower' level -> value at the
    nearest element and mean within `radius` cm; 'cand' level -> value of the merged PR candidate;
    'event' level -> broadcast.  NaN where undefined.  `spec` is the cache's feature spec
    (build_pool_table fills it from the dataset); None falls back to sbndvtx.features.FEATURE_SPEC."""

    __name__ = "cache_feature_aggregates"

    def __init__(self, spec=None, radius=3.0):
        self.spec = None if spec is None else [tuple(s[:3]) for s in spec]
        self.radius = radius

    def with_spec(self, spec):
        return CacheFeatureAggregates(spec, self.radius)

    def __call__(self, ev, pool):
        spec = self.spec
        if spec is None:
            from . import features as fmod
            spec = [tuple(s[:3]) for s in fmod.FEATURE_SPEC]
        n = len(pool)
        geo = {"point": ev.points_xyz, "fit": ev.fit_xyz, "shower": ev.shower_xyz}
        cols, names, trees = [], [], {}
        for name, dim, level in spec:
            if name not in ev.features:
                continue
            a = np.asarray(ev.features[name], np.float64)
            if level == "event":
                a = a.reshape(-1)
                for d in range(dim):
                    cols.append(np.full(n, a[d]))
                    names.append(f"cf_event_{name}_{d}")
                continue
            a = a.reshape(-1, dim)
            if level == "cand":
                ok = pool.cand_index >= 0
                for d in range(dim):
                    v = np.full(n, np.nan)
                    v[ok] = a[pool.cand_index[ok], d]
                    cols.append(v)
                    names.append(f"cf_cand_{name}_{d}")
                continue
            xyz = geo[level]
            if len(xyz) == 0 or n == 0:
                for d in range(dim):
                    cols += [np.full(n, np.nan), np.full(n, np.nan)]
                    names += [f"cf_{level}_{name}_{d}_nn", f"cf_{level}_{name}_{d}_mean"]
                continue
            if level not in trees:
                trees[level] = cKDTree(np.asarray(xyz, np.float64))
            t = trees[level]
            _, jn = t.query(pool.xyz)
            i, j, _ = _pairs(cKDTree(pool.xyz), t, self.radius)
            c = _bincount(i, None, n)
            for d in range(dim):
                s = _bincount(i, a[j, d], n)
                cols += [a[jn, d], np.where(c > 0, s / np.maximum(c, 1), np.nan)]
                names += [f"cf_{level}_{name}_{d}_nn", f"cf_{level}_{name}_{d}_mean"]
        if not cols:
            return np.zeros((n, 0), np.float32), []
        return np.stack(cols, 1).astype(np.float32), names


cache_feature_aggregates = CacheFeatureAggregates()
DEFAULT_FEATURE_FNS = (basic_features, physics_features, cache_feature_aggregates)


def pool_features(ev, pool, feature_fns=DEFAULT_FEATURE_FNS):
    """Concatenate feature functions -> (X (n, F) float32, names)."""
    Xs, names = [], []
    for fn in feature_fns:
        X, nm = fn(ev, pool)
        X = np.asarray(X, np.float32)
        if X.shape != (len(pool), len(nm)):
            raise ValueError(f"{getattr(fn, '__name__', fn)}: shape {X.shape} vs ({len(pool)}, {len(nm)} names)")
        Xs.append(X)
        names += list(nm)
    if len(set(names)) != len(names):
        raise ValueError("duplicate feature names across feature functions")
    return (np.concatenate(Xs, 1) if Xs else np.zeros((len(pool), 0), np.float32)), names
