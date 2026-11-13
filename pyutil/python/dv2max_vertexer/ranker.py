"""Final stage at inference: the ranker's candidates and features for ONE variant (the deployed family), vendored
from dv2max/ranker.py (CNN_FEATURES, cnn_positions, cnn_anchors, cnn_features, c2_fields, VARIANT_FEATURES,
variant_pool, event_union; bodies unchanged except that event_union returns a dict of arrays instead of a pandas
DataFrame and computes no labels).  The LambdaMART model and the isotonic calibrator are applied in vertexer.py.

Candidates of one variant, merged exactly as the offline union (pool_ext.expand_pool, 0.05 cm):
  base pool   baselines.build_pool over ranker.pool_sources (r7_cloud: PR-graph vertex rows 'cand', the traditional
              vertex 'tag', cluster ends 'clusend')
  imaging     cloud points within ranker.dense_radius of the variant's top `anchors` CNN positions
  CNN         the pass-1 vote modes and peaks, pass-2 modes, and their nearest imaging points
Features: the shared pool features (baselines.DEFAULT_FEATURE_FNS minus n_pool) + VARIANT_FEATURES.
"""
from contextlib import nullcontext

import numpy as np
from scipy.spatial import cKDTree

try:                                   # optional: small per-event linear algebra is slower with many BLAS threads
    from threadpoolctl import threadpool_limits
except ImportError:                    # pragma: no cover
    threadpool_limits = None

from . import baselines, decode, pool_ext


def _threads1():
    return threadpool_limits(1) if threadpool_limits is not None else nullcontext()


CNN_FEATURES = (
    "p1_d_mode", "p1_mode_score", "p1_mode_rank", "p1_mode_mass1", "p1_mode_spread", "p1_mode_pri", "p1_top_score",
    "p1_d_top", "p1_dens", "p1_pdens", "p1_heat", "p1_any", "p1_pt", "p1_host", "p1_host_rank", "p1_event",
    "c3_score", "c2_shell", "c2_cone",
    "p2_d_mode", "p2_mode_score", "p2_mode_pcrop", "p2_mode_sp", "p2_mode_mass1", "p2_mode_spread", "p2_mode_rank",
    "p2_dens", "p2_pdens", "p2_heat", "p2_any", "p2_pt", "p2_in_crops", "p2_pcrop_max", "p2_best_sp", "p2_d_best",
    "p2_c3_score", "p2_c2_shell", "p2_c2_cone", "cnn_is_mode")
FAR = 100.0
CNN_SRC = 1          # CNN positions enter expand_pool as extra positions of source index 1 ('img' slot)


# ---------------------------------------------------------------------------- CNN candidates / features
def cnn_positions(rec1, rec2):
    """(xyz (m,3), score (m,)) of all CNN hypotheses of a variant (modes + snapped points of both passes)."""
    xs, ss = [], []
    for rec, sc_fn in ((rec1, lambda M: M[:, 3]), (rec2, None)):
        if rec is None or not len(rec["modes"]):
            continue
        M = rec["modes"]
        if sc_fn is None:
            s = M[:, 3] * rec["crops"][M[:, -1].astype(int), 9]
        else:
            s = sc_fn(M)
        xs.append(M[:, :3])
        ss.append(s)
        ok = np.isfinite(M[:, 10])
        xs.append(M[ok, 10:13])
        ss.append(s[ok])
    if not xs:
        return np.zeros((0, 3)), np.zeros(0)
    return np.concatenate(xs), np.concatenate(ss)


def cnn_anchors(rec1, rec2, n):
    """Top-n CNN positions: best pass-2 mode of each crop first (by score x P(crop)), then pass-1 vote modes;
    positions within 1 cm of an earlier anchor are skipped."""
    cand = []
    if rec2 is not None and len(rec2["modes"]):
        M = rec2["modes"]
        sp = M[:, 3] * rec2["crops"][M[:, -1].astype(int), 9]
        for c in np.unique(M[:, -1].astype(int)):
            m = np.flatnonzero(M[:, -1].astype(int) == c)
            j = m[np.argmax(sp[m])]
            cand.append((sp[j] + 1.0, M[j, :3]))
    if rec1 is not None and len(rec1["modes"]):
        for r in rec1["modes"]:
            cand.append((float(r[3]), r[:3]))
    cand.sort(key=lambda t: -t[0])
    out = []
    for _, x in cand:
        if out and np.min(np.linalg.norm(np.asarray(out) - x, axis=1)) < 1.0:
            continue
        out.append(np.asarray(x, np.float64))
        if len(out) >= n:
            break
    return np.asarray(out).reshape(-1, 3)


def _near_max(tree_xyz, vals, Q, r):
    out = np.zeros(len(Q))
    if tree_xyz is None or len(vals) == 0:
        return out
    for i, nb in enumerate(tree_xyz.query_ball_point(Q, r)):
        if nb:
            out[i] = float(np.max(vals[nb]))
    return out


def _density(cells, mass, Q, sigma):
    if len(cells) == 0:
        return np.zeros(len(Q))
    t = cKDTree(np.asarray(cells, np.float64))
    out = np.zeros(len(Q))
    tot = max(float(np.sum(mass)), 1e-12)
    for i, nb in enumerate(t.query_ball_point(Q, 3 * sigma)):
        if nb:
            d2 = ((np.asarray(cells, np.float64)[nb] - Q[i]) ** 2).sum(1)
            out[i] = float((np.asarray(mass, np.float64)[nb] * np.exp(-0.5 * d2 / sigma ** 2)).sum()) / tot
    return out


def _nearest_mode(M, Q, score, order_score):
    """distance, index of the nearest mode row and its rank by order_score."""
    if len(M) == 0:
        return np.full(len(Q), FAR), np.full(len(Q), -1), np.full(len(Q), np.nan)
    d, j = cKDTree(M[:, :3]).query(Q)
    rank = np.empty(len(M))
    rank[np.argsort(-order_score, kind="stable")] = np.arange(len(M))
    return np.minimum(d, FAR), j, rank[j]


def c2_fields(rec, Q, edges, radius=16.0):
    """C2 shell / cone field scores at positions Q from the hot voxels' shell (cumulative ordinal) and cone
    (von Mises-Fisher) predictions; NaN if the record has no C2 outputs."""
    if rec is None or "hot_shell" not in rec or len(rec["hot_xyz"]) == 0:
        return np.full(len(Q), np.nan), np.full(len(Q), np.nan)
    V = np.asarray(rec["hot_xyz"], np.float64)
    S = np.asarray(rec["hot_shell"], np.float64)
    Cn = np.asarray(rec["hot_cone"], np.float64)
    w = np.maximum(np.asarray(rec["hot_any"], np.float64), 1e-6)
    e = np.asarray(edges, np.float64)
    t = cKDTree(V)
    sh, co = np.full(len(Q), np.nan), np.full(len(Q), np.nan)
    for i, nb in enumerate(t.query_ball_point(Q, radius)):
        if not nb:
            continue
        nb = np.asarray(nb)
        vec = Q[i] - V[nb]
        d = np.linalg.norm(vec, axis=1)
        j = np.searchsorted(e, d, side="right")                 # bin index 0..K
        cum = np.concatenate([np.zeros((len(nb), 1)), S[nb], np.ones((len(nb), 1))], 1)
        p = cum[np.arange(len(nb)), j + 1] - cum[np.arange(len(nb)), j]
        sh[i] = float(np.average(np.log(np.clip(p, 1e-4, 1.0)), weights=w[nb]))
        ok = d > 0.3
        if ok.any():
            u = vec[ok] / d[ok, None]
            mu = Cn[nb][ok, :3]
            mu = mu / np.maximum(np.linalg.norm(mu, axis=1, keepdims=True), 1e-6)
            k = np.exp(np.clip(Cn[nb][ok, 3], -4, 8))
            logsinh = k + np.log1p(-np.exp(-2 * k).clip(max=1 - 1e-7)) - np.log(2.0)
            ll = k * (mu * u).sum(1) + np.log(k) - np.log(4 * np.pi) - logsinh
            co[i] = float(np.average(ll, weights=w[nb][ok]))
    return sh, co


def cnn_features(Q, clus, rec1, rec2, edges, edges2=None):
    """(n, len(CNN_FEATURES)) float32 CNN features at positions Q (clus: cluster id of each position's nearest
    imaging point, -1 if none).  edges / edges2: C2 shell edges of pass 1 / pass 2."""
    n = len(Q)
    F = np.full((n, len(CNN_FEATURES)), np.nan, np.float32)
    col = {c: i for i, c in enumerate(CNN_FEATURES)}
    if n == 0:
        return F
    if rec1 is not None:
        M = rec1["modes"]
        vm = M[M[:, 9] == decode.SRC_VOTE] if len(M) else M
        if len(vm):
            d, j, rk = _nearest_mode(vm, Q, vm[:, 3], vm[:, 3])
            F[:, col["p1_d_mode"]] = d
            F[:, col["p1_mode_score"]] = vm[j, 3]
            F[:, col["p1_mode_rank"]] = rk
            F[:, col["p1_mode_mass1"]] = vm[j, 4]
            F[:, col["p1_mode_spread"]] = vm[j, 6]
            F[:, col["p1_mode_pri"]] = vm[j, 7]
            F[:, col["p1_top_score"]] = vm[0, 3]
            F[:, col["p1_d_top"]] = np.minimum(np.linalg.norm(Q - vm[0, :3], axis=1), FAR)
        F[:, col["p1_dens"]] = _density(rec1["cells"], rec1["cell_mass"], Q, 1.0)
        F[:, col["p1_pdens"]] = _density(rec1["cells"], rec1["cell_pmass"], Q, 1.0)
        th = cKDTree(np.asarray(rec1["hot_xyz"], np.float64)) if len(rec1["hot_xyz"]) else None
        F[:, col["p1_heat"]] = _near_max(th, np.asarray(rec1["hot_heat"]), Q, 1.5)
        F[:, col["p1_any"]] = _near_max(th, np.asarray(rec1["hot_any"]), Q, 1.5)
        if "hot_c3" in rec1:
            F[:, col["c3_score"]] = _near_max(th, np.exp(np.asarray(rec1["hot_c3"], np.float64)), Q, 1.5)
        if "hot_pt_xyz" in rec1 and len(rec1["hot_pt_xyz"]):
            F[:, col["p1_pt"]] = _near_max(cKDTree(np.asarray(rec1["hot_pt_xyz"], np.float64)),
                                           np.asarray(rec1["hot_pt_p"]), Q, 1.0)
        ph = rec1.get("p_host", {})
        if ph:
            ids = np.array(list(ph.keys()), np.int64)
            vals = np.array(list(ph.values()), np.float64)
            rank = np.empty(len(vals))
            rank[np.argsort(-vals, kind="stable")] = np.arange(len(vals))
            lut = dict(zip(ids.tolist(), zip(vals.tolist(), rank.tolist())))
            hv = np.array([lut.get(int(c), (np.nan, np.nan)) for c in clus], np.float64).reshape(-1, 2)
            F[:, col["p1_host"]], F[:, col["p1_host_rank"]] = hv[:, 0], hv[:, 1]
        F[:, col["p1_event"]] = rec1.get("event_p", np.nan)
        F[:, col["c2_shell"]], F[:, col["c2_cone"]] = c2_fields(rec1, Q, edges)
    if rec2 is not None and len(rec2["crops"]):
        M = rec2["modes"]
        pc = rec2["crops"][:, 9]
        if len(M):
            sp = M[:, 3] * pc[M[:, -1].astype(int)]
            d, j, rk = _nearest_mode(M, Q, M[:, 3], sp)
            F[:, col["p2_d_mode"]] = d
            F[:, col["p2_mode_score"]] = M[j, 3]
            F[:, col["p2_mode_pcrop"]] = pc[M[j, -1].astype(int)]
            F[:, col["p2_mode_sp"]] = sp[j]
            F[:, col["p2_mode_mass1"]] = M[j, 4]
            F[:, col["p2_mode_spread"]] = M[j, 6]
            F[:, col["p2_mode_rank"]] = rk
            b = int(np.argmax(sp))
            F[:, col["p2_best_sp"]] = sp[b]
            F[:, col["p2_d_best"]] = np.minimum(np.linalg.norm(Q - M[b, :3], axis=1), FAR)
        F[:, col["p2_dens"]] = _density(rec2["cells"], rec2["cell_mass"], Q, 0.5)
        F[:, col["p2_pdens"]] = _density(rec2["cells"], rec2["cell_pmass"], Q, 0.5)
        th = cKDTree(np.asarray(rec2["hot_xyz"], np.float64)) if len(rec2["hot_xyz"]) else None
        F[:, col["p2_heat"]] = _near_max(th, np.asarray(rec2["hot_heat"]), Q, 0.5)
        F[:, col["p2_any"]] = _near_max(th, np.asarray(rec2["hot_any"]), Q, 0.5)
        if len(rec2["hot_pt_xyz"]):
            F[:, col["p2_pt"]] = _near_max(cKDTree(np.asarray(rec2["hot_pt_xyz"], np.float64)),
                                           np.asarray(rec2["hot_pt_p"]), Q, 0.5)
        if "hot_c3" in rec2:
            F[:, col["p2_c3_score"]] = _near_max(th, np.exp(np.asarray(rec2["hot_c3"], np.float64)), Q, 0.5)
        F[:, col["p2_c2_shell"]], F[:, col["p2_c2_cone"]] = c2_fields(rec2, Q, edges if edges2 is None else edges2)
        C = rec2["crops"]
        inside = np.all((Q[:, None, :] >= C[None, :, 3:6]) & (Q[:, None, :] < C[None, :, 6:9]), 2)
        F[:, col["p2_in_crops"]] = inside.sum(1)
        F[:, col["p2_pcrop_max"]] = np.where(inside.any(1), np.max(np.where(inside, pc[None], -1.0), 1), 0.0)
    P, _ = cnn_positions(rec1, rec2)
    if len(P):
        F[:, col["cnn_is_mode"]] = (cKDTree(P).query(Q)[0] < 0.05).astype(np.float32)
    else:
        F[:, col["cnn_is_mode"]] = 0.0
    return F


# ---------------------------------------------------------------------------- union rows of one event
VARIANT_FEATURES = ("x_is_base", "x_n_fitpt", "x_n_img", "n_pool") + CNN_FEATURES
SHARED_DROP = ("n_pool",)        # union-size dependent: replaced by the per-variant n_pool


def variant_pool(evin, pool, anchors, cnnpos, rcfg, fit_xyz):
    """Candidate set of ONE variant, merged exactly as a single-model run would (sbndvtx pool_ext.expand_pool):
    base pool + fit points + imaging points within dense_radius of the variant's anchors + its CNN positions.
    Returns (ExtPool, x_n_img corrected for CNN positions (n,), cnn row flags (n,))."""
    xs, ss = [], []
    if len(fit_xyz):
        xs.append(fit_xyz)
        ss.append(np.zeros(len(fit_xyz), np.int64))
    if len(anchors) and evin.n_points:
        img, _ = pool_ext.extension_points(evin, anchors, rcfg["dense_radius"], include_fit=False,
                                           include_dense=True)
        xs.append(img)
        ss.append(np.ones(len(img), np.int64))
    n_ext = sum(len(x) for x in xs)
    if len(cnnpos):
        xs.append(np.asarray(cnnpos, np.float64))
        ss.append(np.full(len(cnnpos), CNN_SRC, np.int64))
    allx = np.concatenate(xs) if xs else np.zeros((0, 3))
    alls = np.concatenate(ss) if ss else np.zeros(0, np.int64)
    ext = pool_ext.expand_pool(pool, allx, alls)
    n = len(ext)
    nimg = ext.ext_counts[:, 1].astype(np.float64).copy()
    cnn_rows = np.zeros(n, bool)
    if len(cnnpos):
        a = ext.assign[n_ext:]
        nimg -= np.bincount(a, minlength=n)[:n]
        cnn_rows[a] = True
    return ext, nimg, cnn_rows


def _edges_pair(edges):
    """(pass-1 shell edges, pass-2 shell edges) from either one list or a pair of lists."""
    if len(edges) == 2 and np.ndim(edges[0]) == 1:
        return list(edges[0]), list(edges[1])
    return list(edges), list(edges)


def event_union(ev, recs, variants, rcfg, edges, feature_fns):
    """Rows of one event for all its variants.  Each variant's candidates are merged on their own (no
    cross-variant dependence: a variant's rows, positions and features are exactly those of a run with that
    family alone); the pool features are computed once on the union of exact positions.
    Returns (ids DataFrame of union rows, shared feature matrix, shared names, rows_v list (V) of union row
    indices, block_v list (V) of (len(rows_v), len(VARIANT_FEATURES)) per-variant features)."""
    evin = ev.strip_truth()
    with _threads1():
        src = tuple(rcfg["pool_sources"]) if rcfg.get("pool_sources") else baselines.DEFAULT_SOURCES
        pool = baselines.build_pool(evin, sources=src)
        nb = len(pool)
        fit = np.asarray(evin.fit_xyz, np.float64).reshape(-1, 3) if rcfg["expand_fit"] else np.zeros((0, 3))
        per = []
        for fam in variants:
            r1, r2 = recs.get(fam, (None, None))
            per.append(variant_pool(evin, pool, cnn_anchors(r1, r2, rcfg["anchors"]), cnn_positions(r1, r2)[0],
                                    rcfg, fit))
        # union of exact positions: base rows first (identical in every variant), then the new rows
        new = [e.pool.xyz[e.n_base:] for e, _, _ in per]
        allnew = np.concatenate(new) if new else np.zeros((0, 3))
        uniq = np.unique(allnew, axis=0) if len(allnew) else np.zeros((0, 3))
        U = np.concatenate([pool.xyz.reshape(-1, 3), uniq])
        k = len(uniq)
        upool = baselines.Pool(U, np.concatenate([pool.source, np.full(k, len(baselines.SOURCES), np.int64)]),
                               np.concatenate([pool.counts.reshape(-1, len(baselines.SOURCES)),
                                               np.zeros((k, len(baselines.SOURCES)), np.int64)]),
                               np.concatenate([pool.cand_q, np.full(k, np.nan)]),
                               np.concatenate([pool.cand_index, np.full(k, -1, np.int64)]),
                               np.concatenate([pool.wcp_track, np.zeros(k, np.int64)]),
                               np.concatenate([pool.wcp_shower, np.zeros(k, np.int64)]),
                               np.concatenate([pool.wcp_emax, np.zeros(k)]), pool.raw_xyz)
        X, names = baselines.pool_features(evin, upool, feature_fns)
    keep = [i for i, n in enumerate(names) if n not in SHARED_DROP]
    X = np.asarray(X, np.float32)[:, keep]
    names = [names[i] for i in keep]
    n = len(U)
    P = np.asarray(evin.points_xyz, np.float64)
    clus = np.full(n, -1, np.int64)
    if len(P) and n:
        _, j = cKDTree(P).query(U)
        clus = np.asarray(evin.points_cluster_id, np.int64)[j]
    source = np.asarray(baselines.SOURCES + ("none",), dtype=object)[upool.source].copy()
    rows_v, block_v = [], []
    lut = {tuple(x): nb + i for i, x in enumerate(uniq)}
    for fam, (ext, nimg, cnn_rows) in zip(variants, per):
        rv = np.concatenate([np.arange(ext.n_base), np.array([lut[tuple(x)] for x in ext.pool.xyz[ext.n_base:]],
                                                              np.int64)]).astype(np.int64)
        for i in np.flatnonzero(np.arange(len(ext)) >= ext.n_base):
            if source[rv[i]] == "none":
                source[rv[i]] = "fitpt" if ext.ext_counts[i, 0] > 0 else ("cnn" if cnn_rows[i] else "img")
        r1, r2 = recs.get(fam, (None, None))
        e1, e2 = _edges_pair(edges)
        C = cnn_features(U[rv], clus[rv], r1, r2, e1, e2) if (r1 is not None or r2 is not None) else \
            np.full((len(rv), len(CNN_FEATURES)), np.nan, np.float32)
        blk = np.column_stack([np.arange(len(ext)) < ext.n_base, ext.ext_counts[:, 0], nimg,
                               np.full(len(ext), len(ext))]).astype(np.float32)
        rows_v.append(rv)
        block_v.append(np.concatenate([blk, C], 1).astype(np.float32))
    # dv2max_vertexer: a dict of arrays instead of the pandas DataFrame, no label columns (truth-free)
    df = {"x": U[:, 0], "y": U[:, 1], "z": U[:, 2], "source": source.astype(str), "is_base": np.arange(n) < nb,
          "cluster": clus}
    return df, X, names, rows_v, block_v


