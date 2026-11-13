"""Expanded hypothesis pool: vendored from sbndvtx/eval/pool_ext.py (inference subset: EXT_SOURCES, ExtPool,
extension_points, expand_pool; bodies unchanged).

expand_pool(pool, ext_xyz, ext_src, merge_cm=MERGE_CM) -> ExtPool
    Base Pool (baselines.build_pool) plus extra positions (EXT_SOURCES: 'fitpt' = track_fit points, 'img' =
    clustering_pr points).  An extra position within merge_cm of a kept hypothesis is merged into it (counted in
    ext_counts); the others become new hypotheses (Pool.source = len(SOURCES)).
extension_points(ev, anchors, radius=3.0, include_fit=True, include_dense=True) -> (xyz, src)
"""
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from .baselines import MERGE_CM, SOURCES, Pool

EXT_SOURCES = ("fitpt", "img")
EXT_NAMES = ("x_is_base", "x_n_fitpt", "x_n_img")
RADIUS_CM = 3.0


@dataclass
class ExtPool:
    pool: Pool                 # merged hypotheses (base first, then new ones)
    ext_counts: np.ndarray     # (n, 2) extra positions merged into each hypothesis, per EXT_SOURCES
    n_base: int                # hypotheses [0, n_base) are the base pool, unchanged and in order
    assign: np.ndarray         # (m,) hypothesis index of every extra position given to expand_pool

    def __len__(self):
        return len(self.pool)

    def label(self):
        """(n,) source names: base SOURCES name, else 'fitpt' / 'img' (the first extra source merged)."""
        names = np.asarray(SOURCES + ("none",), dtype=object)[self.pool.source].copy()
        new = np.arange(len(self.pool)) >= self.n_base
        names[new] = np.where(self.ext_counts[new, 0] > 0, EXT_SOURCES[0], EXT_SOURCES[1])
        return names


def extension_points(ev, anchors, radius=RADIUS_CM, include_fit=True, include_dense=True):
    """(xyz (m,3) float64, src (m,) int index into EXT_SOURCES): all track_fit points (if include_fit), then all
    clustering_pr points within `radius` cm of any anchor (if include_dense)."""
    xs, ss = [], []
    if include_fit and ev.n_fit:
        f = np.asarray(ev.fit_xyz, np.float64).reshape(-1, 3)
        xs.append(f)
        ss.append(np.zeros(len(f), np.int64))
    A = np.asarray(anchors, np.float64).reshape(-1, 3)
    A = A[np.all(np.isfinite(A), 1)]
    if include_dense and ev.n_points and len(A):
        P = np.asarray(ev.points_xyz, np.float64).reshape(-1, 3)
        idx = cKDTree(P).query_ball_point(A, radius)
        j = np.unique(np.concatenate([np.asarray(i, np.int64) for i in idx]))
        xs.append(P[j])
        ss.append(np.ones(len(j), np.int64))
    if not xs:
        return np.zeros((0, 3)), np.zeros(0, np.int64)
    return np.concatenate(xs), np.concatenate(ss)


def expand_pool(pool, ext_xyz, ext_src, merge_cm=MERGE_CM):
    """Merge extra positions into a base Pool (see module docstring).  Greedy in input order: an extra position
    within merge_cm of the nearest base hypothesis joins it; otherwise it joins the first earlier kept extra
    position within merge_cm, or becomes a new hypothesis."""
    X = np.asarray(ext_xyz, np.float64).reshape(-1, 3)
    S = np.asarray(ext_src, np.int64).reshape(-1)
    nb = len(pool)
    m = len(X)
    assign = np.full(m, -1, np.int64)
    if m and nb:
        d, j = cKDTree(pool.xyz).query(X)
        hit = d < merge_cm
        assign[hit] = j[hit]
    rest = np.flatnonzero(assign < 0)
    new_xyz = []
    if len(rest):
        R = X[rest]
        nbrs = cKDTree(R).query_ball_point(R, merge_cm * (1 - 1e-9))
        local = np.full(len(rest), -1, np.int64)
        for a in range(len(rest)):
            if local[a] >= 0:
                continue
            local[a] = nb + len(new_xyz)
            new_xyz.append(R[a])
            for b in nbrs[a]:
                if b > a and local[b] < 0:
                    local[b] = local[a]
        assign[rest] = local
    k = len(new_xyz)
    n = nb + k
    counts = np.zeros((n, len(EXT_SOURCES)), np.int64)
    if m:
        np.add.at(counts, (assign, S), 1)
    p = pool
    xyz = np.concatenate([p.xyz.reshape(-1, 3), np.asarray(new_xyz, np.float64).reshape(-1, 3)])
    ext = Pool(xyz, np.concatenate([p.source, np.full(k, len(SOURCES), np.int64)]),
               np.concatenate([p.counts.reshape(-1, len(SOURCES)), np.zeros((k, len(SOURCES)), np.int64)]),
               np.concatenate([p.cand_q, np.full(k, np.nan)]), np.concatenate([p.cand_index, np.full(k, -1, np.int64)]),
               np.concatenate([p.wcp_track, np.zeros(k, np.int64)]), np.concatenate([p.wcp_shower, np.zeros(k, np.int64)]),
               np.concatenate([p.wcp_emax, np.zeros(k)]), p.raw_xyz)
    return ExtPool(ext, counts, nb, assign)
