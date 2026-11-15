"""Pass-2 crop framing (Pandora-DL style two-pass framing, arXiv:2502.06637).

frame(mode, pts_xyz, pts_q, pts_cluster, p_host, crop_cm, upstream_frac, radius) -> (lo, hi, center, axis)
  The predicted host mass is the imaging charge within `radius` of the pass-1 mode, each point weighted by its
  charge times the pass-1 P(host) of its cluster (1 if no cluster head output).  Its charge-weighted principal
  axis, oriented towards the mass (downstream of the vertex), places an axis-aligned cube of side crop_cm so
  that the mode sits `upstream_frac` of the side from the upstream face along that axis:
      center = mode + axis * crop_cm * (0.5 - upstream_frac).
  Fewer than 3 weighted points: the cube is centred on the mode.  Truth-free.
"""
import numpy as np
from scipy.spatial import cKDTree


def host_weights(pts_cluster, p_host):
    if not p_host:
        return np.ones(len(pts_cluster))
    ids = np.fromiter(p_host.keys(), np.int64)
    vals = np.fromiter(p_host.values(), np.float64)
    order = np.argsort(ids)
    ids, vals = ids[order], vals[order]
    c = np.asarray(pts_cluster, np.int64)
    pos = np.clip(np.searchsorted(ids, c), 0, max(len(ids) - 1, 0))
    hit = ids[pos] == c if len(ids) else np.zeros(len(c), bool)
    return np.where(hit, vals[pos], 0.0)


def frame(mode, pts_xyz, pts_q, pts_cluster, p_host, crop_cm, upstream_frac, radius, tree=None):
    mode = np.asarray(mode, np.float64).reshape(3)
    X = np.asarray(pts_xyz, np.float64).reshape(-1, 3)
    axis = np.zeros(3)
    center = mode.copy()
    if len(X):
        tree = tree if tree is not None else cKDTree(X)
        idx = np.asarray(tree.query_ball_point(mode, radius), np.int64)
        if len(idx) >= 3:
            w = np.maximum(np.asarray(pts_q, np.float64)[idx], 0) * host_weights(np.asarray(pts_cluster)[idx],
                                                                                   p_host)
            if w.sum() > 0 and (w > 0).sum() >= 3:
                P = X[idx]
                mu = (w[:, None] * P).sum(0) / w.sum()
                R = P - mu
                Cov = (w[:, None, None] * R[:, :, None] * R[:, None, :]).sum(0) / w.sum()
                vals, vecs = np.linalg.eigh(Cov)
                a = vecs[:, -1]
                s = (w * ((P - mode) @ a)).sum()
                if s < 0:
                    a = -a
                if not np.all(np.isfinite(a)):
                    a = np.zeros(3)
                axis = a
                center = mode + axis * crop_cm * (0.5 - upstream_frac)
    half = crop_cm / 2.0
    return center - half, center + half, center, axis


def inside(xyz, lo, hi):
    X = np.asarray(xyz, np.float64).reshape(-1, 3)
    return np.all((X >= lo) & (X < hi), 1)
