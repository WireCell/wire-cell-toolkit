"""Point-level (clustering_pr) neighbourhood and cluster features, vectorized with cKDTree.

All functions take float64 arrays and return float64; callers cast.  Radii in cm.
"""
import numpy as np
from scipy.spatial import cKDTree

from .geometry import X_ANODE, Y_HALF, Z_HI, Z_LO, wall_distances

BALL_RADII = (0.5, 1.0, 3.0)
PCA_RADII = (1.0, 3.0)
VOXEL_CM = 1.0
FAR_RADIUS = 10.0
Q_UNIT = 1e4               # clustering_pr q: median ~6e3 per point
FACE_CM = 3.0              # a cluster "touches" a face if a point is this close
CHUNK = 4096               # query points per neighbourhood chunk (bounds memory)


def shower_tags(points_xyz, tree, shower_xyz, shower_q):
    """Per point: 0 = not in the PR shower_track layer, 1 = track-associated, 2 = shower-associated.
    shower_track points are exact copies of clustering_pr points (100% on the dev cache)."""
    tag = np.zeros(len(points_xyz), np.int8)
    if len(shower_xyz) == 0 or len(points_xyz) == 0:
        return tag
    d, j = tree.query(np.asarray(shower_xyz, np.float64), distance_upper_bound=1e-2)
    ok = np.isfinite(d)
    tq = np.asarray(shower_q)[ok]
    jj = j[ok]
    tag[jj[tq <= 0]] = 1
    tag[jj[tq > 0]] = 2          # shower wins over track
    return tag


def _pca_from_moments(n, m1, m2):
    """Batch PCA from counts n (k,), first moments m1 (k,3) and second moments m2 (k,6) of RELATIVE coords.
    Returns eigenvalues descending (k,3) and principal vector (k,3)."""
    k = len(n)
    nn = np.maximum(n, 1)[:, None]
    mu = m1 / nn
    c = m2 / nn
    C = np.empty((k, 3, 3))
    C[:, 0, 0] = c[:, 0] - mu[:, 0] ** 2
    C[:, 1, 1] = c[:, 1] - mu[:, 1] ** 2
    C[:, 2, 2] = c[:, 2] - mu[:, 2] ** 2
    C[:, 0, 1] = C[:, 1, 0] = c[:, 3] - mu[:, 0] * mu[:, 1]
    C[:, 0, 2] = C[:, 2, 0] = c[:, 4] - mu[:, 0] * mu[:, 2]
    C[:, 1, 2] = C[:, 2, 1] = c[:, 5] - mu[:, 1] * mu[:, 2]
    w, v = np.linalg.eigh(C)
    w = np.maximum(w[:, ::-1], 0.0)
    return w, v[:, :, 2]


def _moments(i, d, n_out, weights=None):
    """Bincount first/second moments of relative vectors d grouped by i."""
    w = np.ones(len(i)) if weights is None else weights
    m1 = np.stack([np.bincount(i, w * d[:, a], n_out) for a in range(3)], 1)
    prods = (d[:, 0] ** 2, d[:, 1] ** 2, d[:, 2] ** 2, d[:, 0] * d[:, 1], d[:, 0] * d[:, 2], d[:, 1] * d[:, 2])
    m2 = np.stack([np.bincount(i, w * p, n_out) for p in prods], 1)
    return m1, m2


def eigen_shape(w):
    """(k,3) descending eigenvalues -> linearity, planarity, sphericity (each in [0,1])."""
    l1 = np.maximum(w[:, 0], 1e-12)
    return (w[:, 0] - w[:, 1]) / l1, (w[:, 1] - w[:, 2]) / l1, w[:, 2] / l1


def sign_free(v):
    """(k,3) unit vectors -> (k,6) outer-product components xx, yy, zz, xy, xz, yz (invariant to v -> -v)."""
    return np.stack([v[:, 0] ** 2, v[:, 1] ** 2, v[:, 2] ** 2, v[:, 0] * v[:, 1], v[:, 0] * v[:, 2],
                     v[:, 1] * v[:, 2]], 1)


def point_neighbourhood(X, Q, tree, chunk=CHUNK):
    """Ball counts/charges at BALL_RADII, PCA at PCA_RADII, charge-weighted centroid shift at 3 cm.
    Query points are processed in chunks to bound memory.  Returns dict of arrays (n, .)."""
    n = len(X)
    out = {"count": np.zeros((n, len(BALL_RADII))), "charge": np.zeros((n, len(BALL_RADII))),
           "shape": np.zeros((n, 5)), "dir": np.zeros((n, 6)), "shift3": np.zeros((n, 3))}
    rmax = max(BALL_RADII + PCA_RADII)
    for c0 in range(0, n, chunk):
        sl = slice(c0, min(c0 + chunk, n))
        m_ = sl.stop - sl.start
        M = cKDTree(X[sl]).sparse_distance_matrix(tree, rmax, output_type="ndarray")
        i, j, dist = M["i"], M["j"], M["v"]
        del M
        for k, r in enumerate(BALL_RADII):
            m = dist <= r
            out["count"][sl, k] = np.bincount(i[m], minlength=m_)
            out["charge"][sl, k] = np.bincount(i[m], Q[j[m]], m_)
        d = X[j] - X[sl][i]
        shape = []
        for r in PCA_RADII:
            m = dist <= r
            cnt = np.bincount(i[m], minlength=m_)
            m1, m2 = _moments(i[m], d[m], m_)
            w, v = _pca_from_moments(cnt, m1, m2)
            lin, pla, _ = eigen_shape(w)
            few = cnt < 3
            shape += [np.where(few, 0.0, lin), np.where(few, 0.0, pla)]
            if r == PCA_RADII[-1]:
                shape.append(np.where(few, 0.0, np.sqrt(w.sum(1)) / r))
                out["dir"][sl] = sign_free(np.where(few[:, None], 0.0, v))
        out["shape"][sl] = np.stack(shape, 1)
        m = dist <= 3.0
        qw = Q[j[m]]
        qs = np.bincount(i[m], qw, m_)
        dm = d[m]
        sh = np.stack([np.bincount(i[m], qw * dm[:, a], m_) for a in range(3)], 1)
        out["shift3"][sl] = sh / np.maximum(qs, 1e-9)[:, None] / 3.0
    return out


def voxelize(X, Q, size=VOXEL_CM):
    """-> (voxel index per point, voxel charge-weighted centre (k,3), voxel charge (k,))."""
    key = np.floor(X / size).astype(np.int64)
    _, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    k = inv.max() + 1 if len(inv) else 0
    qv = np.bincount(inv, Q, k)
    qq = np.maximum(qv, 1e-9)
    cen = np.stack([np.bincount(inv, Q * X[:, a], k) for a in range(3)], 1) / qq[:, None]
    # voxels with zero total charge: use the plain mean
    zero = qv <= 0
    if np.any(zero):
        cnt = np.bincount(inv, minlength=k)
        plain = np.stack([np.bincount(inv, X[:, a], k) for a in range(3)], 1) / np.maximum(cnt, 1)[:, None]
        cen[zero] = plain[zero]
    return inv, cen, qv


def far_shift(X, Q, radius=FAR_RADIUS, chunk=CHUNK):
    """Charge-weighted centroid shift (/radius) and charge within `radius`, on a 1 cm voxel grid.
    One-sided charge (a track end, a vertex with forward prongs) gives |shift| ~ 0.5."""
    n = len(X)
    if n == 0:
        return np.zeros((0, 3)), np.zeros(0)
    inv, cen, qv = voxelize(X, Q)
    vt = cKDTree(cen)
    k = len(cen)
    qs = np.zeros(k)
    sh = np.zeros((k, 3))
    for c0 in range(0, k, chunk):
        sl = slice(c0, min(c0 + chunk, k))
        m_ = sl.stop - sl.start
        M = cKDTree(cen[sl]).sparse_distance_matrix(vt, radius, output_type="ndarray")
        i, j = M["i"], M["j"]
        qs[sl] = np.bincount(i, qv[j], m_)
        sh[sl] = np.stack([np.bincount(i, qv[j] * (cen[j, a] - cen[sl][i, a]), m_) for a in range(3)], 1)
    shift_v = sh / np.maximum(qs, 1e-9)[:, None] / radius
    # per point: move the reference from the voxel centre to the point itself
    shift = shift_v[inv] + (cen[inv] - X) / radius
    return shift, qs[inv]


def cluster_stats(X, Q, cid, wc_clusters=()):
    """Per-cluster statistics broadcast to points.  Returns (per-point dict, per-cluster dict, inverse index)."""
    n = len(X)
    if n == 0:
        z = np.zeros(0)
        return {k: z for k in ("size", "qtot", "ext1", "ext2", "ext3", "lin", "is_wc", "faces", "frac_z",
                               "rank")}, {}, np.zeros(0, np.int64)
    uc, inv = np.unique(cid, return_inverse=True)
    inv = inv.reshape(-1)
    k = len(uc)
    size = np.bincount(inv, minlength=k).astype(np.float64)
    qtot = np.bincount(inv, Q, k)
    cen = np.stack([np.bincount(inv, X[:, a], k) for a in range(3)], 1) / size[:, None]
    d = X - cen[inv]
    m1, m2 = _moments(inv, d, k)
    w, e1 = _pca_from_moments(size, m1, m2)
    lin, _, _ = eigen_shape(w)
    # orient e1 toward +z (downstream) so frac_z = 0 at the upstream end of the cluster
    e1 = np.where((e1[:, 2] < 0)[:, None], -e1, e1)
    proj = np.einsum("ij,ij->i", d, e1[inv])
    pmin = np.full(k, np.inf)
    pmax = np.full(k, -np.inf)
    np.minimum.at(pmin, inv, proj)
    np.maximum.at(pmax, inv, proj)
    ext1 = pmax - pmin
    frac_z = (proj - pmin[inv]) / np.maximum(ext1[inv], 1e-6)
    ext = 2.0 * np.sqrt(w)                     # ~ +-1 sigma extents along the 3 axes
    # TPC faces touched: -x anode, +x anode, -y, +y, upstream z, downstream z
    lo = np.full((k, 3), np.inf)
    hi = np.full((k, 3), -np.inf)
    for a in range(3):
        np.minimum.at(lo[:, a], inv, X[:, a])
        np.maximum.at(hi[:, a], inv, X[:, a])
    faces = ((lo[:, 0] < -X_ANODE + FACE_CM).astype(int) + (hi[:, 0] > X_ANODE - FACE_CM)
             + (lo[:, 1] < -Y_HALF + FACE_CM) + (hi[:, 1] > Y_HALF - FACE_CM)
             + (lo[:, 2] < Z_LO + FACE_CM) + (hi[:, 2] > Z_HI - FACE_CM)).astype(np.float64)
    is_wc = np.isin(uc, np.asarray(list(wc_clusters), dtype=uc.dtype)).astype(np.float64)
    rank = np.empty(k)
    rank[np.argsort(-qtot, kind="stable")] = np.arange(k)
    st = np.concatenate([[0], np.cumsum(size.astype(np.int64))[:-1]])
    o = np.lexsort((proj, inv))
    end_lo, end_hi = o[st], o[st + size.astype(np.int64) - 1]
    oz = np.lexsort((X[:, 2], inv))
    z_lo = oz[st]
    per_cluster = {"id": uc, "size": size, "qtot": qtot, "ext1": ext1, "ext": ext, "lin": lin, "faces": faces,
                   "is_wc": is_wc, "rank": rank, "e1": e1, "center": cen, "zmin": lo[:, 2], "lo": lo, "hi": hi,
                   "end_lo": end_lo, "end_hi": end_hi, "z_lo": z_lo, "pmin": pmin}
    per_point = {"size": size[inv], "qtot": qtot[inv], "ext1": ext1[inv], "ext2": ext[inv, 1], "ext3": ext[inv, 2],
                 "lin": lin[inv], "is_wc": is_wc[inv], "faces": faces[inv], "frac_z": frac_z, "rank": rank[inv]}
    return per_point, per_cluster, inv


def wall_features(X):
    """(n,5): log1p of anode / cathode / y-face / z-face distances, and the TPC side sign(x)."""
    d = wall_distances(X).astype(np.float64)
    return np.concatenate([np.log1p(d), np.sign(X[:, :1])], 1)
