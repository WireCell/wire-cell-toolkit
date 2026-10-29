"""Wire-Cell PR graph from the track_fit layer: segments, nodes, dQ/dx, residual range, Bragg cues, orientation.

Layer format (WCT MultiAlgBlobClustering.cxx fill_bee_points_from_pr_graph; checked on the dev cache):
  * every PR segment dumps its fitted points in trajectory order.  The FIRST and LAST point of a
    segment are the positions of its two PR vertices, shared exactly by all segments on that vertex;
  * then every PR vertex once (include_vertex_points) as a single point;
  * q = 0.1*dQ - 1000, clipped at 0 (dQ in electrons per fit point; 6.6% of points are clipped).
    Interior dx = half the distance to each neighbour (TrackFitting.cxx dQ_dx_multi_fit);
  * blocks repeat when there are several reco-neutrino bundles.
Parsing: split the ordered list at jumps > 1.5 cm, cluster changes, zero steps (shared vertex) and at
"twinned" points (a position that occurs twice in the layer = a PR vertex).  Runs of >= 3 points (or 2
points not both twinned) are segments; all other points are vertex points.  Nodes are the unique
vertex positions (segment ends + vertex points, merged within 0.05 cm).  Node degree = number of
segment ends on the node.  Nothing here uses truth.

Per segment end e (0 = first stored point, 1 = last) the END_FEATURES table holds Bragg, charge,
connectivity, isolation and geometry cues; `orient` turns them into P(end 0 is the start), i.e. the
end the particle came from (the vertex side), with the far end the stopping / exiting end.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .geometry import boundary_distance

JUMP_CM = 1.5            # larger steps separate segments
DUP_CM = 1e-3            # a zero step = two segments sharing a vertex point
TWIN_CM = 0.01           # a point with an exact copy elsewhere in the layer is a PR vertex position
NODE_CM = 0.05           # vertex positions closer than this are one node
DQ_SCALE, DQ_OFFSET = 0.1, -1000.0
DQ_CLIPPED = 5000.0      # dQ used for clipped points (true dQ in [0, 1e4) e)
MIP_DQDX = 47500.0       # e/cm: mode of interior fit dQ/dx on the dev cache (45-50k bin)
MIN_DX = 0.1             # cm, floor on dx
SMOOTH_HALF = 2          # running median over 2*2+1 fit points (~3 cm)
BRAGG_NEAR_CM = 3.0      # "end region" for mean/max dQ/dx
BRAGG_WIN_CM = 12.0      # window for the log-log slope and template LLR
BRAGG_RMIN_CM = 0.3
BRAGG_BETA = -0.42       # Bragg-Kleeman dE/dx ~ R^-0.42
SIGMA_LOG = 0.3          # per-point scatter of log dQ/dx (MIP tracks), for the LLR
DIR_CM = 3.0             # end direction: end point -> point this far along the segment
NO_TAG_CM = 1000.0       # distance used when tag_nu is missing

END_FEATURES = (
    "bragg_mean3",   # log(mean smoothed dQ/dx within 3 cm of the end / segment median)
    "bragg_max3",    # log(max smoothed dQ/dx within 3 cm of the end / segment median)
    "bragg_slope",   # d log(dQ/dx) / d log(R), R in [0.3, 12] cm from this end (< 0: rises toward the end)
    "bragg_llr",     # log-likelihood ratio Bragg template (R^-0.42) vs flat, clipped to [-20, 20]
    "end_dqdx",      # log(dQ/dx / MIP) of the vertex point at this end (shared by all segments there)
    "degree",        # other segment ends on this end's PR node
    "n_ends2",       # other segments with an end within 2 cm
    "n_cand2",       # PR vertex candidates within 2 cm
    "d_tag",         # log1p(distance to tag_nu / cm); log1p(1000) if no tag_nu
    "q_out3",        # log1p(clustering_pr charge within 3 cm BEYOND the end / 1e4)
    "f_out3",        # fraction of the 3 cm ball charge that lies beyond the end
    "shw_out5",      # log1p(shower-tagged charge within 5 cm beyond the end / 1e4)
    "d_wall",        # log1p(distance to nearest crossable TPC boundary / cm)
    "upstream",      # (z_other_end - z_this_end) / L: +1 if this end is the upstream (beam) end
    "higher",        # (y_this_end - y_other_end) / L: +1 if this end is on top
    "chain_len",     # log1p(total length of the OTHER segments on this end's node / cm)
)
EF = {n: k for k, n in enumerate(END_FEATURES)}

# Orientation: start-likeness score per end, P(end0 = start) = sigmoid(S0 - S1), S_e = sum_k w_k f_k(end e).
# Weights: logistic regression (no intercept, flip-augmented) on END_FEATURE differences of 2,643 fit
# segments (L > 5 cm) in the PR-graph component holding the true vertex, 547 dev-cache events, label =
# end nearer the vertex along the graph (scripts/validate_features.py --fit-orientation).  5-fold CV
# accuracy 0.939 (0.935 on segments attached to the vertex); out of sample on 5,760 other r3 events: 0.932
# (19,108 segments), 0.917 on segments with an end within 2 cm of the vertex (protons 0.979, muons 0.935,
# pions 0.921 by range match).  A sign-constrained refit is as accurate (0.932) but keeps only
# bragg_mean3 among the Bragg cues: alone, Bragg orients stopping protons / muons 85% / 83% of the time,
# showers and broken tracks 62%; connectivity and the beam direction carry the rest.  tag_nu / PR-candidate
# cues are left out
# on purpose so that orientation is independent of Wire-Cell's vertex choice.  These 14 numbers are the
# only MC-derived constants in the features.
ORIENT_WEIGHTS = {
    "bragg_mean3": -1.50, "bragg_max3": 0.80, "bragg_slope": -0.58, "bragg_llr": -0.01, "end_dqdx": 0.30,
    "degree": 0.10, "n_ends2": 0.22, "q_out3": 0.42, "f_out3": -0.73, "shw_out5": -0.37, "d_wall": 0.02,
    "upstream": 1.92, "higher": -0.04, "chain_len": 1.16,
}


def decode_dq(q):
    """Bee track_fit q -> (dQ electrons, clipped flag)."""
    q = np.asarray(q, np.float64)
    clipped = q <= 0
    return np.where(clipped, DQ_CLIPPED, (q - DQ_OFFSET) / DQ_SCALE), clipped


def _group_median(vals, grp, n):
    """Median of vals per group id in [0, n); NaN for empty groups."""
    out = np.full(n, np.nan)
    if len(vals) == 0:
        return out
    o = np.lexsort((vals, grp))
    g, v = grp[o], vals[o]
    cnt = np.bincount(g, minlength=n)
    st = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    ok = cnt > 0
    lo = st[ok] + (cnt[ok] - 1) // 2
    hi = st[ok] + cnt[ok] // 2
    out[ok] = 0.5 * (v[lo] + v[hi])
    return out


def _bsum(g, w, n):
    return np.bincount(g, weights=w, minlength=n)


def _unit(v):
    nrm = np.linalg.norm(v, axis=-1, keepdims=True)
    return np.where(nrm > 1e-9, v / np.maximum(nrm, 1e-9), 0.0)


@dataclass
class FitGraph:
    """Parsed track_fit layer.  Per-point arrays have length n_fit; per-segment S; per-node V."""
    xyz: np.ndarray                      # (n,3) float64
    seg: np.ndarray                      # (n,) segment index, -1 for vertex points
    node: np.ndarray                     # (n,) node index for vertex points and segment end points, else -1
    is_end: np.ndarray                   # (n,) bool: first/last point of a segment
    s: np.ndarray                        # (n,) arc length from the segment's first point (0 for vertex points)
    dq: np.ndarray                       # (n,) decoded dQ (e)
    clipped: np.ndarray                  # (n,) bool
    dx: np.ndarray                       # (n,) cm (segment points; node total for vertex points)
    dqdx: np.ndarray                     # (n,) e/cm, raw
    dqdx_s: np.ndarray                   # (n,) e/cm, running median along the segment (raw for vertex points)
    tangent: np.ndarray                  # (n,3) unit, along storage order (0 for vertex points)
    seg_a: np.ndarray                    # (S,) first point index
    seg_b: np.ndarray                    # (S,) last point index (inclusive)
    seg_cluster: np.ndarray              # (S,) clustering_pr cluster id
    length: np.ndarray                   # (S,) cm
    straight: np.ndarray                 # (S,) |end1 - end0| / length
    med_dqdx: np.ndarray                 # (S,) median interior dQ/dx (e/cm)
    end_xyz: np.ndarray                  # (S,2,3)
    end_dir: np.ndarray                  # (S,2,3) unit vector from the end INTO the segment
    end_node: np.ndarray                 # (S,2)
    end_feat: np.ndarray                 # (S,2,len(END_FEATURES)) float64
    node_xyz: np.ndarray                 # (V,3)
    node_degree: np.ndarray              # (V,) segment ends on the node
    node_dqdx: np.ndarray                # (V,) e/cm of the vertex point (NaN if the node has no vertex point)
    node_cluster: np.ndarray             # (V,)
    p_start: np.ndarray = field(default=None)   # (S,) P(end 0 is the start), set by orient()

    @property
    def n_seg(self):
        return len(self.seg_a)

    @property
    def n_node(self):
        return len(self.node_xyz)

    def start_end(self):
        """(S,) int: index (0/1) of the more likely start end."""
        return (self.p_start < 0.5).astype(np.int64)


def _parse_runs(x, cid):
    """-> (seg (n,) segment index or -1, seg_a, seg_b) for the ordered track_fit points.

    A point with an exact copy LATER in the layer is a segment end in the segment block; a point whose
    copies are all EARLIER is a PR vertex point (the vertex block).  Among the remaining points the
    twinned ones alternate start/end of segments, so consecutive twinned points are split by parity."""
    n = len(x)
    later = np.zeros(n, bool)
    earlier = np.zeros(n, bool)
    if n > 1:
        pr = cKDTree(x).query_pairs(TWIN_CM, output_type="ndarray")
        if len(pr):
            lo, hi = pr.min(1), pr.max(1)
            later[lo] = True
            earlier[hi] = True
    twin = later | earlier
    h = np.flatnonzero(~(earlier & ~later))          # drop vertex-block points
    m = len(h)
    seg = np.full(n, -1, np.int64)
    if m == 0:
        return seg, np.zeros(0, np.int64), np.zeros(0, np.int64)
    xh, twh = x[h], twin[h]
    steph = np.linalg.norm(np.diff(xh, axis=0), axis=1)
    bb = np.ones(m, bool)
    bb[1:] = (cid[h][1:] != cid[h][:-1]) | (steph > JUMP_CM) | (steph < DUP_CM) | (np.diff(h) != 1)
    prev_tw = np.zeros(m, bool)
    prev_tw[1:] = twh[:-1] & ~bb[1:]
    st = twh & ~prev_tw                              # first point of a stretch of consecutive twinned points
    st_idx = np.flatnonzero(st)
    sid = np.maximum(np.cumsum(st) - 1, 0)
    pos = np.arange(m) - (st_idx[sid] if len(st_idx) else 0)
    offset = bb[st_idx].astype(np.int64)[sid] if len(st_idx) else np.zeros(m, np.int64)
    end_tw = twh & ((pos - offset) % 2 == 0)         # stretch at a run start begins with a START
    brk = bb.copy()
    brk[1:] |= end_tw[:-1]
    run = np.cumsum(brk) - 1
    rcnt = np.bincount(run)
    ra = np.concatenate([[0], np.cumsum(rcnt)[:-1]])
    keep = rcnt >= 2
    rank = np.cumsum(keep) - 1
    seg[h] = np.where(keep[run], rank[run], -1)
    return seg, h[ra[keep]], h[ra[keep] + rcnt[keep] - 1]


def parse_segments(fit_xyz, fit_cluster_id):
    """Physical PR-graph segments of the ordered track_fit layer (the parse used by build_fit_graph).
    -> (segment_id (n,) int32: 0..S-1 in storage order, -1 for PR vertex points;
        index_along (n,) int32: 0 at the segment's first stored point, -1 for vertex points).
    Segments are contiguous in storage order; the first/last point of each is a PR-graph node."""
    x = np.asarray(fit_xyz, np.float64).reshape(-1, 3)
    n = len(x)
    if n == 0:
        return np.zeros(0, np.int32), np.zeros(0, np.int32)
    seg, seg_a, _ = _parse_runs(x, np.asarray(fit_cluster_id, np.int64).reshape(-1))
    idx = np.full(n, -1, np.int64)
    on = seg >= 0
    if len(seg_a):
        idx[on] = np.flatnonzero(on) - seg_a[seg[on]]
    return seg.astype(np.int32), idx.astype(np.int32)


def build_fit_graph(fit_xyz, fit_q, fit_cluster_id):
    """Parse the track_fit layer into a FitGraph with per-point dQ/dx and per-end Bragg cues
    (connectivity/isolation cues are filled by add_context)."""
    x = np.asarray(fit_xyz, np.float64).reshape(-1, 3)
    cid = np.asarray(fit_cluster_id, np.int64)
    n = len(x)
    dq, clipped = decode_dq(fit_q)
    if n == 0:
        z = np.zeros(0)
        zi = np.zeros(0, np.int64)
        return FitGraph(x, zi, zi, np.zeros(0, bool), z, z, np.zeros(0, bool), z, z, z, np.zeros((0, 3)),
                        zi, zi, zi, z, z, z, np.zeros((0, 2, 3)), np.zeros((0, 2, 3)), np.zeros((0, 2), np.int64),
                        np.zeros((0, 2, len(END_FEATURES))), np.zeros((0, 3)), zi, z, zi, np.zeros(0))
    seg, seg_a, seg_b = _parse_runs(x, cid)
    S = len(seg_a)
    step = np.linalg.norm(np.diff(x, axis=0), axis=1)
    on_seg = seg >= 0
    is_end = np.zeros(n, bool)
    is_end[seg_a] = True
    is_end[seg_b] = True

    # arc length, dx, dQ/dx
    ds = np.zeros(n)
    ds[1:] = step
    ds[~on_seg] = 0.0
    ds[seg_a] = 0.0
    cs = np.cumsum(ds)
    first = seg_a[np.maximum(seg, 0)] if S else np.zeros(n, np.int64)
    s = np.where(on_seg, cs - cs[np.where(on_seg, first, 0)], 0.0)
    nxt = np.zeros(n)
    nxt[:-1] = ds[1:]
    dx = 0.5 * (ds + nxt)
    length = s[seg_b] if S else np.zeros(0)

    # nodes: vertex points + segment end points, merged within NODE_CM
    cand_idx = np.flatnonzero(~on_seg | is_end)
    node = np.full(n, -1, np.int64)
    if len(cand_idx):
        P = x[cand_idx]
        pairs = cKDTree(P).query_pairs(NODE_CM, output_type="ndarray")
        m = len(P)
        A = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(m, m)) if len(pairs) else \
            coo_matrix((m, m))
        V, lab = connected_components(A, directed=False)
        node[cand_idx] = lab
        node_xyz = np.zeros((V, 3))
        np.add.at(node_xyz, lab, P)
        node_xyz /= np.bincount(lab, minlength=V)[:, None]
    else:
        V, node_xyz = 0, np.zeros((0, 3))
    end_node = np.stack([node[seg_a], node[seg_b]], 1) if S else np.zeros((0, 2), np.int64)
    node_degree = np.bincount(end_node.ravel(), minlength=V) if S else np.zeros(V, np.int64)
    # vertex dQ/dx: vertex-point dQ over the half steps of every segment end on the node
    half = np.zeros(V)
    if S:
        np.add.at(half, end_node[:, 0], dx[seg_a])
        np.add.at(half, end_node[:, 1], dx[seg_b])
    vp = np.flatnonzero(~on_seg)
    node_dq = np.full(V, np.nan)
    node_dq[node[vp]] = dq[vp]
    has_ends = half > 0
    node_dqdx = np.where(has_ends, node_dq / np.maximum(half, MIN_DX), np.nan)
    node_cluster = np.zeros(V, np.int64)
    if len(cand_idx):
        node_cluster[node[cand_idx]] = cid[cand_idx]

    dqdx = dq / np.maximum(dx, MIN_DX)
    dx = np.where(on_seg, dx, np.where(node >= 0, half[np.maximum(node, 0)], 0.0))
    dqdx = np.where(on_seg, dqdx, np.where(node >= 0, np.nan_to_num(node_dqdx[np.maximum(node, 0)], nan=0.0),
                                           0.0))
    # smoothed dQ/dx: running median over interior points of the segment (endpoints are shared vertices)
    dqdx_s = dqdx.copy()
    tangent = np.zeros((n, 3))
    med = np.zeros(S)
    if S:
        idx = np.flatnonzero(on_seg)
        sg = seg[idx]
        a, b = seg_a[sg], seg_b[sg]
        short = (b - a) < 2
        lo = np.where(short, a, a + 1)
        hi = np.where(short, b, b - 1)
        win = np.clip(idx[:, None] + np.arange(-SMOOTH_HALF, SMOOTH_HALF + 1)[None], lo[:, None], hi[:, None])
        dqdx_s[idx] = np.median(dqdx[win], axis=1)
        interior = ~is_end[idx]
        med = _group_median(dqdx[idx][interior], sg[interior], S)
        med = np.where(np.isfinite(med), med, _group_median(dqdx[idx], sg, S))
        med = np.where(np.isfinite(med) & (med > 0), med, MIP_DQDX)
        # tangent (storage direction)
        ip = np.clip(idx - 1, a, b)
        inx = np.clip(idx + 1, a, b)
        tangent[idx] = _unit(x[inx] - x[ip])

    end_xyz = np.stack([x[seg_a], x[seg_b]], 1) if S else np.zeros((0, 2, 3))
    straight = np.where(length > 0, np.linalg.norm(end_xyz[:, 1] - end_xyz[:, 0], axis=1) / np.maximum(length, 1e-6),
                        0.0) if S else np.zeros(0)
    end_dir = np.zeros((S, 2, 3))
    if S:
        # direction from the end to the point DIR_CM (or half the length) along the segment
        key = np.where(on_seg, seg * 1e7 + s, np.inf)
        order_ok = np.flatnonzero(on_seg)
        keys = key[order_ok]
        reach = np.minimum(DIR_CM, 0.5 * length)
        j0 = order_ok[np.clip(np.searchsorted(keys, np.arange(S) * 1e7 + reach), 0, len(keys) - 1)]
        j1 = order_ok[np.clip(np.searchsorted(keys, np.arange(S) * 1e7 + length - reach, side="right") - 1,
                              0, len(keys) - 1)]
        j0 = np.clip(j0, seg_a, seg_b)
        j1 = np.clip(j1, seg_a, seg_b)
        j0 = np.where(j0 == seg_a, np.minimum(seg_a + 1, seg_b), j0)
        j1 = np.where(j1 == seg_b, np.maximum(seg_b - 1, seg_a), j1)
        end_dir[:, 0] = _unit(x[j0] - x[seg_a])
        end_dir[:, 1] = _unit(x[j1] - x[seg_b])

    g = FitGraph(x, seg, node, is_end, s, dq, clipped, dx, dqdx, dqdx_s, tangent, seg_a, seg_b,
                 cid[seg_a] if S else np.zeros(0, np.int64), length, straight, med, end_xyz, end_dir, end_node,
                 np.zeros((S, 2, len(END_FEATURES))), node_xyz, node_degree, node_dqdx, node_cluster)
    _bragg_cues(g)
    return g


def _bragg_cues(g):
    """Fill the dQ/dx-based END_FEATURES (bragg_*, end_dqdx) and the geometric ones not needing context."""
    S = g.n_seg
    if S == 0:
        return
    idx = np.flatnonzero((g.seg >= 0) & ~g.is_end)
    sg = g.seg[idx]
    L = g.length[sg]
    med = g.med_dqdx[sg]
    ef = g.end_feat
    for e in (0, 1):
        R = g.s[idx] if e == 0 else L - g.s[idx]
        # near-end mean / max of the smoothed profile, relative to the segment median
        near = (R <= BRAGG_NEAR_CM) & (R > 0)
        gi = sg[near]
        v = g.dqdx_s[idx][near] / med[near]
        cnt = np.bincount(gi, minlength=S)
        mean = _bsum(gi, v, S) / np.maximum(cnt, 1)
        mx = np.zeros(S)
        np.maximum.at(mx, gi, v)
        ef[:, e, EF["bragg_mean3"]] = np.where(cnt > 0, np.log(np.maximum(mean, 1e-3)), 0.0)
        ef[:, e, EF["bragg_max3"]] = np.where(cnt > 0, np.log(np.maximum(mx, 1e-3)), 0.0)
        # log-log slope and Bragg template LLR on unclipped raw points
        w = (R >= BRAGG_RMIN_CM) & (R <= BRAGG_WIN_CM) & ~g.clipped[idx]
        gi = sg[w]
        X = np.log(R[w])
        Y = np.log(np.maximum(g.dqdx[idx][w], 1.0) / MIP_DQDX)
        nn = np.bincount(gi, minlength=S).astype(np.float64)
        sx, sy = _bsum(gi, X, S), _bsum(gi, Y, S)
        sxx, sxy, syy = _bsum(gi, X * X, S), _bsum(gi, X * Y, S), _bsum(gi, Y * Y, S)
        vx = nn * sxx - sx * sx
        ok = (nn >= 4) & (vx > 1e-6 * np.maximum(nn, 1) ** 2)
        ef[:, e, EF["bragg_slope"]] = np.where(ok, (nn * sxy - sx * sy) / np.where(ok, vx, 1.0), 0.0)
        nz = np.maximum(nn, 1)
        rss_f = syy - sy * sy / nz
        aB = (sy - BRAGG_BETA * sx) / nz
        rss_b = syy - 2 * BRAGG_BETA * sxy + BRAGG_BETA ** 2 * sxx - nz * aB * aB
        llr = (rss_f - rss_b) / (2 * SIGMA_LOG ** 2)
        ef[:, e, EF["bragg_llr"]] = np.where(nn >= 3, np.clip(llr, -20, 20), 0.0)
        nd = g.end_node[:, e]
        vq = g.node_dqdx[nd]
        ef[:, e, EF["end_dqdx"]] = np.where(np.isfinite(vq) & (vq > 0), np.log(np.maximum(vq, 1.0) / MIP_DQDX), 0.0)
        ef[:, e, EF["degree"]] = g.node_degree[nd] - 1
        oth = g.end_xyz[:, 1 - e]
        this = g.end_xyz[:, e]
        Ls = np.maximum(g.length, 1e-3)
        ef[:, e, EF["upstream"]] = np.clip((oth[:, 2] - this[:, 2]) / Ls, -1, 1)
        ef[:, e, EF["higher"]] = np.clip((this[:, 1] - oth[:, 1]) / Ls, -1, 1)
        ef[:, e, EF["d_wall"]] = np.log1p(boundary_distance(this))
    # total length of the other segments on each node
    V = g.n_node
    node_len = np.zeros(V)
    np.add.at(node_len, g.end_node[:, 0], g.length)
    np.add.at(node_len, g.end_node[:, 1], g.length)
    for e in (0, 1):
        ef[:, e, EF["chain_len"]] = np.log1p(np.maximum(node_len[g.end_node[:, e]] - g.length, 0.0))


def add_context(g, pts_xyz, pts_q, pts_tree, pts_shower, cand_xyz, tag_nu):
    """Fill connectivity / isolation END_FEATURES using clustering_pr points (tree, q, shower flag),
    PR candidates and tag_nu."""
    S = g.n_seg
    if S == 0:
        return g
    ef = g.end_feat
    E = g.end_xyz.reshape(-1, 3)
    U = g.end_dir.reshape(-1, 3)
    ne = len(E)
    etree = cKDTree(E)
    # other segments with an end within 2 cm
    pr = etree.query_pairs(2.0, output_type="ndarray")
    if len(pr):
        pr = pr[(pr[:, 0] // 2) != (pr[:, 1] // 2)]
    cnt = np.zeros(ne)
    if len(pr):
        # count distinct partner segments per end
        pp = np.concatenate([pr, pr[:, ::-1]])
        key = np.unique(pp[:, 0] * (S + 1) + pp[:, 1] // 2)
        cnt = np.bincount(key // (S + 1), minlength=ne).astype(np.float64)
    ef[:, :, EF["n_ends2"]] = cnt.reshape(S, 2)
    # PR candidates within 2 cm
    c = np.asarray(cand_xyz, np.float64).reshape(-1, 3)
    if len(c):
        nc = cKDTree(c).query_ball_point(E, 2.0, return_length=True)
        ef[:, :, EF["n_cand2"]] = np.asarray(nc, np.float64).reshape(S, 2)
    t = np.asarray(tag_nu, np.float64).reshape(3)
    dt = np.linalg.norm(E - t, axis=1) if np.all(np.isfinite(t)) else np.full(ne, NO_TAG_CM)
    ef[:, :, EF["d_tag"]] = np.log1p(dt).reshape(S, 2)
    # charge beyond each end
    if pts_tree is not None and pts_tree.n > 0:
        M = etree.sparse_distance_matrix(pts_tree, 5.0, output_type="ndarray")
        i, j, d = M["i"], M["j"], M["v"]
        along = np.einsum("ij,ij->i", pts_xyz[j] - E[i], U[i])
        q = pts_q[j]
        out = along < -0.3
        in3 = d <= 3.0
        q_out3 = _bsum(i[in3 & out], q[in3 & out], ne)
        q_all3 = _bsum(i[in3], q[in3], ne)
        shw = out & (pts_shower[j] == 2)
        q_sh5 = _bsum(i[shw], q[shw], ne)
        ef[:, :, EF["q_out3"]] = np.log1p(q_out3 / 1e4).reshape(S, 2)
        ef[:, :, EF["f_out3"]] = (q_out3 / np.maximum(q_all3, 1.0)).reshape(S, 2)
        ef[:, :, EF["shw_out5"]] = np.log1p(q_sh5 / 1e4).reshape(S, 2)
    return g


def orient(g, weights=None):
    """Set g.p_start = P(end 0 is the start end) from a start-likeness score per end."""
    w = ORIENT_WEIGHTS if weights is None else weights
    S = g.n_seg
    if S == 0:
        g.p_start = np.zeros(0)
        return g
    score = np.zeros((S, 2))
    for name, wk in w.items():
        score += wk * g.end_feat[:, :, EF[name]]
    g.p_start = 1.0 / (1.0 + np.exp(-np.clip(score[:, 0] - score[:, 1], -30, 30)))
    return g


def segment_point_features(g):
    """Per fit point, canonical (start -> stop) quantities.  Returns dict of (n,) / (n,3) arrays:
    rr_start (arc length from the likely start end), rr_stop (residual range to the likely stopping end),
    frac (rr_start / L), p (P that the canonical start is right, >= 0.5), dir (unit tangent pointing
    start -> stop), bragg_start / bragg_stop (bragg_llr at each canonical end)."""
    n = len(g.xyz)
    on = g.seg >= 0
    sg = np.maximum(g.seg, 0)
    out = {k: np.zeros(n) for k in ("rr_start", "rr_stop", "frac", "p", "L", "bragg_start", "bragg_stop",
                                    "mean3_start", "mean3_stop")}
    out["dir"] = np.zeros((n, 3))
    if g.n_seg == 0:
        return out
    se = g.start_end()
    flip = se[sg] == 1
    L = g.length[sg]
    rr0 = g.s
    out["rr_start"] = np.where(on, np.where(flip, L - rr0, rr0), 0.0)
    out["rr_stop"] = np.where(on, np.where(flip, rr0, L - rr0), 0.0)
    out["frac"] = np.where(on & (L > 0), out["rr_start"] / np.maximum(L, 1e-6), 0.0)
    ps = np.where(se == 0, g.p_start, 1.0 - g.p_start)
    out["p"] = np.where(on, ps[sg], 0.0)
    out["L"] = np.where(on, L, 0.0)
    llr = g.end_feat[:, :, EF["bragg_llr"]]
    m3 = g.end_feat[:, :, EF["bragg_mean3"]]
    ar = np.arange(g.n_seg)
    out["bragg_start"] = np.where(on, llr[ar, se][sg], 0.0)
    out["bragg_stop"] = np.where(on, llr[ar, 1 - se][sg], 0.0)
    out["mean3_start"] = np.where(on, m3[ar, se][sg], 0.0)
    out["mean3_stop"] = np.where(on, m3[ar, 1 - se][sg], 0.0)
    out["dir"] = np.where(on[:, None], g.tangent * np.where(flip, -1.0, 1.0)[:, None], 0.0)
    return out
