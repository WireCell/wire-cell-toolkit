"""Truth-free network inputs from an sbndvtx Event.  NOTHING here reads event.truth (tests/test_leakage.py).

build_inputs(ev, grid, icfg, ...) -> dict of numpy arrays for one event (or one crop):
  vox_grid   (M,3) int32   voxel indices on `grid` (translation by the grid origin; never given to the net as xyz)
  vox_coord  (M,3) float32 point centroid of each voxel, reco frame (anchor for offsets / decoding, not an input)
  vox_feat   (M,Fv) float32 inputs: mean of the point features, log summed charge, log count, centroid offset in the
                            voxel, rendered Wire-Cell marks, prong-meet geometry, optional pass-2 channels
  vox_cluster(M,)  int64   clustering_pr cluster of the highest-charge point in the voxel
  pt_feat    (P,Fp) float32 per imaging point: cached point features + sub-voxel offset (SPVCNN branch)
  pt2vox     (P,)  int64
  pt_xyz     (P,3) float32 (targets / decoding only)
  ev_feat    (13,) float32 event_summary + log1p(n_clusters) (event context)
No absolute coordinate is an input: positions enter only as wall distances (cached 'walls' columns), sub-voxel
offsets and distances to marks.  Point features 'all' = the 65 cached point columns; 'charge_only' = log_q;
'imaging' = the 40 columns computed from the imaging points and their clusters alone (no Wire-Cell pattern-recognition
product: no track fit, PR graph, PR candidates / tag_nu, shower/track tags or Wire-Cell neutrino-cluster flag); it
needs marks and prong_geometry off, and sets the PR entries of the event context to their 'absent' values.
"""
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from .features import FEATURE_COLUMNS, FEATURE_SPEC

POINT_FEATURES = [n for n, _, lv in FEATURE_SPEC if lv == "point"]
N_POINT_COLS = sum(d for n, d, lv in FEATURE_SPEC if lv == "point")          # 65
MARK_CHANNELS = ("pr_cand", "pr_main", "tag_nu", "wc_vtx", "wcp_track", "wcp_shower", "wcp_other", "wcp_loge")
PRONG_CHANNELS = ("log1p_d_meet", "log1p_n_meet2", "log1p_d_meet_start")
FAR_LOG = float(np.log1p(1000.0))
PRONG_FAR_CM = 30.0
TRACK_PDG = (13, 211, 321, 2212)
SHOWER_PDG = (11, 22)
# point-feature columns that carry the Wire-Cell neutrino decision, and their values when it is absent
# (mark dropout makes an event look like one without a tag_nu / reco neutrino)
MARK_COLUMNS = {("pr_prox", "log1p_d_tag"): FAR_LOG, ("pr_prox", "has_tag"): 0.0, ("cluster", "is_wc_nu"): 0.0}
EVENT_MARK_COLUMNS = {"has_tag": 0.0, "wc_ok": 0.0, "wc_none": 1.0, "log1p_n_wc_reco": 0.0}
# detector box of the whole-event grid (reco frame: x +-201.2, y +-200, z 0.15-501; sbndvtx.features.geometry)
DET_LO = np.array([-205.0, -205.0, -5.0])
DET_HI = np.array([205.0, 205.0, 505.0])


def _col_offsets():
    off, k = {}, 0
    for n in POINT_FEATURES:
        off[n] = k
        k += len(FEATURE_COLUMNS[n])
    return off


COL_OFF = _col_offsets()
# 'imaging' mode (2026-10-06): inputs available before Wire-Cell's pattern recognition
PR_POINT_GROUPS = ("shower_tag", "fit_assoc", "graph_prox", "pr_prox")
IMAGING_COLS = np.array([COL_OFF[n] + i for n in POINT_FEATURES for i, c in enumerate(FEATURE_COLUMNS[n])
                         if n not in PR_POINT_GROUPS and (n, c) != ("cluster", "is_wc_nu")], np.int64)
N_IMAGING_COLS = len(IMAGING_COLS)                                            # 40
PR_EVENT_ABSENT = {"log1p_n_fit": 0.0, "log1p_n_seg": 0.0, "log1p_n_nodes": 0.0, "log1p_n_cand": 0.0, "has_tag": 0.0,
                   "wc_ok": 0.0, "wc_none": 1.0, "log1p_n_wc_reco": 0.0, "log1p_n_seg_gt5": 0.0,
                   "log1p_max_degree": 0.0}


def mark_cols(icfg):
    """Indices into MARK_CHANNELS of the rendered mark channels (inputs.mark_channels; None = all)."""
    ch = icfg.get("mark_channels") or MARK_CHANNELS
    return [MARK_CHANNELS.index(c) for c in ch]


def n_point_cols(mode):
    return {"charge_only": 1, "imaging": N_IMAGING_COLS}.get(mode, N_POINT_COLS)


def col_index(name, col):
    return COL_OFF[name] + FEATURE_COLUMNS[name].index(col)


# ---------------------------------------------------------------------------- grids
@dataclass
class Grid:
    voxel: float
    origin: np.ndarray          # (3,) float64, reco frame position of grid index 0 (lower corner)
    shape: tuple                # (3,) ints, padded to a multiple of 2**(levels-1)

    def index(self, xyz):
        return np.floor((np.asarray(xyz, np.float64) - self.origin) / self.voxel).astype(np.int64)


def _pad(n, m):
    return int(-(-int(n) // m) * m)


def whole_grid(voxel, levels, shift=(0.0, 0.0, 0.0)):
    """Detector-box grid; `shift` (cm, in [0, voxel)) moves the voxel boundaries (sub-voxel shift aug / TTA)."""
    m = 2 ** (levels - 1)
    origin = DET_LO - np.asarray(shift, np.float64)
    shape = tuple(_pad(np.ceil((DET_HI - origin)[k] / voxel) + 1, m) for k in range(3))
    return Grid(float(voxel), origin, shape)


def crop_grid(lo, hi, voxel, levels, shift=(0.0, 0.0, 0.0)):
    m = 2 ** (levels - 1)
    origin = np.asarray(lo, np.float64) - np.asarray(shift, np.float64)
    shape = tuple(_pad(np.ceil((np.asarray(hi, np.float64) - origin)[k] / voxel) + 1, m) for k in range(3))
    return Grid(float(voxel), origin, shape)


# ---------------------------------------------------------------------------- point features
def point_feature_matrix(ev, mode="all"):
    """(P, 65) float32 cached point features (mode 'all') or (P, 1) log_q (mode 'charge_only')."""
    if mode == "charge_only":
        return np.asarray(ev.features["log_q"], np.float32).reshape(-1, 1)
    return np.concatenate([np.asarray(ev.features[n], np.float32).reshape(ev.n_points, -1) for n in POINT_FEATURES],
                          1)


def apply_gain(F, q, gain, mode="all"):
    """Charge-scale augmentation (q -> gain * q) applied consistently to the cached log-charge columns.
    Returns (F', q').  Charge sums are inverted from their log1p columns; dQ/dx vs MIP shifts by log(gain) only
    on points associated with a fit segment (fit_assoc.near), which are the only points where it is non-zero."""
    if gain == 1.0:
        return F, q
    F = F.copy()
    q = np.asarray(q, np.float32) * np.float32(gain)
    if mode == "charge_only":
        F[:, 0] = np.log1p(np.maximum(q, 0) / 1e3)
        return F, q
    F[:, col_index("log_q", "log1p_q_per_1e3")] = np.log1p(np.maximum(q, 0) / 1e3)
    for name, col in (("density", "log1p_q_r0.5"), ("density", "log1p_q_r1"), ("density", "log1p_q_r3"),
                      ("charge_shift", "log1p_q_r10_per_1e5"), ("cluster", "log1p_qtot")):
        j = col_index(name, col)
        F[:, j] = np.log1p(np.float32(gain) * np.expm1(F[:, j]))
    near = F[:, col_index("fit_assoc", "near")] > 0.5
    F[near, col_index("fit_assoc", "log_dqdx_mip")] += np.float32(np.log(gain))
    return F, q


def drop_mark_columns(F, mode="all"):
    """Set the point-feature columns that encode the Wire-Cell neutrino decision to their 'absent' values."""
    if mode == "charge_only":
        return F
    F = F.copy()
    for (name, col), v in MARK_COLUMNS.items():
        F[:, col_index(name, col)] = v
    return F


def event_features(ev, drop_marks=False, imaging=False, absent=()):
    es = np.asarray(ev.features["event_summary"], np.float32).reshape(-1).copy()
    cols = FEATURE_COLUMNS["event_summary"]
    if drop_marks:
        for c, v in EVENT_MARK_COLUMNS.items():
            es[cols.index(c)] = v
    if imaging:
        for c, v in PR_EVENT_ABSENT.items():
            es[cols.index(c)] = v
    for c in absent:
        es[cols.index(c)] = PR_EVENT_ABSENT[c]
    nc = np.log1p(np.asarray(ev.features["n_clusters"], np.float32).reshape(-1))
    return np.concatenate([es, nc]).astype(np.float32)


N_EVENT_FEATURES = len(FEATURE_COLUMNS["event_summary"]) + 1


# ---------------------------------------------------------------------------- marks
def mark_positions(ev, drop_marks=False):
    """{channel: (positions (m,3), values (m,))} of the Wire-Cell marks (truth-free reco quantities)."""
    out = {}
    cand = np.asarray(ev.cand_xyz, np.float64).reshape(-1, 3)
    cq = np.asarray(ev.cand_q, np.float64)
    out["pr_cand"] = (cand, np.ones(len(cand)))
    if drop_marks:
        z = (np.zeros((0, 3)), np.zeros(0))
        for c in MARK_CHANNELS[1:]:
            out[c] = z
        return out
    main = cq >= 15000
    out["pr_main"] = (cand[main], np.ones(int(main.sum())))
    out["tag_nu"] = (np.asarray(ev.tag_nu, np.float64).reshape(1, 3), np.ones(1)) if ev.has_tag_nu else \
        (np.zeros((0, 3)), np.zeros(0))
    vtx = [np.asarray(w.vtx, np.float64) for w in ev.wc_reco if np.all(np.isfinite(w.vtx))]
    out["wc_vtx"] = (np.asarray(vtx).reshape(-1, 3), np.ones(len(vtx)))
    trk, shw, oth, le, lev = [], [], [], [], []
    for w in ev.wc_reco:
        for p in w.particles:
            s = np.asarray(p.start, np.float64)
            if not np.all(np.isfinite(s)):
                continue
            a = abs(int(p.pdg))
            (trk if a in TRACK_PDG else shw if a in SHOWER_PDG else oth).append(s)
            le.append(s)
            lev.append(np.log1p(max(float(p.energy), 0.0) / 100.0))
    out["wcp_track"] = (np.asarray(trk).reshape(-1, 3), np.ones(len(trk)))
    out["wcp_shower"] = (np.asarray(shw).reshape(-1, 3), np.ones(len(shw)))
    out["wcp_other"] = (np.asarray(oth).reshape(-1, 3), np.ones(len(oth)))
    out["wcp_loge"] = (np.asarray(le).reshape(-1, 3), np.asarray(lev, np.float64))
    return out


def render_marks(vox_xyz, marks, sigma):
    """(M, len(MARK_CHANNELS)) max over marks of value * exp(-d^2 / 2 sigma^2) within 3 sigma."""
    M = len(vox_xyz)
    out = np.zeros((M, len(MARK_CHANNELS)), np.float32)
    if M == 0:
        return out
    tree = None
    for k, c in enumerate(MARK_CHANNELS):
        P, v = marks.get(c, (np.zeros((0, 3)), np.zeros(0)))
        ok = np.all(np.isfinite(P), 1) if len(P) else np.zeros(0, bool)
        P, v = P[ok], v[ok]
        if len(P) == 0:
            continue
        if tree is None:
            tree = cKDTree(vox_xyz)
        for p, val in zip(P, v):
            idx = tree.query_ball_point(p, 3 * sigma)
            if not idx:
                continue
            idx = np.asarray(idx)
            d2 = ((vox_xyz[idx] - p) ** 2).sum(1)
            out[idx, k] = np.maximum(out[idx, k], val * np.exp(-0.5 * d2 / sigma ** 2))
    return out


# ---------------------------------------------------------------------------- prong-meet geometry
def ray_closest(Pi, u, Pj, v, max_cos=0.97):
    """Closest approach of the lines Pi + t u and Pj + s v (unit u, v; row-wise).  Returns (t, s, gap, midpoint,
    near_parallel) with t, s the line parameters (cm) of the two closest points."""
    w0 = Pi - Pj
    b = (u * v).sum(1)
    d = (u * w0).sum(1)
    e = (v * w0).sum(1)
    par = np.abs(b) > max_cos
    den = np.where(par, 1.0, 1.0 - b * b)
    t = (b * e - d) / den
    s = (e - b * d) / den
    ci, cj = Pi + t[:, None] * u, Pj + s[:, None] * v
    return t, s, np.linalg.norm(ci - cj, axis=1), 0.5 * (ci + cj), par


def end_tokens(ev):
    """sbndvtx end tokens of the PR-graph segments of a truth-free event (None without a track_fit layer)."""
    from . import tokens
    if ev.n_fit == 0:
        return None
    return tokens.end_tokens(tokens.analyze(ev.strip_truth() if ev.truth is not None else ev))


def prong_meet_points(ev, min_len=3.0, pair_cm=10.0, max_gap=1.5, back_cm=5.0, max_cos=0.97, et=None):
    """Closest-approach points of pairs of back-extended PR-graph segment axes (judge fix: explicit multi-prong
    geometry).  Each segment end with length >= min_len defines a ray from the end away from its segment
    (x(t) = end - t * dir_in, dir_in = unit vector into the segment).  For end pairs of different segments
    closer than pair_cm, the closest approach of the two lines is kept when it lies within [-1, back_cm] cm of
    both ends along their rays, the gap is < max_gap and the lines are not near-parallel.  et: precomputed
    end_tokens(ev).
    Returns (meet (K,3) float64, start_pair (K,) bool = both ends are likely segment starts (P(start) > 0.5))."""
    empty = (np.zeros((0, 3)), np.zeros(0, bool))
    if ev.n_fit == 0:
        return empty
    et = end_tokens(ev) if et is None else et
    if et is None:
        return empty
    keep = np.asarray(et["length"]) >= min_len
    if keep.sum() < 2:
        return empty
    P = np.asarray(et["xyz"], np.float64)[keep]
    D = np.asarray(et["dir"], np.float64)[keep]
    D = D / np.maximum(np.linalg.norm(D, axis=1, keepdims=True), 1e-9)
    seg = np.asarray(et["seg"])[keep]
    ps = np.asarray(et["p_start"], np.float64)[keep]
    pairs = cKDTree(P).query_pairs(pair_cm, output_type="ndarray")
    if len(pairs) == 0:
        return empty
    i, j = pairs[:, 0], pairs[:, 1]
    ok = seg[i] != seg[j]
    i, j = i[ok], j[ok]
    if len(i) == 0:
        return empty
    t, s, gap, mid, par = ray_closest(P[i], -D[i], P[j], -D[j], max_cos)     # rays away from the segments
    good = (~par) & (t >= -1.0) & (t <= back_cm) & (s >= -1.0) & (s <= back_cm) & (gap < max_gap)
    if not good.any():
        return empty
    meet = mid[good]
    start = (ps[i] > 0.5)[good] & (ps[j] > 0.5)[good]
    return meet, start


def render_prong(vox_xyz, meet, start):
    """(M, 3): log1p distance to the nearest meet point (capped), log1p meet points within 2 cm, log1p distance to
    the nearest start-start meet point."""
    M = len(vox_xyz)
    out = np.full((M, len(PRONG_CHANNELS)), np.log1p(PRONG_FAR_CM), np.float32)
    out[:, 1] = 0.0
    if M == 0 or len(meet) == 0:
        return out
    t = cKDTree(meet)
    d, _ = t.query(vox_xyz, distance_upper_bound=PRONG_FAR_CM)
    out[:, 0] = np.log1p(np.minimum(d, PRONG_FAR_CM))
    out[:, 1] = np.log1p(np.asarray(t.query_ball_point(vox_xyz, 2.0, return_length=True), np.float64))
    if start.any():
        d2, _ = cKDTree(meet[start]).query(vox_xyz, distance_upper_bound=PRONG_FAR_CM)
        out[:, 2] = np.log1p(np.minimum(d2, PRONG_FAR_CM))
    return out


# ---------------------------------------------------------------------------- main builder
def n_voxel_features(icfg, pass_cfg):
    base = n_point_cols(icfg["point_features"])
    n = base + 2 + 3
    if icfg["point_features"] != "charge_only":
        n += (len(mark_cols(icfg)) if icfg["marks"] else 0) + (len(PRONG_CHANNELS) if icfg["prong_geometry"] else 0)
    if pass_cfg.get("kind") == "crop":
        n += int(bool(pass_cfg.get("heat_input"))) + int(bool(pass_cfg.get("mode_input")))
    return n


def n_point_features(icfg):
    return n_point_cols(icfg["point_features"]) + 3


def build_inputs(ev, grid, icfg, pass_cfg, keep=None, gain=1.0, drop_marks=False, heat=None, mode=None,
                 prong=None):
    """Inputs of one event (or one crop) on `grid`.  `ev` must be truth-free.
    keep: boolean point mask (crop box, point dropout, vertex hole) applied before voxelisation.
    heat: (xyz (h,3), p (h,)) pass-1 heat voxels (pass-2 heat channel); mode: (3,) pass-1 mode (Gaussian channel).
    prong: precomputed prong_meet_points(ev) (it is independent of the grid)."""
    if ev.truth is not None:
        raise ValueError("build_inputs needs a truth-free event (ev.strip_truth())")
    mode_f = icfg["point_features"]
    charge_only = mode_f == "charge_only"
    F = point_feature_matrix(ev, mode_f)
    q = np.asarray(ev.points_q, np.float32)
    F, q = apply_gain(F, q, gain, mode_f)
    if drop_marks:
        F = drop_mark_columns(F, mode_f)
    if mode_f == "imaging":
        F = F[:, IMAGING_COLS]
    X = np.asarray(ev.points_xyz, np.float64).reshape(-1, 3)
    cid = np.asarray(ev.points_cluster_id, np.int64)
    m = np.ones(len(X), bool) if keep is None else np.asarray(keep, bool).copy()
    g = grid.index(X)
    inside = np.all((g >= 0) & (g < np.asarray(grid.shape)), 1)
    n_outside = int((m & ~inside).sum())
    m &= inside
    if not m.any():
        raise ValueError(f"{ev.key}: no points on the grid")
    X, F, q, cid, g = X[m], F[m], q[m], cid[m], g[m]
    P = len(X)
    s1, s2 = grid.shape[1], grid.shape[2]
    key = (g[:, 0] * s1 + g[:, 1]) * s2 + g[:, 2]
    ukey, inv = np.unique(key, return_inverse=True)
    M = len(ukey)
    cnt = np.bincount(inv, minlength=M).astype(np.float64)
    vg = np.stack([ukey // (s1 * s2), (ukey // s2) % s1, ukey % s2], 1).astype(np.int32)
    # per-voxel sums in one pass: sort points by voxel, reduceat over the voxel starts (every voxel is occupied)
    order = np.argsort(inv, kind="stable")
    starts = np.r_[0, np.flatnonzero(np.diff(inv[order])) + 1]
    sums = np.add.reduceat(np.concatenate([X[order], F[order].astype(np.float64),
                                           np.maximum(q[order], 0).astype(np.float64)[:, None]], 1), starts, axis=0)
    cen = sums[:, :3] / cnt[:, None]
    Fm = sums[:, 3:3 + F.shape[1]] / cnt[:, None]
    qs = sums[:, -1]
    vox_center = grid.origin + (vg.astype(np.float64) + 0.5) * grid.voxel
    parts = [Fm, np.log1p(qs / 1e4)[:, None], np.log1p(cnt)[:, None], (cen - vox_center) / grid.voxel]
    if not charge_only:
        if icfg["marks"]:
            parts.append(render_marks(cen, mark_positions(ev, drop_marks),
                                      max(icfg["mark_sigma_cm"], grid.voxel))[:, mark_cols(icfg)])
        if icfg["prong_geometry"]:
            meet, start = prong if prong is not None else prong_meet_points(ev)
            parts.append(render_prong(cen, meet, start))
    if pass_cfg.get("kind") == "crop":
        if pass_cfg.get("heat_input"):
            h = np.zeros(M, np.float32)
            if heat is not None and len(heat[0]) and M:
                hx, hp = np.asarray(heat[0], np.float64), np.asarray(heat[1], np.float64)
                d, j = cKDTree(hx).query(cen, distance_upper_bound=1.5)
                ok = np.isfinite(d)
                h[ok] = hp[j[ok]]
            parts.append(h[:, None])
        if pass_cfg.get("mode_input"):
            gm = np.zeros(M, np.float32)
            if mode is not None and M:
                d2 = ((cen - np.asarray(mode, np.float64)) ** 2).sum(1)
                gm = np.exp(-0.5 * d2 / pass_cfg["mode_sigma"] ** 2)
            parts.append(np.asarray(gm, np.float32)[:, None])
    vox_feat = np.concatenate([np.asarray(p, np.float32).reshape(M, -1) for p in parts], 1)
    # cluster of the highest-charge point of each voxel
    vcl = np.full(M, -1, np.int64)
    if P:
        o = np.lexsort((-q, inv))
        first = o[np.r_[True, inv[o][1:] != inv[o][:-1]]]
        vcl[inv[first]] = cid[first]
    sub = ((X - (grid.origin + (g.astype(np.float64) + 0.5) * grid.voxel)) / grid.voxel).astype(np.float32)
    return {"vox_grid": vg, "vox_coord": cen.astype(np.float32), "vox_feat": vox_feat, "vox_cluster": vcl,
            "pt_feat": np.concatenate([F, sub], 1).astype(np.float32), "pt2vox": inv.astype(np.int64),
            "pt_xyz": X.astype(np.float32),
            "ev_feat": event_features(ev, drop_marks, mode_f == "imaging", icfg.get("event_absent", ())),
            "grid_shape": np.asarray(grid.shape, np.int64), "n_outside": n_outside}
