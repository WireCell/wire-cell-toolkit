"""Physics features: the contract between feature code, the cache and the models.

CONTRACT (keep this interface):

FEATURE_SPEC : list of (name, dim, level) tuples.
    level is one of
      'point'  -> one row per clustering_pr point   (event.points_*)
      'fit'    -> one row per track_fit point       (event.fit_*)
      'shower' -> one row per shower_track point    (event.shower_*)   (no features at this level now)
      'cand'   -> one row per PR vertex candidate   (event.cand_*)
      'event'  -> one row per event

compute_features(event) -> dict {name: np.ndarray}
    * `event` is an sbndvtx.event.Event WITHOUT truth (event.truth is None): features are truth-free.
      Inputs: points, fit, shower, cand, tag_nu, wc_reco (all x-fixed reco frame).
    * shape (event.level_size(level), dim) for element levels, (dim,) for 'event'.
    * every name in FEATURE_SPEC present, no extra names, all values finite, float32.
    * deterministic (no RNG).

FP32_FEATURES : names kept float32 in the cache.  Other 'point' features are float16, so they are
    O(1)-scaled (logs, fractions, unit vectors).  Coordinates are never features.

COLUMN NAMES of every multi-dim feature: FEATURE_COLUMNS[name] (tuple of column names).

WHAT THE FEATURES MEAN (details in fitgraph.py / local.py / tokens.py):
  * track_fit is the Wire-Cell PR graph (segments + vertices).  q decodes to dQ = 10*(q+1000) e
    (0 = clipped); dQ/dx uses WC's dx; MIP ~ 47.5k e/cm.  Every segment gets per-end Bragg cues
    (dQ/dx rise, log-log slope, Bragg-template LLR), connectivity (node degree, attached ends),
    isolation (charge beyond the end) and geometry (upstream, wall distance) cues, combined into
    P(start) = which end the particle came from.  "start" = vertex side, "stop" = Bragg / exiting end.
  * point features: local density/charge (0.5/1/3 cm), local PCA (1/3 cm), charge-weighted centroid
    shift at 3 and 10 cm (one-sided charge), PR shower/track tag, the nearest fit segment's oriented
    quantities (residual range to the stop end, distance from the start end, Bragg scores,
    orientation vector x confidence), PR-graph proximity, cluster stats, wall distances, PR proximity.
  * cand features: tokens.CAND_FEATURES local physics descriptors (prongs leaving vs ending,
    occupied directions, one-sidedness, cluster, tag_nu distance), plus the hand-made vertex score.
VALIDATION (scripts/validate_features.py; results/features_validation.json = 5,760 r3 events outside the
dev cache, 4,028 with a target; truth used only for scoring):
  * orientation P(start): 0.917 on L > 5 cm segments with an end within 2 cm of the true vertex
    (proton 0.979 / muon 0.935 / pion 0.921 by range match; L 5-10 cm 0.902, > 50 cm 0.958),
    0.932 on all segments of the vertex's PR-graph component.
  * vertex pick within 1 cm: tag_nu 0.510 nue / 0.507 BNB; hand heuristic_vertex 0.510 / 0.499;
    out-of-fold HistGBM over candidate_positions + CAND_FEATURES 0.617 / 0.581 (tag_nu-present subset
    0.642 vs 0.535 / 0.645 vs 0.577); pool oracle 0.758 / 0.721 (0.678 / 0.647 without cluster ends).
  * compute_features: 1.16 s/event mean on the login node (1 BLAS thread) -> 0.11 h for 22k events / 64 cores.
Token-level helpers for models: sbndvtx.features.tokens (analyze, end_tokens, candidate_positions,
cluster_end_positions, describe_positions, trajectory_tokens, wc_particle_tokens, heuristic_vertex).
"""
from contextlib import nullcontext

import numpy as np
from scipy.spatial import cKDTree

try:                                   # small (n,3) linear algebra is much slower with many BLAS threads
    from threadpoolctl import threadpool_limits
except ImportError:                    # pragma: no cover
    threadpool_limits = None

from . import fitgraph, local, tokens
from .fitgraph import END_FEATURES, MIP_DQDX, build_fit_graph
from .tokens import (CAND_FEATURES, KINDS, analyze, candidate_positions, cluster_end_positions,
                     describe_positions, end_tokens, heuristic_score, heuristic_vertex, trajectory_tokens,
                     wc_particle_tokens)

LEVELS = ("point", "fit", "shower", "cand", "event")
FIT_ASSOC_CM = 2.0          # a point is 'on' a fit segment if the nearest segment fit point is this close
FAR_CM = 1000.0             # distance used when a reference object does not exist
BRAGG_SCALE = 10.0          # bragg_llr / BRAGG_SCALE, clipped to +-2, in point/fit features

FEATURE_COLUMNS = {
    "log_q": ("log1p_q_per_1e3",),
    "log_cluster_size": ("log1p_cluster_size",),
    "density": ("log1p_n_r0.5", "log1p_n_r1", "log1p_n_r3", "log1p_q_r0.5", "log1p_q_r1", "log1p_q_r3"),
    "local_shape": ("linearity_r1", "planarity_r1", "linearity_r3", "planarity_r3", "spread_r3"),
    "local_dir": ("dxx", "dyy", "dzz", "dxy", "dxz", "dyz"),
    "charge_shift": ("shift3_x", "shift3_y", "shift3_z", "shift10_x", "shift10_y", "shift10_z",
                     "log1p_q_r10_per_1e5"),
    "shower_tag": ("is_shower", "is_track"),
    "fit_assoc": ("near", "log1p_d_fit", "log_dqdx_mip", "log_dqdx_rel", "log1p_rr_start", "log1p_rr_stop",
                  "frac_along", "log1p_len", "bragg_start", "bragg_stop", "orient_conf", "dir_x", "dir_y",
                  "dir_z", "straightness"),
    "graph_prox": ("log1p_d_node2", "log1p_node2_degree", "log1p_d_start_end", "log1p_d_stop_end"),
    "cluster": ("log1p_size", "log1p_qtot", "log1p_ext1", "log1p_ext2", "log1p_ext3", "linearity", "is_wc_nu",
                "faces", "frac_z", "log1p_rank"),
    "walls": ("log1p_d_anode", "log1p_d_cathode", "log1p_d_yface", "log1p_d_zface", "tpc_side"),
    "pr_prox": ("log1p_d_cand", "log1p_d_tag", "has_tag"),
    "fit_log_q": ("log1p_q_per_1e3",),
    "fit_track": ("is_vertex_point", "is_seg_end", "log1p_degree", "log_dqdx_mip", "log_dqdx_s_mip",
                  "log_dqdx_rel", "clipped", "log1p_rr_start", "log1p_rr_stop", "frac_along", "log1p_len",
                  "straightness", "orient_conf", "dir_x", "dir_y", "dir_z", "bragg_start", "bragg_stop",
                  "mean3_start", "mean3_stop", "log1p_dx", "seg_is_wc"),
    "cand_desc": CAND_FEATURES,
    "cand_misc": ("is_main", "heuristic_score"),
    "n_clusters": ("n_clusters",),
    "event_summary": ("log1p_n_points", "log1p_n_fit", "log1p_n_seg", "log1p_n_nodes", "log1p_n_cand",
                      "has_tag", "wc_ok", "wc_none", "log1p_n_wc_reco", "log1p_qtot_per_1e5", "log1p_n_seg_gt5",
                      "log1p_max_degree"),
}
_LEVEL_OF = {"log_q": "point", "log_cluster_size": "point", "density": "point", "local_shape": "point",
             "local_dir": "point", "charge_shift": "point", "shower_tag": "point", "fit_assoc": "point",
             "graph_prox": "point", "cluster": "point", "walls": "point", "pr_prox": "point",
             "fit_log_q": "fit", "fit_track": "fit", "cand_desc": "cand", "cand_misc": "cand",
             "n_clusters": "event", "event_summary": "event"}
FEATURE_SPEC = [(n, len(c), _LEVEL_OF[n]) for n, c in FEATURE_COLUMNS.items()]
FP32_FEATURES = set()


def _log_q(q):
    return np.log1p(np.maximum(np.asarray(q, np.float64), 0) / 1000.0)[:, None]


def _nearest(tree_pts, P):
    """Distance from each P to the nearest of tree_pts (FAR_CM if empty) and its index (-1)."""
    if len(tree_pts) == 0 or len(P) == 0:
        return np.full(len(P), FAR_CM), np.full(len(P), -1, np.int64)
    d, j = cKDTree(tree_pts).query(P)
    return np.minimum(d, FAR_CM), j


def _fit_level(an):
    g = an.graph
    n = len(g.xyz)
    F = np.zeros((n, len(FEATURE_COLUMNS["fit_track"])))
    sp = fitgraph.segment_point_features(g)
    if n == 0:
        return F, sp
    on = g.seg >= 0
    sg = np.maximum(g.seg, 0)
    med = g.med_dqdx[sg] if g.n_seg else np.ones(n)
    deg = np.where(g.node >= 0, g.node_degree[np.maximum(g.node, 0)] if g.n_node else 0, 0)
    conf = np.where(on, 2 * sp["p"] - 1, 0.0)
    wc = np.isin(g.seg_cluster, np.asarray(an.wc_clusters, np.int64)).astype(np.float64)
    F[:, 0] = ~on
    F[:, 1] = g.is_end
    F[:, 2] = np.log1p(deg)
    F[:, 3] = np.log(np.maximum(g.dqdx, 1.0) / MIP_DQDX)
    F[:, 4] = np.log(np.maximum(g.dqdx_s, 1.0) / MIP_DQDX)
    F[:, 5] = np.where(on, np.log(np.maximum(g.dqdx_s, 1.0) / np.maximum(med, 1.0)), 0.0)
    F[:, 6] = g.clipped
    F[:, 7] = np.log1p(sp["rr_start"])
    F[:, 8] = np.log1p(sp["rr_stop"])
    F[:, 9] = sp["frac"]
    F[:, 10] = np.log1p(sp["L"])
    F[:, 11] = np.where(on, g.straight[sg] if g.n_seg else 0.0, 0.0)
    F[:, 12] = conf
    F[:, 13:16] = sp["dir"]
    F[:, 16] = np.clip(sp["bragg_start"] / BRAGG_SCALE, -2, 2)
    F[:, 17] = np.clip(sp["bragg_stop"] / BRAGG_SCALE, -2, 2)
    F[:, 18] = sp["mean3_start"]
    F[:, 19] = sp["mean3_stop"]
    F[:, 20] = np.log1p(g.dx)
    F[:, 21] = np.where(on, wc[sg] if g.n_seg else 0.0, 0.0)
    return F, sp


def _point_level(an, fitF, sp):
    """All point-level features (dict name -> (n, dim))."""
    X, Q = an.X, an.Q
    n = len(X)
    out = {"log_q": _log_q(Q)}
    cp = an.cl_point
    out["log_cluster_size"] = np.log1p(cp["size"])[:, None]
    nb = local.point_neighbourhood(X, Q, an.tree)
    out["density"] = np.concatenate([np.log1p(nb["count"]), np.log1p(nb["charge"] / local.Q_UNIT)], 1)
    out["local_shape"] = nb["shape"]
    out["local_dir"] = nb["dir"]
    sh10, q10 = local.far_shift(X, Q)
    out["charge_shift"] = np.concatenate([nb["shift3"], sh10, np.log1p(q10 / 1e5)[:, None]], 1)
    out["shower_tag"] = np.stack([an.shower == 2, an.shower == 1], 1).astype(np.float64)

    g = an.graph
    fa = np.zeros((n, len(FEATURE_COLUMNS["fit_assoc"])))
    segpts = np.flatnonzero(g.seg >= 0)
    d, j = _nearest(g.xyz[segpts], X)
    fa[:, 1] = np.log1p(d)
    near = d < FIT_ASSOC_CM
    if np.any(near):
        k = segpts[j[near]]
        conf = 2 * sp["p"][k] - 1
        fa[near, 0] = 1.0
        fa[near, 2] = fitF[k, 4]
        fa[near, 3] = fitF[k, 5]
        fa[near, 4] = fitF[k, 7]
        fa[near, 5] = fitF[k, 8]
        fa[near, 6] = fitF[k, 9]
        fa[near, 7] = fitF[k, 10]
        fa[near, 8] = fitF[k, 16]
        fa[near, 9] = fitF[k, 17]
        fa[near, 10] = conf
        fa[near, 11:14] = sp["dir"][k] * conf[:, None]
        fa[near, 14] = fitF[k, 11]
    out["fit_assoc"] = fa

    gp = np.zeros((n, 4))
    if g.n_node:
        n2 = np.flatnonzero(g.node_degree >= 2)
        dn, jn = _nearest(g.node_xyz[n2], X)
        gp[:, 0] = np.log1p(dn)
        gp[:, 1] = np.where(jn >= 0, np.log1p(g.node_degree[n2][np.maximum(jn, 0)] if len(n2) else 0), 0.0)
    else:
        gp[:, 0] = np.log1p(FAR_CM)
    if g.n_seg:
        et = tokens.end_tokens(an)
        real = et["length"] > 3.0
        st = real & (et["p_start"] >= 0.5)
        sp_ = real & (et["p_start"] < 0.5)
        gp[:, 2] = np.log1p(_nearest(et["xyz"][st], X)[0])
        gp[:, 3] = np.log1p(_nearest(et["xyz"][sp_], X)[0])
    else:
        gp[:, 2:] = np.log1p(FAR_CM)
    out["graph_prox"] = gp

    out["cluster"] = np.stack([np.log1p(cp["size"]), np.log1p(cp["qtot"] / local.Q_UNIT), np.log1p(cp["ext1"]),
                               np.log1p(cp["ext2"]), np.log1p(cp["ext3"]), cp["lin"], cp["is_wc"], cp["faces"],
                               cp["frac_z"], np.log1p(cp["rank"])], 1) if n else np.zeros((0, 10))
    out["walls"] = local.wall_features(X)
    has_tag = bool(np.all(np.isfinite(an.tag)))
    pp = np.zeros((n, 3))
    pp[:, 0] = np.log1p(_nearest(an.cand, X)[0])
    pp[:, 1] = np.log1p(np.linalg.norm(X - an.tag, axis=1)) if has_tag else np.log1p(FAR_CM)
    pp[:, 2] = float(has_tag)
    out["pr_prox"] = pp
    return out


def _cand_level(an):
    ev = an.event
    D = describe_positions(an, an.cand)
    misc = np.stack([(np.asarray(ev.cand_q) >= 15000).astype(np.float64), heuristic_score(D)], 1) \
        if len(D) else np.zeros((0, 2))
    return D, misc


def _event_level(an):
    ev, g = an.event, an.graph
    return np.array([np.log1p(ev.n_points), np.log1p(ev.n_fit), np.log1p(g.n_seg), np.log1p(g.n_node),
                     np.log1p(ev.n_cand), float(ev.has_tag_nu), float(ev.wc_status == "ok"),
                     float(ev.wc_status == "none"), np.log1p(len(ev.wc_reco)), np.log1p(an.Q.sum() / 1e5),
                     np.log1p(int((g.length > 5).sum()) if g.n_seg else 0),
                     np.log1p(int(g.node_degree.max()) if g.n_node else 0)], np.float64)


def compute_features(event):
    """Truth-free physics features for `event` (see module docstring for the contract).
    Runs with BLAS limited to one thread (see ISSUES.md E2)."""
    with (threadpool_limits(1) if threadpool_limits is not None else nullcontext()):
        return _compute_features(event)


def _compute_features(event):
    an = analyze(event)
    fitF, sp = _fit_level(an)
    feats = _point_level(an, fitF, sp)
    feats["fit_log_q"] = _log_q(event.fit_q)
    feats["fit_track"] = fitF
    feats["cand_desc"], feats["cand_misc"] = _cand_level(an)
    feats["n_clusters"] = np.array([len(an.cl.get("id", []))], np.float64)
    feats["event_summary"] = _event_level(an)
    for name, dim, level in FEATURE_SPEC:
        a = np.asarray(feats[name], np.float64)
        feats[name] = (a.reshape(dim) if level == "event" else a.reshape(-1, dim)).astype(np.float32)
    return feats


def validate_features(feats, event, spec=None):
    """Raise ValueError if `feats` breaks the contract for `event`."""
    spec = FEATURE_SPEC if spec is None else spec
    names = {s[0] for s in spec}
    extra = set(feats) - names
    if extra:
        raise ValueError(f"features not in FEATURE_SPEC: {sorted(extra)}")
    for name, dim, level in (tuple(s[:3]) for s in spec):
        if level not in LEVELS:
            raise ValueError(f"{name}: unknown level {level!r}")
        if name not in feats:
            raise ValueError(f"missing feature {name}")
        a = np.asarray(feats[name])
        want = (dim,) if level == "event" else (event.level_size(level), dim)
        if a.shape != want:
            raise ValueError(f"{name}: shape {a.shape} != {want}")
        if a.dtype.kind == "f" and not np.all(np.isfinite(a)):
            raise ValueError(f"{name}: non-finite values")
