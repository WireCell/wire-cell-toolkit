"""Exact y-mirror (and approximate x-mirror) of a truth-free Event: vendored from sbndvtx/augment.py (odd_columns,
MIRROR_ODD_COLUMNS, mirror_features, mirror_event; truth handling reduced to `truth is None`).  The pass-2 'fy' test-time
view uses mirror_event(ev, 'y') with features='flip' (odd feature columns negated, exact in float16).
"""
import copy
import dataclasses

import numpy as np

AXES = {"x": 0, "y": 1, "z": 2}


# ============================================================================ mirrors
def odd_columns(axis, feature_columns=None):
    """{feature name: tuple of column names} that change sign under a mirror of `axis` ('x' or 'y').
    Rule on FEATURE_COLUMNS names: vector components '<stem>_<axis>' of a shift / direction, sign-free
    products 'd<ab>' with the axis exactly once, and 'tpc_side' (= sign(x)) for x."""
    if axis not in ("x", "y"):
        raise ValueError(f"mirror axis must be 'x' or 'y', got {axis!r}")
    if feature_columns is None:
        from .features import FEATURE_COLUMNS as feature_columns
    out = {}
    for name, cols in feature_columns.items():
        odd = []
        for c in cols:
            stem, _, comp = c.rpartition("_")
            if comp == axis and stem in ("shift3", "shift10", "dir"):
                odd.append(c)
            elif len(c) == 3 and c[0] == "d" and set(c[1:]) <= set("xyz") and c[1] != c[2] and c[1:].count(axis) == 1:
                odd.append(c)
            elif c == "tpc_side" and axis == "x":
                odd.append(c)
        if odd:
            out[name] = tuple(odd)
    return out


MIRROR_ODD_COLUMNS = {
    "y": {"local_dir": ("dxy", "dyz"), "charge_shift": ("shift3_y", "shift10_y"), "fit_assoc": ("dir_y",),
          "fit_track": ("dir_y",)},
    "x": {"local_dir": ("dxy", "dxz"), "charge_shift": ("shift3_x", "shift10_x"), "fit_assoc": ("dir_x",),
          "fit_track": ("dir_x",), "walls": ("tpc_side",)},
}


def _axes(axis):
    ax = tuple(axis)
    if not ax or any(a not in ("x", "y") for a in ax) or len(set(ax)) != len(ax):
        raise ValueError(f"axis must be 'x', 'y' or 'xy', got {axis!r}")
    return ax


def _sign(axis):
    s = np.ones(3, np.float32)
    for a in _axes(axis):
        s[AXES[a]] = -1.0
    return s


def mirror_features(features, axis="y", feature_columns=None):
    """New feature dict with the odd columns of `axis` ('x', 'y' or 'xy') negated.  dtypes are kept
    (negation is exact in float16).  Names not in FEATURE_COLUMNS pass through unchanged."""
    if feature_columns is None:
        from .features import FEATURE_COLUMNS as feature_columns
    out = {k: np.array(v, copy=True) for k, v in features.items()}
    for a in _axes(axis):
        for name, cols in odd_columns(a, feature_columns).items():
            if name not in out:
                continue
            arr = out[name]
            idx = [feature_columns[name].index(c) for c in cols]
            if arr.ndim == 1:
                arr[idx] = -arr[idx]
            else:
                arr[:, idx] = -arr[:, idx]
    return out


def _flip(v, s):
    return None if v is None else (np.asarray(v, np.float32) * s).astype(np.float32)


def mirror_truth(truth, axis="y"):
    """Inference events carry no truth (the original mirrors every truth position)."""
    if truth is None:
        return None
    raise ValueError("dv2max_vertexer.mirror: events must be truth-free")


def mirror_event(ev, axis="y", features="flip"):
    """Mirrored copy of an Event (inputs, Wire-Cell reco rows and truth).  axis 'y' (exact), 'x' (approximate
    cathode mirror, see module docstring) or 'xy'.  features: 'flip' (negate odd columns), 'recompute'
    (compute_features on the mirrored truth-free event) or 'drop' ({}).  Point order, cluster ids and
    track_fit segment ids are unchanged; key / index / fold are kept."""
    s = _sign(axis)
    wc = copy.deepcopy(ev.wc_reco)
    for w in wc:
        w.vtx = _flip(w.vtx, s)
        for p in w.particles:
            p.start = _flip(p.start, s)
    out = dataclasses.replace(
        ev, points_xyz=_flip(ev.points_xyz, s).reshape(-1, 3), fit_xyz=_flip(ev.fit_xyz, s).reshape(-1, 3),
        shower_xyz=_flip(ev.shower_xyz, s).reshape(-1, 3), cand_xyz=_flip(ev.cand_xyz, s).reshape(-1, 3),
        tag_nu=_flip(ev.tag_nu, s).reshape(3), wc_reco=wc, truth=mirror_truth(ev.truth, axis), features={})
    if features == "flip":
        out.features = mirror_features(ev.features, axis)
    elif features == "recompute":
        from .features import compute_features
        out.features = compute_features(out.strip_truth())
    elif features != "drop":
        raise ValueError(f"features must be 'flip', 'recompute' or 'drop', got {features!r}")
    return out
