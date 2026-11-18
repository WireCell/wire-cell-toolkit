"""Truth-free inference of one family's ensemble with test-time views and grid shifts: vendored from dv2max/infer.py
(pass1_prepare / pass1_finish / pass2_prepare / pass2_finish and helpers; bodies unchanged except that collate /
to_device come from batch.py and mirror_event from mirror.py).  Dropped: the DataLoader dataset (PrepareDataset).

Views: 'id', 'fy' (y mirror, exact), 'mx' (cathode x mirror, approximate), 'mxfy'; shifts: 1, 2, 4 or 8
half-voxel grid offsets.  All views x shifts of one event run as one batch through every member; their votes
are mapped back to the reco frame and pooled before mode finding (decode.py).

pass-1 record  (decode.compact_record): modes (vote modes then peaks, MODE_FIELDS), binned vote density, hot
               voxels (heat, any), hot points (point head), cluster P(host), event P(visible nu)
pass-2 record  per pass-1 mode of rank < modes_from_pass1: the framed crop (crop.frame), its P(target in
               crop), its modes and vote density; concatenated with a 'crop' index column
"""
import numpy as np
import torch

from . import batch as batchmod
from . import crop as cropmod
from . import decode
from . import inputs as inp
from . import mirror as sa

VIEW_AXIS = {"id": "", "fy": "y", "mx": "x", "mxfy": "xy"}
SHIFTS = {1: [(0, 0, 0)], 2: [(0, 0, 0), (0.5, 0.5, 0.5)],
          4: [(0, 0, 0), (0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5)],
          8: [(a, b, c) for a in (0, 0.5) for b in (0, 0.5) for c in (0, 0.5)]}
CROP_FIELDS = ("mx", "my", "mz", "lox", "loy", "loz", "hix", "hiy", "hiz", "p_crop", "rank1", "score1")


def axis_matrix(axis):
    s = np.ones(3)
    for a in axis:
        s["xyz".index(a)] = -1.0
    return np.diag(s)


def _inference_sample(x):
    x.update({"has_target": np.float32(0), "target": np.zeros(3, np.float32), "visible": np.zeros((0, 3), np.float32),
              "loc_ok": np.float32(1), "host_cluster": np.int64(-1), "crop_label": np.float32(0), "key": ""})
    return x


def forward_members(models, batch, device, autocast_dtype):
    """Outputs of every ensemble member on one collated batch (all views of one event)."""
    b = batchmod.to_device(batch, device)
    outs = []
    with torch.no_grad():
        for m in models:
            if autocast_dtype is not None:
                with torch.autocast(device_type=torch.device(device).type, dtype=autocast_dtype):
                    outs.append(m(b))
            else:
                outs.append(m(b))
    return outs


# ---------------------------------------------------------------------------- pass 1
def pass1_prepare(ev, cfg, pcfg, views=None):
    """CPU part of pass-1 inference (runs in DataLoader workers): the inputs of every view x shift, collated."""
    if ev.truth is not None:
        raise ValueError("inference needs a truth-free event")
    icfg = pcfg["infer"]
    samples, Rs = [], []
    use_prong = cfg["inputs"]["prong_geometry"] and cfg["inputs"]["point_features"] != "charge_only"
    for v in (views or icfg["views"]):
        axis = VIEW_AXIS[v]
        evm = sa.mirror_event(ev, axis) if axis else ev
        prong = inp.prong_meet_points(evm) if use_prong else None
        for sh in SHIFTS[icfg["shifts"]]:
            grid = inp.whole_grid(pcfg["voxel"], pcfg["levels"], np.asarray(sh) * pcfg["voxel"])
            samples.append(_inference_sample(inp.build_inputs(evm, grid, cfg["inputs"], pcfg, prong=prong)))
            Rs.append(axis_matrix(axis))
    return {"key": ev.key, "batch": batchmod.collate(samples), "Rs": Rs,
            "points": np.asarray(ev.points_xyz, np.float64)}


def pass1_finish(outs, prep, icfg):
    """Decode the members' outputs of all views into the pooled pass-1 record."""
    vr = [decode.event_views(o, e, icfg, R) for e, R in enumerate(prep["Rs"]) for o in outs]
    pv = decode.pool_views(vr)
    modes = decode.find_modes(pv, icfg, prep["points"])
    return decode.compact_record(pv, modes, icfg)


def pass1_record(models, ev, cfg, pcfg, device, autocast_dtype=None, views=None):
    """Pooled pass-1 record of one truth-free event (prepare + forward + decode in one call)."""
    prep = pass1_prepare(ev, cfg, pcfg, views)
    return pass1_finish(forward_members(models, prep["batch"], device, autocast_dtype), prep, pcfg["infer"])


def pass1_modes_for_crops(rec1, k):
    """Up to k pass-1 positions to crop at: vote modes first, then peaks (skipping peaks within 1 cm of a mode)."""
    M = rec1["modes"]
    if len(M) == 0:
        return np.zeros((0, 3)), np.zeros(0)
    vote = M[M[:, 9] == decode.SRC_VOTE]
    peak = M[M[:, 9] == decode.SRC_PEAK]
    rows = [r for r in vote[:k]]
    for r in peak:
        if len(rows) >= k:
            break
        if not rows or np.min(np.linalg.norm(np.asarray(rows)[:, :3] - r[:3], axis=1)) > 1.0:
            rows.append(r)
    R = np.asarray(rows).reshape(-1, M.shape[1])
    return R[:, :3], R[:, 3]


# ---------------------------------------------------------------------------- pass 2
def pass2_prepare(ev, rec1, cfg, pcfg):
    """CPU part of pass-2 inference: crops framed at the pass-1 modes of `rec1` (same family), inputs of every
    crop x view x shift in one collated batch.  Returns None if no crop holds imaging points."""
    if ev.truth is not None:
        raise ValueError("inference needs a truth-free event")
    icfg = pcfg["infer"]
    modes, scores = pass1_modes_for_crops(rec1, icfg["modes_from_pass1"])
    use_prong = cfg["inputs"]["prong_geometry"] and cfg["inputs"]["point_features"] != "charge_only"
    prongs, mirrored = {}, {}
    X = np.asarray(ev.points_xyz, np.float64)
    heat = (rec1["hot_xyz"], rec1["hot_heat"])
    crops, samples, Rs, owner = [], [], [], []
    for r, (m, s1) in enumerate(zip(modes, scores)):
        lo, hi, _, _ = cropmod.frame(m, X, ev.points_q, ev.points_cluster_id, rec1.get("p_host", {}),
                                     pcfg["crop_cm"], pcfg["upstream_frac"], pcfg["frame_radius"])
        if cropmod.inside(X, lo, hi).sum() < 3 and len(X):     # mode in empty space: frame on the nearest point
            m = X[int(np.argmin(((X - m) ** 2).sum(1)))]
            lo, hi, _, _ = cropmod.frame(m, X, ev.points_q, ev.points_cluster_id, rec1.get("p_host", {}),
                                         pcfg["crop_cm"], pcfg["upstream_frac"], pcfg["frame_radius"])
        if cropmod.inside(X, lo, hi).sum() < 3:
            continue
        crops.append((np.asarray(m, np.float64), lo, hi, r, s1))
        for v in icfg["views"]:
            axis = VIEW_AXIS[v]
            R = axis_matrix(axis)
            if axis not in mirrored:
                mirrored[axis] = sa.mirror_event(ev, axis) if axis else ev
                prongs[axis] = inp.prong_meet_points(mirrored[axis]) if use_prong else None
            evm = mirrored[axis]
            c1, c2 = lo @ R, hi @ R
            lom, him = np.minimum(c1, c2), np.maximum(c1, c2)
            keep = cropmod.inside(evm.points_xyz, lom, him)
            hm = (np.asarray(heat[0], np.float64) @ R, heat[1])
            for sh in SHIFTS[icfg["shifts"]]:
                grid = inp.crop_grid(lom, him, pcfg["voxel"], pcfg["levels"], np.asarray(sh) * pcfg["voxel"])
                samples.append(_inference_sample(inp.build_inputs(evm, grid, cfg["inputs"], pcfg, keep=keep,
                                                                  heat=hm, mode=m @ R, prong=prongs[axis])))
                Rs.append(R)
                owner.append(len(crops) - 1)
    if not crops:
        return None
    return {"key": ev.key, "batch": batchmod.collate(samples), "Rs": Rs, "owner": owner, "crops": crops, "points": X}


def pass2_finish(outs, prep, icfg):
    """Decode every crop (its views x shifts x members pooled) into the pass-2 record."""
    X = prep["points"]
    crops, crop_modes, cells, masses, pmasses, cell_crop = [], [], [], [], [], []
    hot_xyz, hot_heat, hot_any, hot_pt, hot_pt_p = [], [], [], [], []
    extra = {k: [] for k in EXTRA_HOT}
    owner = np.asarray(prep["owner"])
    for c, (m, lo, hi, r, s1) in enumerate(prep["crops"]):
        vr = [decode.event_views(o, e, icfg, prep["Rs"][e]) for e in np.flatnonzero(owner == c) for o in outs]
        pv = decode.pool_views(vr)
        md = decode.find_modes(pv, icfg, X[cropmod.inside(X, lo, hi)])
        rc = decode.compact_record(pv, md, icfg)
        crops.append([*m, *lo, *hi, rc["event_p"], r, s1])
        crop_modes.append(np.concatenate([md, np.full((len(md), 1), len(crops) - 1)], 1))
        cells.append(rc["cells"])
        masses.append(rc["cell_mass"])
        pmasses.append(rc["cell_pmass"])
        cell_crop.append(np.full(len(rc["cells"]), len(crops) - 1, np.int16))
        hot_xyz.append(rc["hot_xyz"])
        hot_heat.append(rc["hot_heat"])
        hot_any.append(rc["hot_any"])
        if "hot_pt_xyz" in rc:
            hot_pt.append(rc["hot_pt_xyz"])
            hot_pt_p.append(rc["hot_pt_p"])
        for k in EXTRA_HOT:                              # C2 / C3 arms: fields at the crop's hot voxels
            if k in rc:
                extra[k].append(rc[k])
    return pass2_pack(crops, crop_modes, cells, masses, pmasses, cell_crop, hot_xyz, hot_heat, hot_any, hot_pt,
                      hot_pt_p, extra)


EXTRA_HOT = ("hot_c3", "hot_shell", "hot_cone")


def pass2_pack(crops, crop_modes, cells, masses, pmasses, cell_crop, hot_xyz, hot_heat, hot_any, hot_pt, hot_pt_p,
               extra=None):
    nf = len(decode.MODE_FIELDS) + 1
    cat = lambda L, shape, dt: np.concatenate(L).astype(dt) if L else np.zeros(shape, dt)     # noqa: E731
    rec = {"crops": np.asarray(crops, np.float64).reshape(-1, len(CROP_FIELDS)),
           "modes": cat(crop_modes, (0, nf), np.float64), "cells": cat(cells, (0, 3), np.float32),
           "cell_mass": cat(masses, (0,), np.float32), "cell_pmass": cat(pmasses, (0,), np.float32),
           "cell_crop": cat(cell_crop, (0,), np.int16), "hot_xyz": cat(hot_xyz, (0, 3), np.float32),
           "hot_heat": cat(hot_heat, (0,), np.float32), "hot_any": cat(hot_any, (0,), np.float32),
           "hot_pt_xyz": cat(hot_pt, (0, 3), np.float32), "hot_pt_p": cat(hot_pt_p, (0,), np.float32)}
    for k, L in (extra or {}).items():
        if L and len(L) == len(hot_xyz):                 # every crop carries the field (same model config)
            rec[k] = np.concatenate(L).astype(np.float32)
    return rec


def empty_pass2():
    return pass2_pack(*([[]] * 11))


def pass2_record(models, ev, rec1, cfg, pcfg, device, autocast_dtype=None):
    """Pass-2 record of one truth-free event (prepare + forward + decode in one call)."""
    prep = pass2_prepare(ev, rec1, cfg, pcfg)
    if prep is None:
        return empty_pass2()
    return pass2_finish(forward_members(models, prep["batch"], device, autocast_dtype), prep, pcfg["infer"])


def best_pass2(rec2):
    """(xyz, score) of the best pass-2 mode by score x P(target in crop); None if no crop."""
    M = rec2["modes"]
    if len(M) == 0:
        return None, float("nan")
    cp = rec2["crops"][M[:, -1].astype(int), 9]
    s = M[:, 3] * cp
    vote = M[:, 9] == decode.SRC_VOTE
    if vote.any():
        s = np.where(vote, s, -np.inf)
    i = int(np.argmax(s))
    return M[i, :3], float(s[i])


def best_pass1(rec1):
    M = rec1["modes"]
    if len(M) == 0:
        return None, float("nan")
    return M[0, :3], float(M[0, 3])
