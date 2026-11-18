"""Vote decoding.

Every voxel with P(any) >= vote_min_p votes at anchor + offset (fp64 in numpy) with weight P(any); the
n_vote strongest votes are kept per view.  Votes of all test-time views and ensemble members are mapped back to
the reco frame and POOLED before mode finding (weights / number of members), so vote densities are fused inside
each mode and coordinates are never averaged across modes.  Mean-shift (Gaussian kernel, bandwidth
ms_bandwidth) from NMS-selected seeds finds the modes; a mode's score is its primary-weighted vote share
  score = sum_{votes within 1 cm} w * p_pri / sum_all w * p_pri
(p_pri = the pass's vertexness at the voting voxel).  Modes are kept after NMS (nms_cm) up to top_k; the
highest-vertexness voxels are added as 'peak' hypotheses (n_peaks).  With snap, each mode also gets its nearest
imaging point (the ranker / CV decide which to use).
Everything here is truth-free.
"""
import numpy as np
import torch
from scipy.spatial import cKDTree

MODE_FIELDS = ("x", "y", "z", "score", "mass1", "mass3", "spread", "pri", "peak", "src", "sx", "sy", "sz", "snap_d")
SRC_VOTE, SRC_PEAK = 0, 1


def nms(xyz, score, radius, k=None):
    """Greedy NMS: indices of kept rows (highest score first)."""
    xyz = np.asarray(xyz, np.float64).reshape(-1, 3)
    order = np.argsort(-np.asarray(score, np.float64), kind="stable")
    keep = []
    for i in order:
        if keep and np.min(np.linalg.norm(xyz[keep] - xyz[i], axis=1)) < radius:
            continue
        keep.append(int(i))
        if k is not None and len(keep) >= k:
            break
    return np.asarray(keep, np.int64)


def mean_shift(P, w, seeds, bandwidth, iters=15, tol=1e-4):
    """Weighted Gaussian-kernel mean shift of `seeds` over points P with weights w (numpy, fp64)."""
    P = np.asarray(P, np.float64)
    w = np.asarray(w, np.float64)
    x = np.asarray(seeds, np.float64).copy()
    if len(P) == 0 or len(x) == 0:
        return x
    tree = cKDTree(P)
    r = 3.0 * bandwidth
    for _ in range(iters):
        moved = 0.0
        for i in range(len(x)):
            idx = tree.query_ball_point(x[i], r)
            if not idx:
                continue
            idx = np.asarray(idx)
            d2 = ((P[idx] - x[i]) ** 2).sum(1)
            k = w[idx] * np.exp(-0.5 * d2 / bandwidth ** 2)
            s = k.sum()
            if s <= 0:
                continue
            nx = (k[:, None] * P[idx]).sum(0) / s
            moved = max(moved, float(np.linalg.norm(nx - x[i])))
            x[i] = nx
        if moved < tol:
            break
    return x


def event_views(out, e, icfg, R=None):
    """Numpy view record of event e of a forward pass.  R: 3x3 map of the view (x_view = R x_reco); positions
    are mapped back with R^T (mirrors: R = R^T)."""
    m = (out["batch"] == e)
    coord = out["coord"][m].double().cpu().numpy()
    heat = torch.sigmoid(out["heat"][m].float()).cpu().numpy().astype(np.float64)
    p_any = torch.sigmoid(out["any"][m].float()).cpu().numpy().astype(np.float64)
    if "offset" in out:
        off = out["offset"][m].double().cpu().numpy()
        scale = np.exp(out["log_scale"][m].double().cpu().numpy())
    else:
        off = np.zeros_like(coord)
        scale = np.ones_like(coord)
    votes = coord + off
    ok = p_any >= icfg["vote_min_p"]
    idx = np.flatnonzero(ok)
    if len(idx) > icfg["n_vote"]:
        idx = idx[np.argsort(-p_any[idx], kind="stable")[:icfg["n_vote"]]]
    rec = {"votes": votes[idx], "w": p_any[idx], "pw": heat[idx], "vscale": scale[idx].mean(1),
           "vox": coord, "heat": heat, "any": p_any}
    if "cluster_logit" in out:
        cm = out["cluster_batch"] == e
        rec["clusters"] = dict(zip(out["cluster_id"][cm].cpu().numpy().tolist(),
                                   torch.sigmoid(out["cluster_logit"][cm].float()).cpu().numpy().tolist()))
    if "event_logit" in out:
        rec["event_p"] = float(torch.sigmoid(out["event_logit"][e].float()))
    if "point_logit" in out:
        pm = out["pt_batch"] == e
        rec["pt_xyz"] = out["pt_xyz"][pm].double().cpu().numpy()
        rec["pt_p"] = torch.sigmoid(out["point_logit"][pm].float()).cpu().numpy().astype(np.float64)
    if "lw" in out and "cluster_lw" in out:
        # C3 hierarchical score: log P(cluster) + log P(voxel | cluster), per voxel
        cm = out["cluster_batch"] == e
        cl = out["cluster_lw"][cm].float()
        lpc = torch.log_softmax(torch.cat([cl, cl.new_zeros(1)]), 0)[:-1]
        grp = out["vox_group"][m]
        loc = {int(g): i for i, g in enumerate(torch.nonzero(cm).squeeze(1).tolist())}
        gi = torch.tensor([loc[int(g)] for g in grp.tolist()], device=grp.device)
        s = out["lw"][m].float()
        lse = torch.zeros(len(cl), device=s.device).scatter_reduce(0, gi, s, "amax", include_self=False)
        z = torch.zeros(len(cl), device=s.device).index_add_(0, gi, torch.exp(s - lse[gi]))
        lpv = s - (lse[gi] + torch.log(z[gi]))
        rec["c3"] = (lpc[gi] + lpv).cpu().numpy().astype(np.float64)
    if "shell" in out:
        rec["shell"] = torch.sigmoid(out["shell"][m].float()).cpu().numpy()
        rec["cone"] = out["cone"][m].float().cpu().numpy()
    if R is not None:
        R = np.asarray(R, np.float64)
        for k in ("votes", "vox", "pt_xyz"):
            if k in rec:
                rec[k] = rec[k] @ R
        if "cone" in rec:
            rec["cone"] = np.concatenate([rec["cone"][:, :3] @ R, rec["cone"][:, 3:]], 1)
    return rec


def pool_views(views):
    """Pool the records of several views / ensemble members of one event: votes concatenated with weights / n;
    voxel fields averaged at the first view's voxels (other views mapped by nearest position): heat, any, the
    C2 shell (cumulative) probabilities and cone vectors (mean direction, mean log-kappa), the C3 hierarchical
    score (log of the mean probability); cluster and event probabilities averaged."""
    n = len(views)
    cat = lambda k: np.concatenate([v[k] for v in views], 0)          # noqa: E731
    out = {"votes": cat("votes"), "w": cat("w") / n, "pw": cat("pw"), "n_views": n}
    ref = views[0]
    out["vox"] = ref["vox"]
    maps = []                                          # per view: None = same voxels as ref, else nearest rows
    for v in views:
        if len(v["vox"]) == len(ref["vox"]) and np.allclose(v["vox"], ref["vox"]):
            maps.append(None)
        elif len(ref["vox"]) and len(v["vox"]):
            maps.append(cKDTree(v["vox"]).query(ref["vox"])[1])
        else:
            maps.append(False)

    def mean_field(key, fn=lambda a: a):
        acc, k = None, 0
        for v, j in zip(views, maps):
            if j is False or key not in v:
                continue
            a = fn(np.asarray(v[key], np.float64))
            a = a if j is None else a[j]
            acc = a.copy() if acc is None else acc + a
            k += 1
        return acc / max(k, 1) if acc is not None else None

    out["heat"] = mean_field("heat")
    out["any"] = mean_field("any")
    if out["heat"] is None:
        out["heat"], out["any"] = np.zeros(len(ref["vox"])), np.zeros(len(ref["vox"]))
    cl = {}
    for v in views:
        for c, p in v.get("clusters", {}).items():
            cl.setdefault(int(c), []).append(p)
    out["clusters"] = {c: float(np.mean(p)) for c, p in cl.items()}
    ep = [v["event_p"] for v in views if "event_p" in v]
    out["event_p"] = float(np.mean(ep)) if ep else float("nan")
    if "pt_xyz" in ref:
        out["pt_xyz"] = ref["pt_xyz"]
        P = np.zeros(len(ref["pt_xyz"]))
        for v in views:
            if len(v["pt_xyz"]) == len(ref["pt_xyz"]):
                P += v["pt_p"]
            else:
                d, j = cKDTree(v["pt_xyz"]).query(ref["pt_xyz"])
                P += v["pt_p"][j]
        out["pt_p"] = P / n
    if "c3" in ref:
        out["c3"] = np.log(np.maximum(mean_field("c3", np.exp), 1e-30))
    if "shell" in ref:
        out["shell"] = mean_field("shell")
        out["cone"] = mean_field("cone")
    return out


def find_modes(pv, icfg, points_xyz=None):
    """Modes of a pooled view record: (m, len(MODE_FIELDS)) float64 array, best first; vote modes then peaks."""
    rows = []
    P, w, pw = pv["votes"], pv["w"], pv["pw"]
    if len(P):
        seeds = P[nms(P, w, max(icfg["ms_bandwidth"], 0.5), icfg["n_seeds"])]
        ms = mean_shift(P, w, seeds, icfg["ms_bandwidth"], icfg["ms_iters"])
        tree = cKDTree(P)
        wp = w * pw
        tot = max(wp.sum(), 1e-12)
        feats = []
        for x in ms:
            i1 = np.asarray(tree.query_ball_point(x, 1.0), np.int64)
            i3 = np.asarray(tree.query_ball_point(x, 3.0), np.int64)
            m1 = w[i1].sum() if len(i1) else 0.0
            m3 = w[i3].sum() if len(i3) else 0.0
            sc = wp[i1].sum() / tot if len(i1) else 0.0
            pri = wp[i1].sum() / max(m1, 1e-12) if len(i1) else 0.0
            sp = np.sqrt((w[i3] * ((P[i3] - x) ** 2).sum(1)).sum() / max(m3, 1e-12)) if len(i3) else 0.0
            feats.append((sc, m1, m3, sp, pri))
        feats = np.asarray(feats, np.float64).reshape(-1, 5)
        keep = nms(ms, feats[:, 0], icfg["nms_cm"], icfg["top_k"])
        vt = cKDTree(pv["vox"]) if len(pv["vox"]) else None
        for i in keep:
            pk = 0.0
            if vt is not None:
                nb = vt.query_ball_point(ms[i], 1.0)
                pk = float(pv["heat"][nb].max()) if nb else 0.0
            rows.append([*ms[i], *feats[i], pk, SRC_VOTE])
    if len(pv["vox"]) and icfg.get("n_peaks", 0):
        for i in nms(pv["vox"], pv["heat"], max(icfg["nms_cm"], 1.0), icfg["n_peaks"]):
            h = float(pv["heat"][i])
            rows.append([*pv["vox"][i], h, 0.0, 0.0, 0.0, h, h, SRC_PEAK])
    M = np.asarray(rows, np.float64).reshape(-1, 10)
    snap = np.full((len(M), 4), np.nan)
    if len(M) and points_xyz is not None and len(points_xyz):
        d, j = cKDTree(np.asarray(points_xyz, np.float64)).query(M[:, :3])
        snap[:, :3] = np.asarray(points_xyz, np.float64)[j]
        snap[:, 3] = d
    return np.concatenate([M, snap], 1)


def bin_votes(votes, w, pw, bin_cm):
    """Sparse vote-density summary: (cells (c,3) centres, mass (c,), primary mass (c,)) on a bin_cm grid."""
    if len(votes) == 0:
        return np.zeros((0, 3), np.float32), np.zeros(0, np.float32), np.zeros(0, np.float32)
    g = np.floor(np.asarray(votes, np.float64) / bin_cm).astype(np.int64)
    u, inv = np.unique(g, axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    m = np.bincount(inv, w, len(u))
    mp = np.bincount(inv, w * pw, len(u))
    keep = m >= 1e-4 * max(m.sum(), 1e-12)
    return ((u[keep] + 0.5) * bin_cm).astype(np.float32), m[keep].astype(np.float32), mp[keep].astype(np.float32)


def compact_record(pv, modes, icfg):
    """What the inference store keeps per event and pass: modes, binned vote density, hot voxels / points,
    cluster P(host), event probability (and the C2 / C3 fields at hot voxels when present)."""
    hot = np.flatnonzero(np.maximum(pv["heat"], pv["any"]) >= icfg["hot_min_p"])
    if len(hot) > icfg["max_hot"]:
        hot = hot[np.argsort(-pv["heat"][hot], kind="stable")[:icfg["max_hot"]]]
    cells, mass, pmass = bin_votes(pv["votes"], pv["w"], pv["pw"], icfg["vote_bin_cm"])
    rec = {"modes": modes.astype(np.float64), "cells": cells, "cell_mass": mass, "cell_pmass": pmass,
           "hot_xyz": pv["vox"][hot].astype(np.float32), "hot_heat": pv["heat"][hot].astype(np.float32),
           "hot_any": pv["any"][hot].astype(np.float32), "p_host": pv.get("clusters", {}),
           "event_p": pv.get("event_p", float("nan")), "vote_total": float(pv["w"].sum())}
    if "pt_p" in pv:
        hp = np.flatnonzero(pv["pt_p"] >= icfg["hot_min_p"])
        if len(hp) > icfg["max_hot"]:
            hp = hp[np.argsort(-pv["pt_p"][hp], kind="stable")[:icfg["max_hot"]]]
        rec["hot_pt_xyz"] = pv["pt_xyz"][hp].astype(np.float32)
        rec["hot_pt_p"] = pv["pt_p"][hp].astype(np.float32)
    if "c3" in pv:
        rec["hot_c3"] = pv["c3"][hot].astype(np.float32)
    if "shell" in pv:
        rec["hot_shell"] = pv["shell"][hot].astype(np.float32)
        rec["hot_cone"] = pv["cone"][hot].astype(np.float32)
    return rec
