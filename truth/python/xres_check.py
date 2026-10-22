#!/usr/bin/env python3
"""xres_check.py -- is the drift coordinate of the reconstructed points consistent with the truth (pr_ml doc 21)?

The live BlobSampler places a point at x = xorig_f + dirx_f * (t_slice + reco_toff) * speed (clus/src/BlobSampler.cxx
time2drift) and the TruthPointLabeler places every drifted depo at x = xorig_f + dirx_f * (t_depo + depo_toff +
reco_toff') * speed' (Aux::time2drift with the labeler's reco_toff' / speed' read from the detector_volumes metadata).
Same xorig_f (first collection wire centre) and the same dirx_f, so a per-point residual along x is a TIME residual:
    dx = dirx * speed * (t_slice - t_depo - depo_toff)  (+ the (toff - toff') and (speed - speed') terms if they differ).
This script measures that residual on the R2 outputs without any toolkit code: nearest drifted depo of every live
point (3-D, in the labeler frame), the residual vector split into the component along the local true direction and the
perpendicular remainder, and the perpendicular remainder split by the track's angle to the drift axis.  A constant
offset shows up as a signed median of the perpendicular x residual on tracks ACROSS the drift (|dir_x| small) and is
invisible on tracks along it; a speed mismatch shows up as a slope of the x residual versus drift distance.

Inputs: the pctree dump of a clustering job with truth_point_label=true (the "3d" cloud carries x, y, z, t and the
wpid of every point; the "truth" cloud the local direction dir_x/y/z, tid, instance_valid, d_traj), the drifted depos
that job labelled with (sim-depos.tar.bz2: data (N,7) t q x y z L T, info (N,4) id pdg gen child), the wires file for
xorig_f and the face extents, and the drift constants of the job (reco_toff / speed of the sampler; the labeler's when
they differ).

  python3 xres_check.py --det pdhd --event cosmic0 --pctree .../clus_fold2/pctree-evt.tar.gz --depos .../sim-depos.tar.bz2
      --wires protodunehd-wires-larsoft-v1.json.bz2 --depo-toff-us 314 --reco-toff-us 0 --speed-mmus 1.576 --out DIR
Output: DIR/<event>.json (per face + per direction class), DIR/<event>.md row, DIR/<event>_*.png.
"""
import argparse, bz2, io, json, os, sys, tarfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_3d import load_tensors, PCTree  # noqa: E402

US = 1000.0   # WCT internal units: ns = 1, so us = 1000; mm = 1


def unwrap(x):
    return list(x.values())[0] if isinstance(x, dict) and len(x) == 1 else x


def load_faces(wires_file):
    """{(anode_ident, face_index): dict(xorig, ymin, ymax, zmin, zmax)} from a WCT wires file (face index = position
    of the face in the anode's list, as WirePlaneId::face() reports it)."""
    store = json.load(bz2.open(wires_file))["Store"]
    out = {}
    for a in store["anodes"]:
        a = unwrap(a)
        for fidx, fi in enumerate(a["faces"]):
            f = unwrap(store["faces"][fi])
            ys, zs = [], []
            xorig = None
            for pidx, pi in enumerate(f["planes"]):
                p = unwrap(store["planes"][pi])
                for wi in (p["wires"][0], p["wires"][-1]):
                    w = unwrap(store["wires"][wi])
                    for q in (unwrap(store["points"][w["tail"]]), unwrap(store["points"][w["head"]])):
                        ys.append(q["y"]); zs.append(q["z"])
                if pidx == 2:
                    w = unwrap(store["wires"][p["wires"][0]])
                    t, h = unwrap(store["points"][w["tail"]]), unwrap(store["points"][w["head"]])
                    xorig = 0.5 * (t["x"] + h["x"])
            out[(a["ident"], fidx)] = dict(xorig=xorig, ymin=min(ys), ymax=max(ys), zmin=min(zs), zmax=max(zs))
    return out


def load_depos(fn):
    D = I = None
    with tarfile.open(fn) as t:
        for m in t.getmembers():
            if not m.isfile():
                continue
            a = np.load(io.BytesIO(t.extractfile(m).read()), allow_pickle=False)
            if "depo_data" in m.name:
                D = a
            elif "depo_info" in m.name:
                I = a
    return D.astype(float), I


def decode_wpid(wpid):
    """WirePlaneId packing (iface/src/WirePlaneId.cxx): layer mask bits 0-2, face bit 3, apa from bit 4."""
    wpid = wpid.astype(int)
    return wpid >> 4, (wpid >> 3) & 1


def assign_depo_faces(D, faces, margin=30.0):
    """face key of every depo: the face whose (y, z) extent contains it and whose collection plane is nearest in x
    (the drifted depo sits on its face's response plane, a few cm on the drift side of the wires)."""
    keys = list(faces)
    best = np.full(len(D), -1, dtype=int)
    bestd = np.full(len(D), np.inf)
    for ki, k in enumerate(keys):
        f = faces[k]
        m = (D[:, 3] >= f["ymin"] - margin) & (D[:, 3] <= f["ymax"] + margin) & (D[:, 4] >= f["zmin"] - margin) & (D[:, 4] <= f["zmax"] + margin)
        d = np.abs(D[:, 2] - f["xorig"])
        upd = m & (d < bestd)
        best[upd] = ki; bestd[upd] = d[upd]
    return keys, best


def robust_fit(x, y):
    """median-of-slopes line y = a + b x (Theil-Sen on a subsample), returns a, b."""
    n = len(x)
    if n < 20:
        return float("nan"), float("nan")
    rng = np.random.default_rng(1)
    idx = rng.choice(n, size=min(n, 4000), replace=False)
    xs, ys = x[idx], y[idx]
    i = rng.integers(0, len(xs), size=20000); j = rng.integers(0, len(xs), size=20000)
    ok = np.abs(xs[i] - xs[j]) > 50.0
    if ok.sum() < 100:
        return float(np.median(ys)), 0.0
    b = np.median((ys[i][ok] - ys[j][ok]) / (xs[i][ok] - xs[j][ok]))
    a = np.median(ys - b * xs)
    return float(a), float(b)


def stats(v):
    v = np.asarray(v, dtype=float)
    if len(v) == 0:
        return {"n": 0}
    return {"n": int(len(v)), "median": float(np.median(v)), "mean": float(np.mean(v)), "p16": float(np.percentile(v, 16)),
            "p84": float(np.percentile(v, 84)), "rms": float(np.sqrt(np.mean(v * v)))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--det", required=True)
    ap.add_argument("--event", required=True)
    ap.add_argument("--pctree", required=True)
    ap.add_argument("--depos", required=True)
    ap.add_argument("--wires", required=True)
    ap.add_argument("--depo-toff-us", type=float, required=True, help="truth_time_offset of the labeler (us)")
    ap.add_argument("--reco-toff-us", type=float, default=0.0, help="BlobSampler time_offset (us)")
    ap.add_argument("--speed-mmus", type=float, required=True, help="BlobSampler drift_speed (mm/us)")
    ap.add_argument("--label-toff-us", type=float, default=None, help="labeler reco_time_offset if it differs (us)")
    ap.add_argument("--label-speed-mmus", type=float, default=None, help="labeler drift_speed if it differs (mm/us)")
    ap.add_argument("--pitch-mm", type=float, required=True)
    ap.add_argument("--max-match-mm", type=float, default=60.0, help="points farther than this from any depo are unmatched")
    ap.add_argument("--out", required=True)
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    label_toff = a.reco_toff_us if a.label_toff_us is None else a.label_toff_us
    label_speed = a.speed_mmus if a.label_speed_mmus is None else a.label_speed_mmus

    faces = load_faces(a.wires)
    tens = load_tensors(a.pctree)
    tops = sorted(dp for dp, (m, _) in tens.items() if m.get("datatype") == "pctree" and dp.endswith("/live"))
    P = {k: [] for k in ("x", "y", "z", "t", "wpid")}
    T = {k: [] for k in ("tid", "instance_valid", "dir_x", "dir_y", "dir_z", "d_traj", "cls")}
    for dp in tops:
        try:
            tree = PCTree(tens, dp)
        except RuntimeError:
            continue
        if "3d" not in tree.arrays or "truth" not in tree.arrays:
            continue
        p3, tr = tree.arrays["3d"], tree.arrays["truth"]
        if len(p3["x"]) != len(tr["tid"]):
            print("WARNING %s: 3d %d points, truth %d" % (dp, len(p3["x"]), len(tr["tid"])))
            continue
        for k in P:
            P[k].append(np.asarray(p3[k]))
        for k in T:
            T[k].append(np.asarray(tr[k]))
    res = {"det": a.det, "event": a.event, "faces": {}, "classes": {}, "inputs": vars(a)}
    if not P["x"]:
        res["error"] = "no live tree with 3d + truth clouds"
        json.dump(res, open(os.path.join(a.out, a.event + ".json"), "w"), indent=1)
        print(a.event, "no points"); return
    P = {k: np.concatenate(v).astype(float) for k, v in P.items()}
    T = {k: np.concatenate(v).astype(float) for k, v in T.items()}
    apa, face = decode_wpid(P["wpid"])
    n = len(P["x"])
    D, I = load_depos(a.depos)
    # the archive also holds the undrifted priors (gen 1); only the drifted depos (gen 0) are in the labeler frame
    # (pr_ml doc 22: including the priors moved the doc 21 medians by <= 0.02 mm)
    act = I[:, 2] == 0
    D, I = D[act], I[act]
    keys, dface = assign_depo_faces(D, faces)
    from scipy.spatial import cKDTree

    # per-point outputs
    r = np.full((n, 3), np.nan)          # residual point - nearest depo (labeler frame)
    dt_us = np.full(n, np.nan)           # t_slice - t_depo - depo_toff (us)
    ddrift = np.full(n, np.nan)          # drift distance of the point from its wire plane (mm)
    fkey = np.full(n, -1, dtype=int)
    for ki, k in enumerate(keys):
        f = faces[k]
        sel = (apa == k[0]) & (face == k[1])
        if not sel.any():
            continue
        xo = f["xorig"]
        dirx = float(np.sign(np.median(P["x"][sel] - xo))) or 1.0
        # the sampler's own relation, fitted: x = A + B t  ->  B = dirx*speed, A = xo + dirx*toff*speed
        A, B = robust_fit(P["t"][sel] / US, P["x"][sel])
        fit_speed = abs(B); fit_toff = (A - xo) / B if B else float("nan")
        xpred = xo + dirx * (P["t"][sel] / US + a.reco_toff_us) * a.speed_mmus
        samp_dev = P["x"][sel] - xpred
        dsel = dface == ki
        fres = {"npoints": int(sel.sum()), "ndepos": int(dsel.sum()), "xorig_mm": xo, "dirx": dirx,
                "sampler_fit": {"speed_mmus": fit_speed, "toff_us": fit_toff, "max_abs_dev_mm_vs_config": float(np.max(np.abs(samp_dev)))}}
        if dsel.sum() < 10:
            fres["note"] = "no depos on this face"
            res["faces"]["a%df%d" % k] = fres; continue
        xd = xo + dirx * (D[dsel, 0] / US + a.depo_toff_us + label_toff) * label_speed
        pts = np.column_stack([xd, D[dsel, 3], D[dsel, 4]])
        tree = cKDTree(pts)
        q = np.column_stack([P["x"][sel], P["y"][sel], P["z"][sel]])
        dist, idx = tree.query(q, k=1, distance_upper_bound=a.max_match_mm)
        ok = np.isfinite(dist)
        ii = np.where(sel)[0]
        r[ii[ok]] = q[ok] - pts[idx[ok]]
        tdep = D[dsel, 0][idx[ok]] / US
        dt_us[ii[ok]] = P["t"][sel][ok] / US - tdep - a.depo_toff_us
        ddrift[ii[ok]] = dirx * (P["x"][sel][ok] - xo)
        fkey[ii] = ki
        rx = r[ii[ok], 0] * dirx          # positive = point farther from the anode than its depo
        fres.update({"matched_frac": float(ok.mean()), "dx_drift_mm": stats(rx), "dt_us": stats(dt_us[ii[ok]]),
                     "dy_mm": stats(r[ii[ok], 1]), "dz_mm": stats(r[ii[ok], 2])})
        aa, bb = robust_fit(ddrift[ii[ok]], rx)
        fres["dx_vs_drift"] = {"intercept_mm": aa, "slope": bb, "slope_pct": 100 * bb, "drift_range_mm": [float(np.min(ddrift[ii[ok]])), float(np.max(ddrift[ii[ok]]))]}
        res["faces"]["a%df%d" % k] = fres
    # direction classes (all faces): residual perpendicular to the local true direction
    okm = np.isfinite(r[:, 0])
    dirv = np.column_stack([T["dir_x"], T["dir_y"], T["dir_z"]])
    nd = np.linalg.norm(dirv, axis=1)
    hasdir = okm & (nd > 0.5) & (T["tid"] >= 0)
    dirn = np.zeros_like(dirv); dirn[nd > 0] = dirv[nd > 0] / nd[nd > 0, None]
    rpar = np.sum(r * dirn, axis=1)
    rperp = r - rpar[:, None] * dirn
    dirx_pt = np.array([faces[keys[k]]["xorig"] if k >= 0 else 0 for k in fkey])  # placeholder to build the sign
    sign = np.array([np.sign(P["x"][i] - faces[keys[fkey[i]]]["xorig"]) if fkey[i] >= 0 else 1.0 for i in range(n)])
    rperp_x = rperp[:, 0] * sign
    rperp_yz = np.hypot(rperp[:, 1], rperp[:, 2])
    ax = np.abs(dirn[:, 0])
    classes = {"across_drift (|dir_x|<0.3)": hasdir & (ax < 0.3), "oblique (0.3-0.8)": hasdir & (ax >= 0.3) & (ax < 0.8),
               "along_drift (|dir_x|>=0.8)": hasdir & (ax >= 0.8), "all_with_dir": hasdir}
    for name, m in classes.items():
        res["classes"][name] = {"n": int(m.sum()), "rperp_x_mm": stats(rperp_x[m]), "rperp_yz_mm": stats(rperp_yz[m]),
                                "rperp_mm": stats(np.linalg.norm(rperp[m], axis=1)), "d_traj_mm": stats(T["d_traj"][m]),
                                "dx_drift_mm": stats((r[:, 0] * sign)[m]), "dt_us": stats(dt_us[m])}
    res["summary"] = {"npoints": int(n), "matched_frac": float(okm.mean()), "with_dir_frac": float(hasdir.mean()),
                      "dx_drift_mm_all": stats((r[:, 0] * sign)[okm]), "dt_us_all": stats(dt_us[okm]),
                      "pitch_mm": a.pitch_mm, "label_toff_us": label_toff, "label_speed_mmus": label_speed}
    json.dump(res, open(os.path.join(a.out, a.event + ".json"), "w"), indent=1)
    c = res["classes"]
    def fmt(s, k="median"):
        return "%.2f" % s[k] if s.get("n", 0) else "-"
    row = "| %s | %s | %d | %.2f | %s (%d) | %s (%d) | %s (%d) | %s | %s | %s |" % (
        a.det, a.event, n, okm.mean(),
        fmt(c["across_drift (|dir_x|<0.3)"]["rperp_x_mm"]), c["across_drift (|dir_x|<0.3)"]["n"],
        fmt(c["oblique (0.3-0.8)"]["rperp_x_mm"]), c["oblique (0.3-0.8)"]["n"],
        fmt(c["along_drift (|dir_x|>=0.8)"]["rperp_yz_mm"]), c["along_drift (|dir_x|>=0.8)"]["n"],
        fmt(res["summary"]["dt_us_all"]), fmt(res["summary"]["dx_drift_mm_all"], "p16") + " / " + fmt(res["summary"]["dx_drift_mm_all"], "p84"),
        "; ".join("%s %.1f%%" % (k, v["dx_vs_drift"]["slope_pct"]) for k, v in res["faces"].items() if "dx_vs_drift" in v and v["npoints"] > 200))
    open(os.path.join(a.out, a.event + ".md"), "w").write(row + "\n")
    print(row)
    if a.no_plots:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    lim = 6 * a.pitch_mm
    for name, m in classes.items():
        if name == "all_with_dir" or m.sum() < 5:
            continue
        axs[0].hist(np.clip(rperp_x[m], -lim, lim), bins=80, range=(-lim, lim), histtype="step", label="%s n=%d" % (name, m.sum()))
        axs[1].hist(np.clip(rperp_yz[m], 0, lim), bins=60, range=(0, lim), histtype="step", label=name)
    axs[0].set_xlabel("perpendicular residual, drift component (mm; + = away from anode)"); axs[0].legend(fontsize=7); axs[0].axvline(0, color="k", lw=0.5)
    axs[1].set_xlabel("perpendicular residual, transverse (yz) magnitude (mm)"); axs[1].legend(fontsize=7)
    mm = okm & np.isfinite(ddrift)
    if mm.sum() > 10:
        axs[2].hist2d(ddrift[mm], np.clip((r[:, 0] * sign)[mm], -lim, lim), bins=[60, 60], cmap="viridis")
    axs[2].set_xlabel("drift distance of the point (mm)"); axs[2].set_ylabel("x residual point - nearest depo (mm)")
    fig.suptitle("%s %s: live points vs drifted depos (labeler frame)" % (a.det, a.event))
    fig.tight_layout(); fig.savefig(os.path.join(a.out, a.event + "_xres.png"), dpi=110); plt.close(fig)


if __name__ == "__main__":
    main()
