#!/usr/bin/env python3
"""validate_2d.py -- E0.5: signal processing against the true depos in 2-D (pr_ml doc 19 sec 6).

Per anode of an event, with the SP frame archive (FrameFileSink: frame_<tag>_<ident>.npy, channels_, tickinfo_) and the
DRIFTED depo file of the same Stage B job (DepoFileSink: depo_data_0/depo_info_0, the priors carry the G4 identity; in
Stage A index mode the depo id is the SimEnergyDeposit index and truth.npz's sed_ids maps it to the track id), builds the
NATIVE truth image of every plane: each depo projected onto (channel, tick) with its charge, keyed by track id, so every
pixel holds q_true, tid_argmax and n_tid.  The projection uses the wires file (wire direction and pitch per plane, wire
-> channel), the frame's tickinfo (t0, tick) and the imaging time offset (the BlobDepoFill time_offset of the detector's
truth tier; the residual timing lag measured here is the sec 6.3 self-check of that offset).  Charge is spread with the
depo's longitudinal / transverse extents (DepoBlobAssoc's Gaussians, nsigma) over ticks and wires when --smear is on.

Metrics per plane (from fdvd_sim/stageB/validate_truth.py, extended per track and class):
  closure = sum(reco) / sum(truth) over the frame and inside the dilated truth region (+-PAD ticks, +-1 channel);
  truth segments (per channel, contiguous ticks with truth > THR): efficiency (any reco charge in [start-PAD, end+PAD]),
  reco/true charge vs segment charge, per class of the segment's top track (track / shower / other by pdg);
  fake ROIs: reco runs with no truth in +-FPAD ticks, +-1 channel; timing: reco - truth charge centroid (ticks) for
  segments with Q > QTIME (median must be 0 +- 1 tick when the offset is right); ambiguity floor: fraction of reco
  charge on pixels with n_tid >= 2 (bounds any 3-D label purity); per-track efficiency (tracks with >= MINQ true charge).
usage:
  validate_2d.py --out <dir> --det <name> --wires <wires.json.bz2> --time-offset-us <us> [--tags gauss,dnnsp]
      [--smear] [--nsigma 3] [--speed-mmus <drift speed>] --event <tag> --frames <archive>[,<archive>...] --depos <sim-depos>
      [--truth truth.npz] [--anodes 0,1,...] [--event ... (repeat the group per event)]
"""
import argparse, bz2, io, json, os, re, sys, tarfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import truth_table as T  # noqa: E402

THR = 10.0       # e per pixel: truth threshold for segments
PAD = 10         # ticks around a truth segment
FPAD = 20        # ticks of truth-free margin for a fake ROI
QTIME = 5000.0   # e: minimum segment charge for the timing measurement
MINQ = 2e4       # e: minimum true charge of a track for the per-track efficiency
TRACK_PDGS = {13, 211, 321, 2212, 3222, 3112}


def unwrap(lst):
    return lst


def load_wires(fn):
    """{anode ident: {face index: {plane index: dict(chan, wire_pos (N,3) centers, pitch_dir, pitch, wdir, x)}}}"""
    d = json.load(bz2.open(fn))["Store"]
    anodes, faces, planes, wires = d["anodes"], d["faces"], d["planes"], d["wires"]
    pts = np.array([[p["Point"]["x"], p["Point"]["y"], p["Point"]["z"]] for p in d["points"]])
    out = {}
    for a in anodes:
        A = a["Anode"]
        out[A["ident"]] = {}
        for fi, fidx in enumerate(A["faces"]):
            F = faces[fidx]["Face"]
            out[A["ident"]][fi] = {}
            for pi, pidx in enumerate(F["planes"]):
                P = planes[pidx]["Plane"]
                ws = [wires[w]["Wire"] for w in P["wires"]]
                tails = pts[[w["tail"] for w in ws]]
                heads = pts[[w["head"] for w in ws]]
                cen = 0.5 * (tails + heads)
                wdir = heads - tails
                wdir = (wdir.T / np.linalg.norm(wdir, axis=1)).T
                # pitch direction in the y-z plane, perpendicular to the wire
                pdir = np.stack([np.zeros(len(ws)), -wdir[:, 2], wdir[:, 1]], axis=1)
                pdir = (pdir.T / np.linalg.norm(pdir, axis=1)).T
                chan = np.array([w["channel"] for w in ws])
                if len(ws) > 1:
                    s = np.einsum("ij,j->i", cen - cen[0], pdir[0])
                    order = np.argsort(s, kind="stable")
                    s, chan, cen = s[order], chan[order], cen[order]
                    pitch = float(np.median(np.diff(s)))
                else:
                    pitch, s = 1.0, np.zeros(1)
                # wire index along the pitch direction = round((p - cen0).pdir / pitch): cen0 is the first wire
                out[A["ident"]][fi][pi] = dict(chan=chan, s=s, pdir=pdir[0], pitch=pitch,
                                               cen0=cen[0], x=float(np.median(cen[:, 0])),
                                               ymin=float(min(tails[:, 1].min(), heads[:, 1].min())), ymax=float(max(tails[:, 1].max(), heads[:, 1].max())),
                                               zmin=float(min(tails[:, 2].min(), heads[:, 2].min())), zmax=float(max(tails[:, 2].max(), heads[:, 2].max())))
    return out


def load_frames(fn, tags):
    out = {}
    with tarfile.open(fn) as tf:
        mem = {m.name: m for m in tf.getmembers()}
        for n, m in mem.items():
            mm = re.match(r"frame_([A-Za-z]+?)(\d*)_(\d+)\.npy", n)
            if not mm:
                continue
            tag, ev = mm.group(1), mm.group(3)
            ident = int(mm.group(2)) if mm.group(2) else -1   # -1: one frame for every anode (SBND dnnsp_<event>)
            if tags and tag not in tags:
                continue
            arr = np.load(io.BytesIO(tf.extractfile(m).read()), allow_pickle=False).astype(np.float64)
            ch = np.load(io.BytesIO(tf.extractfile(mem[n.replace("frame_", "channels_")]).read()), allow_pickle=False)
            ti = np.load(io.BytesIO(tf.extractfile(mem[n.replace("frame_", "tickinfo_")]).read()), allow_pickle=False)
            out[(tag, ident)] = (arr, ch, ti)
    return out


def load_depos(fn):
    with tarfile.open(fn) as tf:
        mem = {m.name: m for m in tf.getmembers()}
        data = np.load(io.BytesIO(tf.extractfile(mem["depo_data_0.npy"]).read()), allow_pickle=False).astype(np.float64)
        info = np.load(io.BytesIO(tf.extractfile(mem["depo_info_0.npy"]).read()), allow_pickle=False).astype(np.int64)
    # active depos (gen 0) with the key id of their prior chain root (legacy files: child == own row => prior = next row)
    gen = info[:, 2]
    nprior = int((gen > 0).sum())
    legacy = nprior > 0 and bool((info[gen > 0, 3] == np.arange(len(info))[gen > 0]).all())
    key_id = info[:, 0].copy()
    if nprior:
        child = info[:, 3].copy()
        if legacy:
            child[gen > 0] = np.arange(len(info))[gen > 0] - 1
        # walk: for each prior row, assign its id to its child, repeat so deeper priors win
        for g in sorted(np.unique(gen[gen > 0]), reverse=True):
            rows = np.where(gen == g)[0]
            key_id[child[rows]] = key_id[rows]
    act = gen == 0
    return data[act], info[act], key_id[act]


def runs(mask_row):
    d = np.diff(np.concatenate(([0], mask_row.view(np.int8), [0])))
    return list(zip(np.where(d == 1)[0], np.where(d == -1)[0]))


def assign_faces(data, geo, xmax=500.0):
    """(anode, face) of every drifted depo: the face whose wire extent holds its (y, z) and whose plane x is nearest
    (the drifted depos sit on the response plane of the face they drifted to); -1 when none within xmax."""
    best_d = np.full(len(data), np.inf); best_a = np.full(len(data), -1); best_f = np.full(len(data), -1)
    for a, faces in geo.items():
        for fi, fg in faces.items():
            g0 = fg[0]
            inside = (data[:, 3] >= g0["ymin"] - 50) & (data[:, 3] <= g0["ymax"] + 50) & (data[:, 4] >= g0["zmin"] - 50) & (data[:, 4] <= g0["zmax"] + 50)
            dx = np.abs(data[:, 2] - g0["x"])
            better = inside & (dx < xmax) & (dx < best_d)
            best_d[better] = dx[better]; best_a[better] = a; best_f[better] = fi
    return best_a, best_f


def project(data, key, face_geo, t0, tick, nticks, toff, speed, smear, nsigma, sel):
    """per plane: (q_true grid (nwires, nticks), tid_argmax, n_tid, per-track q) of the depos assigned to the face"""
    d, k = data[sel], key[sel]
    t = d[:, 0] + toff
    tk = (t - t0) / tick
    out = {}
    for pi, g in face_geo.items():
        s = np.einsum("ij,j->i", d[:, 2:5] - g["cen0"], g["pdir"]) / g["pitch"]
        w = np.rint(s).astype(int)
        nw = len(g["chan"])
        q = np.abs(d[:, 1])
        if smear:
            sig_w = np.maximum(d[:, 6] / g["pitch"], 1e-3)
            sig_t = np.maximum(d[:, 5] / speed / tick, 1e-3)
            # spread each depo over a grid of +-nsigma in wires and ticks (truncated at 5x5 cells for cost)
            dw = np.arange(-2, 3); dt = np.arange(-2, 3)
            acc = {}
            W = {}
            for i_w in dw:
                ww = w + i_w
                fw = np.exp(-0.5 * ((ww - s) / sig_w) ** 2)
                for i_t in dt:
                    tt = np.rint(tk).astype(int) + i_t
                    ft = np.exp(-0.5 * ((tt - tk) / sig_t) ** 2)
                    W[(i_w, i_t)] = (ww, tt, fw * ft)
            norm = sum(v[2] for v in W.values())
            cols, rows, vals, keys = [], [], [], []
            for (i_w, i_t), (ww, tt, f) in W.items():
                ok = (ww >= 0) & (ww < nw) & (tt >= 0) & (tt < nticks)
                cols.append(ww[ok]); rows.append(tt[ok]); vals.append(q[ok] * f[ok] / norm[ok]); keys.append(k[ok])
            ww = np.concatenate(cols); tt = np.concatenate(rows); qq = np.concatenate(vals); kk = np.concatenate(keys)
        else:
            tt = np.rint(tk).astype(int)
            ok = (w >= 0) & (w < nw) & (tt >= 0) & (tt < nticks)
            ww, tt, qq, kk = w[ok], tt[ok], q[ok], k[ok]
        grid = np.zeros((nw, nticks))
        if len(ww) == 0:
            out[pi] = (grid, np.full((nw, nticks), -1, dtype=np.int64), np.zeros((nw, nticks), dtype=np.int32), {})
            continue
        np.add.at(grid, (ww, tt), qq)
        # per (pixel, tid) sums -> argmax and multiplicity
        pix = ww.astype(np.int64) * nticks + tt
        order = np.lexsort((kk, pix))
        pix_o, k_o, q_o = pix[order], kk[order], qq[order]
        grp = np.concatenate(([True], (pix_o[1:] != pix_o[:-1]) | (k_o[1:] != k_o[:-1])))
        gidx = np.cumsum(grp) - 1
        qsum = np.bincount(gidx, weights=q_o)
        gp = pix_o[grp]; gk = k_o[grp]
        ntid = np.zeros(nw * nticks, dtype=np.int32)
        np.add.at(ntid, gp, 1)
        best = np.full(nw * nticks, -1, dtype=np.int64); bestq = np.zeros(nw * nticks)
        for p_, kk_, qq_ in zip(gp, gk, qsum):
            if qq_ > bestq[p_]:
                bestq[p_] = qq_; best[p_] = kk_
        qtrack = {}
        for kk_, qq_ in zip(gk, qsum):
            qtrack[int(kk_)] = qtrack.get(int(kk_), 0.0) + float(qq_)
        out[pi] = (grid, best.reshape(nw, nticks), ntid.reshape(nw, nticks), qtrack)
    return out


def analyse(frames, geo, data, key, tid_of_key, cls_of_tid, anode, toff, speed, smear, nsigma, tags, depo_anode, depo_face):
    res = {}
    for (tag, ident), (arr, ch, ti) in frames.items():
        if ident != anode and ident != -1:
            continue
        # The imaging slices and BlobDepoFill measure time from the frame's first tick (tbin0), not from the frame's
        # absolute time stamp ti[0]: the detector's scanned time_offset already contains the readout-window start
        # (PDHD: 314 us = 250 us pre-trigger window + 64 us response-plane transit).  So tick = (t + toff)/tick - tbin0.
        frame_t0, tick, tbin0 = float(ti[0]), float(ti[1]), float(ti[2])
        t0 = tbin0 * tick
        nticks = arr.shape[1]
        chan_row = {int(c): i for i, c in enumerate(ch)}
        nrow = len(ch)
        planes = {}
        for fi, face_geo in geo[anode].items():
            sel = (depo_anode == anode) & (depo_face == fi)
            if not sel.any():
                continue
            P = project(data, key, face_geo, t0, tick, nticks, toff, speed, smear, nsigma, sel)
            for pi, (grid, best, ntid, qtrack) in P.items():
                g = face_geo[pi]
                rows = np.array([chan_row.get(int(c), -1) for c in g["chan"]])
                ok = rows >= 0
                if not ok.any():
                    continue
                if pi not in planes:
                    planes[pi] = dict(tru=np.zeros((nrow, nticks)), bst=np.full((nrow, nticks), -1, dtype=np.int64), bstq=np.zeros((nrow, nticks)),
                                      nt=np.zeros((nrow, nticks), dtype=np.int32), rows=set(), qtrack={}, faces=[])
                Pl = planes[pi]
                Pl["tru"][rows[ok]] += grid[ok]
                Pl["nt"][rows[ok]] += ntid[ok]
                # the dominant track per pixel across faces: keep the larger contribution
                upd = grid[ok] > Pl["bstq"][rows[ok]]
                sub = Pl["bst"][rows[ok]]; sub[upd] = best[ok][upd]; Pl["bst"][rows[ok]] = sub
                subq = Pl["bstq"][rows[ok]]; subq[upd] = grid[ok][upd]; Pl["bstq"][rows[ok]] = subq
                Pl["rows"].update(int(r) for r in rows[ok]); Pl["faces"].append(fi)
                for kk_, qq_ in qtrack.items():
                    Pl["qtrack"][kk_] = Pl["qtrack"].get(kk_, 0.0) + qq_
        for pi, Pl in planes.items():
            if True:
                rows_ = np.array(sorted(Pl["rows"]))
                reco = arr[rows_]
                tru = Pl["tru"][rows_]; bst = Pl["bst"][rows_]; nt = Pl["nt"][rows_]; qtrack = Pl["qtrack"]
                fi = "+".join(str(f) for f in Pl["faces"])
                A = {"tag": tag, "face": fi, "plane": pi, "nch": int(len(rows_)), "frame_t0_us": frame_t0 / 1000.0}
                A["q_true"] = float(tru.sum()); A["q_reco"] = float(reco.sum())
                # whole-plane time-profile cross-correlation: the lag (ticks) of truth onto reco, 0 when the offset is right
                pt, pr = tru.sum(axis=0), np.clip(reco, 0, None).sum(axis=0)
                if pt.sum() > 0 and pr.sum() > 0:
                    best = (0, -2.0)
                    for lag in range(-400, 401, 2):
                        sh = np.roll(pt, lag)
                        if lag > 0: sh[:lag] = 0
                        elif lag < 0: sh[lag:] = 0
                        c = np.corrcoef(sh, pr)[0, 1]
                        if c > best[1]: best = (lag, c)
                    A["profile_lag_ticks"] = best[0]; A["profile_corr"] = float(best[1])
                else:
                    A["profile_lag_ticks"] = None; A["profile_corr"] = None
                A["closure"] = float(reco.sum() / tru.sum()) if tru.sum() > 0 else None
                tmask = tru > THR
                reg = tmask.copy()
                for s_ in range(1, PAD + 1):
                    reg[:, s_:] |= tmask[:, :-s_]; reg[:, :-s_] |= tmask[:, s_:]
                reg2 = reg.copy(); reg2[1:] |= reg[:-1]; reg2[:-1] |= reg[1:]
                A["closure_region"] = float(reco[reg2].sum() / tru[reg2].sum()) if tru[reg2].sum() > 0 else None
                # ambiguity floor
                pos = np.clip(reco, 0, None)
                A["reco_frac_on_multi_tid"] = float(pos[(nt >= 2)].sum() / pos.sum()) if pos.sum() > 0 else None
                A["true_frac_on_multi_tid"] = float(tru[(nt >= 2)].sum() / tru.sum()) if tru.sum() > 0 else None
                # segments
                seg_q, seg_r, dts, cls_eff = [], [], [], {}
                for r in np.where(tmask.any(axis=1))[0]:
                    for s0, s1 in runs(tmask[r]):
                        q = tru[r, s0:s1].sum()
                        w0, w1 = max(0, s0 - PAD), min(nticks, s1 + PAD)
                        qr = reco[r, w0:w1].sum()
                        seg_q.append(q); seg_r.append(qr)
                        # class of the segment's dominant track
                        kk = bst[r, s0:s1]; kk = kk[kk >= 0]
                        if len(kk):
                            top = np.bincount(kk - kk.min()).argmax() + kk.min()
                            c = cls_of_tid(tid_of_key(int(top)))
                            e = cls_eff.setdefault(c, [0, 0, 0.0, 0.0]); e[0] += 1; e[1] += int(qr > 0); e[2] += q; e[3] += qr
                        if q > QTIME:
                            tt = np.arange(s0, s1); ct = (tru[r, s0:s1] * tt).sum() / q
                            ww = np.arange(w0, w1); rr = np.clip(reco[r, w0:w1], 0, None)
                            if rr.sum() > 0:
                                dts.append((rr * ww).sum() / rr.sum() - ct)
                seg_q = np.array(seg_q); seg_r = np.array(seg_r)
                A["nseg"] = int(len(seg_q)); A["seg_eff"] = float((seg_r > 0).mean()) if len(seg_q) else None
                A["seg_eff_q_gt_1000"] = float((seg_r[seg_q > 1000] > 0).mean()) if (seg_q > 1000).any() else None
                A["timing_median_ticks"] = float(np.median(dts)) if dts else None
                A["timing_n"] = len(dts)
                A["class_eff"] = {k: {"nseg": v[0], "eff": v[1] / v[0], "closure": v[3] / v[2] if v[2] > 0 else None} for k, v in cls_eff.items()}
                # fakes
                tfree = tru > 0
                tf2 = tfree.copy(); tf2[1:] |= tfree[:-1]; tf2[:-1] |= tfree[1:]
                cum = np.concatenate([np.zeros((len(tf2), 1)), np.cumsum(tf2, axis=1)], axis=1)
                nfake, qfake = 0, 0.0
                nz = reco != 0
                for r in np.where(nz.any(axis=1))[0]:
                    for s0, s1 in runs(nz[r]):
                        w0, w1 = max(0, s0 - FPAD), min(nticks, s1 + FPAD)
                        if cum[r, w1] - cum[r, w0] == 0:
                            nfake += 1; qfake += float(reco[r, s0:s1].sum())
                A["fake_rois"] = nfake; A["fake_rois_per_channel"] = nfake / max(ok.sum(), 1); A["fake_q_frac"] = qfake / float(reco.sum()) if reco.sum() > 0 else None
                # per-track efficiency: a track counts as found when the reco charge in its pixels' windows exceeds 10 %
                found = tot = 0
                for kk_, qq_ in qtrack.items():
                    if qq_ < MINQ:
                        continue
                    tot += 1
                    m = (bst == kk_) & tmask
                    if m.any():
                        rr = reco[m]
                        if rr.sum() > 0.1 * qq_:
                            found += 1
                A["ntracks_minq"] = tot; A["track_eff"] = found / tot if tot else None
                res["%s_f%s_p%d" % (tag, fi, pi)] = A
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True); ap.add_argument("--det", required=True); ap.add_argument("--wires", required=True)
    ap.add_argument("--time-offset-us", type=float, required=True); ap.add_argument("--tags", default="gauss,dnnsp")
    ap.add_argument("--smear", action="store_true"); ap.add_argument("--nsigma", type=float, default=3.0)
    ap.add_argument("--speed-mmus", type=float, default=1.6)
    ap.add_argument("--event", action="append", required=True); ap.add_argument("--frames", action="append", default=[])
    ap.add_argument("--depos", action="append", default=[]); ap.add_argument("--truth", action="append", default=[])
    ap.add_argument("--anodes", action="append", default=[])
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    geo = load_wires(a.wires)
    tags = set(a.tags.split(",")) if a.tags else None
    toff = a.time_offset_us * 1000.0
    per_event = {}
    for i, ev in enumerate(a.event):
        frames = {}
        for fn in a.frames[i].split(","):
            frames.update(load_frames(fn, tags))
        data, info, key = load_depos(a.depos[i])
        tab = T.read_truth_npz(a.truth[i]) if i < len(a.truth) and a.truth[i] else None
        if tab is not None and len(tab["sed_ids"]):
            sed = tab["sed_ids"][:, 0]
            def tid_of_key(k, sed=sed):
                return int(sed[k]) if 0 <= k < len(sed) else int(k)
            pf = tab["truth_pf"]; tid2pdg = dict(zip(pf[:, T.pf_col("trackid")].astype(int), pf[:, T.pf_col("pdg")].astype(int)))
        else:
            tid_of_key = lambda k: int(k)  # noqa: E731
            tid2pdg = dict(zip(info[:, 0].tolist(), info[:, 1].tolist()))
        def cls_of_tid(t):
            p = abs(tid2pdg.get(abs(t), 0))
            return "track" if (p in TRACK_PDGS or p > 1000000000) else ("shower" if p in (11, 22) else "other")
        anodes = [int(x) for x in a.anodes[i].split(",")] if i < len(a.anodes) and a.anodes[i] else sorted({ident for (_, ident) in frames if ident >= 0} or set(geo))
        depo_anode, depo_face = assign_faces(data, geo)
        r = {}
        for an in anodes:
            if an not in geo:
                continue
            r.update({"a%d_%s" % (an, k): v for k, v in analyse(frames, geo, data, key, tid_of_key, cls_of_tid, an, toff, a.speed_mmus / 1000.0 * 1000.0, a.smear, a.nsigma, tags, depo_anode, depo_face).items()})
        per_event[ev] = r
        for k, v in r.items():
            print(ev, k, "closure %.3f region %.3f seg_eff %s timing %s n=%s lag %s corr %s fake/ch %.3f multi_tid %s" % (
                v["closure"] or 0, v["closure_region"] or 0, None if v["seg_eff"] is None else round(v["seg_eff"], 3), None if v["timing_median_ticks"] is None else round(v["timing_median_ticks"], 2),
                v["timing_n"], v["profile_lag_ticks"], None if v["profile_corr"] is None else round(v["profile_corr"], 3), v["fake_rois_per_channel"], None if v["reco_frac_on_multi_tid"] is None else round(v["reco_frac_on_multi_tid"], 3)))
    json.dump(per_event, open(os.path.join(a.out, "per_event.json"), "w"), indent=1)
    # summary per (tag, plane) over events and anodes
    summ = {}
    for ev, r in per_event.items():
        for k, v in r.items():
            kk = "%s_p%d" % (v["tag"], v["plane"])
            s = summ.setdefault(kk, {"closure": [], "closure_region": [], "seg_eff": [], "timing": [], "lag": [], "fake_per_ch": [], "multi_tid": [], "track_eff": []})
            for src, dst in (("closure", "closure"), ("closure_region", "closure_region"), ("seg_eff", "seg_eff"), ("timing_median_ticks", "timing"), ("profile_lag_ticks", "lag"),
                             ("fake_rois_per_channel", "fake_per_ch"), ("reco_frac_on_multi_tid", "multi_tid"), ("track_eff", "track_eff")):
                if v.get(src) is not None:
                    s[dst].append(v[src])
    summary = {k: {m: ({"mean": float(np.mean(x)), "min": float(np.min(x)), "max": float(np.max(x)), "n": len(x)} if x else None) for m, x in v.items()} for k, v in summ.items()}
    summary["_args"] = vars(a)
    json.dump(summary, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    with open(os.path.join(a.out, "summary.md"), "w") as f:
        f.write("| tag_plane | closure | closure(region) | seg eff | timing (ticks) | profile lag (ticks) | fake ROI / ch | reco on n_tid>=2 | track eff | n |\n|---|---|---|---|---|---|---|---|---|---|\n")
        for k, v in summary.items():
            if k.startswith("_"):
                continue
            cell = lambda m: ("%.3f" % v[m]["mean"]) if v.get(m) else "-"  # noqa: E731
            f.write("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n" % (k, cell("closure"), cell("closure_region"), cell("seg_eff"), cell("timing"), cell("lag"), cell("fake_per_ch"), cell("multi_tid"), cell("track_eff"), v["closure"]["n"] if v.get("closure") else 0))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
