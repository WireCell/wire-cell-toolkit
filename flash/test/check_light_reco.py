#!/usr/bin/env python3
"""Stage-3 validation of the WCT-native PDHD light reconstruction
against the LArSoft reference products carried in the same work dir
(see flash/docs/validation.md).

Compares, per event:
 1. OpDecon "decon" frames vs the in-file LArSoft "deconv" frames
    (per-channel correlation, peak ratio, peak shift);
 2. OpHitFinder hits vs converted PerOpHitTree hits (efficiency,
    surplus rate, PE ratio of matched hits);
 3. OpFlashFinder flashes vs converted LArSoft flashes (time-matched
    pairs: dt, PE ratio on common channels, dy/dz).

Writes PNG plots to <workdir>/light_validation/ and prints a summary.

Usage:
  python3 check_light_reco.py <workdir> [<workdir> ...]
"""
import io
import os
import sys
import tarfile

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TICK = 16.0  # ns


def load_npy(path):
    out = {}
    with tarfile.open(path) as tf:
        for m in tf.getmembers():
            if m.name.endswith(".npy"):
                out[m.name] = np.load(io.BytesIO(tf.extractfile(m).read()))
    return out


def get_tensor(path, idx):
    with tarfile.open(path) as tf:
        names = sorted(m.name for m in tf.getmembers() if m.name.endswith("_array.npy"))
        return np.load(io.BytesIO(tf.extractfile(names[idx]).read()))


def snippet_starts(raw_row):
    nz = np.nonzero(raw_row)[0]
    if not len(nz):
        return []
    starts = [nz[0]]
    for j in nz[1:]:
        if j >= starts[-1] + 1024:
            starts.append(j)
    return starts


def compare_decon(wct, ref, outdir):
    evt = [k for k in wct if k.startswith("frame_decon_")][0].split("_")[-1][:-4]
    o_all, r_all = wct[f"frame_decon_{evt}.npy"], ref[f"frame_deconv_{evt}.npy"]
    raw = ref[f"frame_raw_{evt}.npy"]
    chans = ref[f"channels_deconv_{evt}.npy"]
    corrs, ratios, shifts = [], [], []
    for i in range(len(chans)):
        for s in snippet_starts(raw[i]):
            r, o = r_all[i, s:s + 1024], o_all[i, s:s + 1024]
            if np.abs(r).max() < 1.0:
                continue
            cc = [np.dot(r, np.roll(o, k)) for k in range(-8, 9)]
            k = int(np.argmax(cc)) - 8
            osh = np.roll(o, k)
            corrs.append(np.dot(r, osh) / (np.linalg.norm(r) * np.linalg.norm(osh) + 1e-30))
            ratios.append(osh.max() / r.max())
            shifts.append(-k)
    corrs, ratios, shifts = map(np.array, (corrs, ratios, shifts))
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.5))
    axes[0].hist(corrs, bins=40, range=(0.9, 1.0))
    axes[0].set_xlabel("correlation (shift-aligned)")
    axes[1].hist(ratios, bins=40, range=(0.5, 1.5))
    axes[1].set_xlabel("peak ratio wct/ref")
    axes[2].hist(shifts, bins=17, range=(-8.5, 8.5))
    axes[2].set_xlabel("peak shift [ticks]")
    fig.suptitle(f"OpDecon vs LArSoft deconv, evt {evt} ({len(corrs)} snippets)")
    fig.tight_layout()
    fig.savefig(f"{outdir}/decon_compare.png", dpi=110)
    plt.close(fig)
    return dict(n=len(corrs), corr_med=np.median(corrs), ratio_med=np.median(ratios),
                shift_med=np.median(shifts))


def compare_hits(h_wct, h_ref, outdir):
    eff = tot = 0
    pe_ratios, dts = [], []
    for r in h_ref:
        ch, t, pe = int(r[0]), r[1], r[5]
        cand = h_wct[h_wct[:, 0] == ch]
        if not len(cand):
            continue
        tot += 1
        j = np.argmin(np.abs(cand[:, 1] - t))
        dt = cand[j, 1] - t
        # 16-tick double-ULP quantization, or inside the hit width.
        if abs(dt) <= max(16 * TICK, 0.5 * r[2]):
            eff += 1
            dts.append(dt / TICK)
            if pe > 0:
                pe_ratios.append(cand[j, 5] / pe)
    surplus = 0
    for w in h_wct:
        cand = h_ref[h_ref[:, 0] == int(w[0])]
        if not len(cand) or np.abs(cand[:, 1] - w[1]).min() > max(16 * TICK, 0.5 * w[2]):
            surplus += 1
    pe_ratios, dts = np.array(pe_ratios), np.array(dts)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
    axes[0].hist(dts, bins=33, range=(-16.5, 16.5))
    axes[0].set_xlabel("matched-hit dt [ticks]")
    axes[1].hist(np.log10(np.clip(pe_ratios, 1e-3, 1e3)), bins=40, range=(-1.5, 2.5))
    axes[1].set_xlabel("log10(PE wct/ref), matched hits")
    fig.suptitle(f"OpHitFinder vs PerOpHitTree (eff {eff}/{tot}, surplus {surplus}/{len(h_wct)})")
    fig.tight_layout()
    fig.savefig(f"{outdir}/hit_compare.png", dpi=110)
    plt.close(fig)
    return dict(eff=eff, tot=tot, surplus=surplus, nwct=len(h_wct),
                pe_ratio_med=np.median(pe_ratios) if len(pe_ratios) else np.nan)


def compare_flashes(m_wct, s_wct, m_ref, s_ref, outdir):
    # Match by time, nearest within 10 us, greedy by wct flash size.
    order = np.argsort(-m_wct[:, 1:].sum(axis=1))
    used = set()
    rows = []
    for f in order:
        t = m_wct[f, 0]
        d = np.abs(m_ref[:, 0] - t)
        for j in np.argsort(d):
            if j in used or d[j] > 10e3:
                break
            used.add(j)
            common = np.intersect1d(np.nonzero(m_wct[f, 1:])[0], np.nonzero(m_ref[j, 1:])[0])
            ratio = (m_wct[f, 1 + common].sum() / m_ref[j, 1 + common].sum()) if len(common) else np.nan
            rows.append((m_wct[f, 0] - m_ref[j, 0], ratio,
                         s_wct[f, 2] - s_ref[j, 2], s_wct[f, 3] - s_ref[j, 3],
                         m_ref[j, 1:].sum()))
            break
    rows = np.array(rows) if rows else np.zeros((0, 5))
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.5))
    axes[0].hist(rows[:, 0] / 1e3, bins=40, range=(-10, 10))
    axes[0].set_xlabel("flash dt wct-ref [us]")
    axes[1].hist(np.clip(rows[:, 1], 0, 3), bins=30, range=(0, 3))
    axes[1].set_xlabel("PE ratio (common channels)")
    axes[2].scatter(rows[:, 2] / 10, rows[:, 3] / 10, s=10)
    axes[2].set_xlabel("dy [cm]")
    axes[2].set_ylabel("dz [cm]")
    fig.suptitle(f"OpFlashFinder vs LArSoft flashes ({len(rows)} matched, "
                 f"{m_wct.shape[0]} wct / {m_ref.shape[0]} ref)")
    fig.tight_layout()
    fig.savefig(f"{outdir}/flash_compare.png", dpi=110)
    plt.close(fig)
    return dict(nmatch=len(rows), nwct=m_wct.shape[0], nref=m_ref.shape[0],
                dt_med_us=np.median(rows[:, 0]) / 1e3 if len(rows) else np.nan,
                pe_ratio_med=np.nanmedian(rows[:, 1]) if len(rows) else np.nan)


def main(workdirs):
    for wd in workdirs:
        outdir = os.path.join(wd, "light_validation")
        os.makedirs(outdir, exist_ok=True)
        wct = load_npy(f"{wd}/light-frames-wct.tar.bz2")
        ref = load_npy(f"{wd}/light-frames.tar.bz2")
        d = compare_decon(wct, ref, outdir)
        h = compare_hits(get_tensor(f"{wd}/opflash_pdhd-wct.tar.gz", 2),
                         get_tensor(f"{wd}/opflash_pdhd.tar.gz", 2), outdir)
        f = compare_flashes(get_tensor(f"{wd}/opflash_pdhd-wct.tar.gz", 0),
                            get_tensor(f"{wd}/opflash_pdhd-wct.tar.gz", 1),
                            get_tensor(f"{wd}/opflash_pdhd.tar.gz", 0),
                            get_tensor(f"{wd}/opflash_pdhd.tar.gz", 1), outdir)
        print(f"== {wd}")
        print(f"  decon : {d['n']} snippets, corr med {d['corr_med']:.4f}, "
              f"peak ratio med {d['ratio_med']:.3f}, shift med {d['shift_med']:+.0f} ticks")
        print(f"  hits  : eff {h['eff']}/{h['tot']}, surplus {h['surplus']}/{h['nwct']}, "
              f"PE ratio med {h['pe_ratio_med']:.3f}")
        print(f"  flash : matched {f['nmatch']} (wct {f['nwct']} / ref {f['nref']}), "
              f"dt med {f['dt_med_us']:+.2f} us, PE ratio med {f['pe_ratio_med']:.3f}")


if __name__ == "__main__":
    main(sys.argv[1:])
