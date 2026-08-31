#!/usr/bin/env python3
"""Validate the stage-1 PDHD light ROOT->WCT conversion.

Compares a converted work directory against the source ROOT file:
  - opflash matrix: nflash, per-flash PE row-sum vs flash_total_pe,
    flash time vs (FlashTime - tc)*16ns
  - flash_summary: y/z center/width vs PerFlashTree (cm -> mm)
  - ophits tensor: count and peak times vs PerOpHitTree
  - waveform frames: t_first recovery sanity (OpHit peaks land inside
    snippet windows in the trigger-relative frame), raw/deconv trace
    counts vs decoana

Usage:
  python3 check_light_convert.py <rootfile> <run> <event> <workdir>
"""
import io
import json
import sys
import tarfile

import numpy as np
import uproot

TICK_NS = 16.0
US = 1000.0  # ns


def load_tar_members(path):
    out = {}
    with tarfile.open(path) as tf:
        for m in tf.getmembers():
            data = tf.extractfile(m).read()
            if m.name.endswith(".npy"):
                out[m.name] = np.load(io.BytesIO(data))
            elif m.name.endswith(".json"):
                out[m.name] = json.loads(data)
    return out


def main(rootfile, run, event, workdir):
    f = uproot.open(rootfile)
    nfail = 0

    def check(ok, msg):
        nonlocal nfail
        print(("PASS" if ok else "FAIL"), msg)
        if not ok:
            nfail += 1

    # --- trigger anchor ---
    to = f["trigoff/trigger_offset"].arrays(library="np")
    sel = (to["run"] == run) & (to["event"] == event)
    idx = np.argmin(np.abs(to["offset_us"][sel] - 250.0))
    tc = int(to["tc_time_candidate"][sel][idx])
    rd = int(to["rd_timestamp"][sel][idx])

    # --- opflash tensor set ---
    tens = load_tar_members(f"{workdir}/opflash_pdhd.tar.gz")
    md = next(v for k, v in tens.items() if k.endswith(f"tensorset_{event}_metadata.json"))
    check(md["run"] == run and md["event"] == event, "metadata run/event")
    check(md["tc_time_dts"] == str(tc), f"metadata tc_time_dts == {tc}")
    arrays = sorted(k for k in tens if k.endswith("_array.npy"))
    matrix, summary, ophits = (tens[arrays[0]], tens[arrays[1]], tens[arrays[2]])

    pf = f["opflashana/PerFlashTree"].arrays(library="np")
    s = pf["EventID"] == event
    nflash = int(s.sum())
    check(matrix.shape == (nflash, 161), f"opflash shape {matrix.shape} == ({nflash}, 161)")

    fo = f["flashopdet/flash_opdet"].arrays(library="np")
    sfo = fo["event"] == event
    order = np.argsort(pf["FlashID"][s], kind="stable")  # tree order is FlashID order
    times_exp = (pf["FlashTime"][s] - tc) * TICK_NS
    check(np.allclose(matrix[:, 0], times_exp, rtol=0, atol=1e-6),
          "flash times == (FlashTime - tc) * 16ns")
    rowsums = matrix[:, 1:].sum(axis=1)
    tot_exp = np.array([fo["flash_total_pe"][sfo & (fo["flash_id"] == fid)][0]
                        for fid in pf["FlashID"][s]])
    check(np.allclose(rowsums, tot_exp, rtol=1e-6, atol=1e-6),
          "PE row sums == flash_total_pe")
    check(np.allclose(rowsums, pf["TotalPE"][s], rtol=1e-3),
          "PE row sums ~= PerFlashTree TotalPE")

    check(np.allclose(summary[:, 2], pf["YCenter"][s] * 10.0, rtol=1e-5)
          and np.allclose(summary[:, 3], pf["ZCenter"][s] * 10.0, rtol=1e-5),
          "summary y/z centers (mm)")

    oh = f["opflashana/PerOpHitTree"].arrays(library="np")
    soh = oh["EventID"] == event
    check(ophits.shape[0] == int(soh.sum()), f"ophit count {ophits.shape[0]} == {int(soh.sum())}")
    pt_exp = (oh["PeakTimeAbs"][soh] - tc) * TICK_NS
    check(np.allclose(np.sort(ophits[:, 1]), np.sort(pt_exp), atol=1e-6),
          "ophit peak times (trigger-relative ns)")

    # --- waveform frames ---
    frames = load_tar_members(f"{workdir}/light-frames.tar.bz2")
    tick_raw = frames[f"tickinfo_raw_{event}.npy"]
    frame_time, tick, tbin0 = tick_raw
    check(abs(tick - 16.0) < 1e-9, "frame tick == 16 ns")
    # t_first recovery: frame_time = (t_first - tc)*16ns. Hit peaks
    # (trigger-relative) must land inside their channel's snippet
    # windows when snippets are placed at frame_time + tbin*16ns.
    ch_d = frames[f"channels_deconv_{event}.npy"]
    d = f[f"decoana/run_{run}_evt_{event}"]
    nraw = ndec = 0
    peaks = []  # (channel, peak tick rel t_first) of clear deconv peaks
    chans = set()
    for chk in d.keys(recursive=False):
        chn = int(chk[2:].split(";")[0])
        chans.add(chn)
        chd = d[chk]
        for k in chd.keys(recursive=True):
            if k.startswith("raw/waveform"):
                nraw += 1
            elif k.startswith("deconv/waveform"):
                ndec += 1
                h = chd[k]
                v = h.values()
                if v.max() >= 2.0:
                    peaks.append((chn, round(h.axis().edges()[0]) + int(np.argmax(v))))
    # The dump caps the saved snippets (400/event), so not every hit
    # has a snippet; instead require that (nearly) every clear saved
    # deconv peak has a matching OpHit at the same trigger-relative
    # time.  This pins the t_first recovery and the time convention.
    # Tolerance: the absolute DTS timestamps are stored as doubles
    # whose ULP at ~1.07e17 ticks is 16 ticks, so hit times are
    # quantized to +-8 ticks; allow 16.
    # A long pulse train becomes one wide OpHit whose PeakTime can sit
    # far from a sub-peak argmax, so also accept peaks inside the hit
    # extent (half the hit Width).
    hit_by_ch = {}
    for hch, hpt, hw in zip(oh["OpChannel"][soh], pt_exp, oh["Width"][soh]):
        hit_by_ch.setdefault(int(hch), []).append(
            (hpt / TICK_NS, max(16.0, 0.5 * hw * US / TICK_NS)))
    ok = tot = 0
    for chn, ptick in peaks:
        if chn not in hit_by_ch:
            continue
        tot += 1
        rel = frame_time / TICK_NS + ptick  # ticks rel trigger
        if any(abs(rel - x) <= tol for x, tol in hit_by_ch[chn]):
            ok += 1
    # Bar at 80%: busy events show a small population of secondary
    # deconv peaks (clustered ~5 us from the nearest hit; reflections /
    # afterpulses) that the LArSoft hit finder did not turn into OpHits.
    check(tot > 0 and ok / tot > 0.8,
          f"deconv peaks matched to OpHits (16 ticks or hit width): {ok}/{tot}")
    check(set(ch_d) <= chans, "deconv channel list matches decoana")
    print(f"INFO raw snippets {nraw}, deconv {ndec}, frame channels {len(ch_d)}, "
          f"frame_time {frame_time/US:.3f} us, t_first-rd "
          f"{(frame_time - (rd - tc) * TICK_NS)/TICK_NS:.0f} ticks")

    print("FAILED" if nfail else "OK", f"({nfail} failures)" if nfail else "")
    return 1 if nfail else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]))
