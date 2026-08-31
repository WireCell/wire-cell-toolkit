#!/usr/bin/env python3
'''
spcmp comparison plotter.

Companion to the "spcmp" snakemake workflow.  Given the per-plane (splat, osp,
spng) signal frames produced by that workflow, make a single figure with three
plots (one per U, V, W plane).  Each plot overlays the splat, osp and spng
waveforms taken from a common channel that lies roughly at the middle of the
ideal track (which is perpendicular to that plane's wires and centered in the
transverse plane).

The signal files are passed as positional arguments.  Their (plane, kind) is
parsed from the file name, which the workflow forms as:

    <detector>-signals-<kind>-p<plane>.npz     kind in {splat, osp, spng}

Frame loading (and the associated channel-id ordering / gap handling) is reused
from wirecell.test.ssss so this behaves like the sibling `wcpy test plot-ssss`.
'''

import re
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wirecell import units
from wirecell.test import ssss

# The kinds of signal we compare, in a fixed plotting order with fixed styling.
KINDS = ("splat", "osp", "spng")
STYLES = {
    "splat": dict(color="black",   lw=1.6, ls="-",  zorder=3),
    "osp":   dict(color="tab:blue", lw=1.2, ls="-", zorder=2),
    "spng":  dict(color="tab:red",  lw=1.2, ls="--", zorder=2),
}
PLANE_LETTERS = "UVW"

# Per-detector channel-id ranges bounding the three planes.  Mirrors the
# detection logic in wirecell.test (ssss_args).  Values are channel-id
# boundaries so plane p spans [ranges[p], ranges[p+1]).
DETECTOR_RANGES = {
    "pdhd": [0, 800, 1600, 2560],   # U=0-799, V=800-1599, W=1600-2559
    "pdvd": [0, 476, 1428, 2488],   # U=0-475, V=952-1427, X=1904-2487 (with gaps)
}
NCHAN_RANGES = {
    2560: [0, 800, 1600, 2560],
    2488: [0, 476, 1428, 2488],
    8256: [0, 2400, 4800, 8256],
}


def parse_inputs(files):
    '''Return {(plane:int, kind:str): filename}.'''
    pat = re.compile(r"-signals-(splat|osp|spng)-p(\d+)\.npz$")
    out = {}
    for fname in files:
        m = pat.search(fname)
        if not m:
            raise ValueError(f"cannot parse kind/plane from file name: {fname}")
        kind = m.group(1)
        plane = int(m.group(2))
        out[(plane, kind)] = fname
    return out


def channel_ranges(detector, nchan):
    '''Return channel-id boundary list for the detector (name first, nchan fallback).'''
    if detector in DETECTOR_RANGES:
        return DETECTOR_RANGES[detector]
    if nchan in NCHAN_RANGES:
        return NCHAN_RANGES[nchan]
    raise ValueError(f"unknown channel layout for {detector=} {nchan=}")


def wf_by_id(frame, chid):
    '''Waveform (1D array) for absolute channel id chid from an ssss Frame.'''
    c0 = int(frame.extent[2])
    c1 = int(frame.extent[3])
    if c0 <= chid < c1:
        return frame.frame[chid - c0]
    return np.zeros(frame.frame.shape[1], dtype=frame.frame.dtype)


def middle_channel_id(frame, id_lo, id_hi, thresh_frac=0.05):
    '''Channel id near the middle of the splat activity within [id_lo, id_hi).

    The ideal track is perpendicular to the plane's wires and centered in the
    transverse plane, so it deposits on a contiguous band of channels.  We take
    the (charge-weighted) centroid of that band as "the middle of the track".
    '''
    c0 = int(frame.extent[2])
    ids = np.arange(id_lo, id_hi)
    rows = ids - c0
    valid = (rows >= 0) & (rows < frame.frame.shape[0])
    ids = ids[valid]
    rows = rows[valid]
    qch = np.abs(frame.frame[rows]).sum(axis=1)
    if qch.max() <= 0:
        # No activity: fall back to the geometric middle of the range.
        return int((id_lo + id_hi) // 2)
    active = qch > thresh_frac * qch.max()
    # Charge-weighted centroid over the active band.
    centroid = np.sum(ids[active] * qch[active]) / np.sum(qch[active])
    return int(round(centroid))


def time_axis(frame):
    '''Return times in microseconds for each sample of frame.'''
    t0, tf, _, _ = frame.extent
    nt = frame.frame.shape[1]
    return (t0 + np.arange(nt) * frame.tick) / units.us


def signal_window(times_us, waveform, pad_us=8.0, thresh_frac=0.02):
    '''Return (lo, hi) time window (us) bounding significant |waveform|, padded.'''
    amp = np.abs(waveform)
    if amp.max() <= 0:
        return times_us[0], times_us[-1]
    hits = np.where(amp > thresh_frac * amp.max())[0]
    return times_us[hits[0]] - pad_us, times_us[hits[-1]] + pad_us


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--detector", required=True, help="detector name (pdhd, pdvd)")
    ap.add_argument("-o", "--output", required=True, help="output figure file")
    ap.add_argument("--title", default="", help="extra title text")
    ap.add_argument("files", nargs="+", help="the <det>-signals-<kind>-p<plane>.npz files")
    args = ap.parse_args()

    inputs = parse_inputs(args.files)
    planes = sorted({p for (p, _k) in inputs})

    fig, axes = plt.subplots(len(planes), 1, figsize=(9, 3.0 * len(planes)),
                             squeeze=False)
    axes = axes[:, 0]

    for ax, plane in zip(axes, planes):
        # Load the three signal frames for this plane's WCT job.
        frames = {}
        for kind in KINDS:
            key = (plane, kind)
            if key not in inputs:
                raise ValueError(f"missing input for plane {plane} kind {kind}")
            frames[kind] = ssss.load_frame(inputs[key])

        splat = frames["splat"]
        nchan = splat.frame.shape[0]
        ranges = channel_ranges(args.detector, nchan)
        id_lo, id_hi = ranges[plane], ranges[plane + 1]

        chid = middle_channel_id(splat, id_lo, id_hi)

        # Common time axis (all three share the same tick; splat sets the grid).
        times_us = time_axis(splat)

        for kind in KINDS:
            fr = frames[kind]
            wf = wf_by_id(fr, chid)
            t = time_axis(fr)
            ax.plot(t, wf, label=kind, **STYLES[kind])

        # Zoom the time axis to the splat signal region for a detailed view.
        lo, hi = signal_window(times_us, wf_by_id(splat, chid))
        ax.set_xlim(lo, hi)

        letter = PLANE_LETTERS[plane] if plane < len(PLANE_LETTERS) else str(plane)
        ax.set_title(f"{letter}-plane  (channel {chid})")
        ax.set_ylabel("signal")
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color="0.6", lw=0.6, zorder=1)
        ax.legend(loc="upper right", fontsize="small")

    axes[-1].set_xlabel("time [us]")
    sup = f"{args.detector.upper()} spcmp: splat / osp / spng"
    if args.title:
        sup += f"  {args.title}"
    fig.suptitle(sup)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(args.output)


if __name__ == "__main__":
    main()
