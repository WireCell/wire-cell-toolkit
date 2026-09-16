#!/usr/bin/env python3
'''
morse "spcmp" comparison plotter.

For each wire plane, one ideal line source (perpendicular to the plane's wires,
parallel to the wire planes) is processed in two cases -- 1 m in front of the
collection plane (diffusion) and at the response plane (no diffusion) -- giving a
3 (plane: U,V,W) x 2 (case) grid.  Each panel overlays the splat "true signal"
and the sim+OSP signal for the channel nearest the middle of that plane's line.

Zoom: each COLUMN (case) has one shared time window, +/- nsigma (default 5)
about the W (collection) plane's splat peak.  W is the cleanest signal, so this
keeps the U/V/W rows of a column on a common, noise-free time axis (the two
columns get their own windows).  All six panels share a single y-axis scale.

splat is shifted later by --splat-tick-offset ticks (default 1) to align with
sim+OSP (OSP applies a 129-tick response roll vs splat's 128).

Input files: the per (plane, case) splat/osp frames named by the workflow as
    <detector>-spcmp-<kind>-p<plane>-<case>.npz   kind in {splat, osp}
'''

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from morse_spcmp_lib import (
    US, CASES, CASE_TITLES, PLANE_LETTERS, Frame, parse_inputs,
    channel_ranges, middle_channels, peak_time_and_sigma)

KINDS = ("splat", "osp")
STYLES = {
    "splat": dict(color="black",    lw=1.6, ls="-", zorder=3, label="splat"),
    "osp":   dict(color="tab:blue", lw=1.2, ls="-", zorder=2, label="sim+OSP"),
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--detector", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--nsigma", type=float, default=5.0,
                    help="half-window in sigma of the 1m peak (default 5)")
    ap.add_argument("--splat-tick-offset", type=float, default=1.0,
                    help="shift splat waveforms by this many ticks to align with "
                         "OSP (default 1)")
    ap.add_argument("--title", default="")
    ap.add_argument("files", nargs="+")
    args = ap.parse_args()

    inp = parse_inputs(args.files)
    planes = sorted({p for (p, _c, _k) in inp})
    frames = {key: Frame(f) for key, f in inp.items()}

    for (p, c, k), fr in frames.items():
        if k == "splat":
            fr.toffset_us = args.splat_tick_offset * fr.tick / US
    splat_label = ("splat (+%g tick)" % args.splat_tick_offset
                   if args.splat_tick_offset else "splat")

    chan = middle_channels(args.detector, frames, planes)

    # One x-window per column (case), +/- nsigma about the W-plane splat peak.
    wplane = max(planes)
    xwin_col = {}
    for c in CASES:
        wsplat = frames[(wplane, c, "splat")]
        tpk, sig = peak_time_and_sigma(wsplat.times_us(), wsplat.wf(chan[wplane]))
        xwin_col[c] = (tpk - args.nsigma * sig, tpk + args.nsigma * sig)

    # Global y-range over all in-window samples (shared y-axis).
    gymin, gymax = 0.0, 0.0
    for p in planes:
        for c in CASES:
            lo, hi = xwin_col[c]
            for k in KINDS:
                key = (p, c, k)
                if key not in frames:
                    continue
                fk = frames[key]
                t, w = fk.times_us(), fk.wf(chan[p])
                m = (t >= lo) & (t <= hi)
                if np.any(m):
                    gymin = min(gymin, float(w[m].min()))
                    gymax = max(gymax, float(w[m].max()))
    pad = 0.05 * ((gymax - gymin) or 1.0)
    gymin, gymax = gymin - pad, gymax + pad

    nrow, ncol = len(planes), len(CASES)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 2.6 * nrow),
                             squeeze=False, sharey=True)
    for ir, p in enumerate(planes):
        letter = PLANE_LETTERS[p] if p < len(PLANE_LETTERS) else str(p)
        for ic, c in enumerate(CASES):
            ax = axes[ir][ic]
            for k in KINDS:
                key = (p, c, k)
                if key not in frames:
                    continue
                fk = frames[key]
                style = dict(STYLES[k])
                if k == "splat":
                    style["label"] = splat_label
                ax.plot(fk.times_us(), fk.wf(chan[p]), **style)
            ax.set_xlim(*xwin_col[c])
            ax.set_ylim(gymin, gymax)
            ax.grid(True, alpha=0.3)
            ax.axhline(0, color="0.6", lw=0.6, zorder=1)
            if ir == 0:
                ax.set_title(CASE_TITLES[c])
            if ic == 0:
                ax.set_ylabel(f"{letter}-plane\nsignal")
            ax.text(0.02, 0.95, f"ch {chan[p]}", transform=ax.transAxes,
                    va="top", ha="left", fontsize="small", color="0.4")
            if ir == 0 and ic == ncol - 1:
                ax.legend(loc="upper right", fontsize="small")
    for ic in range(ncol):
        axes[-1][ic].set_xlabel("time [us]")

    sup = f"{args.detector.upper()} morse spcmp: splat vs sim+OSP"
    if args.title:
        sup += f"  {args.title}"
    fig.suptitle(sup)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(args.output)


if __name__ == "__main__":
    main()
