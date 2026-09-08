#!/usr/bin/env python3
'''
morse "spcmp" comparison plotter.

Companion to the morse workflow's spcmp subgraph.  For each wire plane, one ideal
line source (perpendicular to the plane's wires, parallel to the wire planes) is
processed in two cases:

  drift1m : placed 1 m in front of the collection plane (depos diffuse).
  resp    : placed at the response plane (no diffusion).

This makes a 3 (plane: U,V,W) x 2 (case: 1m, response) grid of panels.  Each
panel overlays the splat "true signal" and the sim+OSP signal for the channel
nearest the middle of that plane's line source.

Zoom: each COLUMN (case) has one shared time window, +/- nsigma (default 5)
about the W (collection) plane's splat peak.  W is the cleanest signal, so this
keeps the U/V/W rows of a column on a common, noise-free time axis (the two
columns get their own windows).  All six panels share a single y-axis scale.

Input files are the per (plane, case) splat/osp frames named by the workflow as:

    <detector>-spcmp-<kind>-p<plane>-<case>.npz   kind in {splat, osp}

Frames are read directly from the .npz (frame_/channels_/tickinfo_ arrays) so
the row<->channel-id mapping is correct even for detectors whose channel idents
are non-contiguous (e.g. pdvd).
'''

import re
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# WCT system of units: time base is ns; 1 us = 1000.
US = 1000.0

KINDS = ("splat", "osp")
CASES = ("drift1m", "resp")
CASE_TITLES = {"drift1m": "1 m drift (diffusion)",
               "resp": "response plane (no diffusion)"}
STYLES = {
    "splat": dict(color="black",    lw=1.6, ls="-", zorder=3, label="splat"),
    "osp":   dict(color="tab:blue", lw=1.2, ls="-", zorder=2, label="sim+OSP"),
}
PLANE_LETTERS = "UVW"

# Channel-id boundaries bounding the three planes (plane p spans [r[p], r[p+1])).
DETECTOR_RANGES = {
    "pdhd": [0, 800, 1600, 2560],
    "pdsp": [0, 800, 1600, 2560],
    "pdvd": [0, 476, 1428, 2488],   # idents in 3 blocks (with gaps)
}
NCHAN_RANGES = {
    2560: [0, 800, 1600, 2560],
    2488: [0, 476, 1428, 2488],
    8256: [0, 2400, 4800, 8256],
}


class Frame:
    '''A frame loaded straight from a WCT .npz, with an id->row map.'''
    def __init__(self, fname):
        d = np.load(fname)
        fkey = next(k for k in d.files if k.startswith("frame"))
        ckey = next(k for k in d.files if k.startswith("channels"))
        tkey = next(k for k in d.files if k.startswith("tickinfo"))
        self.mat = d[fkey]                      # (nrows, nticks)
        self.chans = np.asarray(d[ckey]).astype(int)
        ti = d[tkey]
        self.t0 = float(ti[0])
        self.tick = float(ti[1])
        self.id2row = {int(c): i for i, c in enumerate(self.chans)}
        # Optional plotting time shift (us), e.g. to align splat with OSP.
        self.toffset_us = 0.0

    def wf(self, chid):
        row = self.id2row.get(int(chid))
        if row is None:
            return np.zeros(self.mat.shape[1], dtype=self.mat.dtype)
        return self.mat[row]

    def times_us(self):
        nt = self.mat.shape[1]
        return (self.t0 + np.arange(nt) * self.tick) / US + self.toffset_us


def parse_inputs(files):
    pat = re.compile(r"-spcmp-(splat|osp)-p(\d+)-(drift1m|resp)\.npz$")
    out = {}
    for fname in files:
        m = pat.search(fname)
        if not m:
            raise ValueError(f"cannot parse kind/plane/case from: {fname}")
        out[(int(m.group(2)), m.group(3), m.group(1))] = fname
    return out


def channel_ranges(detector, chan_max):
    if detector in DETECTOR_RANGES:
        return DETECTOR_RANGES[detector]
    nchan = chan_max + 1
    if nchan in NCHAN_RANGES:
        return NCHAN_RANGES[nchan]
    raise ValueError(f"unknown channel layout for {detector=} (max id {chan_max})")


def middle_channel_id(frame, id_lo, id_hi, thresh_frac=0.05):
    '''Charge-weighted centroid channel id of activity within [id_lo, id_hi).'''
    ids = np.array([c for c in frame.chans if id_lo <= c < id_hi], dtype=int)
    if len(ids) == 0:
        return int((id_lo + id_hi) // 2)
    rows = np.array([frame.id2row[c] for c in ids])
    qch = np.abs(frame.mat[rows]).sum(axis=1)
    if qch.max() <= 0:
        return int((id_lo + id_hi) // 2)
    active = qch > thresh_frac * qch.max()
    return int(round(np.sum(ids[active] * qch[active]) / np.sum(qch[active])))


def peak_time_and_sigma(times_us, wf, thresh_frac=0.05):
    '''Return (peak_time_us, sigma_us) of the dominant peak by |wf| moments.'''
    amp = np.abs(wf)
    if amp.max() <= 0:
        return 0.5 * (times_us[0] + times_us[-1]), (times_us[-1] - times_us[0]) / 10.0
    sel = amp > thresh_frac * amp.max()
    w, t = amp[sel], times_us[sel]
    tbar = np.sum(w * t) / np.sum(w)
    var = np.sum(w * (t - tbar) ** 2) / np.sum(w)
    return tbar, float(np.sqrt(max(var, 1e-12)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--detector", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--nsigma", type=float, default=5.0,
                    help="half-window in sigma of the 1m peak (default 5)")
    ap.add_argument("--splat-tick-offset", type=float, default=1.0,
                    help="shift splat waveforms by this many ticks to align with "
                         "OSP (OSP applies a 129-tick response roll vs splat's "
                         "128; default 1)")
    ap.add_argument("--title", default="")
    ap.add_argument("files", nargs="+")
    args = ap.parse_args()

    inp = parse_inputs(args.files)
    planes = sorted({p for (p, _c, _k) in inp})
    frames = {key: Frame(f) for key, f in inp.items()}

    # Shift splat in time by the requested number of ticks so it aligns with the
    # OSP response roll (applied to windowing and curves alike).
    for (p, c, k), fr in frames.items():
        if k == "splat":
            fr.toffset_us = args.splat_tick_offset * fr.tick / US
    splat_label = ("splat (+%g tick)" % args.splat_tick_offset
                   if args.splat_tick_offset else "splat")

    # Per-plane middle channel (from the 1 m splat).
    chan = {}
    for p in planes:
        splat_1m = frames[(p, "drift1m", "splat")]
        ranges = channel_ranges(args.detector, int(splat_1m.chans.max()))
        chan[p] = middle_channel_id(splat_1m, ranges[p], ranges[p + 1])

    # One x-window PER COLUMN (case), shared by all three rows.  The window is
    # +/- nsigma about the W (collection) plane's splat peak -- W is the cleanest
    # signal, so this avoids the broadly-distributed noise hits that otherwise
    # widen the U/V windows.
    wplane = max(planes)
    xwin_col = {}
    for c in CASES:
        wsplat = frames[(wplane, c, "splat")]
        tpk, sig = peak_time_and_sigma(wsplat.times_us(), wsplat.wf(chan[wplane]))
        xwin_col[c] = (tpk - args.nsigma * sig, tpk + args.nsigma * sig)

    # Global y-range over all in-window samples (shared y-axis).
    xwin = {}
    gymin, gymax = 0.0, 0.0
    for p in planes:
        for c in CASES:
            lo, hi = xwin_col[c]
            xwin[(p, c)] = (lo, hi)
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
            ax.set_xlim(*xwin[(p, c)])
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
