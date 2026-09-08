#!/usr/bin/env python3
'''
morse "spcmp" Fourier / downsampling report.

Post-processes the per (plane, case) central splat and sim+OSP waveforms from the
spcmp subgraph and writes a multi-page PDF:

  Page 1 -- Fourier amplitude spectra (0 -> Nyquist) of the central splat and
            sim+OSP waveforms, for each plane (U,V,W) x each case (1 m drift,
            response plane).

  Pages 2.. -- for each case, the effect of two 4x downsamplings on the central
            splat and sim+OSP waveforms, shown in BOTH the time and Fourier
            domains:
              * "rebin"    : sum every N consecutive samples into one.
              * "resample" : low-pass in the Fourier domain (drop the frequency
                             bins near the original Nyquist) to 1/N the samples.

Input files: the per (plane, case) splat/osp frames named by the workflow as
    <detector>-spcmp-<kind>-p<plane>-<case>.npz   kind in {splat, osp}
'''

import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from morse_spcmp_lib import (
    US, CASES, CASE_TITLES, PLANE_LETTERS, Frame, parse_inputs,
    middle_channels, peak_time_and_sigma, rebin, resample, amp_spectrum)

KINDS = ("splat", "osp")
KIND_TITLE = {"splat": "splat", "osp": "sim+OSP"}
SPEC_STYLE = {
    "splat": dict(color="black",    lw=1.4, label="splat"),
    "osp":   dict(color="tab:blue", lw=1.2, label="sim+OSP"),
}
# original / rebin / resample styling for the downsampling page.
DS_STYLE = {
    "orig":     dict(color="0.5",      lw=1.2, ls="-",  label="original"),
    "rebin":    dict(color="tab:orange", lw=1.0, ls="-", marker="o", ms=3),
    "resample": dict(color="tab:green",  lw=1.0, ls="-", marker="s", ms=3),
}


def page_spectra(pdf, detector, frames, chan, planes):
    '''Page 1: |FFT| of splat and OSP central waveforms, 0->Nyquist.'''
    nrow, ncol = len(planes), len(CASES)
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 2.6 * nrow),
                             squeeze=False, sharex=True)
    for ir, p in enumerate(planes):
        letter = PLANE_LETTERS[p] if p < len(PLANE_LETTERS) else str(p)
        for ic, c in enumerate(CASES):
            ax = axes[ir][ic]
            for k in KINDS:
                key = (p, c, k)
                if key not in frames:
                    continue
                fr = frames[key]
                f, X = amp_spectrum(fr.wf(chan[p]), fr.tick / US)
                ax.semilogy(f, X + 1e-12, **SPEC_STYLE[k])
            ax.grid(True, which="both", alpha=0.3)
            if ir == 0:
                ax.set_title(CASE_TITLES[c])
            if ic == 0:
                ax.set_ylabel(f"{letter}-plane\n|amplitude|")
            if ir == 0 and ic == ncol - 1:
                ax.legend(loc="upper right", fontsize="small")
    for ic in range(ncol):
        axes[-1][ic].set_xlabel("frequency [MHz]")
    fig.suptitle(f"{detector.upper()} spcmp: Fourier amplitude spectra "
                 f"(central waveforms, 0->Nyquist)")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    pdf.savefig(fig)
    plt.close(fig)


def page_downsample(pdf, detector, frames, chan, planes, case, n):
    '''One page per case: rebin-N vs resample-N, in time and frequency, for the
    central splat and sim+OSP waveforms of each plane.'''
    cols = [("splat", "time"), ("splat", "freq"), ("osp", "time"), ("osp", "freq")]
    nrow, ncol = len(planes), len(cols)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.0 * ncol, 2.6 * nrow),
                             squeeze=False)
    for ir, p in enumerate(planes):
        letter = PLANE_LETTERS[p] if p < len(PLANE_LETTERS) else str(p)
        # Time-domain zoom window from the (clean) splat peak; reused for the OSP
        # time panel so the noisy OSP waveform cannot blow up the window.
        twin = None
        skey = (p, case, "splat")
        if skey in frames:
            sfr = frames[skey]
            st = sfr.t0 / US + np.arange(sfr.mat.shape[1]) * (sfr.tick / US)
            tpk, sig = peak_time_and_sigma(st, sfr.wf(chan[p]))
            twin = (tpk - 8 * sig, tpk + 8 * sig)
        for ic, (kind, dom) in enumerate(cols):
            ax = axes[ir][ic]
            key = (p, case, kind)
            if key not in frames:
                continue
            fr = frames[key]
            x = fr.wf(chan[p])
            tick_us = fr.tick / US
            t0 = fr.t0 / US
            xr = rebin(x, n)
            xs = resample(x, n)

            if dom == "time":
                t = t0 + np.arange(len(x)) * tick_us
                # bin-center times for the downsampled series
                tr = t0 + (np.arange(len(xr)) * n + (n - 1) / 2.0) * tick_us
                ts = t0 + (np.arange(len(xs)) * n + (n - 1) / 2.0) * tick_us
                ax.plot(t, x, **DS_STYLE["orig"])
                ax.plot(tr, xr, label="rebin-%d" % n, **DS_STYLE["rebin"])
                ax.plot(ts, xs, label="resample-%d" % n, **DS_STYLE["resample"])
                if twin:
                    ax.set_xlim(*twin)
                if ir == nrow - 1:
                    ax.set_xlabel("time [us]")
            else:  # freq
                f, X = amp_spectrum(x, tick_us)
                fr_, Xr = amp_spectrum(xr, tick_us * n)
                fs_, Xs = amp_spectrum(xs, tick_us * n)
                ax.semilogy(f, X + 1e-12, **DS_STYLE["orig"])
                ax.semilogy(fr_, Xr + 1e-12, label="rebin-%d" % n, **DS_STYLE["rebin"])
                ax.semilogy(fs_, Xs + 1e-12, label="resample-%d" % n, **DS_STYLE["resample"])
                if ir == nrow - 1:
                    ax.set_xlabel("frequency [MHz]")
            ax.grid(True, which="both", alpha=0.3)
            if ir == 0:
                ax.set_title("%s  %s" % (KIND_TITLE[kind],
                                         "time" if dom == "time" else "Fourier"))
            if ic == 0:
                ax.set_ylabel(f"{letter}-plane")
            if ir == 0 and ic == 0:
                ax.legend(loc="upper right", fontsize="x-small")
    fig.suptitle("%s spcmp: 4x downsampling (rebin vs resample) -- %s"
                 % (detector.upper(), CASE_TITLES[case]))
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    pdf.savefig(fig)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--detector", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("-n", "--downsample", type=int, default=4,
                    help="downsampling factor for rebin/resample (default 4)")
    ap.add_argument("files", nargs="+")
    args = ap.parse_args()

    inp = parse_inputs(args.files)
    planes = sorted({p for (p, _c, _k) in inp})
    frames = {key: Frame(f) for key, f in inp.items()}
    chan = middle_channels(args.detector, frames, planes)

    with PdfPages(args.output) as pdf:
        page_spectra(pdf, args.detector, frames, chan, planes)
        for c in CASES:
            page_downsample(pdf, args.detector, frames, chan, planes, c,
                            args.downsample)


if __name__ == "__main__":
    main()
