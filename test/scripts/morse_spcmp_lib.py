#!/usr/bin/env python3
'''
Shared helpers for the morse "spcmp" plotters (morse-spcmp-plot.py,
morse-spcmp-fft.py).

Frames are read straight from the WCT .npz (frame_/channels_/tickinfo_ arrays)
so the row<->channel-id map is correct even for detectors whose channel idents
are non-contiguous.
'''

import re
import numpy as np

# WCT system of units: time base is ns; 1 us = 1000.
US = 1000.0

CASES = ("drift1m", "resp")
CASE_TITLES = {"drift1m": "1 m drift (diffusion)",
               "resp": "response plane (no diffusion)"}
PLANE_LETTERS = "UVW"

# Channel-id boundaries bounding the three planes (plane p spans [r[p], r[p+1])).
DETECTOR_RANGES = {
    "pdhd": [0, 800, 1600, 2560],
    "pdsp": [0, 800, 1600, 2560],
    "pdvd": [0, 476, 1428, 2488],
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
        self.mat = d[fkey]                       # (nrows, nticks)
        self.chans = np.asarray(d[ckey]).astype(int)
        ti = d[tkey]
        self.t0 = float(ti[0])
        self.tick = float(ti[1])                 # ns
        self.id2row = {int(c): i for i, c in enumerate(self.chans)}
        self.toffset_us = 0.0                     # optional plotting time shift

    def wf(self, chid):
        row = self.id2row.get(int(chid))
        if row is None:
            return np.zeros(self.mat.shape[1], dtype=self.mat.dtype)
        return self.mat[row]

    def times_us(self):
        nt = self.mat.shape[1]
        return (self.t0 + np.arange(nt) * self.tick) / US + self.toffset_us


def parse_inputs(files):
    '''Return {(plane:int, case:str, kind:str): filename}.'''
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


def middle_channels(detector, frames, planes):
    '''Per-plane middle channel id from the 1 m splat frame.'''
    chan = {}
    for p in planes:
        splat_1m = frames[(p, "drift1m", "splat")]
        ranges = channel_ranges(detector, int(splat_1m.chans.max()))
        chan[p] = middle_channel_id(splat_1m, ranges[p], ranges[p + 1])
    return chan


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


# ---- downsampling operations -------------------------------------------------

def rebin(x, n):
    '''WCT "rebin": sum every n consecutive samples into one (drops remainder).'''
    m = len(x) // n
    return x[:m * n].reshape(m, n).sum(axis=1)


def resample(x, n):
    '''WCT "resample": low-pass in the Fourier domain to 1/n the samples.

    FFT the waveform, keep only the lowest 1/n of the frequency bins (around DC),
    drop the rest (those near the original Nyquist), inverse-FFT to len(x)//n
    samples.  Amplitude (not integral) is preserved.
    '''
    L = len(x)
    M = L // n
    X = np.fft.rfft(x)
    Xk = X[:M // 2 + 1]
    y = np.fft.irfft(Xk, n=M)
    return y * (float(M) / float(L))


def amp_spectrum(x, tick_us):
    '''Return (freq_MHz, |rfft|) up to Nyquist for samples spaced tick_us.'''
    X = np.abs(np.fft.rfft(x))
    f = np.fft.rfftfreq(len(x), d=tick_us)   # cycles per us == MHz
    return f, X
