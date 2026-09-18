#!/usr/bin/env python
'''
Convert the trio .npz DepoFluxSplat writes into HDF5 on the DNN training grid.

A "trio" is one channel in each plane at one tick, known to share a cause
because a single deposition produced all three.  DepoFluxSplat emits them
run-length encoded on its own tick grid; the xvunet dataset wants them as a
flat gather list on the 1500-tick training grid, keyed the same way as the
rec/tru frame files.

Two grids are involved and they are not the same length.  DepoFluxSplat's
acceptance window is the *ductor* window (window_start / window_duration come
from tpc.ductor, not tpc.adc), which for pdhd is the 3 ms readout plus one
response-plane transit, floor(3062.281/0.5) = 6124 bins.  So trio ticks run
0..6123 while the ADC readout is 6000 ticks.  The Reframer downstream sets only
nticks -- tbin and toffset both default to 0 -- so it crops to splat tbins
0..5999 with no shift, and ticks at or past 6000 simply have no pixel.  This
script drops them, then rebins by the SPNGResampler ratio to reach 1500.

Output, per frame ident, matching the /{sample}/{dataset} layout of the frame
files:

    /{ident}/trio_uvwt   (K,4) int16    row_U, row_V, row_W, tick
    /{ident}/trio_q      (K,)  float32  charge

The rows index the concatenated (U,V,W) channel axis the dataset builds, so a
trio is three direct gathers into the model's (1, 2560, 1500) tensor.  Getting
there from the global WCT channel ident is NOT the identity: see CHANNEL_MAPS
below.  Ticks are on the output grid, 0..nticks/rebin-1.

Usage:

    uv run python trios_npz_to_h5.py trios.npz -o trios.h5

(h5py lives in the wire-cell-python venv, so run it under uv.)
'''

import argparse
import re
import sys

import h5py
import numpy as np


# APA-local channel (global ident minus the APA base of 2560*tpcid) -> row on
# the concatenated (U,V,W) axis.  Taking the local channel lets one map serve
# any TPC.
#
# WHICH MAP APPLIES IS A PROPERTY OF THE PIPELINE THAT WROTE THE FRAMES, not of
# the trios, and it is not guessable.  'pdhd-spng' was measured against the
# per-view frame files in rebin_fix_input by scanning every rotation and
# reflection and taking the one that put trios on lit pixels: a single sharp
# peak at exactly 1.0000 over 1.1e6 trios, next candidate 0.93.  The tensor
# schema of dnnroi-training-trios.jsonnet stacks its truth as U, V, W1, W2,
# which need not agree.  Use --scan in trios_check_tru.py against the frames
# you actually have rather than assuming either map.
#
#   'apa-local'  everything straight through, W1 then W2
#   'pdhd-spng'  V reflected and rotated by half a plane; the two collection
#                faces in the opposite order, so local 2080..2559 are rows
#                0..479 of the W view
#
# CAVEAT on pdhd-spng: those depos occupy one drift volume, so only W channels
# 2080..2559 appear.  The half covering 1600..2079 is UNVERIFIED; the modulo is
# what a face swap implies, not something the data could test.
CHANNEL_MAPS = {
    'apa-local': (
        lambda c: c,
        lambda c: c,
        lambda c: c,
    ),
    'pdhd-spng': (
        lambda c: c,
        lambda c: 800 + (1999 - c) % 800,
        lambda c: 1600 + (c - 1120) % 960,
    ),
}


def apa_base(*cols):
    '''
    The APA base channel, 2560*tpcid, inferred from the channels present.

    A trio's three channels are all in one APA by construction, so the lowest
    channel seen, rounded down to a multiple of 2560, identifies it.  This is
    inferred rather than taken from the config because the trio file is the
    only thing at hand and it already says which APA it describes.
    '''
    lo = min(int(c.min()) for c in cols)
    return (lo // 2560) * 2560


def expand_runs(runs, value):
    '''
    Expand run-length encoded trios to one entry per tick.

    runs is (R,5) -- chan_U, chan_V, chan_W, first_tick, n_ticks -- and value
    is the per-tick charge in run order, so sum(n_ticks) == len(value).

    Returns (u, v, w, tick, q), each of length len(value).
    '''
    runs = runs.astype(np.int64)
    value = np.asarray(value, dtype=np.float64).ravel()

    n = runs[:, 4]
    total = int(n.sum())
    if total != value.size:
        raise ValueError(f'sum(n_ticks)={total} != len(value)={value.size}; '
                         'the runs and values are not from the same write')

    # Position within each run, as a flat ramp minus each run's start offset.
    starts = np.cumsum(n) - n
    within = np.arange(total, dtype=np.int64) - np.repeat(starts, n)

    tick = np.repeat(runs[:, 3], n) + within
    u = np.repeat(runs[:, 0], n)
    v = np.repeat(runs[:, 1], n)
    w = np.repeat(runs[:, 2], n)
    return u, v, w, tick, value


def to_training_grid(u, v, w, tick, q, nticks, rebin, nchan):
    '''
    Drop ticks outside the readout, rebin the tick axis, re-sum duplicates.

    Rebinning maps several source ticks onto one output tick, so entries that
    were distinct become the same (u,v,w,tick) and their charges must be added
    rather than one of them kept.
    '''
    keep = tick < nticks
    u, v, w, tick, q = u[keep], v[keep], w[keep], tick[keep], q[keep]

    tick = tick // rebin
    ntout = nticks // rebin

    # One int64 key per (u,v,w,tick) so the dedup is a single sort.  The widest
    # value is ~2.5e13 for nchan=2560, ntout=1500, far inside int64.
    key = ((u * nchan + v) * nchan + w) * ntout + tick
    ukey, inverse = np.unique(key, return_inverse=True)
    qsum = np.bincount(inverse, weights=q)

    tick = (ukey % ntout).astype(np.int64)
    rest = ukey // ntout
    w = (rest % nchan).astype(np.int64)
    rest = rest // nchan
    v = (rest % nchan).astype(np.int64)
    u = (rest // nchan).astype(np.int64)
    return u, v, w, tick, qsum


def pack(u, v, w, tick, q):
    '''
    Pack to the output dtypes, checking the int16 range rather than wrapping.
    '''
    uvwt = np.stack([u, v, w, tick], axis=1)
    lo, hi = int(uvwt.min(initial=0)), int(uvwt.max(initial=0))
    if lo < 0 or hi > np.iinfo(np.int16).max:
        raise ValueError(f'trio index range [{lo},{hi}] does not fit in int16')
    return uvwt.astype(np.int16), q.astype(np.float32)


def idents(npz):
    '''
    The frame idents present, in numeric order.
    '''
    found = set()
    for name in npz.files:
        m = re.fullmatch(r'trio_runs_(\d+)', name)
        if m:
            found.add(int(m[1]))
    return sorted(found)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input', help='trio .npz written by DepoFluxSplat')
    ap.add_argument('-o', '--output', required=True, help='output .h5')
    ap.add_argument('--nticks', type=int, default=6000,
                    help='ADC readout length; ticks at or past this are '
                         'dropped (default 6000)')
    ap.add_argument('--rebin', type=int, default=4,
                    help='tick rebin factor, the SPNGResampler ratio '
                         '(default 4, giving 1500)')
    ap.add_argument('--nchan', type=int, default=2560,
                    help='concatenated channel axis length (default 2560)')
    ap.add_argument('--channel-map', default='pdhd-spng',
                    choices=sorted(CHANNEL_MAPS),
                    help='APA-local channel -> concatenated row '
                         '(default pdhd-spng)')
    ap.add_argument('--chan-base', type=int, default=-1,
                    help='APA base channel, 2560*tpcid, subtracted before the '
                         'map.  Default -1 infers it from the channels present.')
    ap.add_argument('--min-charge', type=float, default=0.0,
                    help='drop trios below this charge AFTER rebinning.  '
                         'Default 0: the loss exists to catch the plane that '
                         'missed, and that plane is systematically the '
                         'low-charge one, so a charge cut removes the very '
                         'population being supervised.')
    args = ap.parse_args(argv)

    if args.nticks % args.rebin:
        ap.error(f'--nticks {args.nticks} is not a multiple of '
                 f'--rebin {args.rebin}')

    npz = np.load(args.input)

    offsets = npz.get('trio_tick_offsets')
    if offsets is not None and np.any(np.asarray(offsets)):
        # With nonzero per-plane offsets the three channels of a trio are hot
        # at three different ticks.  One nominal tick plus the offsets is
        # lossless on the source grid, but after rebinning the three planes can
        # land in different output bins, so a single tick column would be
        # wrong.  Refuse rather than write something subtly incorrect.
        ap.error(f'trio_tick_offsets is {np.asarray(offsets).ravel().tolist()}, '
                 'not all zero.  A single output tick column cannot represent '
                 'per-plane offsets across a rebin; this script needs '
                 'extending to per-plane ticks first.')

    ids = idents(npz)
    if not ids:
        ap.error(f'no trio_runs_* arrays in {args.input}')

    ntout = args.nticks // args.rebin
    kept_all = dropped_all = 0
    q_all = q_dropped = 0.0

    with h5py.File(args.output, 'w') as out:
        for ident in ids:
            runs = npz[f'trio_runs_{ident}']
            value = npz[f'trio_value_{ident}']

            u, v, w, tick, q = expand_runs(runs, value)
            nsrc, qsrc = q.size, q.sum()

            base = (apa_base(u, v, w) if args.chan_base < 0
                    else args.chan_base)
            cmap = CHANNEL_MAPS[args.channel_map]
            u, v, w = cmap[0](u - base), cmap[1](v - base), cmap[2](w - base)
            for name, col in (('U', u), ('V', v), ('W', w)):
                if col.min() < 0 or col.max() >= args.nchan:
                    raise ValueError(
                        f'ident {ident}: {name} rows [{col.min()},{col.max()}] '
                        f'outside the {args.nchan}-row concatenated axis under '
                        f'channel map {args.channel_map!r}')

            u, v, w, tick, q = to_training_grid(
                u, v, w, tick, q, args.nticks, args.rebin, args.nchan)

            if args.min_charge > 0:
                keep = q >= args.min_charge
                u, v, w, tick, q = u[keep], v[keep], w[keep], tick[keep], q[keep]

            uvwt, qout = pack(u, v, w, tick, q)

            grp = out.create_group(str(ident))
            grp.create_dataset('trio_uvwt', data=uvwt,
                               compression='gzip', shuffle=True)
            grp.create_dataset('trio_q', data=qout,
                               compression='gzip', shuffle=True)

            kept_all += qout.size
            dropped_all += nsrc
            q_all += float(qout.sum())
            q_dropped += float(qsrc)
            print(f'ident {ident}: {nsrc} source entries -> {qout.size} trios, '
                  f'charge {qsrc:.6g} -> {qout.sum():.6g}')

        out.attrs['source'] = args.input
        out.attrs['nticks_in'] = args.nticks
        out.attrs['rebin'] = args.rebin
        out.attrs['nticks_out'] = ntout
        out.attrs['nchan'] = args.nchan
        out.attrs['min_charge'] = args.min_charge
        out.attrs['channel_map'] = args.channel_map

    print(f'\n{len(ids)} idents, {dropped_all} source entries -> {kept_all} '
          f'trios ({kept_all/dropped_all:.3f}x)')
    print(f'charge {q_dropped:.6g} -> {q_all:.6g} '
          f'({100*(q_dropped-q_all)/q_dropped:.3f}% dropped past tick '
          f'{args.nticks})')
    print(f'wrote {args.output}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
