#!/usr/bin/env python3
'''Analyze check-ctoffset.jsonnet output: drift-distance closure.

usage: check-ctoffset.py [options] {bottom,top} OUT.tar.bz2 [OUT.tar.bz2 ...]

For each collection channel under one of the lines of charge, the peak time
of the "gauss" frame (with parabolic interpolation) times the simulated drift
speed gives the reconstructed distance from the W strips.  Prints the mean
and rms of (reconstructed - true) per line.  With --max-bias and/or
--max-spread, exits 1 if any file fails:

  bias   = mean over lines of the per-line mean            |bias| <= max-bias
  spread = max - min of the per-line means                 spread <= max-spread

The wires file is found via WIRECELL_PATH.  Line placement must match
check-ctoffset.jsonnet.
'''
import argparse
import bz2
import io
import json
import os
import sys
import tarfile

import numpy

# Must match check-ctoffset.jsonnet.
DISTANCES = [0.5, 5, 15, 50, 150, 300, 500, 640]  # cm
SEG_Z0, SEG_DZ, SEG_LEN = 20.0, 180.0, 160.0      # mm
EDGE = 15.0                                        # mm dropped at each line end
WIRES = 'dunevd10kt_3view_30deg_v6_refactored.json.bz2'
IDENTS = {'bottom': 0, 'top': 160}


def find_file(fname):
    for d in os.environ.get('WIRECELL_PATH', '').split(':'):
        path = os.path.join(d, fname)
        if d and os.path.exists(path):
            return path
    raise FileNotFoundError(f'{fname} not found in WIRECELL_PATH')


def unwrap(obj):
    'Wires-file objects may be wrapped as {"Type": {...}}.'
    if len(obj) == 1:
        val = next(iter(obj.values()))
        if isinstance(val, dict):
            return val
    return obj


def collection_z(wires_file, anode_ident):
    'Map channel -> wire center z for face 0, plane 2 of the anode.'
    store = json.load(bz2.open(wires_file))['Store']
    pts = store['points']
    for a in store['anodes']:
        a = unwrap(a)
        if a['ident'] != anode_ident:
            continue
        face = unwrap(store['faces'][a['faces'][0]])
        plane = unwrap(store['planes'][face['planes'][2]])
        ret = {}
        for wi in plane['wires']:
            w = unwrap(store['wires'][wi])
            ret[w['channel']] = 0.5 * (unwrap(pts[w['tail']])['z'] + unwrap(pts[w['head']])['z'])
        return ret
    raise KeyError(f'no anode {anode_ident} in {wires_file}')


def load_frame(fname):
    arrs = {}
    with tarfile.open(fname) as tf:
        for m in tf.getmembers():
            arrs[m.name.split('_')[0]] = numpy.load(io.BytesIO(tf.extractfile(m).read()))
    return arrs['frame'], arrs['channels'], arrs['tickinfo']


def peak_time(wave, t0, tick):
    i = int(numpy.argmax(wave))
    if wave[i] <= 0:
        return None
    di = 0.0
    if 0 < i < len(wave) - 1:
        a, b, c = wave[i - 1], wave[i], wave[i + 1]
        den = a - 2 * b + c
        if den != 0:
            di = 0.5 * (a - c) / den
    return t0 + (i + di) * tick


def analyze(fname, chz, speed):
    'Return [(d_true, nchan, mean, rms)] in cm.  speed in mm/us.'
    frame, chans, (t0, tick, tbin) = load_frame(fname)
    t0 = t0 + tbin * tick           # ns
    rows = []
    for i, d in enumerate(DISTANCES):
        z0 = SEG_Z0 + SEG_DZ * i
        errs = []
        for row, ch in enumerate(chans):
            z = chz.get(int(ch))
            if z is None or not (z0 + EDGE <= z <= z0 + SEG_LEN - EDGE):
                continue
            t = peak_time(frame[row], t0, tick)
            if t is None:
                continue
            errs.append(t / 1000.0 * speed / 10.0 - d)   # ns -> us -> mm -> cm
        errs = numpy.array(errs)
        if len(errs):
            rows.append((d, len(errs), errs.mean(), errs.std()))
        else:
            rows.append((d, 0, numpy.nan, numpy.nan))
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('drift', choices=sorted(IDENTS))
    p.add_argument('files', nargs='+')
    p.add_argument('--speed', type=float, default=1.6056288906533874,
                   help='simulated drift speed in mm/us (as given to the jsonnet)')
    p.add_argument('--max-bias', type=float, default=None, help='cm')
    p.add_argument('--max-spread', type=float, default=None, help='cm')
    args = p.parse_args()

    chz = collection_z(find_file(WIRES), IDENTS[args.drift])
    ok = True
    for fname in args.files:
        rows = analyze(fname, chz, args.speed)
        print(f'== {fname}  ({args.drift}, anode {IDENTS[args.drift]})')
        print('  d_true[cm]  nchan  mean(d_reco-d_true)[cm]  rms[cm]')
        for d, n, m, s in rows:
            print(f'  {d:9.1f}  {n:5d}  {m:+12.3f}  {s:10.3f}')
        means = numpy.array([m for d, n, m, s in rows])
        if numpy.any([n == 0 for d, n, m, s in rows]):
            print('  FAIL: a line has no reconstructed channels')
            ok = False
            continue
        bias, spread = means.mean(), means.max() - means.min()
        verdict = ''
        if args.max_bias is not None and abs(bias) > args.max_bias:
            verdict += f' FAIL: |bias| > {args.max_bias}'
        if args.max_spread is not None and spread > args.max_spread:
            verdict += f' FAIL: spread > {args.max_spread}'
        ok = ok and not verdict
        print(f'  bias {bias:+.3f} cm, spread {spread:.3f} cm{verdict}')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
