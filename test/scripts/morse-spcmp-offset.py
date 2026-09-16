#!/usr/bin/env python3
'''
Print the detlinegen --offset string for one plane/case of the morse "spcmp"
subgraph.

detlinegen places the track center at wp_centers[plane] + offset.  We want the
line, for every plane, at the SAME global transverse point (the collection-plane
transverse center) but at a chosen drift distance in front of the collection
plane along X:

  --case drift1m : a fixed distance (default 1 m) in front of collection, so the
                   depos drift that far and pick up diffusion.
  --case resp    : exactly at the response plane, so the drift distance (to the
                   response plane) is ~0 and the depos have no diffusion.  The
                   response-plane distance is the field-response "origin", read
                   from the detector's FR JSON (differs per detector: ~100 mm for
                   pdsp/pdhd/uboone, ~181 mm for pdvd).

So the offset is (target - wp_centers[plane]) with
  target = collection_center + drift_sign * dist  (dist along X only).

Run with the wire-cell-python environment (needs numpy + wirecell).
'''

import sys
import argparse
import numpy as np

from wirecell import units
from wirecell.gen.linegen import load_wp_spec
from wirecell.util import detectors
from wirecell.sigproc.response import persist as rpersist


def response_origin(detector):
    'Return the FR "origin" (response-plane distance) for detector, in WCT units.'
    field = detectors.load('detectors.jsonnet')[detector]['field']
    return rpersist.load(field).origin


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--detector", required=True)
    ap.add_argument("--apa", type=int, default=0)
    ap.add_argument("--plane", type=int, required=True)
    ap.add_argument("--case", choices=("drift1m", "resp"), required=True)
    ap.add_argument("--drift-mm", type=float, default=1000.0,
                    help="drift1m case: distance in front of collection [mm]")
    args = ap.parse_args()

    try:
        wp_centers, _ = load_wp_spec(args.detector, args.apa)
    except (KeyError, RuntimeError) as err:
        sys.stderr.write(
            f"morse-spcmp-offset: cannot load wire geometry for '{args.detector}'"
            f" (apa {args.apa}): {err!r}\n")
        sys.exit(1)

    coll = wp_centers[2]                                        # collection center
    drift_sign = np.sign(wp_centers[0][0] - wp_centers[2][0])   # toward induction

    if args.case == "drift1m":
        dist = args.drift_mm * units.mm
    else:  # resp: place at the response plane (no diffusion)
        dist = response_origin(args.detector)

    target = np.array([coll[0] + drift_sign * dist, coll[1], coll[2]])
    off = target - wp_centers[args.plane]
    print("{}*mm,{}*mm,{}*mm".format(off[0] / units.mm,
                                     off[1] / units.mm,
                                     off[2] / units.mm))


if __name__ == "__main__":
    main()
