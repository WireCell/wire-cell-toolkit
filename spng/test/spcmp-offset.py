#!/usr/bin/env python3
'''
Print the detlinegen --offset string for one plane of the spcmp workflow.

detlinegen places the track center at wp_centers[plane] + offset.  We want every
plane's track at the SAME global point: the collection-plane transverse center,
moved a fixed drift distance into the drift volume along X.  So the offset is
computed per plane as (target - wp_centers[plane]).

Run with the wire-cell-python environment (needs numpy + wirecell), e.g. via the
same interpreter that provides `wcpy`.
'''

import sys
import argparse
import numpy as np

from wirecell import units
from wirecell.gen.linegen import load_wp_spec


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--detector", required=True)
    ap.add_argument("--apa", type=int, default=0)
    ap.add_argument("--plane", type=int, required=True)
    ap.add_argument("--drift-mm", type=float, default=1000.0,
                    help="track drift distance in front of the collection plane [mm]")
    args = ap.parse_args()

    try:
        wp_centers, _ = load_wp_spec(args.detector, args.apa)
    except (KeyError, RuntimeError) as err:
        sys.stderr.write(
            f"spcmp-offset: cannot load wire geometry for detector "
            f"'{args.detector}' (apa {args.apa}): {err!r}\n"
            f"spcmp-offset: is '{args.detector}' registered in a detectors.jsonnet "
            f"on WIRECELL_PATH?  (pdvd is not yet supported.)\n")
        sys.exit(1)
    coll = wp_centers[2]                                  # collection plane center
    drift_sign = np.sign(wp_centers[0][0] - wp_centers[2][0])  # toward induction
    target = np.array([coll[0] + drift_sign * args.drift_mm * units.mm,
                       coll[1], coll[2]])
    off = target - wp_centers[args.plane]
    print("{}*mm,{}*mm,{}*mm".format(off[0] / units.mm,
                                     off[1] / units.mm,
                                     off[2] / units.mm))


if __name__ == "__main__":
    main()
