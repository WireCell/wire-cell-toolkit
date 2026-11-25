#!/usr/bin/env python3
"""Convert an SBND YZ non-uniformity ROOT map to a YZMap file (json.bz2).

Input: the calibration ROOT file with one TH2 "CzyHist_<plane>_<tpc>" per
wire plane (0=U, 1=V, 2=Y) and TPC (0=East, 1=West): the correction factor C
that makes data dQ/dx uniform in (z, y).

Output: {anode_key: {"<plane>": grid[z][y]}}, as read by the YZMap service and
used by DepoTransformYZ to put the data's YZ non-uniformity into the
simulation.  Pixels without data (C <= 0) get factor 1.

Factor modes (--factor):
  inv-norm  1/C, scaled so its mean over the pixels with data is 1 for each
            plane and TPC (default; keeps the overall charge scale).
  inv       1/C as is.

The YZMap configuration matching the output is printed at the end.

Requires numpy and uproot.
"""
import argparse
import bz2
import json
import sys

import numpy as np
import uproot


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="ROOT file with CzyHist_<plane>_<tpc>")
    ap.add_argument("output", help="output .json.bz2 (or .json)")
    ap.add_argument("--factor", choices=["inv-norm", "inv"], default="inv-norm")
    ap.add_argument("--keys", default="AnodePlane:apa0,AnodePlane:apa1",
                    help="map key of TPC 0 (East) and TPC 1 (West): the anode type:name (default: %(default)s)")
    ap.add_argument("--hist", default="CzyHist", help="histogram name prefix (default: %(default)s)")
    ap.add_argument("--digits", type=int, default=5, help="significant digits kept (default: %(default)s)")
    args = ap.parse_args()

    keys = args.keys.split(",")
    f = uproot.open(args.input)
    out = {}
    binning = None
    for tpc, key in enumerate(keys):
        out[key] = {}
        for plane in range(3):
            h = f[f"{args.hist}_{plane}_{tpc}"]
            c = h.values()                    # shape (nz, ny): x axis is z, y axis is y
            ze, ye = h.axes[0].edges(), h.axes[1].edges()
            this = (len(ze) - 1, len(ye) - 1, ze[1] - ze[0], ye[1] - ye[0], ze[0], ye[0])
            if not (np.allclose(np.diff(ze), ze[1] - ze[0]) and np.allclose(np.diff(ye), ye[1] - ye[0])):
                sys.exit(f"{h.name}: non-uniform binning")
            if binning is None:
                binning = this
            elif not np.allclose(binning, this):
                sys.exit(f"{h.name}: binning {this} differs from {binning}")

            have = c > 0
            g = np.ones_like(c, dtype=float)
            g[have] = 1.0 / c[have]
            if args.factor == "inv-norm":
                g[have] /= g[have].mean()
            out[key][str(plane)] = [[float(f"{v:.{args.digits}g}") for v in row] for row in g]
            print(f"{key} plane {plane}: {have.sum()}/{c.size} pixels with data; factor mean {g[have].mean():.4f} "
                  f"rms {g[have].std():.4f} min {g[have].min():.3f} max {g[have].max():.3f}")

    text = json.dumps(out, separators=(",", ":"))
    opener = bz2.open if args.output.endswith(".bz2") else open
    with opener(args.output, "wt") as fp:
        fp.write(text)

    nz, ny, dz, dy, z0, y0 = binning
    print(f"wrote {args.output}; YZMap configuration (cm):")
    print(f"  nbinsz: {nz}, nbinsy: {ny}, bin_width: {dz:g}*wc.cm, bin_height: {dy:g}*wc.cm, "
          f"zoffset: {0.0 - z0:g}*wc.cm, yoffset: {0.0 - y0:g}*wc.cm")


if __name__ == "__main__":
    main()
