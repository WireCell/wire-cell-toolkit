#!/usr/bin/env python3
"""One-time extraction of the PDHD optical geometry and hardware
channel map from the temporary light-data ROOT files into JSON config
data under cfg/pgrapher/experiment/pdhd/.

Usage:
  python3 extract_pdhd_light_maps.py <onevent_run27305_final.root> <outdir>

Outputs (positions converted cm -> mm, the WCT length unit):
  pdhd-opdet-geom.json : 160 OpDets {opdet, x, y, z} from flashopdet/opdet_geo.
                         For PDHD, offline OpChannel == OpDet (0..159).
  pdhd-opch-map.json   : hardware DAPHNE channel -> OpDet from
                         flashopdet/opch_map (only present in the
                         onevent_run27305 file; 224 of 256 hardware
                         channels, 4 per instrumented OpDet).
"""
import json
import sys

import uproot

CM = 10.0  # mm per cm (WCT length unit is mm)


def main(fname, outdir):
    f = uproot.open(fname)

    g = f["flashopdet/opdet_geo"].arrays(library="np")
    opdets = [
        dict(opdet=int(od), x=float(x) * CM, y=float(y) * CM, z=float(z) * CM)
        for od, x, y, z in zip(g["opdet"], g["x"], g["y"], g["z"])
    ]
    geom = {
        "comment": "PDHD optical detector positions [mm], from flashopdet/opdet_geo "
                   "of the temporary light-data ROOT files. Offline OpChannel == OpDet "
                   "(ChannelsPerOpDet=1), 0..159.",
        "source": fname.split("/")[-1],
        "opdets": opdets,
    }
    with open(f"{outdir}/pdhd-opdet-geom.json", "w") as fp:
        json.dump(geom, fp, indent=1)
    print(f"wrote {outdir}/pdhd-opdet-geom.json ({len(opdets)} opdets)")

    m = f["flashopdet/opch_map"].arrays(library="np")
    rows = [dict(hwch=int(c), opdet=int(od)) for c, od in zip(m["opch"], m["opdet"])]
    chmap = {
        "comment": "PDHD hardware DAPHNE channel -> OpDet map, from flashopdet/opch_map. "
                   "Informational only: all reconstruction-level products use the offline "
                   "OpChannel basis, which equals OpDet. Covers the instrumented channels "
                   "of the 2024 run configuration.",
        "source": fname.split("/")[-1],
        "channels": rows,
    }
    with open(f"{outdir}/pdhd-opch-map.json", "w") as fp:
        json.dump(chmap, fp, indent=1)
    print(f"wrote {outdir}/pdhd-opch-map.json ({len(rows)} hardware channels)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
