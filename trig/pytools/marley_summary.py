#!/usr/bin/env python3
'''
Summarize a MARLEY JSON events file as the per-event kinematics the SNB
pointing proxy needs (see snbpoint.py), in the order LowEDepoSource uses
the events.  Writes a .npz with arrays:

- enu :: neutrino energy [MeV]
- te :: kinetic energy of the leading final electron [MeV]
- cos :: cosine of the angle between that electron and the neutrino
- reaction :: 0 for CC on 40Ar, 1 for elastic scattering on electrons,
  2 for anything else (eg CEvNS, no electron)

Requires only numpy.
'''

import argparse
import json
import sys

import numpy

CC, ES, OTHER = 0, 1, 2
reaction_names = ["CC", "ES", "other"]


def summarize(events):
    '''Return a dict of arrays from a list of MARLEY JSON events.'''
    n = len(events)
    out = dict(enu=numpy.zeros(n), te=numpy.zeros(n), cos=numpy.zeros(n),
               reaction=numpy.full(n, OTHER, dtype=numpy.int8))
    for ind, ev in enumerate(events):
        nu = ev["initial_particles"][0]
        pnu = numpy.array([nu["px"], nu["py"], nu["pz"]])
        out["enu"][ind] = nu["E"]
        finals = ev["final_particles"]
        electrons = [p for p in finals if p["pdg"] == 11]
        if not electrons:
            continue
        e = max(electrons, key=lambda p: p["E"])
        pe = numpy.array([e["px"], e["py"], e["pz"]])
        out["te"][ind] = e["E"] - e["mass"]
        norm = numpy.linalg.norm(pe) * numpy.linalg.norm(pnu)
        out["cos"][ind] = pe @ pnu / norm if norm > 0 else 0.0
        # ES leaves a neutrino in the final state, CC a 40K
        if any(abs(p["pdg"]) in (12, 14, 16) for p in finals):
            out["reaction"][ind] = ES
        elif any(p["pdg"] == 1000190400 for p in finals):
            out["reaction"][ind] = CC
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("marley", help="MARLEY JSON events file")
    p.add_argument("output", help="output .npz")
    args = p.parse_args(argv)
    with open(args.marley) as fp:
        events = json.load(fp)["events"]
    numpy.savez_compressed(args.output, **summarize(events))
    return 0


if __name__ == "__main__":
    sys.exit(main())
