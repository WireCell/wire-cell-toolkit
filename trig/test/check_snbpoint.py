#!/usr/bin/env python3
'''
Checks of the SNB pointing proxy (trig/pytools/snbpoint.py) and the MARLEY
summary (trig/pytools/marley_summary.py) on synthetic inputs.  Run by
test-snbpoint.bats; exits non-zero on failure.
'''

import os
import sys
import types

import numpy

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "pytools"))
import marley_summary  # noqa: E402
import snbpoint  # noqa: E402


def particle(pdg, E, mass, p, charge=0):
    return dict(pdg=pdg, E=E, mass=mass, px=p[0], py=p[1], pz=p[2], charge=charge)


def check_summary():
    nu = particle(12, 20.0, 0, (0, 0, 20.0))
    ar = particle(1000180400, 37224.7, 37224.7, (0, 0, 0))
    me = 0.511
    pe = numpy.sqrt(10.0**2 - me**2)
    cc = dict(initial_particles=[nu, ar],
              final_particles=[particle(11, 10.0, me, (pe, 0, 0), -1),
                               particle(1000190400, 37226.0, 37225.7, (0, 0, 1.0), 1)])
    es = dict(initial_particles=[nu, ar],
              final_particles=[particle(11, 10.0, me, (0, 0, pe), -1), particle(12, 10.5, 0, (0, 0, 10.5))])
    out = marley_summary.summarize([cc, es])
    assert list(out["reaction"]) == [marley_summary.CC, marley_summary.ES], out["reaction"]
    assert numpy.allclose(out["te"], 10.0 - me)
    assert numpy.allclose(out["cos"], [0.0, 1.0])
    assert numpy.allclose(out["enu"], 20.0)


def fake(cos, te, reaction):
    n = len(cos)
    study = types.SimpleNamespace(inwin=numpy.ones(n, dtype=bool))
    kin = types.SimpleNamespace(cos=numpy.asarray(cos, float), te=numpy.asarray(te, float),
                                enu=numpy.full(n, 20.0), reaction=numpy.full(n, reaction))
    return study, kin


def check_information():
    # all electrons at the same angle: mu^2 / s exactly, per event
    c = 0.8
    study, kin = fake([c] * 100, [12.0] * 100, snbpoint.ES)
    p = snbpoint.Pointing(study, kin, snbpoint.Reco())
    found = numpy.ones(100, dtype=bool)
    info = p.information(found, 50, snbpoint.ES)
    assert numpy.isclose(info, 50 * c**2 / ((1 - c**2) / 2)), info
    # half found: half the information
    half = numpy.arange(100) < 50
    assert numpy.isclose(p.information(half, 50, snbpoint.ES), info / 2)
    assert numpy.isclose(p.efficiency(half, snbpoint.ES), 0.5)
    # no CC in the sample
    assert p.information(found, 50, snbpoint.CC) == 0.0

    # isotropic: the unbiased estimate is consistent with no information
    rng = numpy.random.default_rng(1)
    n = 20000
    study, kin = fake(rng.uniform(-1, 1, n), rng.uniform(1, 40, n), snbpoint.CC)
    p = snbpoint.Pointing(study, kin, snbpoint.Reco())
    iso = p.information(numpy.ones(n, dtype=bool), 1000, snbpoint.CC)
    assert iso < 0.5, iso  # per event information of a forward sample is ~1

    # reconstruction smearing and head-tail flips lose information
    study, kin = fake(numpy.full(2000, 0.9), numpy.full(2000, 10.0), snbpoint.ES)
    ideal = snbpoint.Pointing(study, kin, snbpoint.Reco()).information(numpy.ones(2000, bool), 100, snbpoint.ES)
    smeared = snbpoint.Pointing(study, kin, snbpoint.Reco(25, 1, 0.2)).information(
        numpy.ones(2000, bool), 100, snbpoint.ES)
    assert 0 < smeared < ideal, (smeared, ideal)
    c1, _ = snbpoint.Reco(0, 1, 0.5, nrep=2000).measured_moments(numpy.array([1.0]), numpy.array([10.0]))
    assert abs(c1[0]) < 0.1, c1  # half flipped: no mean direction


def check_scales():
    # sigma = 1 rad
    assert numpy.isclose(snbpoint.theta68_deg(1.0), numpy.degrees(snbpoint.theta68_per_sigma))
    assert numpy.isclose(snbpoint.theta68_per_sigma, 1.5096, atol=1e-3)
    assert snbpoint.theta68_deg(0.0) == float("inf")
    # one flavor with half the mean energy has twice the flux
    r = snbpoint.es_per_cc(1.0, {"ve": (0.1, 10, 1), "vx": (0.1, 5, 2)})
    assert numpy.isclose(r, 0.1 + 0.1 * 2 * 2)


if __name__ == "__main__":
    check_summary()
    check_information()
    check_scales()
    print("snbpoint checks pass")
