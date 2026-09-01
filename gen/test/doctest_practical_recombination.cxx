// PracticalBoxRecombination: the Modified Box model with parameters in
// practical units (kV/cm, (kV/cm)(g/cm^2)/MeV, g/cm^3), as every detector
// config in this tree and the WCP prototype express them.
//
// These are PHYSICS ANCHORS, deliberately absolute rather than structural.
// Before doc 88 nothing in the tree would have noticed a factor-of-10 change
// in the quenching term: the three uBooNE-fixture clus doctests assert only
// std::isfinite(), so upstream e6fb7ef3 moved a MIP from 2.10 to 1.37 MeV/cm
// with every test still green.  A bare CHECK on R(MIP) is what catches that.

#include "WireCellGen/PracticalRecombinationModels.h"
#include "WireCellGen/RecombinationModels.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>

using namespace WireCell;

TEST_CASE("practical box recombination reproduces the LArSoft ModBox point")
{
    // LArSoft kModBoxA = 0.93, kModBoxB = 0.212 (kV/cm)(g/cm^2)/MeV,
    // rho = 1.383 g/cm^3, E = 0.5 kV/cm.  A 2.1 MeV/cm MIP must survive at
    // R ~ 0.70 -- the number every LArTPC dQ/dx calibration is built on.
    Gen::PracticalBoxRecombination model(0.5, 0.93, 0.212, 1.383, 23.6e-6);

    const double dEdx = 2.1 * units::MeV / units::cm;
    const double dX = 1.0 * units::cm;
    const double dQ = model(dEdx * dX, dX);
    const double R = dQ * 23.6e-6 / (dEdx * dX);

    // closed form: R = ln(A + B*dEdx/(E*rho)) / (B*dEdx/(E*rho))
    const double xi = 0.212 * 2.1 / (0.5 * 1.383);
    CHECK(R == doctest::Approx(std::log(0.93 + xi) / xi).epsilon(1e-12));
    CHECK(R == doctest::Approx(0.7044).epsilon(1e-3));
    CHECK(R > 0.6);
    CHECK(R < 0.8);
}

TEST_CASE("practical box recombination matches the WCP prototype at the uBooNE point")
{
    // uboone-mabc.jsonnet uBooNE_box_recomb_model, verbatim.
    Gen::PracticalBoxRecombination model(0.273, 1.0, 0.255, 1.38, 23.6e-6);

    const double dX = 1.0 * units::cm;

    // A 2.10 MeV/cm MIP at uBooNE's 0.273 kV/cm gives ~5.5e4 e-/cm.  Pinned
    // to the value the pre-e6fb7ef3 binary produced (doc 88 sec 3).
    CHECK(model(2.1 * units::MeV, dX) == doctest::Approx(55362.11321505277).epsilon(1e-12));

    // prototype (wcp-uboone-bdt/inc/WCPLEEANA/cuts.h lines 394-395):
    //   median_dedx = (exp(dqdx * 23.6e-6*beta/1.38/0.273) - alpha)
    //                 / (beta/1.38/0.273)
    // with dqdx in e-/cm and the result in MeV/cm.
    const double beta = 0.255, alpha = 1.0;
    for (double dqdx : {4.0e4, 5.54e4, 1.0e5, 2.0e5}) {
        const double proto = (std::exp(dqdx * 23.6e-6 * beta / 1.38 / 0.273) - alpha) / (beta / 1.38 / 0.273);
        // dE() returns an energy for the whole dX; divide by dX in cm for MeV/cm.
        const double got = model.dE(dqdx, dX) / units::MeV / (dX / units::cm);
        CHECK(got == doctest::Approx(proto).epsilon(1e-12));
    }

    // forward/inverse round trip
    const double dE_in = 2.1 * units::MeV;
    CHECK(model.dE(model(dE_in, dX), dX) == doctest::Approx(dE_in).epsilon(1e-9));
}

TEST_CASE("the two Box classes are each correct in their own unit convention")
{
    // Same physics point, expressed the two ways.  This is the trap doc 88
    // exists for: feeding practical-unit numbers to the WCT-unit class is
    // wrong by units::cm/units::MeV = 10 in the quenching term.
    const double dX = 1.0 * units::cm, dE = 2.1 * units::MeV;

    Gen::PracticalBoxRecombination practical(0.5, 0.93, 0.212, 1.383, 23.6e-6);
    Gen::BoxRecombination wct(0.5 * units::kilovolt / units::cm, 0.93,
                              0.212 * units::gram / (units::MeV * units::cm2) * (units::kilovolt / units::cm),
                              1.383 * units::gram / units::cm3, 23.6e-6);

    const double R_practical = practical(dE, dX) * 23.6e-6 / dE;
    const double R_wct = wct(dE, dX) * 23.6e-6 / dE;
    CHECK(R_practical == doctest::Approx(R_wct).epsilon(1e-9));
    CHECK(R_wct == doctest::Approx(0.7044).epsilon(1e-3));

    // ... and the mismatched combination is NOT the same, by ~2x.  If this
    // ever starts passing, the two classes have converged and doc 88's
    // uBooNE binding should be revisited.
    Gen::BoxRecombination mismatched(0.5, 0.93, 0.212, 1.383, 23.6e-6);
    const double R_bad = mismatched(dE, dX) * 23.6e-6 / dE;
    CHECK(R_bad != doctest::Approx(R_practical).epsilon(1e-3));
}
