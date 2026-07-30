// PowerBoxRecombination: the SBND free-power Modified Box fit
// (sbnd_xin/docs/55 sec 7g).  Checks the forward against the reference
// python implementation, the bisection inverse round-trip, saturation
// above the monotone branch, and repeat-call determinism.

#include "WireCellGen/RecombinationModels.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

TEST_CASE("powerbox recombination forward matches reference")
{
    Gen::PowerBoxRecombination model;  // defaults = the canonical SBND fit

    // Reference dQ/dx (e/cm) from the independent python implementation of
    //   u = 0.282371*(dEdx/2.1)^1.362179, R = ln(0.93+u)/u,
    //   dQ/dx = 0.855175 * R * dEdx / 23.6e-6
    const struct { double dedx, dqdx; } ref[] = {
        {1.7, 38564.146133},  {2.1, 51897.776028},  {5.0, 121139.693212},
        {10.0, 182657.492846}, {30.0, 251212.560172}, {50.0, 264729.447514},
    };
    for (const auto& r : ref) {
        const double dx = 1 * units::cm;
        const double dQ = model(r.dedx * units::MeV / units::cm * dx, dx);
        CHECK(dQ == doctest::Approx(r.dqdx).epsilon(1e-9));
    }
}

TEST_CASE("powerbox recombination inverse round-trips")
{
    Gen::PowerBoxRecombination model;

    // 0.8 MeV/cm is just above the A<1 zero crossing (~0.754 MeV/cm); the
    // clus/ callers clamp to [0, 50] MeV/cm.
    for (double dedx = 0.8; dedx <= 50.0; dedx += 0.7) {
        const double dx = 0.3 * units::cm;
        const double dE_in = dedx * units::MeV / units::cm * dx;
        const double dQ = model(dE_in, dx);
        const double dE_out = model.dE(dQ, dx);
        CHECK(dE_out == doctest::Approx(dE_in).epsilon(1e-9));
    }
}

TEST_CASE("powerbox recombination inverse edge cases")
{
    Gen::PowerBoxRecombination model;
    const double dx = 1 * units::cm;

    // non-positive charge -> zero energy
    CHECK(model.dE(0.0, dx) == 0.0);
    CHECK(model.dE(-1000.0, dx) == 0.0);

    // dQ/dx above the forward's peak (~77.4 MeV/cm, ~268e3 e/cm) saturates at dedx_max
    const double dE_sat = model.dE(300e3, dx);
    CHECK(dE_sat == doctest::Approx(77.0 * units::MeV).epsilon(1e-12));

    // repeat-call determinism: bit-identical results
    const double dQ = 55e3;
    const double first = model.dE(dQ, dx);
    for (int i = 0; i != 5; ++i) {
        CHECK(model.dE(dQ, dx) == first);
    }
}

TEST_CASE("powerbox recombination configurable")
{
    Gen::PowerBoxRecombination model;
    auto cfg = model.default_configuration();
    CHECK(cfg["A"].asDouble() == doctest::Approx(0.93));
    CHECK(cfg["k"].asDouble() == doctest::Approx(0.282371));
    CHECK(cfg["p"].asDouble() == doctest::Approx(1.362179));
    CHECK(cfg["C"].asDouble() == doctest::Approx(0.855175));
    CHECK(cfg["pivot"].asDouble() == doctest::Approx(2.1));
    CHECK(cfg["Wi"].asDouble() == doctest::Approx(23.6e-6));
    CHECK(cfg["dedx_max"].asDouble() == doctest::Approx(77.0));

    // p=1, k=B/(rho*E)*pivot, C=1 reduces to the plain Modified Box
    Configuration alt;
    const double b_over_re = 0.255 / (1.38 * 0.5);  // sbnd_box_recomb parameters
    alt["p"] = 1.0;
    alt["k"] = b_over_re * 2.1;
    alt["C"] = 1.0;
    alt["A"] = 1.0;
    model.configure(alt);

    Gen::BoxRecombination box(0.5, 1.0, 0.255, 1.38, 23.6e-6);
    const double dx = 1 * units::cm;
    for (double dedx = 1.0; dedx <= 20.0; dedx += 2.5) {
        const double dE_in = dedx * units::MeV / units::cm * dx;
        CHECK(model(dE_in, dx) == doctest::Approx(box(dE_in, dx)).epsilon(1e-9));
    }
}
