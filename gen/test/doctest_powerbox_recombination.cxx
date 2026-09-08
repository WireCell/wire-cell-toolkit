// PowerBoxRecombination: the SBND free-power Modified Box fit
// (sbnd_xin/docs/55 sec 7g).  Checks the forward against the reference
// python implementation, the bisection inverse round-trip, saturation
// above the monotone branch, and repeat-call determinism.

#include "WireCellGen/RecombinationModels.h"
#include "WireCellGen/PracticalRecombinationModels.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>

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
    // PowerBox builds u from dE/dx expressed in MeV/cm, while BoxRecombination's
    // quenching term is the plain WCT-units ratio since upstream e6fb7ef3 fixed a
    // 10x unit error there.  The mapping constant must carry that same factor or
    // the p=1 identity no longer holds.  (PowerBox itself is unchanged, and SBND
    // production uses the fitted A/k/p/C, not this mapping -- doc 87 sec 1.)
    alt["k"] = b_over_re * 2.1 * (units::MeV / units::cm);
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

// doc pdhd/16.  PowerBoxRecombination at p = 1 and k = b'*pivot IS the
// Modified Box at drift field E, with C exposed as an explicit normalization
// the plain Box has no slot for.  That equivalence is what lets
// CheckSTM_Michel carry a CALIBRATED inverse without a new C++ model: the PID
// dQ/dx tables the same configs ship are the Modified Box TIMES 0.85
// (energy_loss/pion_travel/convert_field.C:42,71) and PracticalBoxRecombination
// has nowhere to put that factor, so dE() on a measured dQ/dx lands 1/0.85
// away from the scale its own tables were built on.
//
// Both ProtoDUNE operating points are pinned here.  If this test ever fails,
// the two classes have diverged and every C fitted against one of them is void.
TEST_CASE("powerbox at p=1 reproduces the practical Modified Box")
{
    const double A = 0.93, B = 0.212, rho = 1.38, Wi = 23.6e-6, pivot = 2.1;
    const struct { const char* det; double E; double k; } dets[] = {
        {"pdvd", 0.45,   0.7169082125603865},   // b' = 0.212/(1.38*0.45)
        {"pdhd", 0.4959, 0.6505519170239442},   // b' = 0.212/(1.38*0.4959)
    };
    for (const auto& d : dets) {
        const double bp = B / (rho * d.E);
        CHECK(d.k == doctest::Approx(bp * pivot).epsilon(1e-12));

        Gen::PowerBoxRecombination power(A, d.k, 1.0, 1.0, pivot, Wi, 77.0);
        Gen::PracticalBoxRecombination box(d.E, A, B, rho, Wi);

        for (double dedx = 0.8; dedx <= 50.0; dedx += 0.37) {
            const double dx = 0.6 * units::cm;
            const double dE_in = dedx * units::MeV / units::cm * dx;
            const double q_pow = power(dE_in, dx);
            const double q_box = box(dE_in, dx);
            CHECK(q_pow == doctest::Approx(q_box).epsilon(1e-9));
            // and the inverses agree on the same charge
            CHECK(power.dE(q_box, dx) == doctest::Approx(box.dE(q_box, dx)).epsilon(1e-9));
        }

        // The two inverses are NOT equivalent everywhere, and the range above
        // does not show it.  PowerBox saturates at dedx_max = 77 MeV/cm and
        // returns 0 for dQ/dx <= 0; the Box exponential does neither.  Neither
        // branch is reachable through segment_cal_kine_dQdx, whose own
        // [0, 50 MeV/cm] clamp binds first -- pinned here so that stays true.
        const double dx = 1 * units::cm;
        const double q_at_50 = box(50.0 * units::MeV / units::cm * dx, dx);
        const double q_at_77 = box(77.0 * units::MeV / units::cm * dx, dx);
        CHECK(q_at_50 < q_at_77);
        CHECK(power.dE(q_at_77 * 1.5, dx) == doctest::Approx(77.0 * units::MeV / units::cm * dx));
        CHECK(power.dE(-1.0, dx) == 0.0);

        // Replay segment_cal_kine_dQdx's own guards (PRSegmentFunctions.cxx:2483)
        // so the equivalence claim covers the whole input domain, not just the
        // interpolation range above.
        auto through_seg_cal = [&](WireCell::IRecombinationModel& m, double dqdx) {
            double dE = m.dE(dqdx * dx / units::cm, dx);
            if (dE < 0) dE = 0;
            const double ceil = 50 * units::MeV / units::cm * dx;
            if (dE > ceil) dE = ceil;
            return dE;
        };
        // (a) SATURATION IS UNOBSERVABLE.  dQ/dx >= forward(dedx_max) implies
        // dQ/dx > forward(50), so the Box's answer exceeds 50 MeV/cm and gets
        // clamped to 50 -- exactly where the PowerBox's saturated 77 also lands.
        for (double q : {40000.0, 200000.0,
                         q_at_50 / dx * units::cm * 1.001,
                         q_at_77 / dx * units::cm * 1.5, 1e6, 1e7}) {
            CHECK(through_seg_cal(power, q) == doctest::Approx(through_seg_cal(box, q)).epsilon(1e-9));
        }
        // (b) NON-POSITIVE CHARGE IS A REAL DIVERGENCE, and the PowerBox is the
        // one that is right.  The Modified Box inverse has an offset: at
        // dQ/dx = 0 it returns (1 - A)/beta' MeV/cm, ~0.2, because R = ln(A+u)/u
        // sends dQ/dx to -infinity as dE/dx -> 0 and inverting outside that
        // domain is meaningless.  A fit point the charge solve could not read
        // would therefore contribute ~0.2 MeV/cm of energy under the Box and
        // zero under the PowerBox (doc pdhd/16 sec 6.2 prices this on the
        // arms: the muon chain has no such point, and the effect on the Michel
        // object is under a percent).  The window is dQ/dx in
        // (ln(A)/(beta'*Wi), 0]; below it the Box goes negative and the clamp
        // above puts both at zero again.
        const double bp2 = B / (rho * d.E);
        const double q_zero_cross = std::log(A) / (bp2 * Wi);   // ~ -9000 e/cm
        CHECK(q_zero_cross < 0);
        CHECK(through_seg_cal(power, 0.0) == 0.0);
        CHECK(through_seg_cal(box, 0.0) > 0.15 * units::MeV / units::cm * dx);
        CHECK(through_seg_cal(box, 0.0) < 0.25 * units::MeV / units::cm * dx);
        CHECK(through_seg_cal(power, q_zero_cross * 1.01) == 0.0);
        CHECK(through_seg_cal(box, q_zero_cross * 1.01) == 0.0);   // both zero below the crossing
    }
}

// The calibrated PDVD operating point this round ships, pinned so a config
// typo cannot pass unnoticed: at C = 0.7941 a measured MIP-plateau dQ/dx must
// come back near the muon table's own dE/dx rather than ~15 % below it.
TEST_CASE("powerbox calibrated pdvd operating point")
{
    Gen::PowerBoxRecombination cal(0.93, 0.7169082125603865, 1.0, 0.7941, 2.1, 23.6e-6, 77.0);
    Gen::PracticalBoxRecombination raw(0.45, 0.93, 0.212, 1.38, 23.6e-6);
    const double dx = 1 * units::cm;
    const double dqdx_plateau = 54682.0;                 // e/cm, PDVD data, rr 40-60 cm
    const double dQ = dqdx_plateau * dx / units::cm;
    const double raw_dedx = raw.dE(dQ, dx) / (units::MeV / units::cm * dx);
    const double cal_dedx = cal.dE(dQ, dx) / (units::MeV / units::cm * dx);
    // The table's OWN plateau dE/dx, recovered by inverting convert_field.C on
    // its rr = 59.5 cm entry (54657.7 e/cm / 0.85), is 2.193 MeV/cm.  The
    // uncalibrated inverse lands 17 % below it; the calibrated one lands 8 %
    // ABOVE it, and that 8 % is not an error -- C is fixed by the whole-track
    // energy integral against range, so it also absorbs the charge PDVD's
    // reconstruction does not recover (doc pdvd/50 measures that same deficit
    // differentially as k = 0.86-0.95).  C = 0.85 x k, see doc pdhd/16 sec 5.
    CHECK(raw_dedx == doctest::Approx(1.826).epsilon(2e-3));   // 17 % below 2.193
    CHECK(cal_dedx == doctest::Approx(2.377).epsilon(2e-3));   //  8 % above 2.193
}
