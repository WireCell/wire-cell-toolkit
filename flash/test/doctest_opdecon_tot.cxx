// OpDecon saturation_repair_mode "tot" (time-over-threshold fill).
//
// The reference numbers are from the Python prototype the C++ ports,
// pdvd/docs/qlmatch/scripts/saturation_tot_study.py (kernel / model / Shape /
// methods_pe tot_fill), evaluated on the synthetic template below by
// pdvd/docs/qlmatch/scripts/d31_tot_doctest_ref.py (doc qlmatch/31).

#include "WireCellUtil/doctest.h"

#include "WireCellFlash/OpDecon.h"

#include <cmath>
#include <vector>

using namespace WireCell;
using Flash::OpDecon;

namespace {
    // 120-sample template: exp rise (tau 1.5) to a peak of 100 at sample 10,
    // exp fall (tau 20) after it.  Stored as float, as OpDecon holds it.
    std::vector<float> synth_template()
    {
        std::vector<float> w(120);
        for (int t = 0; t < 120; ++t) {
            w[t] = t < 10 ? (float) (100 * std::exp((t - 10) / 1.5))
                          : (float) (100 * std::exp(-(t - 10) / 20.0));
        }
        return w;
    }
    const double tau_fall = 20.000000062737683;   // Channel.tau_fall of the template
    const std::vector<double> par{0.25, 0.4, 6.0, 90.0, 0.8};

    bool rel_close(double a, double b, double tol = 1e-6)
    {
        return std::abs(a - b) <= tol * std::max(std::abs(b), 1e-300);
    }
}

TEST_CASE("opdecon tot model and width table match the python prototype")
{
    const auto y = OpDecon::tot_model(synth_template(), tau_fall, par);
    REQUIRE(y.size() == 1600);
    CHECK(rel_close(y[0], 6.4087546879110576e-06));
    CHECK(rel_close(y[5], 0.16234818503100848));
    CHECK(rel_close(y[10], 4.6616416534914178));
    CHECK(rel_close(y[13], 28.52858912950618));
    CHECK(rel_close(y[40], 23.894997988354589));
    CHECK(rel_close(y[200], 1.7049952076199291));
    CHECK(rel_close(y[700], 0.0065729265699159605));
    CHECK(rel_close(y[1599], 3.0174456457639008e-07, 1e-5));

    const auto shp = OpDecon::tot_shape(y);
    REQUIRE(shp.lam.size() == 2999);
    CHECK(shp.s[19] == doctest::Approx(1.0));
    CHECK(rel_close(shp.lam[0], 0.001));
    CHECK(rel_close(shp.lam[2998], 0.99769929780441036));
    CHECK(shp.u[0] == 4);
    CHECK(shp.u[1000] == 7);
    CHECK(shp.u[2000] == 10);
    CHECK(shp.u[2998] == 19);
    CHECK(rel_close(shp.w[0], 534.3531120179797));
    CHECK(rel_close(shp.w[1000], 323.92992857026297));
    CHECK(rel_close(shp.w[2000], 117.51216395092177));
    CHECK(rel_close(shp.w[2998], 2.459902036979146));

    CHECK(rel_close(OpDecon::tot_level_for(shp, 5.5), 0.98349139430880206));
    CHECK(rel_close(OpDecon::tot_level_for(shp, 50.3), 0.32175785785273325));
    CHECK(rel_close(OpDecon::tot_level_for(shp, 300.7), 0.012900528940129604));
}

TEST_CASE("opdecon tot level solves width == ToT")
{
    const auto shp = OpDecon::tot_shape(OpDecon::tot_model(synth_template(), tau_fall, par));
    for (double tot : {3.0, 12.0, 40.0, 88.0, 176.0, 400.0}) {
        const double lam = OpDecon::tot_level_for(shp, tot);
        // width of the table at the solved level, interpolated back
        size_t k = 0;
        while (k + 1 < shp.lam.size() && shp.lam[k + 1] <= lam) ++k;
        const double f = (lam - shp.lam[k]) / (shp.lam[k + 1] - shp.lam[k]);
        const double wid = shp.w[k] + f * (shp.w[k + 1] - shp.w[k]);
        CHECK(wid == doctest::Approx(tot).epsilon(1e-3));
    }
}

TEST_CASE("opdecon tot fill: clipped noise-free pulse round trip")
{
    const auto shp = OpDecon::tot_shape(OpDecon::tot_model(synth_template(), tau_fall, par));
    const double ped = 1000.0, rail = 16383.0;
    // python: (d, i, j, filled_sum, full area ratio)
    struct Ref { double d; int i, j; double filled_sum, area_ratio; };
    for (const Ref& r : {Ref{2.0, 113, 146, 776674.17987625906, 0.984988},
                         Ref{6.7, 111, 199, 4162090.673335663, 0.988099},
                         Ref{20.0, 109, 285, 14804279.792909553, 0.988008}}) {
        const double A = r.d * (rail - ped);
        std::vector<double> full(2000, ped);
        for (int k = 0; k < 1600; ++k) full[100 + k] += A * shp.s[k];
        std::vector<float> w(2000);
        double true_area = 0;
        for (int t = 0; t < 2000; ++t) {
            true_area += full[t] - ped;
            w[t] = (float) std::min(full[t], rail);
        }
        int i = 0;
        while (w[i] < rail) ++i;
        int j = i;
        while (j < 2000 && w[j] >= rail) ++j;
        REQUIRE(i == r.i);
        REQUIRE(j == r.j);
        const auto before = w;
        REQUIRE(OpDecon::tot_fill_run(w, i, j, ped, rail, shp));
        double filled = 0, area = 0;
        for (int t = 0; t < 2000; ++t) {
            area += w[t] - ped;
            if (t >= i && t < j) {
                filled += w[t] - ped;
                CHECK(w[t] >= before[t]);          // fill never lowers a sample
            }
            else {
                CHECK(w[t] == before[t]);          // outside the run untouched
            }
        }
        // float storage of the waveform vs python float64: 1e-5 relative
        CHECK(rel_close(filled, r.filled_sum, 1e-5));
        CHECK(area / true_area == doctest::Approx(r.area_ratio).epsilon(1e-4));
    }
}

TEST_CASE("opdecon tot fill refuses a run wider than the shape table")
{
    const auto shp = OpDecon::tot_shape(OpDecon::tot_model(synth_template(), tau_fall, par));
    std::vector<float> w(1000, 16383.0f);
    const auto before = w;
    CHECK_FALSE(OpDecon::tot_fill_run(w, 10, 900, 1000.0, 16383.0, shp));   // 890 > widest 534
    CHECK(w == before);
}
