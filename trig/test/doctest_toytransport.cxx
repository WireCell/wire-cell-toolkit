// Tests of the toy particle transport.

#include "WireCellTrig/ToyTransport.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>

using namespace WireCell;
using namespace WireCell::Trig;

namespace {
    double sum_energy(const std::vector<Deposit>& deps)
    {
        double sum = 0;
        for (const auto& d : deps) sum += d.energy;
        return sum;
    }
}  // namespace

TEST_CASE("trig toy transport attenuation")
{
    ToyTransport tt;
    const double per_cm2g = units::cm2 / units::g * tt.config().density;
    // Klein-Nishina for argon (Z/A = 0.4506) at 1 MeV: water's 0.0707 cm2/g
    // scaled by the electron density per mass, 0.0574 cm2/g.
    CHECK(tt.compton_mu(1 * units::MeV) / per_cm2g == doctest::Approx(0.0574).epsilon(0.01));
    // Compton falls with energy, pair rises above threshold.
    CHECK(tt.compton_mu(5 * units::MeV) < tt.compton_mu(1 * units::MeV));
    CHECK(tt.pair_mu(1 * units::MeV) == 0);
    CHECK(tt.pair_mu(10 * units::MeV) > tt.pair_mu(5 * units::MeV));
    // the attenuation length of MeV gammas is around 10 to 30 cm
    const double lambda1 = 1 / tt.total_mu(1 * units::MeV), lambda5 = 1 / tt.total_mu(5 * units::MeV);
    CHECK(lambda1 > 10 * units::cm);
    CHECK(lambda1 < 15 * units::cm);
    CHECK(lambda5 > 20 * units::cm);
    CHECK(lambda5 < 30 * units::cm);
    // photoelectric dominates well below 100 keV
    CHECK(tt.photo_mu(30 * units::keV) > tt.compton_mu(30 * units::keV));
    // a MIP-like electron loses about 2 MeV/cm
    const double dedx = tt.electron_dedx(2 * units::MeV) / (units::MeV / units::cm);
    CHECK(dedx > 1.8);
    CHECK(dedx < 2.1);
}

TEST_CASE("trig toy transport electron range and conservation")
{
    ToyTransport::Config cfg;
    cfg.multiple_scattering = false;
    cfg.brems_cut = 1e9 * units::MeV;  // no brems
    ToyTransport tt(cfg);
    ToyTransport::rng_t rng(1);
    std::vector<Deposit> deps;
    const double T = 1 * units::MeV;
    auto tally = tt(rng, 11, T, Vector(1, 0, 0), Point(0, 0, 0), deps);
    CHECK(tally.deposited == doctest::Approx(T));
    CHECK(sum_energy(deps) == doctest::Approx(T));
    CHECK(tally.lost == 0);

    // straight: the end point is the path length, about the CSDA range
    double path = 0, xmax = 0;
    for (const auto& d : deps) {
        path += d.length;
        xmax = std::max(xmax, d.pos.x());
    }
    double csda = 0;
    const int n = 10000;
    for (int i = 0; i < n; ++i) {
        const double t = cfg.electron_cut + (T - cfg.electron_cut) * (i + 0.5) / n;
        csda += (T - cfg.electron_cut) / n / tt.electron_dedx(t);
    }
    CHECK(path == doctest::Approx(csda).epsilon(0.02));
    CHECK(path > 3 * units::mm);  // ESTAR argon CSDA at 1 MeV is about 0.5 g/cm2
    CHECK(path < 5 * units::mm);
}

TEST_CASE("trig toy transport showers conserve energy")
{
    ToyTransport tt;
    ToyTransport::rng_t rng(2);
    for (int pdg : {11, -11, 22}) {
        for (double T : {0.5 * units::MeV, 5 * units::MeV, 40 * units::MeV}) {
            std::vector<Deposit> deps;
            auto tally = tt(rng, pdg, T, Vector(0, 0, 1), Point(0, 0, 0), deps);
            const double expected = T + (pdg == -11 ? 2 * 0.51099895 * units::MeV : 0);
            CHECK(tally.deposited == doctest::Approx(expected));
            CHECK(sum_energy(deps) == doctest::Approx(expected));
            for (const auto& d : deps) CHECK(d.energy > 0);
        }
    }
}

TEST_CASE("trig toy transport gamma spread")
{
    ToyTransport tt;
    ToyTransport::rng_t rng(3);
    // energy weighted mean distance of deposits from a 2 MeV gamma's origin
    double sum = 0, wsum = 0;
    for (int i = 0; i < 2000; ++i) {
        std::vector<Deposit> deps;
        tt(rng, 22, 2 * units::MeV, Vector(0, 0, 1), Point(0, 0, 0), deps);
        for (const auto& d : deps) {
            sum += d.energy * d.pos.magnitude();
            wsum += d.energy;
        }
    }
    const double mean = sum / wsum;
    CHECK(mean > 10 * units::cm);
    CHECK(mean < 30 * units::cm);
}

TEST_CASE("trig toy transport ions and neutrals")
{
    ToyTransport tt;
    ToyTransport::rng_t rng(4);
    std::vector<Deposit> deps;
    auto tally = tt(rng, 1000020040, 5.5 * units::MeV, Vector(0, 0, 1), Point(0, 0, 0), deps);
    REQUIRE(deps.size() == 1);
    CHECK(deps[0].energy == doctest::Approx(5.5 * units::MeV));
    CHECK(deps[0].length < 0.1 * units::mm);  // tens of microns
    CHECK(tally.lost == 0);

    deps.clear();
    tally = tt(rng, 2112, 3 * units::MeV, Vector(0, 0, 1), Point(0, 0, 0), deps);
    CHECK(deps.empty());
    CHECK(tally.lost == doctest::Approx(3 * units::MeV));
}
