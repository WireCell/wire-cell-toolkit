// Tests of decay spectra and the MARLEY event reader.

#include "WireCellTrig/Kinematics.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Trig;

TEST_CASE("trig beta spectra")
{
    // mean beta energies: 39Ar 218.9 keV, 85Kr 251.4 keV (NNDC)
    BetaSpectrum ar39(0.565 * units::MeV, 19, BetaSpectrum::Shape::unique_first);
    BetaSpectrum kr85(0.687 * units::MeV, 37, BetaSpectrum::Shape::unique_first);
    CHECK(ar39.mean() == doctest::Approx(0.2189 * units::MeV).epsilon(0.05));
    CHECK(kr85.mean() == doctest::Approx(0.2514 * units::MeV).epsilon(0.05));

    ToyTransport::rng_t rng(1);
    double sum = 0;
    const int n = 100000;
    for (int i = 0; i < n; ++i) {
        const double t = ar39(rng);
        REQUIRE(t >= 0);
        REQUIRE(t <= ar39.q());
        sum += t;
    }
    CHECK(sum / n == doctest::Approx(ar39.mean()).epsilon(0.01));
}

TEST_CASE("trig decay branches")
{
    Configuration cfg;
    cfg["branches"][0]["probability"] = 0.8;
    cfg["branches"][0]["beta"]["q"] = 3.5 * units::MeV;
    cfg["branches"][0]["beta"]["z"] = 20;
    cfg["branches"][1]["probability"] = 0.2;
    cfg["branches"][1]["beta"]["q"] = 2.0 * units::MeV;
    cfg["branches"][1]["beta"]["z"] = 20;
    cfg["branches"][1]["gammas"][0] = 1.5 * units::MeV;
    Decay decay(cfg);
    ToyTransport::rng_t rng(2);
    int ngamma = 0;
    const int n = 20000;
    for (int i = 0; i < n; ++i) {
        auto inter = decay(rng);
        for (const auto& p : inter.primaries) {
            if (p.pdg == 22) {
                ++ngamma;
                CHECK(p.kinetic == 1.5 * units::MeV);
            }
            CHECK(std::abs(p.dir.magnitude() - 1) < 1e-9);
        }
    }
    CHECK(ngamma == doctest::Approx(0.2 * n).epsilon(0.05));

    Configuration alpha;
    alpha["branches"][0]["alpha"] = 5.49 * units::MeV;
    Decay rn(alpha);
    auto inter = rn(rng);
    REQUIRE(inter.primaries.size() == 1);
    CHECK(inter.primaries[0].pdg == 1000020040);
    CHECK(inter.energy == 5.49 * units::MeV);
    CHECK(rn.mean_energy() == 5.49 * units::MeV);

    Configuration bad;
    bad["branches"][0]["beta"]["q"] = 1.0;
    bad["branches"][0]["beta"]["z"] = 1;
    bad["branches"][0]["beta"]["shape"] = "second";
    CHECK_THROWS(Decay(bad));
}

TEST_CASE("trig marley events")
{
    // as written by MARLEY 1.2.1 with format "json"
    const std::string text = R"({"events": [{"Ex": 7.952,
      "final_particles": [
        {"E": 4.056, "charge": -1, "mass": 0.511, "pdg": 11, "px": 1.18, "py": -3.48, "pz": -1.65},
        {"E": 939.713, "charge": 0, "mass": 939.565, "pdg": 2112, "px": -15.47, "py": 1.41, "pz": 6.02}],
      "initial_particles": [
        {"E": 13.0037, "charge": 0, "mass": 0, "pdg": 12, "px": -0.69, "py": 2.47, "pz": -12.75},
        {"E": 37224.7, "charge": 0, "mass": 37224.7, "pdg": 1000180400, "px": 0, "py": 0, "pz": 0}]}]})";
    Configuration doc;
    REQUIRE(Json::Reader().parse(text, doc));
    auto events = MarleyEvents::parse(doc);
    REQUIRE(events.size() == 1);
    CHECK(events[0].energy == doctest::Approx(13.0037 * units::MeV));
    REQUIRE(events[0].primaries.size() == 2);
    CHECK(events[0].primaries[0].pdg == 11);
    CHECK(events[0].primaries[0].kinetic == doctest::Approx(3.545 * units::MeV));
    CHECK(events[0].primaries[1].kinetic == doctest::Approx(0.148 * units::MeV));
    CHECK(events[0].primaries[0].dir.magnitude() == doctest::Approx(1));
}
