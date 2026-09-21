// doc sbnd_xin/119 round 3 -- hoisting the per-(apa,face) constants out of the per-point
// projection.
//
// The round claims a pure CPU saving with byte-identical output.  That claim rests on three
// arithmetic identities, and this file is what makes them evidence rather than assertion:
//
//   * point2wind_cs(p, cos(a), sin(a), pitch, center) reproduces the LEGACY point2wind body
//     exactly -- the same int, for every input, not merely "to within a wire".  A value landing
//     on a .5 boundary is the case that would break it, so the sweep drives one there on purpose.
//   * point2wind_cont_cs likewise reproduces the legacy point2wind_cont body BIT for bit.
//   * drift2time(xsign, xorig, ...) reproduces the legacy IAnodeFace body's expression bit for
//     bit, with xsign/xorig supplied instead of re-derived.
//
// The legacy bodies are written out here literally, rather than called, so that this test keeps
// checking the round's claim even after the production functions are refactored again.  If a
// future change makes the hoisted and the derived form disagree by one ulp, these cases fail
// rather than an A/B gate failing three hours into an arm.
//
// (Grouping::fastgeom's memo of these constants is covered by the output gate, not here: filling
// it needs a real IAnodeFace, and the arithmetic it feeds is what is tested below.)

#include "WireCellUtil/doctest.h"
#include "WireCellClus/Facade_Util.h"
#include "WireCellClus/Graphs.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <set>

using namespace WireCell;
using namespace WireCell::Clus;
using WireCell::Clus::Facade::geo_point_t;

namespace {

    // The bodies exactly as they stood before doc sbnd_xin/119 round 3.
    int legacy_point2wind(const geo_point_t& point, const double angle, const double pitch,
                          const double center)
    {
        double y = cos(angle) * point[2] - sin(angle) * point[1];
        double wind = (y - center) / pitch - 0.5;
        return std::round(wind);
    }
    double legacy_point2wind_cont(const geo_point_t& point, const double angle, const double pitch,
                                  const double center)
    {
        double y = cos(angle) * point[2] - sin(angle) * point[1];
        return (y - center) / pitch - 0.5;
    }
    double legacy_drift2time(const double xsign, const double xorig, const double time_offset,
                             const double drift_speed, const double drift)
    {
        return (drift - xorig) / (xsign * drift_speed) - time_offset;
    }

    // Wire angles and pitches spanning the three detectors this round is measured on: SBND
    // (+-60 deg, 3 mm), PDHD (+-35.7 deg, ~4.7/5.0 mm) and PDVD.  Units are WCT's internal ones,
    // so pitch is in mm.
    struct Plane { double angle, pitch, center; };
    const std::vector<Plane> planes = {
        {  60.0 * M_PI / 180.0, 3.0,    0.0 },
        { -60.0 * M_PI / 180.0, 3.0, -100.0 },
        {   0.0,                3.0,  250.0 },
        {  35.7 * M_PI / 180.0, 4.669,  12.5 },
        { -35.7 * M_PI / 180.0, 4.669, -12.5 },
        {  90.0 * M_PI / 180.0, 4.790, 300.0 },
        {  30.0 * M_PI / 180.0, 5.100, -7.25 },
    };

}  // namespace

TEST_CASE("doc sbnd_xin/119 r3: point2wind_cs reproduces the legacy body exactly")
{
    std::mt19937_64 rng(20260921);
    std::uniform_real_distribution<double> uni(-4000.0, 4000.0);

    size_t n = 0;
    for (const auto& pl : planes) {
        const double ca = cos(pl.angle), sa = sin(pl.angle);
        for (int i = 0; i < 20000; ++i) {
            const geo_point_t p(uni(rng), uni(rng), uni(rng));
            CHECK(Facade::point2wind_cs(p, ca, sa, pl.pitch, pl.center)
                  == legacy_point2wind(p, pl.angle, pl.pitch, pl.center));
            // bit-for-bit, not approximately
            CHECK(Facade::point2wind_cont_cs(p, ca, sa, pl.pitch, pl.center)
                  == legacy_point2wind_cont(p, pl.angle, pl.pitch, pl.center));
            ++n;
        }
        // The angle-taking forms still exist for callers with no memo; they must agree too.
        const geo_point_t q(12.5, -33.25, 807.125);
        CHECK(Facade::point2wind(q, pl.angle, pl.pitch, pl.center)
              == legacy_point2wind(q, pl.angle, pl.pitch, pl.center));
        CHECK(Facade::point2wind_cont(q, pl.angle, pl.pitch, pl.center)
              == legacy_point2wind_cont(q, pl.angle, pl.pitch, pl.center));
    }
    CHECK(n == 20000 * planes.size());
}

TEST_CASE("doc sbnd_xin/119 r3: point2wind_cs agrees at the rounding boundary")
{
    // The one input class where an ulp of disagreement would change the answer: a continuous
    // coordinate sitting exactly on .5, and its immediate neighbours.  Construct points that land
    // there by inverting the projection on the W plane (angle 0 => y == point[2]).
    const double pitch = 3.0, center = 0.0, angle = 0.0;
    const double ca = cos(angle), sa = sin(angle);
    for (int k = -500; k <= 500; ++k) {
        for (double frac : {0.5, 0.5 - 1e-12, 0.5 + 1e-12, 0.0, 0.25, 0.75}) {
            const double z = (k + frac + 0.5) * pitch + center;
            const geo_point_t p(0.0, 0.0, z);
            CHECK(Facade::point2wind_cs(p, ca, sa, pitch, center)
                  == legacy_point2wind(p, angle, pitch, center));
        }
    }
}

TEST_CASE("doc sbnd_xin/119 r3: drift2time with supplied constants is the legacy expression")
{
    std::mt19937_64 rng(20260921);
    std::uniform_real_distribution<double> xd(-4000.0, 4000.0);
    std::uniform_real_distribution<double> od(-5000.0, 5000.0);

    // dirx() is +-1; drift speed and tick offsets in WCT internal units.
    for (double xsign : {1.0, -1.0}) {
        for (double speed : {1.6e-3, 1.076e-3, 0.8e-3}) {
            for (int i = 0; i < 20000; ++i) {
                const double xorig = od(rng), toff = od(rng), drift = xd(rng);
                CHECK(Facade::drift2time(xsign, xorig, toff, speed, drift)
                      == legacy_drift2time(xsign, xorig, toff, speed, drift));
            }
        }
    }
}

TEST_CASE("doc sbnd_xin/119 r3: find_neighbors_nlevel returns the std::set order")
{
    // The BFS used to return std::set<vertex_type>; it now returns a vector, and TrackFitting
    // walks that vector to build a blob set whose ITERATION ORDER reaches the output.  So the
    // vector must be ascending and unique -- the raw BFS discovery order would not be.
    Graphs::Weighted::Graph g(12);
    // A deliberately awkward shape: high-numbered vertices discovered before low-numbered ones,
    // and a cycle, so discovery order and ascending order differ.
    const std::vector<std::pair<size_t, size_t>> edges = {
        {5, 11}, {5, 9}, {5, 2}, {11, 0}, {9, 3}, {2, 7}, {7, 1}, {0, 3}, {1, 4}, {4, 8}, {8, 10},
    };
    for (auto [a, b] : edges) boost::add_edge(a, b, 1.0, g);

    Graphs::Weighted::GraphAlgorithms ga(g);

    for (int nlevel : {0, 1, 2, 3, 4}) {
        for (bool self : {true, false}) {
            auto got = ga.find_neighbors_nlevel(5, nlevel, self);
            // ascending
            CHECK(std::is_sorted(got.begin(), got.end()));
            // unique
            CHECK(std::adjacent_find(got.begin(), got.end()) == got.end());
            // and exactly the membership a std::set would have held
            std::set<size_t> want(got.begin(), got.end());
            CHECK(want.size() == got.size());
            CHECK(std::equal(want.begin(), want.end(), got.begin()));
            // self-inclusion honoured
            const bool has_self = std::find(got.begin(), got.end(), 5u) != got.end();
            if (nlevel >= 0) CHECK(has_self == self);
        }
    }

    // Spot-check the membership itself so a broken BFS cannot pass by being sorted.
    auto one = ga.find_neighbors_nlevel(5, 1, false);
    CHECK(one == std::vector<Graphs::Weighted::vertex_type>{2, 9, 11});
    auto two = ga.find_neighbors_nlevel(5, 2, false);
    CHECK(two == std::vector<Graphs::Weighted::vertex_type>{0, 2, 3, 7, 9, 11});
    // An out-of-range seed is still empty, and a negative nlevel too.
    CHECK(ga.find_neighbors_nlevel(99, 2).empty());
    CHECK(ga.find_neighbors_nlevel(5, -1).empty());
}
