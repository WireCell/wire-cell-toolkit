// doc sbnd_xin/docs/pr/51 -- main-vertex graph audit: unit tests for the
// pure decision geometry (path_overlap_fraction, the op1/op3 gate input).
// Synthetic point paths only -- the full pass edits a PR graph and needs an
// event; its validation record is the pr/51 A/B arms and on-census.

#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <vector>

using namespace WireCell;
using namespace WireCell::Clus::PR;

namespace {

    // A straight path along z at fixed (x, y), n points, step cm apart.
    std::vector<Point> straight(double x, double y, int n, double step = 0.6*units::cm)
    {
        std::vector<Point> pts;
        for (int i = 0; i < n; ++i) {
            pts.emplace_back(x, y, i * step);
        }
        return pts;
    }

    // The pr/51 defaults, restated (C++ member defaults in
    // NeutrinoPatternBase.h).
    const double DUP_TOL = 1.4*units::cm;
    const double DUP_FRAC = 0.7;

}  // namespace

TEST_CASE("pr51 path_overlap_fraction: identical paths overlap fully")
{
    auto a = straight(0, 0, 20);
    CHECK(path_overlap_fraction(a, a, DUP_TOL) == doctest::Approx(1.0));
}

TEST_CASE("pr51 path_overlap_fraction: the 360535 parallel-pair geometry")
{
    // Two parallel ribbons ~1 cm apart over one corridor: 0% at a 0.6 cm
    // tolerance, ~100% at the 1.4 cm default (the measured 360535 pair reads
    // 0/77-80% at 0.6/1.4 cm -- the default must clear the fitter's ~1 cm
    // ribbon separation, which 1.2 cm did not in the first on-census).
    auto a = straight(0, 0, 20);
    auto b = straight(1.0*units::cm, 0, 20);
    CHECK(path_overlap_fraction(a, b, 0.6*units::cm) == doctest::Approx(0.0));
    CHECK(path_overlap_fraction(a, b, DUP_TOL) == doctest::Approx(1.0));
    CHECK(path_overlap_fraction(a, b, DUP_TOL) >= DUP_FRAC);
}

TEST_CASE("pr51 path_overlap_fraction: distant paths do not overlap")
{
    auto a = straight(0, 0, 20);
    auto b = straight(10*units::cm, 0, 20);
    CHECK(path_overlap_fraction(a, b, DUP_TOL) == doctest::Approx(0.0));
}

TEST_CASE("pr51 path_overlap_fraction: partial rider (268067-15003 shape)")
{
    // A short path riding the first half of a longer one, then diverging:
    // the overlapped half counts, the divergent half does not.
    auto corridor = straight(0, 0, 40);
    std::vector<Point> rider;
    for (int i = 0; i < 10; ++i) {
        rider.emplace_back(0.2*units::cm, 0, i * 0.6*units::cm);  // on-corridor
    }
    for (int i = 10; i < 20; ++i) {
        rider.emplace_back(0.2*units::cm + (i - 9) * 1.0*units::cm, 0, i * 0.6*units::cm);  // walks away
    }
    double frac = path_overlap_fraction(rider, corridor, DUP_TOL);
    CHECK(frac > 0.45);
    CHECK(frac < 0.75);
}

TEST_CASE("pr51 path_overlap_fraction: asymmetric by construction")
{
    // Shorter-onto-longer is the op1 convention: the short rider scores
    // high, the long corridor against the short rider scores low.
    auto corridor = straight(0, 0, 40);
    auto rider = straight(0.2*units::cm, 0, 8);
    CHECK(path_overlap_fraction(rider, corridor, DUP_TOL) == doctest::Approx(1.0));
    CHECK(path_overlap_fraction(corridor, rider, DUP_TOL) < 0.5);
}

TEST_CASE("pr51 path_overlap_fraction: degenerate inputs")
{
    std::vector<Point> empty;
    auto a = straight(0, 0, 5);
    CHECK(path_overlap_fraction(empty, a, DUP_TOL) == doctest::Approx(0.0));
    CHECK(path_overlap_fraction(a, empty, DUP_TOL) == doctest::Approx(0.0));
    CHECK(path_overlap_fraction(empty, empty, DUP_TOL) == doctest::Approx(0.0));
    // Zero tolerance: only exact coincidences count.
    CHECK(path_overlap_fraction(a, a, 0.0) == doctest::Approx(1.0));
}
