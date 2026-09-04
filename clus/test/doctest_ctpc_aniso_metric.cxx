// doc pdvd/36: the lattice-normalised (anisotropic) ctpc radius query of
// clus/inc/WireCellClus/CtpcAnisoMetric.h, exercised on a bare NFKDVec tree
// with no detector geometry.
//
// What is pinned here:
//   1. yscale >= 1 is the legacy isotropic query EXACTLY -- same index set,
//      same squared distances, same existence answer -- on a lattice shaped
//      like SBND's (pitch finer than the drift step => s clamps to 1).
//   2. For s < 1 (PDVD U/V shape, 2.96 x 7.65 mm) the returned set equals a
//      brute-force ellipse test over the whole lattice for a grid of query
//      phases, including the corner phases where the isotropic circle and the
//      ellipse disagree most; and exists <=> !radius.empty().
//   3. The trap the design is built around: a point whose Euclidean-nearest
//      lattice neighbour lies OUTSIDE the ellipse while a Euclidean-farther
//      neighbour lies inside it.  A filter on a nearest-only query would say
//      "no"; the circumscribe-and-filter query says "yes".
//   4. ctpc_yscale clamps at 1 and returns 1 (the identity) on 0 / NaN input,
//      so a grouping with no nticks_live_slice metadata can never produce a
//      nonsense metric.
//
// The tree type is exactly Grouping::kd2d_t (NFKDVec::Tree<double,
// IndexStatic>) so the template is instantiated on what production uses.

#include "WireCellClus/CtpcAnisoMetric.h"
#include "WireCellUtil/NFKDVec.h"
#include "WireCellUtil/doctest.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus::Facade;

namespace {

    using kd_t = NFKDVec::Tree<double, NFKDVec::IndexStatic>;
    using query_t = std::array<double, 2>;

    struct lattice_t {
        double dx, dy;               // cell edges: drift step, pitch
        int nx, ny;                  // extent in cells
        std::vector<double> xs, ys;  // flat columns, index = ix*ny + iy
        kd_t::points_type points() const { return {xs, ys}; }
    };

    lattice_t make_lattice(double dx, double dy, int nx, int ny)
    {
        lattice_t L{dx, dy, nx, ny, {}, {}};
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                L.xs.push_back(ix * dx);
                L.ys.push_back(iy * dy);
            }
        }
        return L;
    }

    // Brute-force reference: every lattice point inside the ellipse.
    std::set<size_t> brute_ellipse(const lattice_t& L, const query_t& q, double r, double s)
    {
        std::set<size_t> out;
        for (size_t i = 0; i < L.xs.size(); ++i) {
            const double dx = L.xs[i] - q[0];
            const double dy = (L.ys[i] - q[1]) * s;
            if (dx * dx + dy * dy < r * r) out.insert(i);
        }
        return out;
    }

    std::set<size_t> index_set(const kd_t::results_type& res)
    {
        std::set<size_t> out;
        for (const auto& [ind, d2] : res) out.insert(ind);
        return out;
    }

}  // namespace

TEST_CASE("ctpc aniso metric: yscale derivation clamps at 1 and is the identity on bad input")
{
    // PDVD: drift step 2.9615 mm, pitch 7.65 (U/V) and 5.10 (W)
    CHECK(ctpc_yscale(2.9615, 7.65) == doctest::Approx(2.9615 / 7.65));
    CHECK(ctpc_yscale(2.9615, 5.10) == doctest::Approx(2.9615 / 5.10));
    // SBND: 3.126 mm drift step, 3.000 mm pitch => raw 1.042 => clamped
    CHECK(ctpc_yscale(3.126, 3.000) == 1.0);
    // exactly commensurate
    CHECK(ctpc_yscale(3.0, 3.0) == 1.0);
    // absent metadata (nticks 0), or garbage => identity, never a nonsense scale
    CHECK(ctpc_yscale(0.0, 7.65) == 1.0);
    CHECK(ctpc_yscale(2.96, 0.0) == 1.0);
    CHECK(ctpc_yscale(-2.96, 7.65) == 1.0);
    CHECK(ctpc_yscale(std::numeric_limits<double>::quiet_NaN(), 7.65) == 1.0);
    CHECK(ctpc_yscale(2.96, std::numeric_limits<double>::infinity()) == 1.0);
}

TEST_CASE("ctpc aniso metric: yscale 1 is the legacy isotropic query, exactly")
{
    // SBND-shaped lattice: pitch (3.0) finer than the drift step (3.126)
    const auto L = make_lattice(3.126, 3.000, 12, 12);
    kd_t kd(L.points());
    REQUIRE(kd.npoints() == L.xs.size());

    const double s = ctpc_yscale(L.dx, L.dy);
    REQUIRE(s == 1.0);

    // A grid of query phases across one cell, several radii.
    for (double r : {1.0, 2.0, 2.2, 3.0, 6.0}) {
        for (int i = 0; i < 7; ++i) {
            for (int j = 0; j < 7; ++j) {
                const query_t q{5 * L.dx + (i / 7.0) * L.dx, 5 * L.dy + (j / 7.0) * L.dy};
                const auto legacy = kd.radius<query_t>(r * r, q);
                const auto aniso = ctpc_radius_aniso(kd, q, r, s);
                // identical sequence: same indices, same squared distances, same order
                REQUIRE(aniso.size() == legacy.size());
                for (size_t k = 0; k < legacy.size(); ++k) {
                    CHECK(aniso[k].first == legacy[k].first);
                    CHECK(aniso[k].second == legacy[k].second);
                }
                CHECK(ctpc_exists_within_aniso(kd, q, r, s) == kd.exists_within<query_t>(r * r, q));
            }
        }
    }
}

TEST_CASE("ctpc aniso metric: PDVD-shaped lattice matches the brute-force ellipse")
{
    // PDVD U/V: drift step 2.9615 mm, pitch 7.65 mm => s = 0.387
    const auto L = make_lattice(2.9615, 7.65, 14, 8);
    kd_t kd(L.points());
    const double s = ctpc_yscale(L.dx, L.dy);
    REQUIRE(s < 1.0);
    REQUIRE(s == doctest::Approx(0.3871).epsilon(1e-3));

    int nonempty = 0, differs_from_circle = 0;
    // 1.0 mm is below dx/2, so some phases must come back empty; 2.09 mm is
    // dx/sqrt(2), the full-coverage bound of doc 34 sec 5, and above it every
    // phase must be non-empty.
    for (double r : {1.0, 2.0, 2.09, 3.0, 6.0}) {  // 0.1 cm, 0.2 cm, dx/sqrt(2), 0.3 cm, 0.6 cm
        for (int i = 0; i < 9; ++i) {
            for (int j = 0; j < 9; ++j) {
                // sweep one full cell of phase, centred in the lattice interior
                const query_t q{6 * L.dx + (i / 9.0) * L.dx, 3 * L.dy + (j / 9.0) * L.dy};
                const auto want = brute_ellipse(L, q, r, s);
                const auto got = ctpc_radius_aniso(kd, q, r, s);
                CHECK(index_set(got) == want);
                // the returned distance is the ISOTROPIC squared distance
                for (const auto& [ind, d2] : got) {
                    const double ddx = L.xs[ind] - q[0], ddy = L.ys[ind] - q[1];
                    CHECK(d2 == doctest::Approx(ddx * ddx + ddy * ddy));
                }
                // existence agrees with the enumeration
                CHECK(ctpc_exists_within_aniso(kd, q, r, s) == !want.empty());
                if (!want.empty()) ++nonempty;
                if (want != index_set(kd.radius<query_t>(r * r, q))) ++differs_from_circle;
            }
        }
    }
    // The sweep must actually exercise both the empty and the non-empty
    // answers and must include phases where the ellipse and the circle
    // disagree, or the checks above could pass on a degenerate sweep.
    CHECK(nonempty > 0);
    CHECK(nonempty < 5 * 81);
    CHECK(differs_from_circle > 0);
}

TEST_CASE("ctpc aniso metric: the Euclidean-nearest point can be outside the ellipse while a farther point is inside")
{
    // The trap the design is built around.  s = 0.387 (PDVD U/V), r = 2 mm.
    //   A at (2.1, 0.0): displaced along DRIFT.  Euclidean 2.1 (> r), scaled
    //                    2.1 (> r) => outside the circle AND outside the ellipse.
    //   C at (1.2, 4.0): displaced mostly along PITCH.  Euclidean 4.18 (> r),
    //                    scaled sqrt(1.44 + (0.387*4)^2) = sqrt(3.84) = 1.96
    //                    (< r) => outside the circle, INSIDE the ellipse.
    // A is the Euclidean-nearest.  A "nearest-then-filter" query would test A,
    // reject it, and answer "no"; the circumscribe-and-filter query finds C.
    const double s = 0.387;
    const double r = 2.0;
    kd_t::points_type pts = {{2.1, 1.2}, {0.0, 4.0}};
    kd_t kd(pts);
    const query_t q{0.0, 0.0};

    size_t nearest = 99;
    double nd2 = 0;
    REQUIRE(kd.knn1(q, nearest, nd2));
    CHECK(nearest == 0);                                    // A is the nearest ...
    CHECK(!ctpc_in_ellipse(2.1, 0.0, 0.0, 0.0, r, s));       // ... and outside the ellipse
    CHECK(ctpc_in_ellipse(1.2, 4.0, 0.0, 0.0, r, s));        // C is inside the ellipse
    CHECK(kd.exists_within<query_t>(r * r, q) == false);    // the circle of radius r is empty
    CHECK(kd.radius<query_t>(r * r, q).empty());            // (isotropic answer: nothing)
    CHECK(ctpc_exists_within_aniso(kd, q, r, s) == true);    // the two-level query finds C
    CHECK(index_set(ctpc_radius_aniso(kd, q, r, s)) == std::set<size_t>{1});
}

TEST_CASE("ctpc aniso metric: empty tree and zero radius are safe")
{
    kd_t empty(kd_t::points_type{{}, {}});
    const query_t q{0.0, 0.0};
    CHECK(ctpc_radius_aniso(empty, q, 2.0, 0.387).empty());
    CHECK(ctpc_exists_within_aniso(empty, q, 2.0, 0.387) == false);

    kd_t::points_type pts = {{0.0}, {0.0}};
    kd_t one(pts);
    CHECK(ctpc_radius_aniso(one, q, 0.0, 0.387).empty());
    CHECK(ctpc_exists_within_aniso(one, q, 0.0, 0.387) == false);
    // a point exactly at the query is inside any positive radius
    CHECK(ctpc_exists_within_aniso(one, q, 1e-6, 0.387) == true);
}
