// doc sbnd_xin/docs/pr/51 round 7: unit tests for the robust-vertex-fit
// geometry helpers (mvfit_annulus_pca / mvfit_rout_dyn / mvfit_leg_disagrees)
// behind MyFCN's per-leg dynamic direction-window substitution.  Synthetic
// point paths only -- the full pass needs an event; its validation record is
// the pr/51 round-7 A/B arms.

#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus::PR;

namespace {

    // Straight run from `start` along `dir` (unit), `n` points at `step`.
    std::vector<Point> straight(Point start, Vector dir, int n, double step)
    {
        std::vector<Point> pts;
        for (int i = 0; i < n; i++) pts.push_back(start + dir * (i * step));
        return pts;
    }

    // Elbow at the vertex end: from the origin (the vertex seat) run along
    // dir1 for arm1_len, then turn into dir2 -- the 234638 shape, where the
    // first ~4-5 cm (the UpdateInfo re-seat) points somewhere the far body
    // does not.
    std::vector<Point> elbow(Vector dir1, double arm1_len, Vector dir2, double arm2_len, double step)
    {
        auto pts = straight(Point(0, 0, 0), dir1, static_cast<int>(arm1_len / step) + 1, step);
        const Point corner = pts.back();
        const int n2 = static_cast<int>(arm2_len / step);
        for (int i = 1; i <= n2; i++) pts.push_back(corner + dir2 * (i * step));
        return pts;
    }

    // The pr/51 round-7 C++ defaults (NeutrinoPatternBase.h m_mvfit_*).
    const double RIN = 6.0 * units::cm;    // reseat 4 cm + margin 2 cm
    const double ROUT_FRAC = 0.5;
    const double ROUT_MIN = 9.0 * units::cm;
    const double ROUT_MAX = 18.0 * units::cm;
    const double ANGLE = 20;
    const int MIN_PTS = 5;
    const double MIN_ANISO = 3.0;

    const double STEP = 0.5 * units::cm;
    const Point VTX(0, 0, 0);

    double fold_deg(const Point& a, const Point& b)
    {
        double c = std::abs(a.dot(b)) / (a.magnitude() * b.magnitude());
        if (c > 1.0) c = 1.0;
        return std::acos(c) * 180.0 / M_PI;
    }

}  // namespace

TEST_CASE("mvfit annulus pca reproduces the production window math")
{
    // A straight leg: the (1.5, 6] production shell and any outer shell must
    // both recover the leg axis exactly, with the (0.15 cm)^2 floor on the
    // transverse eigenvalues.
    const Vector dir(0, 0.6, 0.8);
    auto pts = straight(VTX, dir, 61, STEP);  // 30 cm

    const auto inner = mvfit_annulus_pca(pts, VTX, 1.5 * units::cm, 6.0 * units::cm);
    CHECK(inner.npts == 10);  // 1.5..6 cm inclusive at 0.5 cm steps
    // center = kept point closest to the vertex
    CHECK(inner.center.magnitude() == doctest::Approx(1.5 * units::cm));
    CHECK(fold_deg(inner.axes[0], Point(dir)) < 1e-6);
    // floor: a perfect line's transverse eigenvalues are exactly the floor
    CHECK(inner.vals[1] == doctest::Approx(std::pow(0.15 * units::cm, 2)));
    CHECK(inner.vals[2] == doctest::Approx(std::pow(0.15 * units::cm, 2)));

    const auto outer = mvfit_annulus_pca(pts, VTX, RIN, 15.0 * units::cm);
    CHECK(outer.npts == 19);
    CHECK(outer.center.magnitude() == doctest::Approx(6.0 * units::cm));
    CHECK(fold_deg(outer.axes[0], Point(dir)) < 1e-6);
    // longer lever => larger leading eigenvalue => stronger fit weight
    CHECK(outer.vals[0] > inner.vals[0]);
}

TEST_CASE("mvfit straight leg never fires")
{
    auto pts = straight(VTX, Vector(0, 0, 1), 61, STEP);
    const auto inner = mvfit_annulus_pca(pts, VTX, 1.5 * units::cm, 6.0 * units::cm);
    const auto outer = mvfit_annulus_pca(pts, VTX, RIN,
                                         mvfit_rout_dyn((pts.front() - pts.back()).magnitude(),
                                                        ROUT_FRAC, ROUT_MIN, ROUT_MAX));
    CHECK_FALSE(mvfit_leg_disagrees(inner.axes[0], outer, ANGLE, MIN_PTS, MIN_ANISO));
}

TEST_CASE("mvfit 234638-shape elbow fires and recovers the body axis")
{
    // Inner ~5 cm arm bent 30 deg away from the outer body: the fit-visible
    // window reads mostly the bent arm, the outer window reads the body.
    const Vector d1(0, std::sin(30.0 * M_PI / 180.0), std::cos(30.0 * M_PI / 180.0));
    const Vector d2(0, 0, 1);
    auto pts = elbow(d1, 5.0 * units::cm, d2, 25.0 * units::cm, STEP);

    const auto inner = mvfit_annulus_pca(pts, VTX, 1.5 * units::cm, 6.0 * units::cm);
    const double chord = (pts.front() - pts.back()).magnitude();
    const auto outer = mvfit_annulus_pca(pts, VTX, RIN,
                                         mvfit_rout_dyn(chord, ROUT_FRAC, ROUT_MIN, ROUT_MAX));

    CHECK(mvfit_leg_disagrees(inner.axes[0], outer, ANGLE, MIN_PTS, MIN_ANISO));
    // the outer axis is the body axis, not the bent-arm axis
    CHECK(fold_deg(outer.axes[0], Point(d2)) < 10.0);
    CHECK(fold_deg(inner.axes[0], Point(d2)) > 15.0);
}

TEST_CASE("mvfit folded decision immune to point order and axis sign")
{
    const Vector d1(0, std::sin(30.0 * M_PI / 180.0), std::cos(30.0 * M_PI / 180.0));
    auto pts = elbow(d1, 5.0 * units::cm, Vector(0, 0, 1), 25.0 * units::cm, STEP);
    auto rpts = pts;
    std::reverse(rpts.begin(), rpts.end());

    const double chord = (pts.front() - pts.back()).magnitude();
    const double rout = mvfit_rout_dyn(chord, ROUT_FRAC, ROUT_MIN, ROUT_MAX);

    const auto in_a = mvfit_annulus_pca(pts, VTX, 1.5 * units::cm, 6.0 * units::cm);
    const auto out_a = mvfit_annulus_pca(pts, VTX, RIN, rout);
    const auto in_b = mvfit_annulus_pca(rpts, VTX, 1.5 * units::cm, 6.0 * units::cm);
    const auto out_b = mvfit_annulus_pca(rpts, VTX, RIN, rout);

    // identical geometry regardless of traversal order
    CHECK(mvfit_leg_disagrees(in_a.axes[0], out_a, ANGLE, MIN_PTS, MIN_ANISO) ==
          mvfit_leg_disagrees(in_b.axes[0], out_b, ANGLE, MIN_PTS, MIN_ANISO));
    // and an explicit sign flip of the inner axis changes nothing (folded)
    const Point flipped(-in_a.axes[0].x(), -in_a.axes[0].y(), -in_a.axes[0].z());
    CHECK(mvfit_leg_disagrees(flipped, out_a, ANGLE, MIN_PTS, MIN_ANISO) ==
          mvfit_leg_disagrees(in_a.axes[0], out_a, ANGLE, MIN_PTS, MIN_ANISO));
}

TEST_CASE("mvfit shower fan fails the anisotropy gate")
{
    // A cone-like fan: straight trunk to 6 cm then points spraying +-45 deg
    // in y around the z axis -- the outer window is broad, sqrt(l0/l1) small.
    std::mt19937 rng(4242);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    auto pts = straight(VTX, Vector(0, 0, 1), 13, STEP);  // trunk to 6 cm
    for (int i = 0; i < 40; i++) {
        const double r = 6.0 * units::cm + (i % 20) * 0.5 * units::cm;
        const double dy = u(rng) * 0.9 * r;
        pts.push_back(Point(0, dy, std::sqrt(r * r - dy * dy)));
    }
    // make the elbow visible to the angle gate: inner trunk bent 30 deg
    const auto inner = mvfit_annulus_pca(pts, VTX, 1.5 * units::cm, 6.0 * units::cm);
    const auto outer = mvfit_annulus_pca(pts, VTX, RIN, ROUT_MAX);
    CHECK(outer.npts >= MIN_PTS);
    CHECK(std::sqrt(outer.vals[0] / outer.vals[1]) < MIN_ANISO);
    CHECK_FALSE(mvfit_leg_disagrees(inner.axes[0], outer, ANGLE, MIN_PTS, MIN_ANISO));
}

TEST_CASE("mvfit sparse outer window never fires")
{
    // Only 3 points beyond 6 cm: below MIN_PTS regardless of angle.
    const Vector d1(0, 1, 0);
    auto pts = elbow(d1, 5.0 * units::cm, Vector(0, 0, 1), 1.5 * units::cm, STEP);
    const auto inner = mvfit_annulus_pca(pts, VTX, 1.5 * units::cm, 6.0 * units::cm);
    const auto outer = mvfit_annulus_pca(pts, VTX, RIN, ROUT_MAX);
    CHECK(outer.npts < MIN_PTS);
    CHECK_FALSE(mvfit_leg_disagrees(inner.axes[0], outer, ANGLE, MIN_PTS, MIN_ANISO));
}

TEST_CASE("mvfit dynamic outer radius clamps")
{
    CHECK(mvfit_rout_dyn(10 * units::cm, 0.5, ROUT_MIN, ROUT_MAX) == doctest::Approx(ROUT_MIN));
    CHECK(mvfit_rout_dyn(30 * units::cm, 0.5, ROUT_MIN, ROUT_MAX) == doctest::Approx(15 * units::cm));
    CHECK(mvfit_rout_dyn(50 * units::cm, 0.5, ROUT_MIN, ROUT_MAX) == doctest::Approx(ROUT_MAX));
    CHECK(mvfit_rout_dyn(200 * units::cm, 0.5, ROUT_MIN, ROUT_MAX) == doctest::Approx(ROUT_MAX));
}

TEST_CASE("mvfit degenerate inputs are safe")
{
    // empty path, single point, all points inside rin
    std::vector<Point> none;
    CHECK(mvfit_annulus_pca(none, VTX, RIN, ROUT_MAX).npts == 0);
    auto one = straight(VTX, Vector(0, 0, 1), 1, STEP);
    CHECK(mvfit_annulus_pca(one, VTX, RIN, ROUT_MAX).npts == 0);
    auto close = straight(VTX, Vector(0, 0, 1), 8, STEP);  // to 3.5 cm only
    const auto pc = mvfit_annulus_pca(close, VTX, RIN, ROUT_MAX);
    CHECK(pc.npts == 0);
    // an unusable outer window never fires
    const Point zaxis(0, 0, 1);
    CHECK_FALSE(mvfit_leg_disagrees(zaxis, pc, ANGLE, MIN_PTS, MIN_ANISO));
}
