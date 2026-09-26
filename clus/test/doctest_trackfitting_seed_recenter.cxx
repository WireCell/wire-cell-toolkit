// doc pdvd/111 -- the default-OFF seed re-centring of the trajectory fit.
//
// seed_recenter_sigma (default 0 = off) moves each organize_orig_path point
// transversely onto the local charge ridge before the first trajectory_fit
// pass, with TrackFittingUtil::recenter_point_transverse.  These cases pin:
//   * the defaults are the legacy path (sigma 0) and the set/get round trip;
//   * a seed 12 mm off a straight track returns to the ridge, and only its
//     transverse coordinates move;
//   * a five-times brighter branch 40 mm away, inside the 3-sigma kernel, does
//     not capture a seed 12 mm off the main track (it converges to the nearest mode);
//   * the no-op cases return the input bit for bit: no points, zero direction,
//     sigma <= 0, points outside the slab, and a move beyond max_move.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/TrackFitting_Util.h"

#include <cmath>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus;
using WireCell::Clus::TrackFittingUtil::recenter_point_transverse;

namespace {
    // a straight track along z at (x, y) = (0, y0): a 1-mm grid over z in [-60, 60] mm and a
    // transverse charge profile of +-2 mm in y (weights 0.3 / 1 / 0.3)
    void add_track(std::vector<Point>& pts, std::vector<double>& w, double y0, double scale)
    {
        for (int iz = -60; iz <= 60; ++iz) {
            for (int dy = -2; dy <= 2; dy += 2) {
                pts.emplace_back(0.0, y0 + dy, static_cast<double>(iz));
                w.push_back(scale * (dy == 0 ? 1.0 : 0.3));
            }
        }
    }
}

TEST_CASE("doc pdvd/111 seed re-centring: defaults are the legacy path")
{
    TrackFitting tf;
    CHECK(tf.get_parameter("seed_recenter_sigma") == 0.0);
    CHECK(tf.get_parameter("seed_recenter_iter") == 5.0);
    CHECK(tf.get_parameter("seed_recenter_max_move") == 30.0);
    TrackFitting::Parameters p;
    CHECK(p.seed_recenter_sigma == 0.0);
}

TEST_CASE("doc pdvd/111 seed re-centring: set_parameter / get_parameter round trip")
{
    TrackFitting tf;
    tf.set_parameter("seed_recenter_sigma", 10.0);
    tf.set_parameter("seed_recenter_iter", 3.0);
    tf.set_parameter("seed_recenter_max_move", 20.0);
    CHECK(tf.get_parameter("seed_recenter_sigma") == 10.0);
    CHECK(tf.get_parameter("seed_recenter_iter") == 3.0);
    CHECK(tf.get_parameter("seed_recenter_max_move") == 20.0);
    CHECK_THROWS_AS(tf.set_parameter("seed_recenter_radius", 10.0), ValueError);
}

TEST_CASE("doc pdvd/111 recenter_point_transverse: an off-ridge seed returns to the ridge, transversely only")
{
    std::vector<Point> pts;
    std::vector<double> w;
    add_track(pts, w, 0.0, 1.0);
    const Point p0(3.0, 12.0, 7.0);
    const Vector dir(0.0, 0.0, 1.0);
    const Point p = recenter_point_transverse(p0, dir, pts, w, 10.0, 6.0, 20, 30.0);
    CHECK(std::abs(p.y()) < 1.0);
    CHECK(std::abs(p.x()) < 1.0);
    CHECK(p.z() == p0.z());      // the along-track coordinate never moves
    // a direction that is not unit length gives the same answer
    const Point q = recenter_point_transverse(p0, dir * 3.0, pts, w, 10.0, 6.0, 20, 30.0);
    CHECK(q.x() == doctest::Approx(p.x()));
    CHECK(q.y() == doctest::Approx(p.y()));
}

TEST_CASE("doc pdvd/111 recenter_point_transverse: a brighter branch 40 mm away does not capture the seed")
{
    std::vector<Point> pts;
    std::vector<double> w;
    add_track(pts, w, 0.0, 1.0);
    add_track(pts, w, 40.0, 5.0);
    const Point p0(0.0, 12.0, 0.0);
    // sigma 10 mm: the branch (28 mm from the seed) is inside the 3-sigma kernel and five times brighter
    const Point p = recenter_point_transverse(p0, Vector(0, 0, 1), pts, w, 10.0, 6.0, 20, 30.0);
    CHECK(std::abs(p.y()) < 1.0);
}

TEST_CASE("doc pdvd/111 recenter_point_transverse: no-op cases return the input exactly")
{
    std::vector<Point> pts;
    std::vector<double> w;
    const Point p0(1.5, -2.5, 3.5);
    // no points
    CHECK(recenter_point_transverse(p0, Vector(0, 0, 1), pts, w, 10.0, 6.0, 5, 30.0) == p0);
    add_track(pts, w, 0.0, 1.0);
    // zero direction, sigma <= 0
    CHECK(recenter_point_transverse(p0, Vector(0, 0, 0), pts, w, 10.0, 6.0, 5, 30.0) == p0);
    CHECK(recenter_point_transverse(p0, Vector(0, 0, 1), pts, w, 0.0, 6.0, 5, 30.0) == p0);
    // every point outside the slab: the seed sits 440 mm beyond the end of the track along the path direction
    const Point far_z(0.0, 5.0, 500.0);
    CHECK(recenter_point_transverse(far_z, Vector(0, 0, 1), pts, w, 10.0, 6.0, 5, 30.0) == far_z);
    // a move larger than max_move keeps the seed
    const Point off(0.0, 25.0, 0.0);
    CHECK(recenter_point_transverse(off, Vector(0, 0, 1), pts, w, 10.0, 6.0, 20, 5.0) == off);
    // mismatched weights
    std::vector<double> w_short(w.begin(), w.end() - 1);
    CHECK(recenter_point_transverse(p0, Vector(0, 0, 1), pts, w_short, 10.0, 6.0, 5, 30.0) == p0);
}
