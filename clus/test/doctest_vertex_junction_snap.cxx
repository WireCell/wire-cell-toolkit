// doc sbnd_xin/docs/pr/104: unit tests for the pure helpers behind the
// main-vertex junction snap (snap_main_vertex_to_junction):
// vjs_direction_classes (pass-through folding) and vjs_joint_fit (joint
// least-squares line intersection).  Synthetic geometry only -- the full
// pass needs an event; its validation record is the pr/104 A/B arms.

#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus::PR;

TEST_CASE("pr104 vjs_direction_classes: three distinct prongs are three classes")
{
    std::vector<Vector> d = {Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)};
    CHECK(vjs_direction_classes(d, 150) == 3);
}

TEST_CASE("pr104 vjs_direction_classes: a pass-through pair folds to one class")
{
    // 405707 shape at the pre-op0 main vertex: the muon's two halves (170 deg
    // apart) are one track passing through, not two prongs.
    const double a = 170.0 * M_PI / 180.0;
    std::vector<Vector> d = {Vector(0, 0, 1), Vector(std::sin(a), 0, std::cos(a))};
    CHECK(vjs_direction_classes(d, 150) == 1);
    // add a genuine third prong at 90 deg -> two classes
    d.push_back(Vector(1, 0, 0));
    CHECK(vjs_direction_classes(d, 150) == 2);
    // a 120 deg pair is NOT collinear at the 150 deg threshold
    const double b = 120.0 * M_PI / 180.0;
    std::vector<Vector> e = {Vector(0, 0, 1), Vector(std::sin(b), 0, std::cos(b))};
    CHECK(vjs_direction_classes(e, 150) == 2);
}

TEST_CASE("pr104 vjs_direction_classes: empty and null-direction inputs")
{
    std::vector<Vector> none;
    CHECK(vjs_direction_classes(none, 150) == 0);
    std::vector<Vector> z = {Vector(0, 0, 0), Vector(0, 0, 1)};
    CHECK(vjs_direction_classes(z, 150) == 2);  // a null direction never folds
}

TEST_CASE("pr104 vjs_joint_fit: three lines through one point recover it")
{
    const Point J(1.0 * units::cm, -2.0 * units::cm, 3.0 * units::cm);
    std::vector<std::pair<Point, Vector>> lines = {
        {J + Vector(3, 0, 0) * units::cm, Vector(1, 0, 0)},
        {J + Vector(0, 4, 0) * units::cm, Vector(0, 1, 0)},
        {J + Vector(0, 0, 5) * units::cm, Vector(0, 0, 1)},
    };
    Point out;
    double rms = -1;
    REQUIRE(vjs_joint_fit(lines, out, rms));
    CHECK((out - J).magnitude() < 1e-6 * units::cm);
    CHECK(rms < 1e-6 * units::cm);
}

TEST_CASE("pr104 vjs_joint_fit: the 65289 shape -- two vertices 2.4 cm apart, fit lands on the junction")
{
    // M sits 2.4 cm along the muon from the junction J; the muon passes
    // through both, the two protons and the pion all radiate from J.  Lines
    // are anchored where each prong's direction window starts.
    const Point J(0, 0, 0);
    const Vector mu = Vector(0, 0, 1);
    const Point M = J + mu * (2.4 * units::cm);
    std::vector<std::pair<Point, Vector>> lines = {
        {M + mu * (1.0 * units::cm), mu},                                    // muon, anchored at M
        {J + Vector(1, 0, 0) * (1.5 * units::cm), Vector(1, 0, 0)},          // proton 1
        {J + Vector(0, 1, 0) * (1.5 * units::cm), Vector(0, 1, 0)},          // proton 2
        {J + Vector(0.6, -0.8, 0) * (1.5 * units::cm), Vector(0.6, -0.8, 0)} // pion
    };
    Point out;
    double rms = -1;
    REQUIRE(vjs_joint_fit(lines, out, rms));
    CHECK((out - J).magnitude() < 1e-6 * units::cm);
    CHECK((out - M).magnitude() - (out - J).magnitude() > 0.5 * units::cm);  // the tier-B margin
    CHECK(rms < 1.0 * units::cm);
}

TEST_CASE("pr104 vjs_joint_fit: parallel lines are singular")
{
    std::vector<std::pair<Point, Vector>> lines = {
        {Point(0, 0, 0), Vector(0, 0, 1)},
        {Point(1 * units::cm, 0, 0), Vector(0, 0, 1)},
    };
    Point out;
    double rms = -1;
    CHECK_FALSE(vjs_joint_fit(lines, out, rms));
    std::vector<std::pair<Point, Vector>> one = {{Point(0, 0, 0), Vector(1, 0, 0)}};
    CHECK_FALSE(vjs_joint_fit(one, out, rms));
}

TEST_CASE("pr104 vjs_joint_fit: non-intersecting lines report the residual")
{
    // two skew lines 1 cm apart -> solution midway, rms = 0.5 cm
    std::vector<std::pair<Point, Vector>> lines = {
        {Point(0, 0, 0), Vector(1, 0, 0)},
        {Point(0, 0, 1 * units::cm), Vector(0, 1, 0)},
    };
    Point out;
    double rms = -1;
    REQUIRE(vjs_joint_fit(lines, out, rms));
    CHECK(std::abs(out.z() - 0.5 * units::cm) < 1e-6 * units::cm);
    CHECK(std::abs(rms - 0.5 * units::cm) < 1e-6 * units::cm);
}
