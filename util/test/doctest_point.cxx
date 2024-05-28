#include "WireCellUtil/Point.h"
#include "WireCellUtil/BoundingBox.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Math.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <iomanip>
#include <iostream>

using namespace WireCell;

using spdlog::debug;


TEST_CASE("point3d zero is zero")
{
    Point ppz;
    CHECK(ppz.x() == 0);
    CHECK(ppz.y() == 0);
    CHECK(ppz.z() == 0);
}
//// this "feature" is removed 
// TEST_CASE("point invalid is invalid")
// {
//     Point pp1;
//     pp1.invalidate();
//     CHECK(!pp1);
//     Point pp2 = pp1;
//     CHECK(!pp2);
// }
TEST_CASE("point3d can be in vector")
{
    std::vector<Point> pts(3);
    CHECK(pts[2].size() == 3);
    CHECK(pts[2].x() == 0);
    CHECK(pts[2].y() == 0);
    CHECK(pts[2].z() == 0);
    
    std::vector<Point> pts2;
    pts2 = pts;
    CHECK(pts2[2].size() == 3);
    CHECK(pts2[2].x() == 0);
    CHECK(pts2[2].y() == 0);
    CHECK(pts2[2].z() == 0);
}

TEST_CASE("point3d point set")
{
    Point pp1(1, -500, -495);
    Point pp2(1, 500, -495);
    PointSet results;
    results.insert(pp1);
    results.insert(pp2);
    CHECK(2 == results.size());
}

TEST_CASE("point3d dot product")
{
    Point origin(0, 0, 0);
    Vector zdir(0, 0, 1);
    Point pt(0 * units::mm, 3.92772 * units::mm, 5.34001 * units::mm);
    double dot = zdir.dot(pt - origin);
    debug("origin={}, zdir={}, pt={} dot={}",
          origin / units::mm, zdir, pt / units::mm, dot / units::mm);
}

TEST_CASE("point3d rays")
{
    Point p1(1, 2, 3);
    CHECK(p1.x() == 1);
    CHECK(p1.y() == 2);
    CHECK(p1.z() == 3);

    Point p2 = p1;

    CHECK(p1 == p2);

    Point p3(p1);
    CHECK(p1 == p3);

    PointF pf(p1);
    CHECK(Point(pf) == p1);
    CHECK(pf == PointF(p1));

    Point ps = p1 + p2;
    CHECK(ps.x() == 2);

    CHECK(p1.norm().magnitude() == 1.0);

    double eps = (1 - 1e-11);
    Point peps = p1 * eps;
    PointSet pset;
    pset.insert(p1);
    pset.insert(p2);
    pset.insert(p3);
    pset.insert(ps);
    pset.insert(peps);
    std::stringstream ss;
    for (auto pit = pset.begin(); pit != pset.end(); ++pit) {
        ss << *pit << " ";
    }
    debug(ss.str());
    CHECK(pset.size() == 2);

    BoundingBox bb(pset.begin(), pset.end());
    CHECK(bb.inside(p1));

    Point pdiff = p1;
    pdiff.set(3, 2, 1);
    CHECK(p1 != pdiff);

    const Ray r1(Point(3.75, -323.316, -500), Point(3.75, 254.034, 500));
    const Ray r2(Point(2.5, -254.034, 500), Point(2.5, 323.316, -500));
    const Ray c12 = ray_pitch(r1, r2);
    debug("r1={} r2={} rp={}", r1, r2, c12);
}

