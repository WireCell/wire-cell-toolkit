#include "WireCellUtil/Point.h"
#include "WireCellUtil/Intersection.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <vector>

using namespace WireCell;
using namespace std;

TEST_SUITE("intersection3") {

TEST_CASE("z-direction ray, reversed bounds") {
    // note, make this "upside-down" (from max to min) to assure
    // box_intersection() handles mismatch order.
    Ray bounds(Point(10,10,10), Point(0,0,0));

    const Vector dir(0,0,1);
    std::vector<Point> points = {
        Point(5,5,-5),
        Point(5,5,5),
        Point(5,5,15)
    };
    for (const auto& point : points) {
        Ray hits;
        int hm = box_intersection(bounds, point, dir, hits);
        CHECK(hm == 3);
        CHECK(hits.first.x() == point.x());
        CHECK(hits.second.x() == point.x());
        CHECK(hits.first.y() == point.y());
        CHECK(hits.second.y() == point.y());
        CHECK(hits.first.z() == 0);
        CHECK(hits.second.z() == 10);
        auto hdir = hits.second - hits.first;
        CHECK(hdir.dot(dir) > 0);
    }
}

TEST_CASE("pierce corners diagonal") {
    Ray bounds(Point(0,0,0), Point(10,10,10));
    Vector dir(1,1,1);
    dir = dir.norm();
    Point point(5,5,5);
    Ray hits;
    int hm = box_intersection(bounds, point, dir, hits);
    CHECK(hm == 3);
    CHECK(hits.first == bounds.first);
    CHECK(hits.second == bounds.second);
    auto hdir = hits.second - hits.first;
    CHECK(hdir.dot(dir) > 0);
}

TEST_CASE("degenerate cases dir 1 1 -1") {
    Ray bounds(Point(0,0,0), Point(10,10,10));
    Vector dir(1,1,-1);
    dir = dir.norm();

    SUBCASE("point at corner") {
        Point point(0,0,0);
        Ray hits;
        int hm = box_intersection(bounds, point, dir, hits);
        spdlog::debug("hm = {}", hm);
        bool hm_valid = (hm == 1 || hm == 2);
        CHECK(hm_valid);  // FP err could change order
        CHECK(hits.first == point);
        CHECK(hits.first == hits.second);
    }

    SUBCASE("point at bounds.first + dir") {
        Point point = bounds.first + dir;
        Ray hits;
        int hm = box_intersection(bounds, point, dir, hits);
        spdlog::debug("hm = {}", hm);
        CHECK(hm == 1);
        CHECK(hits.second == point);
        auto hdir = hits.second - hits.first;
        CHECK(hdir.dot(dir) > 0);
    }

    SUBCASE("point at bounds.first - dir") {
        Point point = bounds.first - dir;
        Ray hits;
        int hm = box_intersection(bounds, point, dir, hits);
        spdlog::debug("hm = {}", hm);
        CHECK(hm == 2);
        CHECK(hits.first == point);
        auto hdir = hits.second - hits.first;
        CHECK(hdir.dot(dir) > 0);
    }

    SUBCASE("point at bounds.second + dir") {
        Point point = bounds.second + dir;
        Ray hits;
        int hm = box_intersection(bounds, point, dir, hits);
        spdlog::debug("hm = {}", hm);
        CHECK(hm == 1);
        CHECK(hits.second == point);
        auto hdir = hits.second - hits.first;
        CHECK(hdir.dot(dir) > 0);
    }

    SUBCASE("point at bounds.second - dir") {
        Point point = bounds.second - dir;
        Ray hits;
        int hm = box_intersection(bounds, point, dir, hits);
        spdlog::debug("hm = {}", hm);
        CHECK(hm == 2);
        CHECK(hits.first == point);
        auto hdir = hits.second - hits.first;
        CHECK(hdir.dot(dir) > 0);
    }
}

}  // TEST_SUITE("intersection3")
