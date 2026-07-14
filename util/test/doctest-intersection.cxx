#include "WireCellUtil/Intersection.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <random>

using namespace WireCell;
using namespace std;

TEST_SUITE("intersection") {

TEST_CASE("explicit line in bounds") {
    Ray bounds(Point(-5, -500, -500), Point(5, 500, 500));
    Ray line(Point(1, 0, -495), Point(1, 1, -495));
    Vector dir = ray_unit(line);
    Ray hits;
    int hitmask = box_intersection(bounds, line.first, dir, hits);
    spdlog::debug("line: {} bounds:{}", line, bounds);
    spdlog::debug("got: {} {}", hitmask, hits);
    CHECK(hitmask == 3);
}

TEST_CASE("random points") {
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<> dist(-2, 2);

    Ray bounds(Point(-1, -1, -1), Point(1, 1, 1));
    for (int ind = 0; ind < 100; ++ind) {
        Point p1(dist(re), dist(re), dist(re));
        Point p2(dist(re), dist(re), dist(re));
        Vector dir = (p2 - p1).norm();

        spdlog::debug("{} point={} dir={}", ind, p1, dir);

        for (int axis = 0; axis < 3; ++axis) {
            Ray hits;
            int hitmask = box_intersection(axis, bounds, p1, dir, hits);
            spdlog::debug("\t axis={} [{}]hits={}", axis, hitmask, hits);
            CHECK(hitmask >= 0);
        }

        {
            Ray ray(p1, p2), hits;
            int hitmask = box_intersection(bounds, p1, dir, hits);
            spdlog::debug("box: hitmask={} ray={} hits={}", hitmask, ray, hits);
            if (point_contained(p1, bounds) && point_contained(p2, bounds)) {
                CHECK_MESSAGE(hitmask == 3, "Inside box, but not enough hits");
            }
        }
    }
}

}  // TEST_SUITE("intersection")
