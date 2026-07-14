#include "WireCellUtil/Point.h"
#include "WireCellUtil/Intersection.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <cmath>
#include <map>

using namespace WireCell;
using namespace std;

TEST_SUITE("intersection2") {

TEST_CASE("D3Vector operations") {
    D3Vector<int> a(3, 4, 5), b(4, 3, 5), c(-5, -12, -13);

    CHECK(a.dot(b) == 49);

    auto axb = a.cross(b);
    CHECK(axb.x() == 5);
    CHECK(axb.y() == 5);
    CHECK(axb.z() == -7);

    spdlog::debug("a . b : {}", a.dot(b));
    spdlog::debug("a . b x c : {}", a.triplescal(b, c));
}

TEST_CASE("box_intersection grid scan") {
    Ray bounds(Point(0, 0, 0), Point(1, 1, 1));
    Vector direction = Point(1, 1, 1).norm();

    for (double x = -1.1; x <= 1; x += 0.5) {
        for (double y = -1.1; y <= 1; y += 0.5) {
            for (double z = -1.0; z <= 1; z += 0.5) {
                Vector point(x, y, z);
                Ray hits(Point(-111, -111, -111), Point(-222, -222, -222));

                int hitmask = box_intersection(bounds, point, direction, hits);
                spdlog::debug("RESULT: {} p={} hits={}", hitmask, point, hits);
                CHECK(hitmask >= 0);
            }
        }
    }
}

}  // TEST_SUITE("intersection2")
