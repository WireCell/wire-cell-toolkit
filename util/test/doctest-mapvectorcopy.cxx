#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <cmath>
#include <map>
#include <vector>

typedef std::vector<float> sequence;
typedef std::map<int, sequence> enumerated;

using namespace std;

void mutate(enumerated& en)
{
    en[0][0] = 42.0;
    en[1][0] = 6.9;
}

bool equal(double a, double b)
{
    if (a == b) {
        return true;
    }
    return std::abs(a - b) / (a + b) < 0.0001;
}

TEST_SUITE("map vector copy") {

TEST_CASE("map stores vector by value") {
    sequence vec{1.0, 2.0, 3.0};
    enumerated stuff;
    stuff[0] = vec;
    stuff[1] = {10.0, 20.0, 30.0};
    mutate(stuff);

    for (auto iv : stuff) {
        spdlog::debug("{}:", iv.first);
        for (auto x : iv.second) {
            spdlog::debug(" {}", x);
        }
    }

    // map entries were mutated
    CHECK(equal(stuff[0][0], 42.0));
    CHECK(equal(stuff[1][0], 6.9));
    // original vec is unmodified (map stored a copy, not a reference)
    CHECK(equal(vec[0], 1.0));
}

}  // TEST_SUITE("map vector copy")
