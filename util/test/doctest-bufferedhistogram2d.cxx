#include "WireCellUtil/BufferedHistogram2D.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <vector>

using namespace WireCell;
using namespace std;

TEST_SUITE("bufferedhistogram2d") {

TEST_CASE("default construction") {
    BufferedHistogram2D hist;
    CHECK(hist.popx().size() == 0);
    CHECK(false == hist.fill(-1, -1, 0));
    CHECK(hist.xmin() == 0.0);
    CHECK(hist.ymin() == 0.0);
    spdlog::debug("mins: {} {}", hist.xmin(), hist.ymin());
}

TEST_CASE("fill and pop") {
    BufferedHistogram2D hist;
    CHECK(hist.fill(0.5, 3.5));
    vector<double> dat = hist.popx();
    spdlog::debug("dat size {}", dat.size());
    CHECK(dat.size() == 4);
    CHECK(dat[3] == 1.0);
    CHECK(hist.xmin() == 1.0);
    CHECK(hist.ymin() == 0.0);
}

}  // TEST_SUITE("bufferedhistogram2d")
