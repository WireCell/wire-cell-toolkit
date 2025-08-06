#include "WireCellClus/PRSegment.h"

#include "WireCellUtil/Logging.h"

#include "WireCellUtil/doctest.h"

#include <iostream>

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("clus pr segment") {
    PR::Segment seg;

    REQUIRE(! seg.descriptor_valid());

}
