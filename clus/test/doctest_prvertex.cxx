#include "WireCellClus/PRVertex.h"

#include "WireCellUtil/Logging.h"

#include "WireCellUtil/doctest.h"

#include <iostream>

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("clus pr vertex") {
    PR::Vertex vtx;

    REQUIRE(! vtx.fit().valid());

    REQUIRE(! vtx.descriptor_valid());

}
