#include "WireCellClus/PRSegment.h"

#include "WireCellUtil/Logging.h"

#include "WireCellUtil/doctest.h"

#include <iostream>

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("clus pr segment") {
    PR::Segment seg;

    REQUIRE( (seg.ident() == -1) );
    REQUIRE(seg.cluster_id() == -1);

    // Test "chainable setters"
    seg.ident(0).cluster_id(42);
    REQUIRE(seg.ident() == 0);
    REQUIRE(seg.cluster_id() == 42);


    REQUIRE(! seg.descriptor_valid());

}
