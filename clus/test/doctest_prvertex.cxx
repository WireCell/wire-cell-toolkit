#include "WireCellClus/PRVertex.h"

#include "WireCellUtil/Logging.h"

#include "WireCellUtil/doctest.h"

#include <iostream>

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("clus pr vertex") {
    PR::Vertex vtx;

    REQUIRE( (vtx.ident() == -1) );
    REQUIRE(vtx.cluster_id() == -1);

    // Test "chainable setters"
    vtx.ident(0).cluster_id(42);
    REQUIRE(vtx.ident() == 0);
    REQUIRE(vtx.cluster_id() == 42);

    REQUIRE(! vtx.fit().valid());

    REQUIRE(! vtx.descriptor_valid());

}
