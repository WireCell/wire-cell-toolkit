#include "WireCellClus/PRVertex.h"

#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("clus pr vertex basic") {
    PR::Vertex vtx;

    CHECK(! vtx.fit().valid());
    CHECK(! vtx.descriptor_valid());
}

// F8 regression: HasCluster<Vertex> CRTP parameter must be Vertex, not Segment.
// Before the fix, cluster(ptr) called dynamic_cast<Segment*>(this) on a Vertex
// and returned a dangling reference.  After the fix the cast is
// dynamic_cast<Vertex*>(this) which always succeeds, so &ref == &vtx.
TEST_CASE("clus pr vertex HasCluster CRTP") {
    PR::Vertex vtx;

    // cluster(nullptr) stores null and returns *this as Vertex&.
    PR::Vertex& ref = vtx.cluster(nullptr);
    CHECK(&ref == &vtx);   // same object — not a Segment dangling ref
    CHECK(vtx.cluster() == nullptr);

    // Storing a non-null pointer must round-trip through the CRTP setter.
    Facade::Cluster* fake = reinterpret_cast<Facade::Cluster*>(0x1);
    vtx.cluster(fake);
    CHECK(vtx.cluster() == fake);
}
