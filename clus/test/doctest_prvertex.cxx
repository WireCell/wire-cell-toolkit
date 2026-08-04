#include "WireCellClus/PRVertex.h"

#include "WireCellUtil/Units.h"
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

// doc sbnd_xin/docs/pr/32 §11 F4 (was P12) and F1 (was P1).
//
// F4: kMainCandidate is the toolkit's stand-in for the prototype's
// map_cluster_main_candidate_vertices.  It must occupy its own bit -- setting it
// must not disturb kNeutrinoVertex, because MultiAlgBlobClustering and
// PrDisplayDump both branch on kNeutrinoVertex alone and would mis-render every
// candidate as the neutrino vertex if the bits collided.
//
// F1: fit_distance() is the exact size of the divergence -- the gap between the
// continuous fit that the prototype's get_fit_pt() returns and the Steiner-node
// snap that the toolkit reads at eleven expressions.  Pin that it is the
// fit-to-wcpt distance and that an unfitted vertex reports valid()==false, which
// is what makes the wcpt fallback load-bearing.
TEST_CASE("clus pr32 F4 kMainCandidate is independent of kNeutrinoVertex") {
    PR::Vertex vtx;

    CHECK(! vtx.flags_any(PR::VertexFlags::kMainCandidate));
    CHECK(! vtx.flags_any(PR::VertexFlags::kNeutrinoVertex));

    vtx.set_flags(PR::VertexFlags::kMainCandidate);
    CHECK(vtx.flags_any(PR::VertexFlags::kMainCandidate));
    CHECK(! vtx.flags_any(PR::VertexFlags::kNeutrinoVertex));

    // The main vertex is also a candidate; both bits must coexist.
    vtx.set_flags(PR::VertexFlags::kNeutrinoVertex);
    CHECK(vtx.flags_any(PR::VertexFlags::kMainCandidate));
    CHECK(vtx.flags_any(PR::VertexFlags::kNeutrinoVertex));

    vtx.unset_flags(PR::VertexFlags::kMainCandidate);
    CHECK(! vtx.flags_any(PR::VertexFlags::kMainCandidate));
    CHECK(vtx.flags_any(PR::VertexFlags::kNeutrinoVertex));
}

TEST_CASE("clus pr32 F1 fit_distance is the fit-to-wcpt gap") {
    PR::Vertex vtx;
    vtx.wcpt().point = Point(10*units::cm, 0, 0);

    // No fit yet: the fallback in pr32_vtx_pt is what stops fit().point's
    // default (0,0,0) reaching a score ladder.
    CHECK(! vtx.fit().valid());
    CHECK(vtx.fit().point == Point(0, 0, 0));

    // After a fit, fit_distance() is exactly |fit - wcpt| -- the quantisation
    // onto the Steiner lattice that F1 removes.
    vtx.fit().point = Point(10.5*units::cm, 0, 0);
    vtx.fit_index(0);
    CHECK(vtx.fit().valid());
    CHECK(vtx.fit_distance() == doctest::Approx(0.5*units::cm));

    // Fit::reset() clears the index but NOT the point -- the property that makes
    // the third form_map_graph re-validate an already-fitted vertex.
    vtx.fit().reset();
    CHECK(! vtx.fit().valid());
    CHECK(vtx.fit().point == Point(10.5*units::cm, 0, 0));
}
