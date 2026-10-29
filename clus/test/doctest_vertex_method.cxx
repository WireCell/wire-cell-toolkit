// vertex_method (TaggerCheckNeutrino): the names the configuration accepts and
// the placement constant of dv2_split_move.  The placements themselves need an
// event; their validation record is the in-chain fold-0 run of the branch.

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellUtil/doctest.h"

#include <string>

using namespace WireCell::Clus::PR;
using VM = PatternAlgorithms::VertexMethod;

TEST_CASE("vertex_method: every accepted name parses and round-trips")
{
    for (const std::string name : {"original", "dv2_snap", "dv2_split", "dv2_split_move"}) {
        VM m = VM::original;
        REQUIRE(PatternAlgorithms::parse_vertex_method(name, m));
        CHECK(std::string(PatternAlgorithms::vertex_method_name(m)) == name);
    }
}

TEST_CASE("vertex_method: unknown and retired names are rejected and leave the method unchanged")
{
    // "scn" / "dv2max" were the values of the retired dl_vtx_engine knob.
    for (const std::string name : {"", "scn", "dv2max", "dv2", "snap", "split", "split_move", "Original", "dv2_split+move"}) {
        VM m = VM::dv2_split;
        CHECK_FALSE(PatternAlgorithms::parse_vertex_method(name, m));
        CHECK(m == VM::dv2_split);
    }
}

TEST_CASE("vertex_method: the default is the original chain, and only dv2_* methods are DV2-MAX")
{
    PatternAlgorithms pa;
    CHECK(pa.m_vertex_method == VM::original);
    CHECK_FALSE(pa.is_dv2());
    CHECK_FALSE(pa.m_dv2_target_valid);
    for (VM m : {VM::dv2_snap, VM::dv2_split, VM::dv2_split_move}) {
        pa.m_vertex_method = m;
        CHECK(pa.is_dv2());
    }
}

TEST_CASE("vertex_method: dv2_split_move reuses a PR-graph vertex within 1.5 cm (the MyFCN move limit)")
{
    CHECK(PatternAlgorithms::kDv2MoveReuseDistCm == doctest::Approx(1.5));
}
