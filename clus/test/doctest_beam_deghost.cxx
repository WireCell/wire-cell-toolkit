// doc pdvd/125: the geometric link between a deghosted cell and the production
// blobs (same face, overlapping ticks and wires in all three views) and the
// rule that gives a cell to one cluster, on plain integers.

#include "WireCellClus/BeamDeghostFunctions.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell::Clus::PR;

static DeghostBox box(int wpid, int t0, int t1, int u0, int u1, int v0, int v1, int w0, int w1, int cid = -1)
{
    DeghostBox b;
    b.wpid = wpid; b.t0 = t0; b.t1 = t1;
    b.u0 = u0; b.u1 = u1; b.v0 = v0; b.v1 = v1; b.w0 = w0; b.w1 = w1;
    b.cluster_id = cid;
    return b;
}

TEST_CASE("beam_deghost_overlap: half-open ranges, all four must intersect")
{
    const auto a = box(8, 100, 104, 10, 20, 30, 40, 50, 60);
    CHECK(beam_deghost_overlap(a, a) == 4LL * 10 * 10 * 10);
    CHECK(beam_deghost_overlap(a, box(8, 102, 106, 15, 19, 39, 45, 59, 61)) == 2LL * 4 * 1 * 1);
    CHECK(beam_deghost_overlap(a, box(16, 100, 104, 10, 20, 30, 40, 50, 60)) == 0);   // other face
    CHECK(beam_deghost_overlap(a, box(8, 104, 108, 10, 20, 30, 40, 50, 60)) == 0);    // next slice: touching is not overlap
    CHECK(beam_deghost_overlap(a, box(8, 100, 104, 20, 24, 30, 40, 50, 60)) == 0);    // U touching
    CHECK(beam_deghost_overlap(a, box(8, 100, 104, 10, 20, 40, 44, 50, 60)) == 0);    // V touching
    CHECK(beam_deghost_overlap(a, box(8, 100, 104, 10, 20, 30, 40, 60, 64)) == 0);    // W touching
}

TEST_CASE("beam_deghost_assign: largest summed overlap, tie to the smaller cluster id")
{
    std::vector<DeghostBox> owners = {
        box(8, 100, 104, 0, 16, 0, 16, 0, 16, 7),
        box(8, 100, 104, 16, 32, 0, 16, 0, 16, 3),
        box(8, 104, 108, 0, 16, 0, 16, 0, 16, 7),
        box(16, 100, 104, 0, 16, 0, 16, 0, 16, 9),
        box(8, 2000, 2004, 0, 4000, 0, 4000, 0, 4000, 5),   // far in time
    };
    std::vector<DeghostBox> cells = {
        box(8, 100, 104, 4, 8, 4, 8, 4, 8),       // inside cluster 7
        box(8, 100, 104, 14, 18, 0, 4, 0, 4),     // 2 wires in 7, 2 wires in 3: tie -> 3
        box(8, 100, 104, 13, 17, 0, 4, 0, 4),     // 3 wires in 7, 1 in 3 -> 7
        box(8, 100, 104, 40, 44, 0, 4, 0, 4),     // outside in U
        box(8, 108, 112, 0, 4, 0, 4, 0, 4),       // no owner in that slice
        box(16, 100, 104, 0, 4, 0, 4, 0, 4),      // other face -> 9
        box(24, 100, 104, 0, 4, 0, 4, 0, 4),      // face without owners
        box(8, 102, 106, 0, 4, 0, 4, 0, 4),       // spans two slices of cluster 7
    };
    const auto got = beam_deghost_assign(cells, owners);
    REQUIRE(got.size() == cells.size());
    CHECK(got[0] == 7);
    CHECK(got[1] == 3);
    CHECK(got[2] == 7);
    CHECK(got[3] == -1);
    CHECK(got[4] == -1);
    CHECK(got[5] == 9);
    CHECK(got[6] == -1);
    CHECK(got[7] == 7);
    CHECK(beam_deghost_assign({}, owners).empty());
    CHECK(beam_deghost_assign(cells, {}) == std::vector<int>(cells.size(), -1));
}

TEST_CASE("beam_deghost_assign: an owner longer than the others is still found")
{
    // owner 0 starts well before the cell and still covers it; the per-face
    // maximum span bounds the look-back.
    std::vector<DeghostBox> owners = {
        box(8, 0, 400, 0, 16, 0, 16, 0, 16, 1),
        box(8, 200, 204, 100, 116, 0, 16, 0, 16, 2),
        box(8, 396, 400, 100, 116, 0, 16, 0, 16, 2),
    };
    const auto got = beam_deghost_assign({box(8, 396, 400, 0, 4, 0, 4, 0, 4)}, owners);
    CHECK(got[0] == 1);
}

// The visitor's knob defaults: the window is off (low >= high => no bundle is
// picked, nothing is replaced) and a bundle cluster with no cell is dropped.
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

TEST_CASE("clus knob defaults: ClusteringBeamDeghost and PointTreeConcat")
{
    using namespace WireCell;
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("ClusteringBeamDeghost", "doc125_knobdefaults_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    CHECK(cfg["grouping"].asString() == "live");
    CHECK(cfg["deghost_grouping"].asString() == "deghost");
    CHECK(cfg["correction_name"].asString() == "T0Correction");
    CHECK(cfg["beam_window_low"].asDouble() == 0.0);
    CHECK(cfg["beam_window_high"].asDouble() == 0.0);
    CHECK(cfg["drop_empty"].asBool() == true);
    auto jcfg = Factory::lookup<IConfigurable>("PointTreeConcat", "doc125_knobdefaults_probe");
    REQUIRE(jcfg);
    CHECK(jcfg->default_configuration()["multiplicity"].asInt() == 2);
}
