// Tests for the "nu_band_veto_role" per-blob provenance mechanism (doc pr/66):
// ClusteringNeutrino's protect_iso_band_xext veto (record_band_veto knob) can
// stamp a refused pair as per-blob provenance, and every downstream
// merge_clusters() call (i.e. every all-APA merge stage) honors it via
// band_veto_forbids() -- an edge-level veto, dropped before
// boost::connected_components, so a third unmarked cluster can still bridge
// the two.
//
// Covers, in order:
//  * the predicate truth table (band_veto_forbids), standalone;
//  * the carry_singles registry inside merge_clusters(): a member lacking the
//    array contributes zeros, one that has it is carried verbatim;
//  * byte-identity: no writer anywhere => the merged cluster carries no
//    "nu_band_veto_role" key at all;
//  * the edge filter drops the vetoed EDGE only, not the whole connected
//    component -- a third cluster can still transitively bridge the two.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/ClusteringFuncs.h"

#include <algorithm>
#include <unordered_map>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;
using fa_float_t = WireCell::Clus::Facade::float_t;
using fa_int_t = WireCell::Clus::Facade::int_t;

// A cluster with nb BARE blob children -- no "scalar" PC at all.  Safe for
// any test that never triggers the veto's drop-log line (which calls
// get_length()/center_pos(), both of which require a real "scalar" PC on
// every blob) -- i.e. every test below except the last, which uses
// make_geom_cluster() instead.
static Cluster* make_cluster(Grouping& g, int nb)
{
    Cluster& cl = g.make_child();
    for (int i = 0; i < nb; ++i) cl.make_child();
    return &cl;
}

// A cluster whose blobs carry a minimal but COMPLETE "scalar" PC (every key
// Blob::fill_cache() reads), with all wire-index/slice ranges EMPTY
// (min==max) so Cluster::get_length() never reaches into Grouping's
// (here-absent) detector-volume metadata -- it returns 0 safely.  Needed only
// for the test that actually exercises the veto's drop-log line.
static Blob& add_geom_blob(Cluster& cl, double cx, double cy, double cz)
{
    Blob& b = cl.make_child();
    Dataset scalar({
        {"charge", Array(std::vector<fa_float_t>{0})},
        {"center_x", Array(std::vector<fa_float_t>{(fa_float_t) cx})},
        {"center_y", Array(std::vector<fa_float_t>{(fa_float_t) cy})},
        {"center_z", Array(std::vector<fa_float_t>{(fa_float_t) cz})},
        {"wpid", Array(std::vector<fa_int_t>{0})},
        {"npoints", Array(std::vector<fa_int_t>{0})},
        {"slice_index_min", Array(std::vector<fa_int_t>{0})},
        {"slice_index_max", Array(std::vector<fa_int_t>{0})},
        {"u_wire_index_min", Array(std::vector<fa_int_t>{0})},
        {"u_wire_index_max", Array(std::vector<fa_int_t>{0})},
        {"v_wire_index_min", Array(std::vector<fa_int_t>{0})},
        {"v_wire_index_max", Array(std::vector<fa_int_t>{0})},
        {"w_wire_index_min", Array(std::vector<fa_int_t>{0})},
        {"w_wire_index_max", Array(std::vector<fa_int_t>{0})},
        {"max_wire_interval", Array(std::vector<fa_int_t>{0})},
        {"min_wire_interval", Array(std::vector<fa_int_t>{0})},
        {"max_wire_type", Array(std::vector<fa_int_t>{0})},
        {"min_wire_type", Array(std::vector<fa_int_t>{0})},
    });
    b.value().local_pcs()["scalar"] = scalar;
    return b;
}

static Cluster* make_geom_cluster(Grouping& g, int nb)
{
    Cluster& cl = g.make_child();
    for (int i = 0; i < nb; ++i) add_geom_blob(cl, 0, 0, 0);
    return &cl;
}

// Stamp every blob of `c` with `role` in "nu_band_veto_role"/"perblob".
static void stamp(Cluster* c, int role)
{
    c->put_pcarray(std::vector<int>(c->nchildren(), role), "nu_band_veto_role", "perblob");
}

static std::vector<int> perblob_vals(Cluster* c, const std::string& key)
{
    std::vector<int> ret;
    auto& lpcs = c->value().local_pcs();
    auto it = lpcs.find("perblob");
    if (it == lpcs.end()) return ret;
    auto arr = it->second.get(key);
    if (!arr) return ret;
    auto sp = arr->elements<int>();
    ret.assign(sp.begin(), sp.end());
    return ret;
}

static bool has_key(Cluster* c, const std::string& key)
{
    auto& lpcs = c->value().local_pcs();
    auto it = lpcs.find("perblob");
    if (it == lpcs.end()) return false;
    return it->second.get(key) != nullptr;
}

// Build a cluster_connectivity_graph_t over ALL of g's current children (in
// children() order, matching the vertex-index invariant merge_clusters()
// relies on -- see the long comment at the top of ClusteringFuncs.cxx's
// merge_clusters), with an edge for each pair in `edges`.
static cluster_connectivity_graph_t build_graph(
    Grouping& g, const std::vector<std::pair<Cluster*, Cluster*>>& edges)
{
    cluster_connectivity_graph_t graph;
    auto live = g.children();
    std::unordered_map<const Cluster*, int> idx;
    idx.reserve(live.size());
    for (size_t i = 0; i < live.size(); ++i) {
        idx[live[i]] = (int) i;
        boost::add_vertex((int) i, graph);
    }
    for (const auto& [a, b] : edges) {
        boost::add_edge(idx.at(a), idx.at(b), graph);
    }
    return graph;
}

TEST_CASE("nu_band_veto: predicate truth table")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();

    Cluster* band = make_cluster(*g, 2);
    Cluster* nonband = make_cluster(*g, 2);
    Cluster* band2 = make_cluster(*g, 2);
    Cluster* unmarked = make_cluster(*g, 2);
    Cluster* both = make_cluster(*g, 2);   // carries BOTH roles (post-merge mixed)

    stamp(band, band_veto_band);
    stamp(nonband, band_veto_nonband);
    stamp(band2, band_veto_band);
    // unmarked: no array at all.
    both->put_pcarray(std::vector<int>{band_veto_band, band_veto_nonband},
                      "nu_band_veto_role", "perblob");

    // band vs non-band forbids, both orientations.
    CHECK(band_veto_forbids(band, nonband));
    CHECK(band_veto_forbids(nonband, band));
    // Same role never forbids.
    CHECK_FALSE(band_veto_forbids(band, band2));
    // Marked vs unmarked never forbids (fail-open on the absent side).
    CHECK_FALSE(band_veto_forbids(band, unmarked));
    CHECK_FALSE(band_veto_forbids(nonband, unmarked));
    // Neither marked never forbids.
    CHECK_FALSE(band_veto_forbids(unmarked, unmarked));
    // A cluster carrying BOTH roles against something unmarked: the unmarked
    // side has neither role, so still no veto (needs a role on BOTH sides).
    CHECK_FALSE(band_veto_forbids(both, unmarked));
    // ... but against an actual band or non-band, it does forbid (it has the
    // opposite role too).
    CHECK(band_veto_forbids(both, band));      // both has nonband row
    CHECK(band_veto_forbids(both, nonband));   // both has band row
}

TEST_CASE("nu_band_veto: carry_singles zero-fills a member with no array")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();

    Cluster* marked = make_cluster(*g, 2);
    stamp(marked, band_veto_band);
    Cluster* plain = make_cluster(*g, 3);   // no array at all

    // marked--plain does NOT trigger the veto: plain lacks the array, so
    // band_veto_forbids is fail-open false regardless of geometry -- these
    // are bare blobs with no "scalar" PC, safe as long as no edge is dropped.
    CHECK_FALSE(band_veto_forbids(marked, plain));

    auto graph = build_graph(*g, {{marked, plain}});
    auto fresh = merge_clusters(graph, *g);
    REQUIRE(fresh.size() == 1);
    Cluster* merged = fresh.at(0);
    REQUIRE(merged->nchildren() == 5);

    auto roles = perblob_vals(merged, "nu_band_veto_role");
    REQUIRE(roles.size() == 5);
    // 2 rows carried verbatim (band_veto_band), 3 rows zero-filled (plain
    // never wrote the array) -- the array's own documented "unmarked" value,
    // not a sentinel, per the carry_singles comment in ClusteringFuncs.cxx.
    const int n_band = (int) std::count(roles.begin(), roles.end(), band_veto_band);
    const int n_zero = (int) std::count(roles.begin(), roles.end(), band_veto_none);
    CHECK(n_band == 2);
    CHECK(n_zero == 3);
}

TEST_CASE("nu_band_veto: no writer anywhere => no key, byte-identical")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();

    Cluster* a = make_cluster(*g, 2);
    Cluster* b = make_cluster(*g, 2);
    // Neither carries "nu_band_veto_role" -- the knob-off state.

    auto graph = build_graph(*g, {{a, b}});
    auto fresh = merge_clusters(graph, *g);
    REQUIRE(fresh.size() == 1);
    Cluster* merged = fresh.at(0);
    REQUIRE(merged->nchildren() == 4);
    // The carry_singles registry must not CREATE the key when nobody wrote
    // it: an empty-but-present array would still change the "perblob" key
    // set (and so the compiled Bee/tensor output) even with the value all
    // zero.
    CHECK_FALSE(has_key(merged, "nu_band_veto_role"));
}

TEST_CASE("nu_band_veto: edge filter drops the EDGE, not the component")
{
    Points::node_t root;
    Grouping* g = root.value.facade<Grouping>();

    SUBCASE("a direct band<->nonband edge is dropped: no merge happens") {
        Cluster* band = make_geom_cluster(*g, 2);
        Cluster* nonband = make_geom_cluster(*g, 2);
        stamp(band, band_veto_band);
        stamp(nonband, band_veto_nonband);

        auto graph = build_graph(*g, {{band, nonband}});
        auto fresh = merge_clusters(graph, *g);
        // The one edge is dropped before connected_components, so both
        // vertices end up as singleton components -- merge_clusters skips
        // components with fewer than 2 members (ClusteringFuncs.cxx: "if
        // (descs.size() < 2) continue;"), so nothing is merged.
        CHECK(fresh.empty());
        CHECK(band->nchildren() == 2);
        CHECK(nonband->nchildren() == 2);
    }

    SUBCASE("an unmarked bridge cluster still transitively joins them") {
        // doc pr/66 sec2.6's finding, reproduced structurally: the veto is a
        // PAIRWISE edge veto, not a component-level partition.  With edges
        // band-bridge and bridge-nonband (but NOT a direct band-nonband
        // edge), neither edge is individually forbidden (bridge carries no
        // array), so both survive and all three end up in ONE merged
        // cluster despite band and nonband never being tested against each
        // other directly.
        Cluster* band = make_geom_cluster(*g, 2);
        Cluster* bridge = make_geom_cluster(*g, 1);
        Cluster* nonband = make_geom_cluster(*g, 2);
        stamp(band, band_veto_band);
        stamp(nonband, band_veto_nonband);
        // bridge: no array.

        CHECK_FALSE(band_veto_forbids(band, bridge));
        CHECK_FALSE(band_veto_forbids(bridge, nonband));

        auto graph = build_graph(*g, {{band, bridge}, {bridge, nonband}});
        auto fresh = merge_clusters(graph, *g);
        REQUIRE(fresh.size() == 1);
        CHECK(fresh.at(0)->nchildren() == 5);
    }
}
