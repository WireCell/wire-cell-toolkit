// PDVD doc 25 round 12 gate: `TrackFitting::fetch_channel_from_anode` (the
// cold-cache path behind `get_channel_for_wire`) checked only the UPPER wire
// bound (`wire >= wires.size()`), never the lower one.  A negative wire index
// -- routine here, since form_point_association's fallback projection walks
// `fb_wire_{u,v,w} + i` outward in a diamond pattern with no floor at 0 --
// reaches `wires[wire]` with wire cast to a huge size_t and segfaults.
//
// Repro (before the fix): PDVD 039349/32,
// TrackFitting::fetch_channel_from_anode <- get_channel_for_wire <-
// form_point_association (TrackFitting.cxx:2777) <- form_map <-
// do_single_tracking <- find_other_segments <- find_proto_vertex <-
// TaggerCheckNeutrino::run_dual_chain_off_pass, apa=1 face=1 plane=2 wire=-5.
//
// The hot-cache branch of get_channel_for_wire (TrackFitting.cxx:735) already
// gets this right -- `wire >= 0 && wire < size` -- so the fix makes the cold
// path (fetch_channel_from_anode) match it instead of inventing new policy.
// Unknobbed: this is undefined behaviour with an unambiguous intended answer
// (the hot path's own contract), not a documented divergence -- CLAUDE.md
// M13 precedent (see the cold_it->second comment a few lines above the
// call site in get_channel_for_wire).

#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade.h"

#include "WireCellIface/IAnodePlane.h"

#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::PointCloud::Tree;

namespace {

static const char* kConfig = "../clus/test/data/uboone-mabc_config.json";

void configure_components(const Json::Value& configs)
{
    auto log = Log::logger("test");
    for (const auto& comp : configs) {
        std::string type = comp["type"].asString();
        std::string name = comp.get("name", "").asString();
        std::string tn = name.empty() ? type : type + ":" + name;
        try {
            auto icfg = Factory::lookup_tn<IConfigurable>(tn);
            auto cfg = icfg->default_configuration();
            if (comp.isMember("data")) {
                for (const auto& key : comp["data"].getMemberNames()) cfg[key] = comp["data"][key];
            }
            icfg->configure(cfg);
        }
        catch (const std::exception& e) {
            log->trace("Skipping {}: {}", tn, e.what());
        }
    }
}

std::vector<IAnodePlane::pointer> load_anodes()
{
    static bool done = false;
    static std::vector<IAnodePlane::pointer> anodes;
    if (done) return anodes;
    done = true;

    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellClus");
    pm.add("WireCellAux");
    pm.add("WireCellGen");
    pm.add("WireCellSigProc");

    auto full = Persist::load(kConfig);
    Json::Value wcfg(Json::arrayValue);
    if (full.isArray() && !full.empty() && full[0]["type"].asString() == "wire-cell") {
        for (Json::ArrayIndex i = 1; i < full.size(); ++i) wcfg.append(full[i]);
    }
    else {
        wcfg = full;
    }
    configure_components(wcfg);

    for (const auto& comp : wcfg) {
        if (comp["type"].asString() != "AnodePlane") continue;
        std::string nm = comp.get("name", "").asString();
        std::string tn = nm.empty() ? "AnodePlane" : "AnodePlane:" + nm;
        anodes.push_back(Factory::find_tn<IAnodePlane>(tn));
    }
    return anodes;
}

}  // namespace

TEST_CASE("get_channel_for_wire: negative wire index does not crash (doc pdvd/25 sec 13.11, 039349/32)")
{
    auto anodes = load_anodes();
    REQUIRE(!anodes.empty());
    auto anode = anodes.front();
    REQUIRE(anode->nfaces() > 0);

    auto faces = anode->faces();
    int face = 0;
    while (face < static_cast<int>(faces.size()) && !faces[face]) ++face;
    REQUIRE(face < static_cast<int>(faces.size()));

    auto planes = faces[face]->planes();
    REQUIRE(!planes.empty());
    int plane = 0;
    auto wires = planes[plane]->wires();
    REQUIRE(!wires.empty());

    const int apa = anode->ident();
    const int expect_ch0 = wires.front()->channel();
    const int nwires = static_cast<int>(wires.size());

    Points::node_t root;
    auto& grouping = *root.value.facade<Facade::Grouping>();
    grouping.set_anodes(anodes);
    auto& cluster = grouping.make_child();

    // Fresh fitter: every lookup below is a cold-cache miss, i.e. it goes
    // straight through fetch_channel_from_anode -- the path that crashed.
    TrackFitting tf;
    tf.preload_clusters({&cluster});

    SUBCASE("in-range wire returns the real channel (regression, not new policy)") {
        CHECK(tf.get_channel_for_wire(apa, face, plane, 0) == expect_ch0);
    }
    SUBCASE("wire == wires.size() is the pre-existing upper-bound guard") {
        CHECK(tf.get_channel_for_wire(apa, face, plane, nwires) == -1);
    }
    SUBCASE("negative wire (production value -5, 039349/32) returns -1, not a crash") {
        CHECK(tf.get_channel_for_wire(apa, face, plane, -5) == -1);
    }
    SUBCASE("wire == -1 returns -1, not a crash") {
        CHECK(tf.get_channel_for_wire(apa, face, plane, -1) == -1);
    }
}
