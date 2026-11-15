// pr_ml doc 79 (P3 / P4): the outside segment graph, the clump term and their knobs.
//
// Pins (1) the legacy defaults of the new knobs -- TaggerCheckNeutrino pr_fit_only = "" (the full PR chain), pr_external_graph = "",
// TrackFitting dqdx_clump_term = traj_clump_claim = 0 -- so that "knob absent => byte-identical" stays true; (2) the graph-file
// parser's checks (ExternalGraph.h): version, units, vertex ids, segment ends on their vertices, at least two points; (3) the
// clump API: clumps are stored per cluster and report no charges before a fit.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Exceptions.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellClus/ExternalGraph.h"
#include "WireCellClus/TrackFitting.h"

#include <filesystem>
#include <fstream>
#include <string>

using namespace WireCell;
using namespace WireCell::Clus;

namespace {
    // a private instance name: the factory caches one instance per (type, name), and other tests configure the unnamed one
    const char* kProbe = "extgraph_probe";

    std::string write_file(const std::string& name, const std::string& text)
    {
        const std::string path = (std::filesystem::current_path() / name).string();
        std::ofstream(path) << text;
        return path;
    }

    const std::string kGood = R"({"version": 1, "units": "mm",
        "vertices": [{"id": 3, "p": [0, 0, 0]}, {"id": 7, "p": [10, 0, 0]}],
        "segments": [{"id": 0, "v": [3, 7], "pts": [[0, 0, 0], [5, 0.5, 0], [10, 0, 0]], "type": "track"}],
        "clumps": [{"id": 0, "pts": [[20, 1, 1], [21, 1, 1]], "type": "shower"}]})";
}

TEST_CASE("doc 79: TaggerCheckNeutrino pr_fit_only knobs default to the full chain")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("TaggerCheckNeutrino", kProbe);
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE(cfg.isMember("pr_fit_only"));
    CHECK(cfg["pr_fit_only"].asString() == "");
    REQUIRE(cfg.isMember("pr_external_graph"));
    CHECK(cfg["pr_external_graph"].asString() == "");
    REQUIRE(cfg.isMember("pr_external_max_med_dist"));
    CHECK(cfg["pr_external_max_med_dist"].asDouble() == doctest::Approx(6.0));
}

TEST_CASE("doc 79: TrackFitting clump parameters default off and round-trip")
{
    TrackFitting tf;
    CHECK(tf.get_parameter("dqdx_clump_term") == doctest::Approx(0.0));
    CHECK(tf.get_parameter("traj_clump_claim") == doctest::Approx(0.0));
    tf.set_parameter("dqdx_clump_term", 1.0);
    tf.set_parameter("traj_clump_claim", 1.0);
    CHECK(tf.get_parameter("dqdx_clump_term") == doctest::Approx(1.0));
    CHECK(tf.get_parameter("traj_clump_claim") == doctest::Approx(1.0));
}

TEST_CASE("doc 79: read_external_graph parses a version-1 file")
{
    const auto path = write_file("doctest_external_graph_good.json", kGood);
    auto eg = PR::read_external_graph(path);
    REQUIRE(eg.vertices.size() == 2);
    REQUIRE(eg.segments.size() == 1);
    REQUIRE(eg.clumps.size() == 1);
    CHECK(eg.segments[0].v0 == 3);
    CHECK(eg.segments[0].v1 == 7);
    REQUIRE(eg.segments[0].pts.size() == 3);
    CHECK(eg.segments[0].pts[1].y() == doctest::Approx(0.5 * units::mm));
    CHECK(eg.vertices[1].p.x() == doctest::Approx(10 * units::mm));
    CHECK(eg.segments[0].type == "track");
    CHECK(eg.clumps[0].pts.size() == 2);
    std::filesystem::remove(path);
}

TEST_CASE("doc 79: read_external_graph refuses a malformed file")
{
    auto bad = [](const std::string& name, std::string text, const std::string& from, const std::string& to) {
        auto i = text.find(from);
        REQUIRE(i != std::string::npos);
        text.replace(i, from.size(), to);
        const auto path = write_file(name, text);
        CHECK_THROWS_AS(PR::read_external_graph(path), ValueError);
        std::filesystem::remove(path);
    };
    bad("doctest_external_graph_v2.json", kGood, "\"version\": 1", "\"version\": 2");
    bad("doctest_external_graph_cm.json", kGood, "\"units\": \"mm\"", "\"units\": \"cm\"");
    bad("doctest_external_graph_vid.json", kGood, "\"v\": [3, 7]", "\"v\": [3, 9]");
    bad("doctest_external_graph_end.json", kGood, "[[0, 0, 0], [5, 0.5, 0], [10, 0, 0]]", "[[0, 0, 0], [5, 0.5, 0], [10, 0, 1]]");
    bad("doctest_external_graph_one.json", kGood, "[[0, 0, 0], [5, 0.5, 0], [10, 0, 0]]", "[[0, 0, 0]]");
    bad("doctest_external_graph_dup.json", kGood, "{\"id\": 7, \"p\": [10, 0, 0]}", "{\"id\": 3, \"p\": [10, 0, 0]}");
}

TEST_CASE("doc 79: charge clumps are stored per cluster and have no charges before a fit")
{
    TrackFitting tf;
    const Facade::Cluster* key = reinterpret_cast<const Facade::Cluster*>(0x10);   // a key only; never dereferenced
    CHECK(tf.get_clump_charges(key) == nullptr);
    tf.set_charge_clumps(key, {{WireCell::Point(0, 0, 0), WireCell::Point(1, 0, 0)}});
    CHECK(tf.get_clump_charges(key) == nullptr);
    tf.set_charge_clumps(key, {});
    CHECK(tf.get_clump_charges(key) == nullptr);
}
