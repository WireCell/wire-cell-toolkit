// TaggerBeeVisitor (ai-helper issue 33; WireCell/wire-cell-toolkit PR 535 review):
//  * only beam-window main clusters are written, cluster_id = the verdict;
//  * a blob lacking 3d x/y/z or scalar charge is skipped, not dereferenced;
//  * the Bee event index is the one the MABC publishes on the Ensemble
//    (Ensemble::bee_index()), and an Ensemble without one is an error.
#include "WireCellClus/Facade.h"
#include "WireCellClus/IBeeSink.h"
#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/Bee.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/doctest.h"

#include <memory>
#include <string>
#include <vector>

using namespace WireCell;
using namespace WireCell::PointCloud;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus::Facade;

// An IBeeSink that records what it is handed.
class RecordingBeeSink : public WireCell::Clus::IBeeSink {
   public:
    struct Rec { std::string type; size_t index; Configuration js; };
    std::vector<Rec> recs;
    virtual ~RecordingBeeSink() {}
    virtual void acquire() {}
    virtual void release() {}
    virtual size_t write(const Bee::Object& obj, size_t index, int, int, int) {
        auto js = obj.asJson();
        recs.push_back({js["type"].asString(), index, js});
        return index;
    }
};
WIRECELL_FACTORY(RecordingBeeSink, RecordingBeeSink, WireCell::Clus::IBeeSink)

namespace {
    Dataset xyz(std::vector<double> x, std::vector<double> y, std::vector<double> z) {
        return Dataset({{"x", Array(x)}, {"y", Array(y)}, {"z", Array(z)}});
    }
    Dataset cluster_scalar(int is_main, double t0, int stm, int lm_flag) {
        return Dataset({{"flag_main_cluster", Array({is_main})}, {"cluster_t0", Array({t0})},
                        {"flag_STM", Array({stm})}, {"lm_flag", Array({lm_flag})}});
    }
    // Build: live grouping with three clusters.
    //  A: main, t0 in window, STM tagged, lm_flag 2; blobs: good (3 pts, charge 30),
    //     no charge, no z
    //  B: main, t0 outside the window (one good blob)
    //  C: not main, t0 in window (one good blob)
    void build(Points::node_t& root) {
        auto& ens = *root.value.facade<Ensemble>();
        auto& live = ens.make_grouping("live");
        auto* gnode = live.node();
        auto* ca = gnode->insert(Points({{"cluster_scalar", cluster_scalar(1, 1000, 1, 2)}}));
        ca->insert(Points({{"scalar", Dataset({{"charge", Array({30.0})}})},
                           {"3d", xyz({1, 2, 3}, {4, 5, 6}, {7, 8, 9})}}));
        ca->insert(Points({{"scalar", Dataset({{"npoints", Array({2})}})},     // no charge
                           {"3d", xyz({1, 2}, {1, 2}, {1, 2})}}));
        ca->insert(Points({{"scalar", Dataset({{"charge", Array({10.0})}})},
                           {"3d", Dataset({{"x", Array(std::vector<double>{1})},
                                           {"y", Array(std::vector<double>{1})}})}}));  // no z
        auto* cb = gnode->insert(Points({{"cluster_scalar", cluster_scalar(1, 5000, 1, 0)}}));
        cb->insert(Points({{"scalar", Dataset({{"charge", Array({30.0})}})}, {"3d", xyz({1}, {1}, {1})}}));
        auto* cc = gnode->insert(Points({{"cluster_scalar", cluster_scalar(0, 1000, 1, 0)}}));
        cc->insert(Points({{"scalar", Dataset({{"charge", Array({30.0})}})}, {"3d", xyz({1}, {1}, {1})}}));
    }
    std::shared_ptr<RecordingBeeSink> setup(const std::string& name) {
        PluginManager::instance().add("WireCellClus");
        make_RecordingBeeSink_factory();  // the factory lives in this binary, not in a plugin (util doctest precedent)
        auto sink = std::dynamic_pointer_cast<RecordingBeeSink>(
            Factory::lookup<WireCell::Clus::IBeeSink>("RecordingBeeSink", name));
        auto icfg = Factory::lookup<IConfigurable>("TaggerBeeVisitor", name);
        auto cfg = icfg->default_configuration();
        cfg["bee_sink"] = "RecordingBeeSink:" + name;
        cfg["beam_window"][0] = 200.0;
        cfg["beam_window"][1] = 2200.0;
        icfg->configure(cfg);
        return sink;
    }
}

TEST_CASE("tagger bee visitor: candidates, verdicts, blob guards, MABC index")
{
    auto sink = setup("tbv_main");
    auto vis = Factory::find_tn<WireCell::Clus::IEnsembleVisitor>("TaggerBeeVisitor:tbv_main");

    Points::node_t root;
    build(root);
    auto& ens = *root.value.facade<Ensemble>();
    ens.set_bee_index(5);
    REQUIRE_NOTHROW(vis->visit(ens));

    REQUIRE(sink->recs.size() == 4);
    for (const auto& r : sink->recs) {
        CHECK(r.index == 5);
        // only cluster A's good blob: 3 points (B out of window, C not main,
        // A's two broken blobs skipped)
        REQUIRE(r.js["x"].size() == 3);
        CHECK(r.js["q"][0].asDouble() == doctest::Approx(10.0));   // 30 / 3 points
        const int cid = r.js["cluster_id"][0].asInt();
        if (r.type == "tagger_stm" || r.type == "tagger_lm") CHECK(cid == 1);
        else CHECK(cid == 0);                                      // tgm, fc not tagged
    }

    // the next event: whatever index the MABC publishes, not a private counter
    sink->recs.clear();
    ens.set_bee_index(9);
    vis->visit(ens);
    REQUIRE(sink->recs.size() == 4);
    for (const auto& r : sink->recs) CHECK(r.index == 9);
}

TEST_CASE("tagger bee visitor: an Ensemble without a Bee index is an error")
{
    setup("tbv_noindex");
    auto vis = Factory::find_tn<WireCell::Clus::IEnsembleVisitor>("TaggerBeeVisitor:tbv_noindex");
    Points::node_t root;
    build(root);
    auto& ens = *root.value.facade<Ensemble>();
    CHECK_THROWS(vis->visit(ens));
}
