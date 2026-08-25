// doc 76 round 3 -- TrackFitting::reset_for_new_event() must drop the state
// that gives a WRONG ANSWER, not only the state that CRASHES.
//
// Round 2 introduced reset_for_new_event() to fix a use-after-free: a visitor
// that owns a TrackFitting as a member keeps it for the whole process, while
// m_grouping / the cluster and blob sets / global_rb_map / the charge maps all
// point into the per-event Points tree that MultiAlgBlobClustering builds as a
// local.  That half was easy to see, because getting it wrong segfaults.
//
// The other half is silent.  These members survive the event that wrote them,
// carry no "which event am I from" stamp, and are read by the NEXT event:
//
//   * the pattern-recognition results TaggerCheckNeutrino stashes for the
//     downstream consumers (Bee particle flow, PrDisplayDump, the kine and
//     tagger trees).  They are written at the END of visit(), so any event
//     that finishes without writing one -- no vertex found, a knob off, or a
//     selected candidate >0 whose values landed on a FRESH fitter instead of
//     the member one -- leaves its predecessor's copy in place for a consumer
//     to read as its own.  The shower and vertex handles additionally point
//     into the destroyed event tree.
//
//   * m_cluster_xext_cache, the ident-keyed memo behind
//     skip_revert_iso_xext_cut.  Cluster idents REPEAT across events, so it is
//     the one that actually moved a physics number.  It has no public seam to
//     populate from a unit test -- its evidence is the group-vs-per-event A/B
//     in doc 76 sec 10.9, not this file.  Stated here so the gap is on the
//     record rather than implied by its absence.
//
// Reverting any clear below must make the matching case fail: every setter
// used here is public and every getter reports the member directly.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRShower.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"

#include <map>
#include <set>
#include <vector>

using namespace WireCell;
using namespace WireCell::PointCloud::Tree;
using namespace WireCell::Clus;

TEST_CASE("doc76r3 event reset: pattern-recognition results do not cross the boundary")
{
    // Event N.  Stand in for TaggerCheckNeutrino's end-of-visit() stash.
    Points::node_t root;
    auto& grouping = *root.value.facade<Facade::Grouping>();
    auto& cluster = grouping.make_child();

    TrackFitting tf;

    PR::Graph graph;
    auto vtx = std::make_shared<PR::Vertex>();
    auto shower = std::make_shared<PR::Shower>(graph);

    PR::IndexedShowerSet showers;
    showers.insert(shower);
    PR::IndexedShowerSet pi0_showers;
    pi0_showers.insert(shower);
    PR::ShowerIntMap map_shower_pio_id;
    map_shower_pio_id[shower] = 7;
    std::map<int, std::vector<PR::ShowerPtr>> map_pio_id_showers;
    map_pio_id_showers[7] = {shower};
    std::map<int, std::pair<double, int>> map_pio_id_mass;
    map_pio_id_mass[7] = {135.0, 2};

    tf.set_main_vertex(vtx);
    tf.set_showers(showers);
    tf.set_pi0_data(pi0_showers, map_shower_pio_id, map_pio_id_showers, map_pio_id_mass);
    tf.set_dropped_satellite_shower_ids(std::set<int>{3, 4});

    PR::KineInfo kine{};
    kine.kine_reco_Enu = 1234.5;
    tf.set_kine_info(kine);

    PR::VertexScoreboard board{};
    board.filled = true;
    board.route = "dl-accept";
    tf.set_vertex_scoreboard(board);

    // Written, so event N's own consumers see event N's answer.
    CHECK(tf.get_main_vertex() == vtx);
    CHECK(tf.get_showers().size() == 1);
    CHECK(tf.get_kine_info().kine_reco_Enu == doctest::Approx(1234.5));
    CHECK(tf.get_vertex_scoreboard().filled);

    // Event N+1 begins.
    tf.reset_for_new_event();

    // ...and inherits nothing.  Before the round-3 fix every one of these held
    // event N's value, and an event that never re-wrote them served event N's
    // vertex, showers and neutrino energy as its own.
    CHECK(tf.get_main_vertex() == nullptr);
    CHECK(tf.get_showers().empty());
    CHECK(tf.get_pi0_showers().empty());
    CHECK(tf.get_dropped_satellite_shower_ids().empty());
    CHECK(tf.get_kine_info().kine_reco_Enu == doctest::Approx(PR::KineInfo{}.kine_reco_Enu));
    CHECK_FALSE(tf.get_vertex_scoreboard().filled);

    (void)cluster;
}

TEST_CASE("doc76r3 event reset: inert on a fresh fitter, so one-event jobs are unchanged")
{
    // The legacy per-event process calls this once, at the top of the only
    // visit() it will ever run.  It must be a no-op there -- that is the whole
    // reason the round could ship without a knob.
    TrackFitting tf;

    CHECK(tf.get_main_vertex() == nullptr);
    CHECK(tf.get_showers().empty());
    CHECK(tf.get_graph() == nullptr);

    tf.reset_for_new_event();

    CHECK(tf.get_main_vertex() == nullptr);
    CHECK(tf.get_showers().empty());
    CHECK(tf.get_graph() == nullptr);
    CHECK_FALSE(tf.get_vertex_scoreboard().filled);

    // Idempotent: a second call changes nothing either.
    tf.reset_for_new_event();
    CHECK(tf.get_showers().empty());
}
