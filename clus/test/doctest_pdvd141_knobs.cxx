// doc pdvd/141: the default-off knobs of the gap-fill beam chain follow-ups.  Pure predicates and defaults only (the
// full paths need an event; their record is the doc 141 arms):
//   sec 3  fit_face_majority / TrackFitting fit_face_from_cells   (solve face = the face of most cells)
//   sec 4  Grouping ref_bridge_3d_support                          (3D test of the reference-layer bridges)
//   sec 5  other_seg_speck_drop                                    (size rule for a far speck)
//   sec 6  PatternAlgorithms m_first_seg_same_component            (first-segment seed)
//   sec 7  beam_companion_drop                                     (companion containment)
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/BeamParticleFunctions.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/TrackFittingPresets.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus::PR;

TEST_CASE("pdvd141 fit_face_majority: majority, tie to the lowest face, empty")
{
    std::pair<int,int> out{-9, -9};
    std::map<std::pair<int,int>, int> c;
    CHECK_FALSE(fit_face_majority(c, out));
    c[{6, 1}] = 2; c[{6, 0}] = 7;
    CHECK(fit_face_majority(c, out));
    CHECK(out == std::make_pair(6, 0));
    c[{6, 1}] = 7;                       // tie -> lowest (apa, face)
    CHECK(fit_face_majority(c, out));
    CHECK(out == std::make_pair(6, 0));
    c[{2, 1}] = 7;                       // tie over anodes -> lowest apa
    CHECK(fit_face_majority(c, out));
    CHECK(out == std::make_pair(2, 1));
    c[{6, 1}] = 8;
    CHECK(fit_face_majority(c, out));
    CHECK(out == std::make_pair(6, 1));
}

TEST_CASE("pdvd141 fit_face_from_cells: TrackFitting default off, round trip")
{
    Clus::TrackFitting tf;
    CHECK(tf.get_parameter("fit_face_from_cells") == doctest::Approx(0.0));
    tf.set_parameter("fit_face_from_cells", 1.0);
    CHECK(tf.get_parameter("fit_face_from_cells") == doctest::Approx(1.0));
    CHECK(Clus::TrackFittingPresets::create_with_current_values().get_parameters().fit_face_from_cells == doctest::Approx(0.0));
}

TEST_CASE("pdvd141 speck size rule: off never drops, on drops a small far component only")
{
    CHECK_FALSE(other_seg_speck_drop(0, 0, 1, 500 * units::cm));
    CHECK_FALSE(other_seg_speck_drop(10, 0, 1, 500 * units::cm));
    CHECK_FALSE(other_seg_speck_drop(0, 10 * units::cm, 1, 500 * units::cm));
    const int n = 10; const double x = 10 * units::cm;
    CHECK(other_seg_speck_drop(n, x, 9, 10.1 * units::cm));
    CHECK_FALSE(other_seg_speck_drop(n, x, 10, 50 * units::cm));   // not fewer than n terminals
    CHECK_FALSE(other_seg_speck_drop(n, x, 3, 10 * units::cm));    // not longer than x
    WireCell::Clus::PR::PatternAlgorithms pa;
    CHECK(pa.m_other_seg_min_terminals == 0);
    CHECK(pa.m_other_seg_far_chord == 0.0);
    CHECK_FALSE(pa.m_first_seg_same_component);
}

TEST_CASE("pdvd141 companion containment: off never drops, on drops below the cut")
{
    CHECK_FALSE(beam_companion_drop(0, 0.0));
    CHECK_FALSE(beam_companion_drop(-1, 0.0));
    CHECK(beam_companion_drop(0.5, 0.49));
    CHECK_FALSE(beam_companion_drop(0.5, 0.5));
    CHECK_FALSE(beam_companion_drop(0.5, 1.0));
}

TEST_CASE("pdvd141 ref_bridge_3d_support: default off, set/get, not carried by Grouping::from")
{
    using namespace WireCell::PointCloud::Tree;
    using namespace WireCell::Clus::Facade;
    Points::node_t root;
    auto* g = root.value.facade<Grouping>();
    CHECK(g->ref_bridge_3d_support() == 0.0);
    CHECK(g->ref_bridge_3d_exempt_pts() == 0);
    g->set_ref_bridge_3d_support(2 * units::cm, 20);
    auto& cl = g->make_child();
    CHECK(cl.grouping()->ref_bridge_3d_support() == doctest::Approx(2 * units::cm));
    CHECK(cl.grouping()->ref_bridge_3d_exempt_pts() == 20);
    Points::node_t root2;
    auto* g2 = root2.value.facade<Grouping>();
    g2->from(*g);
    CHECK(g2->ref_bridge_3d_support() == 0.0);
    g->set_ref_bridge_3d_support(0, 0);
    CHECK(g->ref_bridge_3d_support() == 0.0);
}
