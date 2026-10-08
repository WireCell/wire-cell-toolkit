// wcp-porting-img icarus/docs/04: the knobs and the fallback added to run the
// SBND pattern-recognition chain on ICARUS, one job per cryostat.
//
//  * TrackFitting electron_lifetime: fitted dQ *= exp(t_drift / tau), off (= 0)
//    by default so no other detector moves.
//  * Grouping::fallback_apa_face (G1): an uncontained vertex fell back to
//    (apa, face) = (0, 0), which does not exist in the ICARUS west cryostat
//    (anodes {2,3}), so wire_angles(0,0) threw.  Unchanged wherever (0,0) exists.
//  * TaggerCheckNeutrino cosmic_y_mid / cosmic_vtx_z_origin (G2/G3): default 0,
//    the legacy y=0 mid-plane and absolute-z BDT feature.
//  * Grouping dead_region (sec 8): a 3-D region (IFiducial) that is dead on
//    every plane; none by default.

#include "WireCellUtil/doctest.h"

#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Units.h"

#include <cmath>
#include <map>
#include <set>

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("icarus/04 electron_lifetime: off by default, factor 1 when off")
{
    TrackFitting tf;
    CHECK(tf.get_parameter("electron_lifetime") == 0.0);
    const double v = 1.5714 * units::mm / units::us;
    CHECK(TrackFitting::electron_lifetime_factor(100 * units::cm, v, 0.0) == 1.0);
    CHECK(TrackFitting::electron_lifetime_factor(100 * units::cm, v, -1.0) == 1.0);
    CHECK(TrackFitting::electron_lifetime_factor(100 * units::cm, 0.0, 3 * units::ms) == 1.0);
}

TEST_CASE("icarus/04 electron_lifetime: exp(t_drift/tau), drift measured from the anode either side")
{
    const double v = 1.5714 * units::mm / units::us;
    const double tau = 3 * units::ms;
    // 148.2 cm = the ICARUS full drift: t = 943 us, factor exp(0.3144) = 1.369.
    const double d = 148.2 * units::cm;
    const double f = TrackFitting::electron_lifetime_factor(d, v, tau);
    CHECK(f == doctest::Approx(std::exp(d / v / tau)));
    CHECK(f == doctest::Approx(1.3694).epsilon(1e-3));
    // x - x_anode is negative on the other drift side: same factor.
    CHECK(TrackFitting::electron_lifetime_factor(-d, v, tau) == doctest::Approx(f));
    // At the anode nothing is corrected.
    CHECK(TrackFitting::electron_lifetime_factor(0.0, v, tau) == 1.0);
}

TEST_CASE("icarus/04 electron_lifetime: set_parameter round trip")
{
    TrackFitting tf;
    tf.set_parameter("electron_lifetime", 3 * units::ms);
    CHECK(tf.get_parameter("electron_lifetime") == doctest::Approx(3 * units::ms));
}

TEST_CASE("icarus/06 electron_lifetime_apa<N>: off by default, round trip, bad names rejected")
{
    TrackFitting tf;
    CHECK(tf.get_parameters().electron_lifetime_apa.empty());
    CHECK(tf.get_parameter("electron_lifetime_apa2") == 0.0);
    tf.set_parameter("electron_lifetime_apa2", 8.6 * units::ms);
    CHECK(tf.get_parameter("electron_lifetime_apa2") == doctest::Approx(8.6 * units::ms));
    CHECK(tf.get_parameter("electron_lifetime_apa3") == 0.0);
    // The scalar is not touched by a per-anode value.
    CHECK(tf.get_parameter("electron_lifetime") == 0.0);
    CHECK_THROWS(tf.set_parameter("electron_lifetime_apa", 1.0));
    CHECK_THROWS(tf.set_parameter("electron_lifetime_apaX", 1.0));
    CHECK_THROWS(tf.get_parameter("electron_lifetime_apa-1"));
}

TEST_CASE("icarus/06 electron_lifetime_for: the anode's value when set and positive, else the scalar")
{
    const double scalar = 4.1 * units::ms;
    const std::map<int, double> none;
    CHECK(TrackFitting::electron_lifetime_for(none, 0, scalar) == scalar);
    const std::map<int, double> by{{0, 4.0 * units::ms}, {1, 4.2 * units::ms}, {2, 0.0}};
    CHECK(TrackFitting::electron_lifetime_for(by, 0, scalar) == 4.0 * units::ms);
    CHECK(TrackFitting::electron_lifetime_for(by, 1, scalar) == 4.2 * units::ms);
    CHECK(TrackFitting::electron_lifetime_for(by, 2, scalar) == scalar);   // 0 = not set
    CHECK(TrackFitting::electron_lifetime_for(by, 3, scalar) == scalar);   // absent
    CHECK(TrackFitting::electron_lifetime_for(by, -1, scalar) == scalar);  // no anode
}

TEST_CASE("icarus/04 G1: fallback_apa_face keeps (0,0) wherever it exists")
{
    using Facade::Grouping;
    // No detector volumes: legacy.
    CHECK(Grouping::fallback_apa_face({}) == std::make_pair(0, 0));
    // SBND-like (two one-sided anodes) and PDHD-like faces: (0,0) present.
    std::set<WirePlaneId> sbnd{WirePlaneId(kAllLayers, 0, 0), WirePlaneId(kAllLayers, 0, 1)};
    CHECK(Grouping::fallback_apa_face(sbnd) == std::make_pair(0, 0));
    std::set<WirePlaneId> east{WirePlaneId(kAllLayers, 0, 0), WirePlaneId(kAllLayers, 1, 0),
                               WirePlaneId(kAllLayers, 0, 1), WirePlaneId(kAllLayers, 1, 1)};
    CHECK(Grouping::fallback_apa_face(east) == std::make_pair(0, 0));
}

TEST_CASE("icarus/04 G1: ICARUS west (anodes {2,3}) falls back to its own lowest face")
{
    using Facade::Grouping;
    std::set<WirePlaneId> west{WirePlaneId(kAllLayers, 1, 3), WirePlaneId(kAllLayers, 0, 3),
                               WirePlaneId(kAllLayers, 1, 2), WirePlaneId(kAllLayers, 0, 2)};
    CHECK(Grouping::fallback_apa_face(west) == std::make_pair(2, 0));
    // A detector without face 0 on its lowest anode.
    std::set<WirePlaneId> odd{WirePlaneId(kAllLayers, 1, 0), WirePlaneId(kAllLayers, 0, 1)};
    CHECK(Grouping::fallback_apa_face(odd) == std::make_pair(0, 1));
}

TEST_CASE("icarus/04 G2/G3: cosmic_y_mid and cosmic_vtx_z_origin default to 0")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("TaggerCheckNeutrino", "icarus04_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("cosmic_y_mid"), "missing knob: cosmic_y_mid");
    REQUIRE_MESSAGE(cfg.isMember("cosmic_vtx_z_origin"), "missing knob: cosmic_vtx_z_origin");
    CHECK(cfg["cosmic_y_mid"].asDouble() == 0.0);
    CHECK(cfg["cosmic_vtx_z_origin"].asDouble() == 0.0);
}

// doc icarus/04 sec 8: the 3-D dead region.  A slab |z| < w stands in for the
// ICARUS face seam; a point inside it is dead on every plane, answered before
// any wire lookup (so no geometry is needed here).
namespace {
    struct SlabZ : public IFiducial {
        double lo, hi;
        SlabZ(double l, double h) : lo(l), hi(h) {}
        bool contained(const Point& p) const override { return p.z() >= lo && p.z() <= hi; }
    };
}

TEST_CASE("icarus/04 dead_region: none by default")
{
    PointCloud::Tree::Points::node_t root;
    Facade::Grouping* g = root.value.facade<Facade::Grouping>();
    REQUIRE(g);
    CHECK_FALSE(g->in_dead_region(Point(0, 0, 0)));
    g->set_dead_region(nullptr);
    CHECK_FALSE(g->in_dead_region(Point(0, 0, 0)));
}

TEST_CASE("icarus/04 dead_region: inside the slab every plane is dead")
{
    PointCloud::Tree::Points::node_t root;
    Facade::Grouping* g = root.value.facade<Facade::Grouping>();
    REQUIRE(g);
    g->set_dead_region(std::make_shared<SlabZ>(-2.0 * units::cm, 2.5 * units::cm));
    const Point in(-100 * units::cm, 30 * units::cm, 1.0 * units::cm);
    CHECK(g->in_dead_region(in));
    CHECK(g->in_dead_region(Point(0, 0, -2.0 * units::cm)));
    CHECK_FALSE(g->in_dead_region(Point(0, 0, 2.6 * units::cm)));
    CHECK_FALSE(g->in_dead_region(Point(0, 0, -300 * units::cm)));

    // The good-point tests answer "dead on all planes" (apa/face are not consulted).
    int slots[6];
    g->test_good_point(in, 0, 1, slots, 0.6 * units::cm, 1);
    CHECK(slots[0] == 0); CHECK(slots[1] == 0); CHECK(slots[2] == 0);
    CHECK(slots[3] == 1); CHECK(slots[4] == 1); CHECK(slots[5] == 1);
    CHECK(g->is_good_point(in, 0, 1));
    CHECK(g->is_good_point_wc(in, 0, 1));
    for (int pind = 0; pind < 3; ++pind) {
        CHECK(g->get_closest_dead_chs(in, 1, 0, 1, pind));
    }
}

TEST_CASE("icarus/04 dead_region: MultiAlgBlobClustering round-trips the key, empty by default")
{
    PluginManager::instance().add("WireCellClus");
    auto icfg = Factory::lookup<IConfigurable>("MultiAlgBlobClustering", "icarus04_dr_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE_MESSAGE(cfg.isMember("dead_region"), "missing knob: dead_region");
    CHECK(cfg["dead_region"].asString() == "");
}
