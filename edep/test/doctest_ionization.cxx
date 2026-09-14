// Tests for the edep/ Ionization node (IPrimaryVertexSet -> ITrackSegmentSet).
//
// The "config" case is hermetic (no Geant4 run) and always runs.  The "run"
// case drives a real edep-sim event and so needs a detector GDML and the
// Geant4 dataset environment; it is enabled only when WIRECELL_EDEP_GDML points
// at a GDML file, and otherwise self-skips.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include "WireCellIface/IPrimaryTracker.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellAux/SimplePrimaryParticle.h"
#include "WireCellAux/SimplePrimaryVertex.h"
#include "WireCellAux/SimplePrimaryVertexSet.h"

#include <cmath>
#include <cstdlib>

using namespace WireCell;

TEST_CASE("edep ionization config")
{
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellEdep");

    auto icfg = Factory::lookup<IConfigurable>("Ionization");
    REQUIRE(icfg);

    auto cfg = icfg->default_configuration();
    CHECK(cfg.isMember("gdml"));
    CHECK(cfg.isMember("physics_list"));
    CHECK(cfg.isMember("macro"));
    CHECK(cfg.isMember("w_quanta"));

    cfg["gdml"] = "some/geometry.gdml";
    icfg->configure(cfg);  // does not start Geant4 (that is lazy on first call)
}

TEST_CASE("edep ionization run")
{
    const char* gdml = std::getenv("WIRECELL_EDEP_GDML");
    if (!gdml) {
        MESSAGE("set WIRECELL_EDEP_GDML to a GDML file to run the live case");
        return;
    }

    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellEdep");

    auto icfg = Factory::lookup<IConfigurable>("Ionization");
    auto cfg = icfg->default_configuration();
    cfg["gdml"] = gdml;
    icfg->configure(cfg);

    auto node = Factory::find<IPrimaryTracker>("Ionization");

    // One forward muon at the origin.
    const double p = 1000.0 * units::MeV;
    const double mass = 105.6583745 * units::MeV;
    const double energy = std::sqrt(p * p + mass * mass);
    auto mu = std::make_shared<Aux::SimplePrimaryParticle>(13, Vector(0, 0, p), energy);
    IPrimaryParticle::vector parts{mu};
    auto vtx = std::make_shared<Aux::SimplePrimaryVertex>(Point(0, 0, 0), 0.0, parts);
    IPrimaryVertex::vector vtxs{vtx};
    IPrimaryVertexSet::pointer pvs = std::make_shared<Aux::SimplePrimaryVertexSet>(1, vtxs);

    IPrimaryTracker::output_pointer out;
    bool ok = (*node)(pvs, out);
    CHECK(ok);
    REQUIRE(out);
    CHECK(out->ident() == 1);
    REQUIRE(out->segments());
    MESSAGE("produced " << out->segments()->size() << " track segments");
    CHECK(out->segments()->size() > 0);
}
