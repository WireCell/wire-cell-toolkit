// Tests for the edep ParticleTracking node (IPrimaryVertexSet -> ISimTruth).
//
// The "config" case is hermetic and checks the component on/off defaults (only
// segments on).  The "run" case drives a real edep-sim event and needs a
// detector GDML and the Geant4 dataset environment; it is enabled only when
// WIRECELL_EDEP_GDML points at a GDML file, and otherwise self-skips.

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

TEST_CASE("edep particletracking config")
{
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellEdep");

    auto icfg = Factory::lookup<IConfigurable>("ParticleTracking");
    REQUIRE(icfg);

    auto cfg = icfg->default_configuration();
    CHECK(cfg["segments"].asBool() == true);       // only segments on by default
    CHECK(cfg["primaries"].asBool() == false);
    CHECK(cfg["trajectories"].asBool() == false);
    CHECK(cfg["photons"].asBool() == false);

    icfg->configure(cfg);  // does not start Geant4 (lazy on first call)
}

TEST_CASE("edep particletracking run")
{
    const char* gdml = std::getenv("WIRECELL_EDEP_GDML");
    if (!gdml) {
        MESSAGE("set WIRECELL_EDEP_GDML to a GDML file to run the live case");
        return;
    }

    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellEdep");

    auto icfg = Factory::lookup<IConfigurable>("ParticleTracking");
    auto cfg = icfg->default_configuration();
    cfg["gdml"] = gdml;
    cfg["trajectories"] = true;  // also exercise trajectory-tree building
    icfg->configure(cfg);

    auto node = Factory::find<IPrimaryTracker>("ParticleTracking");

    const double p = 1000.0 * units::MeV;
    const double mass = 105.6583745 * units::MeV;
    const double energy = std::sqrt(p * p + mass * mass);
    auto mu = std::make_shared<Aux::SimplePrimaryParticle>(13, Vector(0, 0, p), energy);
    IPrimaryParticle::vector parts{mu};
    auto vtx = std::make_shared<Aux::SimplePrimaryVertex>(Point(0, 0, 0), 0.0, parts);
    IPrimaryVertex::vector vtxs{vtx};
    IPrimaryVertexSet::pointer pvs = std::make_shared<Aux::SimplePrimaryVertexSet>(1, vtxs);

    IPrimaryTracker::output_pointer truth;
    const bool ok = (*node)(pvs, truth);
    CHECK(ok);
    REQUIRE(truth);
    CHECK(truth->ident() == 1);

    REQUIRE(truth->segments());
    MESSAGE("segments=" << truth->segments()->segments()->size());
    CHECK(truth->segments()->segments()->size() > 0);

    REQUIRE(truth->trajectories());
    MESSAGE("trajectories=" << truth->trajectories()->trajectories()->size());
    CHECK(truth->trajectories()->trajectories()->size() > 0);

    CHECK(!truth->primaries());  // off by default
    CHECK(!truth->photons());    // off by default
}
