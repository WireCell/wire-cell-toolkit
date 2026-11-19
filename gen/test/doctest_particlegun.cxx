// Hermetic tests for Gen::ParticleGun (no Geant4, no randomness needed for the
// fixed-value cases).

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include "WireCellIface/IPrimaryVertexSetSource.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellAux/ParticleInfo.h"

#include <cmath>

using namespace WireCell;

TEST_CASE("gen particlegun name to pdg")
{
    CHECK(Aux::ParticleInfo::name_to_pdg("muon") == 13);
    CHECK(Aux::ParticleInfo::name_to_pdg("antimuon") == -13);
    CHECK(Aux::ParticleInfo::name_to_pdg("electron") == 11);
    CHECK(Aux::ParticleInfo::name_to_pdg("positron") == -11);
    CHECK(Aux::ParticleInfo::name_to_pdg("proton") == 2212);
    CHECK(Aux::ParticleInfo::name_to_pdg("gamma") == 22);
    CHECK(Aux::ParticleInfo::name_to_pdg("13") == 13);  // numeric string
}

TEST_CASE("gen particlegun single muon")
{
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellGen");

    auto icfg = Factory::lookup<IConfigurable>("ParticleGun");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    // default: muon, 1 GeV kinetic, +z, count 1
    icfg->configure(cfg);

    auto gun = Factory::find<IPrimaryVertexSetSource>("ParticleGun");
    REQUIRE(gun);

    IPrimaryVertexSetSource::output_pointer out;
    REQUIRE((*gun)(out));
    REQUIRE(out);
    CHECK(out->ident() == 1);
    REQUIRE(out->vertices());
    REQUIRE(out->vertices()->size() == 1);
    auto vtx = out->vertices()->at(0);
    REQUIRE(vtx->particles());
    REQUIRE(vtx->particles()->size() == 1);
    auto part = vtx->particles()->at(0);

    CHECK(part->pdg() == 13);

    const double mass = Aux::ParticleInfo::pdg_to_mass(13);
    const double total = 1000.0 * units::MeV + mass;         // KE + mass
    const double pmag = std::sqrt(total * total - mass * mass);
    CHECK(part->energy() == doctest::Approx(total));
    const Vector mom = part->momentum();
    CHECK(mom.x() == doctest::Approx(0.0));
    CHECK(mom.y() == doctest::Approx(0.0));
    CHECK(mom.z() == doctest::Approx(pmag));  // +z direction

    // count==1: next call is EOS, then false.
    IPrimaryVertexSetSource::output_pointer eos;
    REQUIRE((*gun)(eos));
    CHECK(!eos);
    CHECK(!(*gun)(eos));
}
