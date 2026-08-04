#include "WireCellUtil/doctest.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Units.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellIface/ITrackSegmentSampler.h"
#include "WireCellAux/SimpleTrackSegment.h"
#include "WireCellAux/SimpleTrackSegmentSet.h"

#include <cmath>
#include <memory>

using namespace WireCell;
using WireCell::Aux::SimpleTrackSegment;
using WireCell::Aux::SimpleTrackSegmentSet;

namespace {

// Instantiate + configure a named component with its defaults, as a WCT job's
// config sequence would before any component that references it.
void preconfigure(const std::string& tn)
{
    auto icfg = Factory::lookup_tn<IConfigurable>(tn);
    icfg->configure(icfg->default_configuration());
}

ITrackSegmentSampler::pointer make_sampler(const std::string& name, Configuration extra)
{
    auto icfg = Factory::lookup_tn<IConfigurable>("TrackSegmentSampler:" + name);
    auto cfg = icfg->default_configuration();
    for (const auto& key : extra.getMemberNames()) {
        cfg[key] = extra[key];
    }
    icfg->configure(cfg);
    return Factory::lookup_tn<ITrackSegmentSampler>("TrackSegmentSampler:" + name);
}

double total_charge(const IDepoSet::pointer& ds)
{
    double tot = 0;
    for (const auto& depo : *ds->depos()) {
        tot += depo->charge();
    }
    return tot;
}

}  // namespace

TEST_CASE("gen tracksegmentsampler recombination MIP") {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellGen");

    preconfigure("BoxRecombination");
    Configuration extra;
    extra["ionization"] = "recombination";
    extra["recombination"] = "BoxRecombination";
    auto sampler = make_sampler("mip", extra);

    // A MIP-like 1 cm segment: dE/dx ~ 2.1 MeV/cm at the default 500 V/cm.
    const double energy = 2.1 * units::MeV;
    const double length = 1 * units::cm;
    const Point start(10 * units::cm, 0, 0);
    const Point stop(10 * units::cm, 0, length);
    ITrackSegment::vector segs{std::make_shared<SimpleTrackSegment>(
        start, stop, 0.0, 0.05 * units::us, energy, 0.0, -1.0, length, 42, 13)};

    ITrackSegmentSampler::output_pointer out;
    REQUIRE((*sampler)(std::make_shared<SimpleTrackSegmentSet>(3, segs), out));
    REQUIRE(out != nullptr);
    CHECK(out->ident() == 3);

    auto depos = out->depos();
    REQUIRE(depos != nullptr);
    // 10 mm at the default 1 mm step.
    REQUIRE(depos->size() == 10);

    // Total charge matches the model called directly, and is NEGATIVE.
    auto model = Factory::lookup_tn<IRecombinationModel>("BoxRecombination");
    const double expected = (*model)(energy, length);
    const double got = total_charge(out);
    CHECK(got < 0);
    CHECK(got == doctest::Approx(expected).epsilon(1e-9));

    // Physics anchor: recombination factor ~0.7 at 500 V/cm, so the electron
    // count for 2.1 MeV/cm is ~6e4/cm (23.6 eV ionization W-value).
    const double nele = std::abs(got);
    CHECK(nele > 5.5e4);
    CHECK(nele < 7.0e4);
    const double rfactor = nele * 23.6 * units::eV / energy;
    CHECK(rfactor > 0.6);
    CHECK(rfactor < 0.8);

    // The model's inverse recovers the deposited energy.
    CHECK(model->dE(got, length) == doctest::Approx(energy).epsilon(1e-6));

    // Geometry/time interpolation: samples at the sub-interval centers.
    const auto& first = depos->at(0);
    const auto& last = depos->at(9);
    CHECK(first->pos().z() == doctest::Approx(0.5 * units::mm));
    CHECK(last->pos().z() == doctest::Approx(9.5 * units::mm));
    CHECK(first->pos().x() == doctest::Approx(10 * units::cm));
    CHECK(first->time() == doctest::Approx(0.05 * 0.05 * units::us));
    CHECK(last->time() == doctest::Approx(0.95 * 0.05 * units::us));
    // Truth keys and per-sample energy conservation.
    CHECK(first->id() == 42);
    CHECK(first->pdg() == 13);
    CHECK(first->energy() == doctest::Approx(energy / 10));
    CHECK(first->extent_long() == 0.0);
    CHECK(first->extent_tran() == 0.0);
}

TEST_CASE("gen tracksegmentsampler quanta mode") {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellGen");

    Configuration extra;
    extra["ionization"] = "quanta";
    auto sampler = make_sampler("quanta", extra);

    const Point start(0, 0, 0);
    const Point stop(0, 0, 2 * units::mm);
    const double nele = 1.2e4;

    ITrackSegment::vector segs{std::make_shared<SimpleTrackSegment>(
        start, stop, 0.0, 1 * units::ns, 0.42 * units::MeV, 0.14 * units::MeV, nele,
        2 * units::mm, 7, 11)};
    ITrackSegmentSampler::output_pointer out;
    REQUIRE((*sampler)(std::make_shared<SimpleTrackSegmentSet>(1, segs), out));
    REQUIRE(out->depos()->size() == 2);
    CHECK(total_charge(out) == doctest::Approx(-nele));

    // A segment lacking n_electrons is an error in quanta mode.
    ITrackSegment::vector bad{std::make_shared<SimpleTrackSegment>(
        start, stop, 0.0, 1 * units::ns, 0.42 * units::MeV)};
    ITrackSegmentSampler::output_pointer bad_out;
    CHECK_THROWS((*sampler)(std::make_shared<SimpleTrackSegmentSet>(2, bad), bad_out));
}

TEST_CASE("gen tracksegmentsampler edge cases") {
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellGen");

    preconfigure("BoxRecombination");
    auto sampler = make_sampler("edges", Configuration());

    // Zero-length segment: one depo at the (degenerate) midpoint.
    const Point at(1 * units::mm, 2 * units::mm, 3 * units::mm);
    ITrackSegment::vector segs{std::make_shared<SimpleTrackSegment>(
        at, at, 5 * units::ns, 5 * units::ns, 0.1 * units::MeV)};
    ITrackSegmentSampler::output_pointer out;
    REQUIRE((*sampler)(std::make_shared<SimpleTrackSegmentSet>(1, segs), out));
    REQUIRE(out->depos()->size() == 1);
    CHECK(out->depos()->at(0)->pos() == at);
    CHECK(out->depos()->at(0)->time() == 5 * units::ns);
    CHECK(out->depos()->at(0)->charge() < 0);

    // EOS passes through.
    ITrackSegmentSampler::output_pointer eos_out;
    REQUIRE((*sampler)(nullptr, eos_out));
    CHECK(eos_out == nullptr);

    // Empty set -> empty depo set, same ident.
    ITrackSegmentSampler::output_pointer empty_out;
    REQUIRE((*sampler)(std::make_shared<SimpleTrackSegmentSet>(9, ITrackSegment::vector{}),
                       empty_out));
    REQUIRE(empty_out != nullptr);
    CHECK(empty_out->ident() == 9);
    CHECK(empty_out->depos()->size() == 0);
}
