// Test the edep Ionization stripper: ISimTruth -> ITrackSegmentSet.  Hermetic
// (no Geant4): feed a hand-built ISimTruth and check the segments come out
// as-is.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/ISimTruthSegments.h"

#include "WireCellAux/SimpleTrackSegment.h"
#include "WireCellAux/SimpleTrackSegmentSet.h"
#include "WireCellAux/SimpleSimTruth.h"

using namespace WireCell;

TEST_CASE("edep ionization stripper")
{
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellEdep");

    auto node = Factory::find<ISimTruthSegments>("Ionization");
    REQUIRE(node);

    // A sim truth carrying two segments (and nothing else).
    ITrackSegment::vector segs;
    segs.push_back(std::make_shared<Aux::SimpleTrackSegment>(Point(0, 0, 0), Point(0, 0, 10), 0, 1, 5.0));
    segs.push_back(std::make_shared<Aux::SimpleTrackSegment>(Point(0, 0, 10), Point(0, 0, 20), 1, 2, 6.0));
    auto tss = std::make_shared<Aux::SimpleTrackSegmentSet>(42, segs);
    auto truth = std::make_shared<Aux::SimpleSimTruth>(42, nullptr, nullptr, tss, nullptr);

    ISimTruthSegments::output_pointer out;
    const bool ok = (*node)(truth, out);
    CHECK(ok);
    REQUIRE(out);
    CHECK(out->ident() == 42);
    REQUIRE(out->segments());
    CHECK(out->segments()->size() == 2);

    // The stripper re-emits the same object, no copy.
    CHECK(out.get() == tss.get());
}
