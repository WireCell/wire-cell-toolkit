#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Units.h"
#include "WireCellAux/SimpleTrackSegment.h"
#include "WireCellAux/SimpleTrackSegmentSet.h"

using namespace WireCell;
using namespace WireCell::Aux;

TEST_CASE("aux simpletracksegment accessors") {
    const Point start(1 * units::mm, 2 * units::mm, 3 * units::mm);
    const Point stop(1 * units::mm, 2 * units::mm, 8 * units::mm);
    const double t0 = 10 * units::ns;
    const double t1 = 11 * units::ns;
    const double energy = 1.05 * units::MeV;
    const double secondary = 0.35 * units::MeV;
    const double nele = 3.0e4;
    const double tlen = 6 * units::mm;

    SimpleTrackSegment ts(start, stop, t0, t1, energy, secondary, nele, tlen, 42, 13);

    CHECK(ts.start() == start);
    CHECK(ts.stop() == stop);
    CHECK(ts.start_time() == t0);
    CHECK(ts.stop_time() == t1);
    CHECK(ts.energy() == energy);
    CHECK(ts.secondary() == secondary);
    CHECK(ts.n_electrons() == nele);
    CHECK(ts.track_length() == tlen);
    CHECK(ts.id() == 42);
    CHECK(ts.pdg() == 13);
}

TEST_CASE("aux simpletracksegment defaults") {
    const Point start(0, 0, 0);
    const Point stop(3 * units::mm, 4 * units::mm, 0);

    SimpleTrackSegment ts(start, stop, 0, 0, 2.1 * units::MeV);

    // default track_length falls back to the start->stop distance
    CHECK(ts.track_length() == doctest::Approx(5 * units::mm));
    // unknown partition/count defaults
    CHECK(ts.secondary() == 0.0);
    CHECK(ts.n_electrons() < 0);
    CHECK(ts.id() == 0);
    CHECK(ts.pdg() == 0);

    // a producer-given track length may exceed the straight-line distance
    SimpleTrackSegment folded(start, stop, 0, 0, 2.1 * units::MeV, 0, -1, 7 * units::mm);
    CHECK(folded.track_length() == 7 * units::mm);
}

TEST_CASE("aux simpletracksegmentset") {
    ITrackSegment::vector segs;
    for (int ind = 0; ind < 3; ++ind) {
        const Point start(0, 0, ind * units::mm);
        const Point stop(0, 0, (ind + 1) * units::mm);
        segs.push_back(std::make_shared<SimpleTrackSegment>(start, stop, ind * units::ns, (ind + 1) * units::ns,
                                                            0.21 * units::MeV, 0.07 * units::MeV, 6.0e3,
                                                            1 * units::mm, ind, 11));
    }

    SimpleTrackSegmentSet tss(5, segs);
    CHECK(tss.ident() == 5);
    auto got = tss.segments();
    REQUIRE(got != nullptr);
    REQUIRE(got->size() == 3);

    // access polymorphically through the interface
    ITrackSegment::pointer one = got->at(1);
    CHECK(one->start().z() == 1 * units::mm);
    CHECK(one->stop().z() == 2 * units::mm);
    CHECK(one->id() == 1);
    CHECK(one->pdg() == 11);
    CHECK(one->n_electrons() == doctest::Approx(6.0e3));
}
