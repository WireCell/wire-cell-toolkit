#include "WireCellUtil/Pimpos.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <sstream>

using namespace std;
using namespace WireCell;

TEST_SUITE("pimpos") {

TEST_CASE("basic geometry") {
    const double pitch_dist = 3 * units::mm;
    const int nwires = 2001;
    const double halfwireextent = pitch_dist * 0.5 * (nwires - 1);
    spdlog::debug("Wires at extremes +/- {}mm", halfwireextent / units::mm);
    Pimpos pimpos(nwires, -halfwireextent, halfwireextent);

    for (int ind = 0; ind < 3; ++ind) {
        std::ostringstream oss;
        oss << pimpos.axis(ind);
        spdlog::debug("axis{}: {}", ind, oss.str());
    }

    const Point zero(0, 0, 0);

    {
        auto val = pimpos.relative(zero);
        std::ostringstream oss;
        oss << val;
        spdlog::debug("relative: {}", oss.str());
        CHECK(val == -1.0 * pimpos.origin());
    }

    {
        auto val = pimpos.distance(zero, 0);
        spdlog::debug("distance 0: {}", val);
        CHECK(val == 0.0);
    }
}

TEST_CASE("region binning") {
    const double pitch_dist = 3 * units::mm;
    const int nwires = 2001;
    const double halfwireextent = pitch_dist * 0.5 * (nwires - 1);
    Pimpos pimpos(nwires, -halfwireextent, halfwireextent);

    auto rb = pimpos.region_binning();
    spdlog::debug("[{},{}] [{},{}] binsize:{}mm",
                  rb.range().first, rb.range().second,
                  rb.irange().first, rb.irange().second,
                  rb.binsize() / units::mm);

    CHECK(rb.nbins() == nwires);
    CHECK(rb.min() == -3001.5 * units::mm);
    CHECK(rb.max() == 3001.5 * units::mm);
    CHECK(rb.binsize() == pitch_dist);

    CHECK(rb.bin(0) == nwires / 2);
    CHECK(rb.center(rb.bin(0)) == 0.0);

    CHECK(rb.inside(0.0));
    const double outside = halfwireextent + pitch_dist;
    CHECK(!rb.inside(outside));
}

TEST_CASE("closest wire and reflect") {
    const double pitch_dist = 3 * units::mm;
    const int nwires = 2001;
    const double halfwireextent = pitch_dist * 0.5 * (nwires - 1);
    Pimpos pimpos(nwires, -halfwireextent, halfwireextent);

    const int center_wire = nwires / 2;
    auto center_wi = pimpos.closest(0.0);
    CHECK(center_wi.first == center_wire);
    CHECK(center_wi.second == 0);

    const int center_imp = pimpos.wire_impact(center_wire);
    auto ref1 = pimpos.reflect(center_wire, center_imp + 4);
    CHECK(ref1 == center_imp - 4);
}

}  // TEST_SUITE("pimpos")
