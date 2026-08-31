#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include <WireCellUtil/DetectorWires.h>
#include <vector>

using namespace WireCell;
using namespace WireCell::WireSchema;

TEST_SUITE("detectorwires uboone") {

TEST_CASE("wire plane generation and store validation") {
    std::vector<size_t> expected_nwires = {2400, 2400, 3456};
    StoreDB storedb;
    for (size_t ind = 0; ind < 3; ++ind) {
        auto& plane = get_append(storedb, ind, 0, 0, 0);
        size_t nwires = DetectorWires::plane(storedb, plane, DetectorWires::uboone[ind]);
        CHECK(nwires == expected_nwires[ind]);
    }
    DetectorWires::flat_channels(storedb);
    Store store(std::make_shared<StoreDB>(storedb));
    REQUIRE_NOTHROW(validate(store));
    spdlog::debug("valid");
}

}  // TEST_SUITE("detectorwires uboone")
