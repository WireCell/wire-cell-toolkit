#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::WireSchema;

TEST_SUITE("wireschema generate") {

TEST_CASE("generate wires and validate") {
    StoreDB storedb;
    auto& plane = get_append(storedb, 0, 0, 0, 0);

    Ray pitch(Point(0,0,0), Point(0,0,1));
    Ray bounds(Point(0,0,0), Point(0,10,10));

    spdlog::debug("pitch=({},{},{}) -> ({},{},{})",
                  pitch.first[0], pitch.first[1], pitch.first[2],
                  pitch.second[0], pitch.second[1], pitch.second[2]);
    spdlog::debug("bounds=({},{},{}) -> ({},{},{})",
                  bounds.first[0], bounds.first[1], bounds.first[2],
                  bounds.second[0], bounds.second[1], bounds.second[2]);

    int nwires = generate(storedb, plane, pitch, bounds);
    spdlog::debug("nwires={}", nwires);

    for (const auto& wire : storedb.wires) {
        spdlog::debug("{},{},{} ({},{},{}) -> ({},{},{})",
                      wire.ident, wire.channel, wire.segment,
                      wire.tail[0], wire.tail[1], wire.tail[2],
                      wire.head[0], wire.head[1], wire.head[2]);
    }

    CHECK(nwires > 0);

    Store store(std::make_shared<StoreDB>(storedb));
    bool valid = true;
    try {
        validate(store);
    }
    catch (ValueError&) {
        valid = false;
    }
    CHECK(valid);
}

}  // TEST_SUITE("wireschema generate")
