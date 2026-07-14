#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <set>

using namespace std;
using namespace WireCell;

TEST_SUITE("wireschema") {

TEST_CASE("load twice shares db") {
    /// BIG FAT NOTE: don't directly load() from user code.  Use
    /// IWireSchema component named WireSchemaFile instead.
    const char* filename = "microboone-celltree-wires-v2.1.json.bz2";
    auto store1 = WireSchema::load(filename);
    auto store2 = WireSchema::load(filename);

    CHECK(store1.db().get() == store2.db().get());

    spdlog::debug("Total store has:");
    spdlog::debug("\t{} detectors", store1.detectors().size());
    spdlog::debug("\t{} anodes", store1.anodes().size());
    spdlog::debug("\t{} faces", store1.faces().size());
    spdlog::debug("\t{} planes", store1.planes().size());
    spdlog::debug("\t{} wires", store1.wires().size());

    CHECK(store1.wires().size() > 0);

    auto& a0 = store1.anode(0);
    spdlog::debug("Got anode 0 with {} faces", store1.faces(a0).size());
}

TEST_CASE("bogus anode id throws") {
    const char* filename = "microboone-celltree-wires-v2.1.json.bz2";
    auto store1 = WireSchema::load(filename);

    const int bogusid = 0xdeadbeaf;
    bool caught = false;
    try {
        store1.anode(bogusid);
    }
    catch (WireCell::Exception& e) {
        caught = true;
        spdlog::debug("Correctly caught exception with unlikely anode id {}:\n{}", bogusid, e.what());
    }
    CHECK(caught);
}

TEST_CASE("plane wire stats") {
    const char* filename = "microboone-celltree-wires-v2.1.json.bz2";
    auto store1 = WireSchema::load(filename);

    for (const auto& anode : store1.anodes()) {
        for (const auto& face : store1.faces(anode)) {
            for (const auto& plane : store1.planes(face)) {
                auto bb = store1.bounding_box(plane);
                auto wp = store1.wire_pitch(plane);
                auto ex = wp.first.cross(wp.second);

                auto wires = store1.wires(plane);
                std::set<int> channels;
                for (const auto& w : wires) {
                    channels.insert(w.channel);
                }
                auto mm = std::minmax_element(channels.begin(), channels.end());

                spdlog::debug("anode:{} face:{} plane:{}:", anode.ident, face.ident, plane.ident);
                spdlog::debug("\tdrift: {}", ex);
                spdlog::debug("\twire:  {}", wp.first);
                spdlog::debug("\tpitch: {}", wp.second);
                spdlog::debug("\tbb: {}cm", bb.bounds() / units::cm);
                spdlog::debug("\tchids: {} in [{}...{}] reading {} wires",
                              channels.size(), *mm.first, *mm.second, wires.size());

                CHECK(wires.size() > 0);
                CHECK(channels.size() > 0);
            }
            auto wplane = store1.planes(face).back();
            auto wires = store1.wires(wplane);

            const double zmax = wires.back().head.z();
            const double zmin = wires.front().head.z();
            const double dd = zmax - zmin;
            const int n = wires.size();
            spdlog::debug("\tW plane: {} wires, [{},{}]m, extent={}cm, pitch={}mm",
                          n, zmin / units::m, zmax / units::m,
                          dd / units::cm, dd / (n - 1) / units::mm);

            CHECK(n > 1);
        }
    }
}

}  // TEST_SUITE("wireschema")
