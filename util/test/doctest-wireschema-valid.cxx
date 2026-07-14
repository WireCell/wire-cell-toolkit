#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

TEST_SUITE("wireschema valid") {

TEST_CASE("wireschema valid") {
    const std::string fname = "microboone-celltree-wires-v2.1.json.bz2";
    WireSchema::Correction cor = WireSchema::Correction::pitch;
    spdlog::debug("loading {} with correction {}.", fname, (int)cor);
    auto store = WireSchema::load(fname.c_str(), cor);
    CHECK_NOTHROW(WireSchema::validate(store));
}

}  // TEST_SUITE("wireschema valid")
