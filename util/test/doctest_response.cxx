#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Response.h"

using namespace WireCell;

TEST_CASE("util field response")
{
    std::string fr_file = "garfield-11impacts-dune.json.bz2";

    auto fr = Response::Schema::load(fr_file.c_str());
    auto fravg = Response::wire_region_average(fr);
}
