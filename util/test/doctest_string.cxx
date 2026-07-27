#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"

using namespace WireCell::String;

TEST_CASE("string format")
{
    SUBCASE("ordered positional") {
        CHECK("abc" == format("a%sc", "b"));
    }
    SUBCASE("numbered positional ") {
        CHECK("cba" == format("%2%b%1%", "a", "c"));
    }
}
