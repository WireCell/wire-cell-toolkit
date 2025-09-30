#include "WireCellSpng/Testing.h"
#include "WireCellSpng/Util.h"

using namespace WireCell::SPNG;

TEST_CASE("spng util vshape")
{
    std::vector<int64_t> store = {2,3,4};
    torch::IntArrayRef tshape(store);
    auto vec = vshape(tshape);
    REQUIRE(vec.size() == 3);
    CHECK(2 == vec[0]);
    CHECK(3 == vec[1]);
    CHECK(4 == vec[2]);
}

