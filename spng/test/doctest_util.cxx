#include "WireCellSpng/Testing.h"
#include "WireCellSpng/Util.h"

using namespace WireCell::Torch;

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

TEST_CASE("spng util linear shape")
{
    CHECK(4 == linear_shape(2,3));
    std::vector<int64_t> vshape1 = {2,2};
    std::vector<int64_t> vshape2 = {3,3};

    torch::IntArrayRef tshape1(vshape1);
    torch::IntArrayRef tshape2(vshape2);


    // both give vectors back
    auto vshapemin = linear_shape(vshape1, vshape2);
    std::vector<int64_t> tshapemin = linear_shape(tshape1, tshape2);

    for (size_t ind=0; ind<2; ++ind) {
        CHECK(4 == vshapemin[ind]);
        CHECK(4 == tshapemin[ind]);
    }
}
