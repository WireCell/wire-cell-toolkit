#include "WireCellSpng/Testing.h"
#include "WireCellSpng/Convo.h"

using namespace WireCell::SPNG;

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
