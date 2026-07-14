#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include "WireCellUtil/MultiArray.h"

using namespace std;

TEST_SUITE("boost multi array") {

TEST_CASE("2D array create fill and copy") {
    typedef boost::multi_array<double, 2> array_type;

    size_t l_size = 10;
    size_t t_size = 3;
    array_type ar(boost::extents[l_size][t_size]);
    spdlog::debug("Dimensions: {}", (int)ar.dimensionality);
    spdlog::debug("Dimension 0 is size {}", ar.shape()[0]);
    spdlog::debug("Dimension 1 is size {}", ar.shape()[1]);

    CHECK(l_size == ar.shape()[0]);
    CHECK(t_size == ar.shape()[1]);

    for (size_t l_ind = 0; l_ind < l_size; ++l_ind) {
        for (size_t t_ind = 0; t_ind < t_size; ++t_ind) {
            ar[l_ind][t_ind] = l_ind * t_ind;
        }
    }
    for (size_t t_ind = 0; t_ind < t_size; ++t_ind) {
        for (size_t l_ind = 0; l_ind < l_size; ++l_ind) {
            spdlog::debug("\t[{}][{}] = {}", l_ind, t_ind, ar[l_ind][t_ind]);
        }
    }

    SUBCASE("copy assignment") {
        array_type ar2(boost::extents[l_size][t_size]);
        ar2 = ar;
        CHECK(l_size == ar2.shape()[0]);
        CHECK(t_size == ar2.shape()[1]);
        for (size_t l_ind = 0; l_ind < l_size; ++l_ind) {
            for (size_t t_ind = 0; t_ind < t_size; ++t_ind) {
                CHECK(ar2[l_ind][t_ind] == ar[l_ind][t_ind]);
            }
        }
    }

    SUBCASE("copy initialization") {
        array_type ar3 = ar;
        CHECK(l_size == ar3.shape()[0]);
        CHECK(t_size == ar3.shape()[1]);
        for (size_t l_ind = 0; l_ind < l_size; ++l_ind) {
            for (size_t t_ind = 0; t_ind < t_size; ++t_ind) {
                CHECK(ar3[l_ind][t_ind] == ar[l_ind][t_ind]);
            }
        }
    }
}

TEST_CASE("3D array creation and fill") {
    int nx = 3, ny = 4, nz = 5;
    typedef boost::multi_array<double, 3> array_type3;

    array_type3 B(boost::extents[nx][ny][nz]);
    CHECK(B.shape()[0] == (size_t)nx);
    CHECK(B.shape()[1] == (size_t)ny);
    CHECK(B.shape()[2] == (size_t)nz);

    typedef array_type3::index index3;
    for (index3 i = 0; i < nx; ++i) {
        for (index3 j = 0; j < ny; ++j) {
            for (index3 k = 0; k < nz; ++k) {
                B[i][j][k] = i * j * k;
            }
        }
    }
    CHECK(B[0][0][0] == 0.0);
    CHECK(B[1][1][1] == 1.0);
    CHECK(B[2][3][4] == 24.0);
}

TEST_CASE("multi array from byte pointer ref") {
    std::vector<std::size_t> shape = {23746, 6};
    std::vector<double> dvec(shape[0] * shape[1]);

    const std::byte* bytes = (std::byte*)dvec.data();
    const double* data = reinterpret_cast<const double*>(bytes);

    auto ref = boost::const_multi_array_ref<double, 2>(data, shape);
    auto ma = boost::multi_array<double, 2>(ref);
    CHECK(ma.shape()[0] == shape[0]);
    CHECK(ma.shape()[1] == shape[1]);
}

}  // TEST_SUITE("boost multi array")
