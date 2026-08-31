#include "WireCellUtil/Eigen.h"
#include "WireCellUtil/doctest.h"

using Scalar = int;
using COLM = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using ROWM = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

TEST_SUITE("eigen rowcol2") {

TEST_CASE("matrix dimensions") {
    const int data[8] = {0,1,2,3,4,5,6,7};
    const int shape[2] = {2,4};

    COLM c2c = Eigen::Map<const COLM>(data, shape[0], shape[1]);
    COLM r2c = Eigen::Map<const ROWM>(data, shape[0], shape[1]);
    ROWM c2r = Eigen::Map<const COLM>(data, shape[0], shape[1]);
    ROWM r2r = Eigen::Map<const ROWM>(data, shape[0], shape[1]);

    CHECK(c2c.rows() == shape[0]);
    CHECK(c2r.rows() == shape[0]);
    CHECK(r2c.rows() == shape[0]);
    CHECK(r2r.rows() == shape[0]);

    CHECK(c2c.cols() == shape[1]);
    CHECK(c2r.cols() == shape[1]);
    CHECK(r2c.cols() == shape[1]);
    CHECK(r2r.cols() == shape[1]);
}

TEST_CASE("element access respects memory layout") {
    const int data[8] = {0,1,2,3,4,5,6,7};
    const int shape[2] = {2,4};

    COLM c2c = Eigen::Map<const COLM>(data, shape[0], shape[1]);
    COLM r2c = Eigen::Map<const ROWM>(data, shape[0], shape[1]);
    ROWM c2r = Eigen::Map<const COLM>(data, shape[0], shape[1]);
    ROWM r2r = Eigen::Map<const ROWM>(data, shape[0], shape[1]);

    // origin element is 0 for all layouts
    CHECK(c2c(0,0) == 0);
    CHECK(r2c(0,0) == 0);
    CHECK(c2r(0,0) == 0);
    CHECK(r2r(0,0) == 0);

    // col-major: data[1] is (1,0); row-major: data[1] is (0,1)
    CHECK(c2c(1,0) == 1);
    CHECK(c2r(1,0) == 1);
    CHECK(r2c(0,1) == 1);
    CHECK(r2r(0,1) == 1);

    CHECK(c2c(1,1) == 3);
    CHECK(c2r(1,1) == 3);
    CHECK(r2c(0,3) == 3);
    CHECK(r2r(0,3) == 3);
}

}  // TEST_SUITE("eigen rowcol2")
