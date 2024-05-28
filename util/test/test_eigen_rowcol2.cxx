#include "WireCellUtil/Eigen.h"
#include <cassert>

using Scalar = int;
using COLM = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using ROWM = Eigen::Array<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

int main()
{
    const int data[8] = {0,1,2,3,4,5,6,7};
    const int shape[2] = {2,4};

    COLM c2c = Eigen::Map<const COLM>(data, shape[0], shape[1]);
    COLM r2c = Eigen::Map<const ROWM>(data, shape[0], shape[1]);
    ROWM c2r = Eigen::Map<const COLM>(data, shape[0], shape[1]);
    ROWM r2r = Eigen::Map<const ROWM>(data, shape[0], shape[1]);

    assert(c2c.rows() == shape[0]);
    assert(c2r.rows() == shape[0]);
    assert(r2c.rows() == shape[0]);
    assert(r2r.rows() == shape[0]);

    assert(c2c.cols() == shape[1]);
    assert(c2r.cols() == shape[1]);
    assert(r2c.cols() == shape[1]);
    assert(r2r.cols() == shape[1]);

    assert(c2c(0,0) == 0);
    assert(r2c(0,0) == 0);
    assert(c2r(0,0) == 0);
    assert(r2r(0,0) == 0);

    assert(c2c(1,0) == 1);
    assert(c2r(1,0) == 1);
    assert(r2c(0,1) == 1);
    assert(r2r(0,1) == 1);

    assert(c2c(1,1) == 3);
    assert(c2r(1,1) == 3);
    assert(r2c(0,3) == 3);
    assert(r2r(0,3) == 3);

    return 0;
}
