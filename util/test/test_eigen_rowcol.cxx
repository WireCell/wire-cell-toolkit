#include "WireCellUtil/Eigen.h"
#include <iostream>

using DEFM = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic>; // should be ColMajor
using COLM = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using ROWM = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

COLM get_mapped_cw()
{
    std::vector<float> col_major{11,21,12,22,13,23};
    Eigen::Map<COLM> ret(col_major.data(), 2,3);
    return ret;
}
ROWM get_mapped_rw()
{
    std::vector<float> row_major{11, 12, 13, 21, 22, 23};
    Eigen::Map<ROWM> ret(row_major.data(), 2,3);
    return ret;
}
    

COLM get_colwise()
{
    COLM ret(2,3);
    for (int major=0; major<2; ++major) {
        for (int minor=0; minor<3; ++minor) {
            ret(major,minor) = (major+1)*10 + minor+1;
        }
    }
    return ret;
}
ROWM get_rowwise()
{
    ROWM ret(2,3);
    for (int major=0; major<2; ++major) {
        for (int minor=0; minor<3; ++minor) {
            ret(major,minor) = (major+1)*10 + minor+1;
        }
    }
    return ret;
}

void dump_def(DEFM arr)
{
    std::cout << "DEFM" << "("<<arr.rows()<<","<<arr.cols()<<")\n" << arr << "\n";
    for (int major=0; major<2; ++major) {
        for (int minor=0; minor<3; ++minor) {
            std::cout << arr(major,minor)
                      << "=" << arr.data()[minor*2+major] << " ";
        }
    }
    std::cout << "\n";
}
void dump_cw(COLM arr)
{
    std::cout << "COLM" << "("<<arr.rows()<<","<<arr.cols()<<")\n" << arr << "\n";
    for (int major=0; major<2; ++major) {
        for (int minor=0; minor<3; ++minor) {
            std::cout << arr(major,minor)
                      << "=" << arr.data()[minor*2+major] << " ";
        }
    }
    std::cout << "\n";
}
void dump_rw(ROWM arr)
{
    std::cout << "ROWM" << "("<<arr.rows()<<","<<arr.cols()<<")\n" << arr << "\n";
    for (int major=0; major<2; ++major) {
        for (int minor=0; minor<3; ++minor) {
            std::cout << arr(major,minor)
                      << "=" << arr.data()[major*3+minor] << " ";
        }
    }
    std::cout << "\n";
}

void dump_data(const COLM& arr)
{
    std::cout << "COLM:";
    const int size = arr.rows() * arr.cols();
    for (int ind=0; ind<size; ++ind) {
        std::cout << arr(ind) << " ";
    }
    std::cout << std::endl;
}
void dump_data(float* data, int size=6)
{
    std::cout << "data:";
    for (int ind=0; ind<size; ++ind) {
        std::cout << data[ind] << " ";
    }
    std::cout << std::endl;
}

int main()
{
    std::cout << "mapped-cw:\n";
    dump_def(get_mapped_cw());
    dump_cw(get_mapped_cw());
    dump_rw(get_mapped_cw());
    std::cout << "mapped-rw:\n";
    dump_def(get_mapped_rw());
    dump_cw(get_mapped_rw());
    dump_rw(get_mapped_rw());

    std::cout << "column-wise:\n";
    dump_def(get_colwise());
    dump_cw(get_colwise());
    dump_rw(get_colwise());
    std::cout << "row-wise:\n";
    dump_def(get_rowwise());
    dump_cw(get_rowwise());
    dump_rw(get_rowwise());

    auto cw = get_colwise();
    std::cout << "cw is row-major:" << cw.IsRowMajor << "\n";
    dump_data(cw);
    dump_data(cw.data());

    auto cat = cw.transpose();
    std::cout << "cat is row-major:" << cat.IsRowMajor << "\n";
    dump_data(cat);
    dump_data(cat.data());

    COLM cct = cw.transpose();
    std::cout << "cct is row-major:" << cct.IsRowMajor << "\n";
    dump_data(cct);
    dump_data(cct.data());


}
/*
  conclusions: 

  - default is indeed column-wise storage order

  - a copy between ROWM and COLM transforms the underlying data.

  - indices given via the callable operator are invariant.

  - transpose:

    - a transpose() returns same array with IsRowMajor toggled.

    - a IsRowMajor COLM copied to a COLM applies the transpose to the
      data.

    - access via indicies is always as expected but .data() order
      interpretation must check C++ type against IsRowMajor!
 */
