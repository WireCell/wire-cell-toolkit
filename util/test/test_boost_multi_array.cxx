#include "WireCellUtil/Testing.h"

#include "WireCellUtil/MultiArray.h"

#include <iostream>

using namespace std;

int main()
{
    typedef boost::multi_array<double, 2> array_type;
    // typedef array_type::index index;

    size_t l_size = 10;
    size_t t_size = 3;
    array_type ar(boost::extents[l_size][t_size]);
    cerr << "Dimensions: " << ar.dimensionality << endl;
    cerr << "Dimension 0 is size " << ar.shape()[0] << endl;
    cerr << "Dimension 1 is size " << ar.shape()[1] << endl;

    Assert(l_size == ar.shape()[0]);
    Assert(t_size == ar.shape()[1]);
    for (size_t l_ind = 0; l_ind < l_size; ++l_ind) {
        for (size_t t_ind = 0; t_ind < t_size; ++t_ind) {
            ar[l_ind][t_ind] = l_ind * t_ind;
        }
    }
    for (size_t t_ind = 0; t_ind < t_size; ++t_ind) {
        for (size_t l_ind = 0; l_ind < l_size; ++l_ind) {
            cerr << "\t[" << l_ind << "][" << t_ind << "] = " << ar[l_ind][t_ind];
        }
        cerr << endl;
    }

    array_type ar2(boost::extents[l_size][t_size]);
    ar2 = ar;
    Assert(l_size == ar2.shape()[0]);
    Assert(t_size == ar2.shape()[1]);

    array_type ar3 = ar;
    Assert(l_size == ar3.shape()[0]);
    Assert(t_size == ar3.shape()[1]);

    int nx = 3, ny = 4, nz = 5;

    typedef boost::multi_array<double, 3> array_type3;

    array_type3 B(boost::extents[nx][ny][nz]);
    typedef array_type3::index index3;
    for (index3 i = 0; i < nx; ++i) {
        for (index3 j = 0; j < ny; ++j) {
            for (index3 k = 0; k < nz; ++k) {
                const double foo = i * j * k;

                B[i][j][k] = foo;
            }
        }
    }

    {
        std::vector<std::size_t> shape = {23746,6};
        std::vector<double> dvec(shape[0]*shape[1]);
        
        const std::byte* bytes = (std::byte*)dvec.data();
        const double* data = reinterpret_cast<const double*>(bytes);

        auto ref = boost::const_multi_array_ref<double,2>(data, shape);
        auto ma = boost::multi_array<double, 2> (ref);
        Assert(ma.shape()[0] == shape[0]);
        Assert(ma.shape()[1] == shape[1]);
    }


    return 0;
}
