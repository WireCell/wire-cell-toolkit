#include "WireCellUtil/MultiArray.h"

#include <iostream>

using namespace std;

template<typename ArrayType>
void resize(ArrayType& arr, const std::vector<size_t>& sv)
{
    arr.resize(sv);
}
// template<typename ArrayType, typename ...Sizes>
// void resize(ArrayType& arr, Sizes... sizes)
// {
//     std::vector<int> shape = {sizes...};
//     arr.resize(shape);
// }

int main()
{
    const int ndim = 3;

    typedef boost::multi_array<double, ndim> array_type;
    // typedef array_type::index index;

    array_type ar3(boost::extents[3][4][5]);
    resize(ar3,{3,4,5});
    //resize(ar3,3,4,5);

    boost::general_storage_order<ndim> gso = ar3.storage_order();

    cout << "i sh st so as\n";
    for (int ind = 0; ind < ndim; ++ind) {
        cout << ind << " " << ar3.shape()[ind] << " " << ar3.strides()[ind] << " " << gso.ordering(ind) << " "
             << gso.ascending(ind) << endl;
    }

    cout << "\n" << ar3.size() << " " << ar3.num_elements() << " " << ar3.num_dimensions() << endl;
    return 0;
}
