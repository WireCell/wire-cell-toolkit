#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/Util.h"
#include "WireCellUtil/Testing.h"

#include <boost/multi_array.hpp>
#include <iostream>

using namespace WireCell;

int main()
{
    const std::vector<size_t> shape = {1, 2, 9, 9};
    size_t nele = 1;
    for (auto s : shape) {
        nele *= s;
    }
    auto st = std::make_shared<Aux::SimpleTensor >(shape, (float*)nullptr);
    Assert(st->size() == sizeof(float) * nele);
    auto& d = st->store();
    Assert(d.size() == st->size());

    boost::multi_array_ref<float, 4> mar((float*) d.data(), shape);
    Assert(mar.num_elements() == nele);

    for (size_t i = 0; i < shape[0]; ++i) {
        for (size_t j = 0; j < shape[1]; ++j) {
            for (size_t k = 0; k < shape[2]; ++k) {
                for (size_t l = 0; l < shape[3]; ++l) {
                    mar[i][j][k][l] = (i + 1) * 1000 + (j + 1) * 100 + (k + 1) * 10 + l + 1;
                }
            }
        }
    }

    std::cout << "print 4d tensor\n";
    std::cout << Aux::dump<float>(st);

    std::cout << "\nspecify limit 4/dim (default 10): \n";
    std::cout << Aux::dump<float>(st, 4) << "\n";

    std::cout << "\npartial print {0, 1, : , :} (not fully polished yet)\n";
    std::cout << Aux::dump<float>(st, {0, 1}) << "\n";

    return 0;
}
