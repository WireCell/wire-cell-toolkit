#include "WireCellSig/Decon2DResponse.h"
#include "WireCellAux/Util.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Logging.h"

#include <iostream>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

using namespace WireCell;

namespace {
    ITensorSet::pointer gen_tensorset(const unsigned int nrow = 6000, const unsigned int ncol = 800)
    {
        Array::array_xxf iarr = Array::array_xxf::Random(nrow, ncol);
        auto oiten = Aux::eigen_array_to_itensor<std::complex<float>>(iarr);
        ITensor::vector *itv = new ITensor::vector;
        itv->push_back(ITensor::pointer(oiten));
        int seqno = 0;
        Configuration md;
        return std::make_shared<Aux::SimpleTensorSet>(seqno, md, ITensor::shared_vector(itv));
    }

    // void print(const ITensorSet::pointer itens)
    // {
    //     Assert(itens);
    //     Assert(itens->tensors()->size() == 1);
    //     auto iten = itens->tensors()->front();
    //     auto oarr = Aux::itensor_to_eigen_array<std::complex<float>>(iten);
    //     std::cout << oarr << "\n";
    // }
}  // namespace

int main()
{
    Log::add_stdout(true, "debug");
    const int nloop = 1;

    auto iitens = gen_tensorset();
    ITensorSet::pointer oitens;

    auto node = std::make_shared<Sig::Decon2DResponse>();
    auto cfg = node->default_configuration();

    auto t1 = Clock::now();
    for (int i = 0; i < nloop; ++i) {
        (*node)(iitens, oitens);
    }
    auto t2 = Clock::now();
    std::cout << "Decon2DResponse: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / (float) nloop << " ms/loop"
              << std::endl;

    return 0;
}
