// Example for FFT resampling with zero-padding tricks
#include "WireCellUtil/Waveform.h"

#include "WireCellAux/DftTools.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

#include <algorithm>

// for FFT

#include "WireCellUtil/EigenFFT.h"

#include <complex>
#include <iostream>

using namespace WireCell;

int main()
{
    PluginManager& pm = PluginManager::instance();
    pm.add("WireCellAux");
    auto idft = Factory::lookup_tn<IDFT>("FftwDFT");

    std::vector<float> a = {1, 2, 3, 2, 1};
    // can be sampled to 10 ticks: 1 , 1.35279 , 2 , 2.69443 , 3 , 2.69443 , 2 , 1.35279 , 1 , 0.905573

    auto tran = Aux::DftTools::fwd_r2c(idft, a);

    std::cout << " tran = " << std::endl;
    std::cout << tran.size() << std::endl;
    auto mag = WireCell::Waveform::magnitude(tran);
    for (auto x : mag) {
        std::cout << x << "  ";
    }
    std::cout << std::endl;

    // zero-padding
    int inSmps = tran.size();
    tran.resize(inSmps * 2);
    if (inSmps % 2 == 0) {
        // even number of samples, eg, inSmps=6, 012345 -(zero-padding)-> 012345,000000 -> 012000,000345
        std::rotate(tran.begin() + inSmps / 2, tran.begin() + inSmps, tran.end());
    }
    else {
        // odd number, eg, inSmps=5, 01234,00000 -> 01200,00034
        std::rotate(tran.begin() + (inSmps + 1) / 2, tran.begin() + inSmps, tran.end());
    }

    //
    std::cout << " tran = " << std::endl;
    std::cout << tran.size() << std::endl;
    for (auto x : tran) {
        std::cout << x.real() << "  ";
    }
    std::cout << std::endl;

    // inverse FFT
    auto b = Aux::DftTools::inv_c2r(idft, tran);
    float scale = tran.size() / inSmps;
    //
    std::cout << " b = " << std::endl;
    std::cout << b.size() << std::endl;
    for (auto x : b) {
        std::cout << x * scale << "  ";
    }
    std::cout << std::endl;

    return 0;
}
