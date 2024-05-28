

#include "WireCellUtil/EigenFFT.h"

#include <iostream>
using namespace std;

int main()
{
    Eigen::FFT<float> fft;
    std::vector<float> timevec = {1, 2, 3, 2, 1};
    std::vector<std::complex<float> > freqvec;

    fft.fwd(freqvec, timevec);
    // manipulate freqvec
    fft.inv(timevec, freqvec);

    for (auto x : timevec) {
        cerr << x << " ";
    }
    cerr << endl;
    for (auto x : freqvec) {
        cerr << x << " ";
    }
    cerr << endl;

    // stored "plans" get destroyed with fft destructor

    return 0;
}
