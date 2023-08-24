#include "WireCellUtil/ExecMon.h"
#include "WireCellUtil/Waveform.h"

#include <random>
#include <iostream>
using namespace std;
using namespace WireCell::Waveform;

int main()
{
    WireCell::ExecMon em("test median starts");
    std::random_device rd;
    std::default_random_engine re(rd());
    std::uniform_real_distribution<> dist(0, 1000);

    em("generate randoms");

    const int nmaxbins = 10000;
    std::vector<float> data(nmaxbins);
    for (int ind = 0; ind < nmaxbins; ++ind) {
        data[ind] = dist(re);
    }

    const int ntimes = 10000;
    em("start scaling tests");
    int sizes[] = {100, 1000, 10000, 0};
    float dummy = 0;  // don't optimize away
    for (int ind = 0; sizes[ind]; ++ind) {
        auto size = sizes[ind];
        std::vector<float> chunk(data.begin(), data.begin() + size);

        em("start median");
        float med1 = median(chunk);
        cerr << ind << ": median of " << size << " = " << med1 << endl;
        for (int tmp = 0; tmp < ntimes; ++tmp) {
            dummy += median(chunk);
        }
        cerr << em("end medium (many times)") << endl;

        em("start median_binned");
        float med2 = median_binned(chunk);
        cerr << ind << ": median_binned of " << size << " = " << med2 << endl;
        for (int tmp = 0; tmp < ntimes; ++tmp) {
            dummy += median_binned(chunk);
        }
        cerr << em("end medium_binned (many times)") << endl;
        cerr << "Median difference: " << med1 - med2 << endl;

        em("done with chunked test");
    }
    cerr << "dummy value: " << dummy << endl; // to avoid warning of used var

    cerr << em.summary() << endl;
    return 0;
}
