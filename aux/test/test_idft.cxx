// Test IDFT implementations.
#include "WireCellUtil/Waveform.h"

#include "aux_test_dft_helpers.h"

#include <chrono>
#include <vector>
#include <thread>
#include <numeric>
#include <iostream>

using namespace WireCell;
using namespace WireCell::Aux::Test;

template<typename ValueType>
void dump(ValueType* data, int nrows, int ncols, std::string msg="")
{
    std::cerr << msg << "("<<nrows<<","<<ncols<<")\n";
    for (int irow=0; irow<nrows; ++irow) {
        for (int icol=0; icol<ncols; ++icol) {
            std::cerr << data[irow*ncols + icol] << " ";
        }
        std::cerr << "\n";
    }
}

static
void test_1d_zero(IDFT::pointer dft, int size = 1024)
{
    std::vector<IDFT::complex_t> inter(size,0), freq(size,0);

    dft->fwd1d(inter.data(), freq.data(), inter.size());
    assert_flat_value(freq, czero);
    dft->inv1d(freq.data(), inter.data(), freq.size());
    assert_flat_value(inter, czero); 
}
static
void test_1d_impulse(IDFT::pointer dft, int size=1024)
{
    std::vector<IDFT::complex_t> inter(size,0), freq(size,0), back(size,0);
    inter[0] = 1.0;

    dft->fwd1d(inter.data(), freq.data(), freq.size());
    assert_flat_value(freq);

    dft->inv1d(freq.data(), back.data(), back.size());
    assert_impulse_at_index(back);
}



static
void test_2d_zero(IDFT::pointer dft, int size = 1024)
{
    int stride=size, nstrides=size;
    std::vector<IDFT::complex_t> inter(stride*nstrides,0);
    std::vector<IDFT::complex_t> freq(stride*nstrides,0);

    dft->fwd2d(inter.data(), freq.data(), nstrides, stride);
    assert_flat_value(inter, 0);
    dft->inv2d(freq.data(), inter.data(), nstrides, stride);
    assert_flat_value(freq, 0);
}
static
void test_2d_impulse(IDFT::pointer dft, int nrows=128, int ncols=128)
{
    const int size = nrows*ncols;
    std::vector<IDFT::complex_t> inter(size,0), freq(size,0), back(size,0);
    inter[0] = 1.0;
    dft->fwd2d(inter.data(), freq.data(), nrows, ncols);
    assert_flat_value(freq);

    dft->inv2d(freq.data(), back.data(), nrows, ncols);
    assert_impulse_at_index(back);

}


static void assert_on_axis(const std::vector<IDFT::complex_t>& freq,
                           int axis, int nrows=128, int ncols=128)
{
    for (int irow=0; irow<nrows; ++irow) {
        for (int icol=0; icol<ncols; ++icol) {
            int ind = irow*ncols + icol;
            auto val = std::abs(freq[ind]);
            if (axis) {
                if (irow==0) {
                    assert_small(std::abs(val - 1.0));
                }
                else {
                    assert_small(val);
                }
            }
            else {
                if (icol==0) {
                    assert_small(std::abs(val - 1.0));
                }
                else {
                    assert_small(val);
                }
            }
        }
    }
}

void test_1b_impulse(IDFT::pointer dft, int axis, int nrows=128, int ncols=128)
{
    const int size = nrows*ncols;
    std::cerr << "1b impulse freq axis="<<axis << " nrows="<<nrows<<" ncols="<<ncols<<"\n";

    std::vector<IDFT::complex_t> inter(size,0), freq(size,0), back(size,0);
    inter[0] = 1.0;

    dft->fwd1b(inter.data(), freq.data(), nrows, ncols, axis);
    dump(freq.data(), nrows, ncols, "freq");
    assert_on_axis(freq, axis, nrows, ncols);

    dft->inv1b(freq.data(), back.data(), nrows, ncols, axis);
    dump(back.data(), nrows, ncols, "back");
    assert_impulse_at_index(back, 0);


    std::vector<IDFT::complex_t> inplace(size,0);
    inplace[0] = 1.0;
    dft->fwd1b(inplace.data(), inplace.data(), nrows, ncols, axis);
    assert_on_axis(inplace, axis, nrows, ncols);

    std::vector<IDFT::complex_t> inback(inplace.begin(), inplace.end());
    dft->inv1b(inback.data(), inback.data(), nrows, ncols, axis);
    assert_impulse_at_index(inback, 0);
}

void fwdrev(IDFT::pointer dft, int id, int ntimes, int size)
{
    int stride=size, nstrides=size;
    std::vector<IDFT::complex_t> inter(stride*nstrides,0);
    std::vector<IDFT::complex_t> freq(stride*nstrides,0);

    // std::cerr << "running " << id << std::endl;

    while (ntimes) {
        //std::cerr << ntimes << "\n";
        dft->fwd2d(inter.data(), freq.data(), nstrides, stride);
        dft->inv2d(freq.data(), inter.data(), nstrides, stride);

        --ntimes;        
        auto tot [[maybe_unused]] = Waveform::sum(inter);
        assert(std::real(tot) == 0);
    }
    //std::cerr << "finished " << id << std::endl;
}

static
void test_2d_threads(IDFT::pointer dft, int nthreads, int nloops, int size = 1024)
{
    using namespace std::chrono;

    steady_clock::time_point t1 = steady_clock::now();

    std::vector<std::thread> workers;

    //std::cerr << "Starting workers\n";
    for (int ind=0; ind<nthreads; ++ind) {
        workers.emplace_back(fwdrev, dft, ind, nloops, size);
    }
    //std::cerr << "Waiting for workers\n";
    for (auto& worker : workers) {
        if (worker.joinable()) {
            worker.join();
        }
    }
    steady_clock::time_point t2 = steady_clock::now();
    duration<double> dt1 = duration_cast<duration<double>>(t2 - t1);
    std::cerr << "ndfts: " << nthreads*nloops
              << " " << nthreads << " " << nloops
              << " " << dt1.count() << std::endl;
}


template<typename ValueType>
void test_2d_transpose(IDFT::pointer dft, int nrows, int ncols)
{
    std::vector<ValueType> arr(nrows*ncols);
    std::iota(arr.begin(), arr.end(), 0);

    std::vector<ValueType> arr2(nrows*ncols, 0);
    std::vector<ValueType> arr3(arr.begin(), arr.end());

    dft->transpose(arr.data(), arr2.data(), nrows, ncols);
    dft->transpose(arr3.data(), arr3.data(), nrows, ncols);

    for (int irow=0; irow<nrows; ++irow) {
        for (int icol=0; icol<ncols; ++icol) {
            assert( arr[irow*ncols + icol] == arr2[icol*nrows + irow]);
            assert( arr[irow*ncols + icol] == arr3[icol*nrows + irow]);
        }
    }
    dump(arr.data(), nrows, ncols, "original");
    dump(arr2.data(), ncols, nrows, "copy");
    dump(arr3.data(), ncols, nrows, "inplace");

}


int main(int argc, char* argv[])
{
    DftArgs args;
    int rc = make_dft_args(args, argc, argv);
    if (rc) { return rc; }
    auto idft = make_dft(args.tn, args.pi, args.cfg);

    test_1d_zero(idft);
    test_1d_impulse(idft);
    test_2d_zero(idft);
    test_2d_impulse(idft);

    test_1b_impulse(idft, 0, 2, 8);
    test_1b_impulse(idft, 1, 2, 8);
    test_1b_impulse(idft, 0, 8, 2);
    test_1b_impulse(idft, 1, 8, 2);

    test_2d_transpose<IDFT::scalar_t>(idft, 2, 8);
    test_2d_transpose<IDFT::scalar_t>(idft, 8, 2);
    test_2d_transpose<IDFT::complex_t>(idft, 2, 8);
    test_2d_transpose<IDFT::complex_t>(idft, 8, 2);

    std::vector<int> sizes = {128,256,512,1024};
    for (auto size : sizes) {
        int ndouble=3, ntot=2*16384/size;
        while (ndouble) {
            int nthread = 1<<ndouble;
            int nloop = ntot/nthread;
            --ndouble;
            std::cerr << "size=" << size << " nthread=" << nthread << " nloop=" << nloop << "\n";
            test_2d_threads(idft, nthread, nloop, size);
        }
    }

    return 0;
}
