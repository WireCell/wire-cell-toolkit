#include "WireCellUtil/Units.h"
#include "WireCellUtil/TimeKeeper.h"

#include "WireCellSpng/Testing.h"

#include <iostream>

using namespace WireCell;

TEST_CASE("spng torch basics") {
    at::Tensor a = at::ones({2, 2}, at::kInt);
    at::Tensor b = at::randn({2, 2});
    auto c = a + b.to(at::kInt);
    
}

using namespace torch::indexing;

static bool quiet = false;

template<typename T>
void dump(const T& ten, const std::string& msg="")
{
    if (quiet) return;
    if (msg.size()) {
        std::cerr << msg << "\n";
    }
    if (ten.dtype() == torch::kComplexFloat) {
        std::cerr << "real\n";
        std::cerr << torch::real(ten) << "\n";
        std::cerr << "imag\n";
        std::cerr << torch::imag(ten) << "\n";
    }
    else{
        std::cerr << ten << "\n";
    }    
}

struct Tens {

    std::vector<long int> shape;    
    torch::TensorOptions topt;
    torch::Tensor s;            // a "signal"
    torch::Tensor f0;           // a 1D filter on dim 0
    torch::Tensor f1;           // a 1D filter on dim 1
    torch::Tensor f;            // their outer product
    torch::Tensor F0;           // f0 in Fourier domain
    torch::Tensor F1;           // f1 in Fourier domain
    torch::Tensor F;            // f in Fourier domain
    Tens(const std::vector<long int>& shape_, const torch::Device& device)
        : shape(shape_)
        , topt(torch::TensorOptions().dtype(torch::kFloat32).device(device).requires_grad(false))
        , s(torch::rand(shape, topt))
        , f0(torch::rand(shape[0], topt))
        , f1(torch::rand(shape[1], topt))
        , f(torch::outer(f0, f1))
        , F0(torch::fft::fft(f0).reshape({-1,1}))
        , F1(torch::fft::fft(f1))
        , F(torch::fft::fft2(f))
        { }
};


static
void test_spng_torch_convo_onestep(const Tens& tens)
{
    torch::fft::ifft2(torch::fft::fft2(tens.s)*tens.F);
}
// static
// void test_spng_torch_convo_twostep(const Tens& tens)
// {
//     auto S1 = torch::fft::fft2(tens.s, torch::nullopt, {1});
//     auto SF1 = S1*tens.F1;

//     auto SF01 = torch::fft::fft2(SF1, torch::nullopt, {0});

//     auto M01 = SF01 * tens.F0;
//     auto m01 = torch::fft::ifft2(M01);
// }


static
void test_spng_torch_convo(const Tens& tens, bool smallness_checks=false)
{
    // Cyclically convolve a 2D uniform random array with two 1D uniform random
    // arrays independently spanning the two dimensions in two ways: via
    // one-shot 2D DFT and via 2 1D, per-dimension DFTs.

    // In Python/Numpy:
    /*
      import numpy
      shape = (5, 10)
      s = numpy.random.uniform(0, 1, size=shape)
      f0 = numpy.random.uniform(0, 1, size=(shape[0],))
      f1 = numpy.random.uniform(0, 1, size=(shape[1],))
      f = numpy.outer(f0, f1)
      S = numpy.fft.fft2(s)
      F = numpy.fft.fft2(f)
      M = S*F
      m = numpy.fft.ifft2(M)
      S1 = numpy.fft.fft2(s, axes=(1,))
      F1 = numpy.fft.fft(f1)
      SF1 = S1*F1
      SF01 = numpy.fft.fft2(SF1, axes=(0,))
      F0 = numpy.fft.fft(f0).reshape(shape[0], 1)
      M01 = SF01*F0
      m01 = numpy.fft.ifft2(M01)
      dr = numpy.abs(numpy.real(m) - numpy.real(m01))
      assert numpy.all(dr < 1e-14)
     */

    auto f = torch::outer(tens.f0, tens.f1);

    {                           // peek at array shape
        auto sizes = f.sizes();
        // std::cerr << sizes << "\n";
        size_t size = sizes.size();
        REQUIRE(size == 2);
        REQUIRE(sizes[0] == tens.shape[0]);
        REQUIRE(sizes[1] == tens.shape[1]);
    }

    // One shot 2D convo
    auto S = torch::fft::fft2(tens.s);
    auto F = torch::fft::fft2(f);
    auto M = S*F;
    auto m = torch::fft::ifft2(M); // complex

    // Per dimension convo
    auto S1 = torch::fft::fft2(tens.s, torch::nullopt, {1});
    auto F1 = torch::fft::fft(tens.f1);
    dump(S1, "S1");
    dump(F1, "F1");

    // auto SF1 = torch::zeros(shape);
    auto SF1 = S1*F1;
    dump(SF1, "SF1");

    auto SF01 = torch::fft::fft2(SF1, torch::nullopt, {0});
    dump(SF01, "SF01");

    auto F0 = torch::fft::fft(tens.f0).reshape({-1,1});
    auto M01 = SF01 * F0;
    auto m01 = torch::fft::ifft2(M01);

    auto m_r = torch::real(m);
    auto m_i = torch::imag(m);
    dump(m, "m");
    auto m01_r = torch::real(m01);
    auto m01_i = torch::imag(m01);
    dump(m01, "m01");

    if (smallness_checks) {
        auto im = torch::max(torch::abs(m_i)).item<float>();
        auto im_01 = torch::max(torch::abs(m01_i)).item<float>();
        CHECK(im < 1e-8);
        CHECK(im_01 < 1e-12);
        auto dr = torch::abs(m - m01);
        auto dr_max = torch::max(dr).item<float>();
        CHECK(dr_max < 1e-6);

        dump(dr, "dr");
        auto small = torch::all( dr < 1e-6 );
        dump(small, "small");

        CHECK(small.item<bool>());
    }
}

static void small(const torch::Device& device)
{
    quiet = false;
    std::vector<long int> shape = {2,5};
    Tens tens(shape, device);
    test_spng_torch_convo(tens, true);
}

TEST_CASE("spng torch convo small cpu")
{
    small(torch::kCPU);
}
TEST_CASE("spng torch convo small gpu")
{
    small(torch::kCUDA);
}

static void perf(const std::string& msg,
                 const torch::Device& device,
                 const std::vector<long int>& shape = {1024,8192}, 
                 const size_t tries = 100)
{
    quiet = true;
    TimeKeeper tk("spng torch convo speed test");
    Tens tens(shape, device);
    tk("start");
        
    // for (size_t ind=0; ind<tries; ++ind) {
    //     test_spng_torch_convo(tens);
    // }
    // tk(msg);

    for (size_t ind=0; ind<tries; ++ind) {
        test_spng_torch_convo_onestep(tens);
    }
    tk(msg + " one-step");

    // for (size_t ind=0; ind<tries; ++ind) {
    //     test_spng_torch_convo_twostep(tens);
    // }
    // tk(msg + " two-step");

    for (size_t ind=0; ind<tries; ++ind) {
        tens.s.quantile(0.5, 1);
    }
    tk(msg + " quantile");

    for (size_t ind=0; ind<tries; ++ind) {
        tens.s.median(1);
    }
    tk(msg + " median");
    // Fixme: add Waveform:: versions

    std::cerr << tk.summary() << "\n";
}

const std::vector<long int> power2 = {1024,8192};
const std::vector<long int> duneind = {800,6000};
const std::vector<long int> dunecol = {960,6000};

TEST_CASE("spng torch convo perf gpu")
{
    const size_t tries=1000;
    perf("GPU 1024x8192 x1000", torch::kCUDA, power2, tries);
    // perf("GPU 960x6000 1000", torch::kCUDA, dunecol, tries);
    // perf("GPU 800x6000 1000", torch::kCUDA, duneind, tries);
}
TEST_CASE("spng torch convo perf cpu")
{
    const size_t tries=100;
    perf("CPU 1024x8192 x100", torch::kCPU, power2, tries);
    // perf("CPU 960x6000 x100", torch::kCPU, dunecol, tries);
    // perf("CPU 800x6000 x100", torch::kCPU, duneind, tries);

}
