/*

  Apples-to-apples performance comparison between ArrayFire and LibTorch for
  bottleneck operations identified from OmnibusSigProc.

  It is built in a stand-alone manner with a script:

    ./spng/test/build-afvslt

  Run all tests like:

    jsonnet spng/test/afvslt.jsonnet | OMP_NUM_THREADS=1 ./afvslt

  Run narrowed tests:

    jsonnet spng/test/afvslt.jsonnet -A tests=convo -A devices=gpu   -A techs=af,lt | OMP_NUM_THREADS=1 ./afvslt

  To add a new benchmark:

  0) Add Perf::run_<name>() ABC.
  1) Add generic templated implementation PerfT::run_<name>().
  2) If needed, extend the generic API in namespace "vs" and AF+LT imp.
  3) Extend the factory in perf() to call run_<name>() given a "kind" of "name".
  4) Extend afvslt.jsonnet to include config for the new "name"

 */

// This main program is really a couple libraries all smashed together.  Each
// array technology is kept in an #ifdef/#endif and the testing scaffolding is
// type independent.
#define USE_ARRAYFIRE
#define USE_LIBTORCH
#define USE_EIGEN               // eigen to hold array but algs are fftw3 or std::

#include <vector>
#include <sstream>

using size_type = long long;
using shape_type = std::vector<size_type>;

std::string shape_str(const shape_type& shape)
{
    std::stringstream ss;
    ss << "(";
    std::string comma = "";
    for (const auto dim : shape) {
        ss << comma << dim;
        comma = " ";
    }
    ss << ")";
    return ss.str();
}

namespace vs {
    // A small subset of technology-independent array API which are specialize
    // for a given array technology below.
    //
    // General limitations unless commented otherwise:
    //
    // - The array dtype and shape of args and return value are same.
    //
    // - An array may be 2, 1 or 0 (scalar) dimensions of real or complex
    // single-precision float.

    // Return the "shape" of array
    template<typename T> shape_type shape(const T& arr);

    // Return the real value of the zeroth element of array 
    template<typename T> float to_float(const T& arr);

    // Element-wise arith. return arr1 {+,*,/} arr2.
    template<typename T> T add(const T& arr1, const T& arr2);
    template<typename T> T mul(const T& arr1, const T& arr2);
    template<typename T> T mul(const T& arr1, float val);

    // Some kind of norm of the array.
    template<typename T> float norm(const T& arr);
    template<typename T> T real(const T& arr);
    template<typename T> T imag(const T& arr);
    template<typename T> T fft2(const T& arr);
    template<typename T> T ifft2(const T& arr);
    template<typename T> T median(const T& arr, int dim=0);
    template<typename T> T sort(const T& arr, int dim=0);

    // Call to finish calculation (for lazy libs like AF).
    template<typename T> T& eval(T& arr);
    // Finish all lazy calculations
    template<typename T> void sync();
}

#include <unordered_map>
#include <string>
#include <fstream>
#include <random>
#include <iostream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

using time_point = std::chrono::high_resolution_clock::time_point;

time_point now()
{
    return std::chrono::high_resolution_clock::now();
}
double tdiffms(time_point t)
{
    return std::chrono::duration<double, std::milli>(now() - t).count();
}

// Get a value of a key in an object, setting it to default if missing.
template<typename T> T getdef(json& jobj, const std::string& key, const T& def) {
    auto ret = jobj.value<T>(key, def);
    jobj[key] = ret;
    return ret;
}
    
size_type shape_size(const shape_type& shape)
{
    size_t nele = 1;
    for (auto& dim : shape) {
        nele *= dim;
    }
    return nele;
}

shape_type to_shape(const json& jshape)
{
    shape_type shape;
    for (const auto& ele : jshape) {
        size_t size = ele;
        shape.push_back(size);
    }
    return shape;
}

// Partially abstract base class for use by templated mid class with
// implementations by a concrete base class.
struct Perf {

    std::default_random_engine rng;
    std::uniform_real_distribution<float> runiform;

    Perf(long unsigned int seed = 12345)
        : rng{seed}
        , runiform(0,1)
        {}

    // Template class can call for one-time init prior to a run_*().  
    // Must return cfg, possibly after modification.
    virtual json prerun(json tst) { return tst; };


    // Tech-independent adding of a named, random array.  
    virtual void add(const std::string& name, const shape_type& shape) = 0;
    virtual void add(const std::string& name, const shape_type& shape, const float* data) = 0;

    // The tests
    virtual json run_arith(json tst) = 0;
    virtual json run_convo(json tst) = 0;
    virtual json run_median(json tst) = 0;
    virtual json run_sort(json tst) = 0;

    std::unique_ptr<float[]> randu(size_t n) {
        auto dat = std::make_unique<float[]>(n);
        for (size_t ind=0; ind<n; ++ind) {
            dat[ind] = runiform(rng);
        }
        return dat;
    }

};


struct TestConfig {
    size_t repeat;
    shape_type shape;
    std::string device;
    time_point start{now()};

    TestConfig(const json& tst)
        : repeat(tst["repeat"])
        , shape(tst["shape"])
        , device(tst["device"])
        {
        }

    void go() {
        start = now();
    }

    json results() {
        json ret;
        double ttot = tdiffms(start);
        ret["time"] = ttot;
        ret["dt"] = ttot / repeat;
        return ret;
    }
        

};


// A "bag of arrays" and test.
template<typename T>
struct PerfT : public Perf {

    std::unordered_map<std::string, T> arrays;

    json run_arith(json tst) {
        TestConfig tc(tst);

        add("A", tc.shape);
        add("B", tc.shape);
        add("C", tc.shape);

        T& A = arrays["A"];
        T& B = arrays["B"];
        T& C = arrays["C"];

        tc.go();

        for (size_t count = 0; count < tc.repeat; ++count) {
            A = vs::add(A, vs::add(B,C));
            A = vs::add(A, vs::mul(B,C));
            A = vs::mul(A, 1.0/vs::norm(A));

            // A += B + C;
            // A += B * C;
            // A /= vs::max(A);
        }
        vs::sync<T>();
        return tc.results();
    }

    json run_convo(json tst) {
        TestConfig tc(tst);

        add("signal", tc.shape);
        add("response", tc.shape);
        add("total", tc.shape);
        const T& s = arrays["signal"];
        const T& r = arrays["response"];

        // meaningless accumulate to assure complete eval.
        T tot = arrays["total"];

        auto R = vs::fft2(r);

        tc.go();

        for (size_t count = 0; count < tc.repeat; ++count) {
            // auto t0 = now();

            // auto t = now();
            auto S = vs::fft2(s);
            // double dt_S = tdiffms(t);

            // t = now();
            auto SR = vs::mul(S,R);
            // double dt_SR = tdiffms(t);

            // t = now();
            auto m = vs::ifft2(SR);
            // double dt_m = tdiffms(t);
            
            // t = now();
            // double dt_sync = tdiffms(t);

            // double dt = tdiffms(t0);
            // std::cerr << count << ": S:" << dt_S << " SR:" << dt_SR << " m:" << dt_m
            //           << " sync:" << dt_sync << " tot:" << dt << "\n";

            tot = vs::add(tot, vs::real(m));
        }
        vs::sync<T>();

        return tc.results();
    }

    json run_median(json tst) {
        TestConfig tc(tst);

        add("signal", tc.shape);
        const T& s = arrays["signal"];

        shape_type shape1d = {tc.shape[0]};
        add("total", shape1d);
        T tot = arrays["total"];

        tc.go();

        for (size_t count = 0; count < tc.repeat; ++count) {
            tot = vs::add(tot, vs::median(s, 1));
        }
        vs::sync<T>();

        return tc.results();
    }

    json run_sort(json tst) {

        TestConfig tc(tst);

        add("signal", tc.shape);
        const T& s = arrays["signal"];

        add("total", tc.shape);
        T tot = arrays["total"];

        tc.go();

        for (size_t count = 0; count < tc.repeat; ++count) {
            tot = vs::add(tot, vs::sort(s, 1));
        }
        vs::sync<T>();

        return tc.results();
    }
    
};


#ifdef USE_ARRAYFIRE
#include <arrayfire.h>

struct PerfAF : public PerfT<af::array> {
    
    virtual json prerun(json tst) {
        const std::string device = getdef<std::string>(tst, "device", "cpu");

        if (device == "gpu" || device == "cuda") {
            af::setBackend(AF_BACKEND_CUDA);
        }
        // I can not get an OpenCL stack build working.
        // else if (device == "opencl") {
        //     af::setBackend(AF_BACKEND_OPENCL);
        // }
        else if (device == "cpu") {
            af::setBackend(AF_BACKEND_CPU);
        }
        else {
            throw (std::runtime_error("unsupported device: " + device));
        }
        tst["tech"] = "af";
        return tst;
    }

    virtual void add(const std::string& name, const shape_type& shape) {
        auto data = randu(shape_size(shape));
        af::dim4 dims(shape.size(), shape.data());
        arrays[name] = af::array(dims, data.get()); // type from float*
    }

    virtual void add(const std::string& name, const shape_type& shape, const float* data) {
        af::dim4 dims(shape.size(), shape.data());
        arrays[name] = af::array(dims, data); // type from float*
    }

};


namespace vs {


template<> shape_type shape(const af::array& arr)
{
    auto dims = arr.dims();
    const size_t n = dims.ndims();
    shape_type ret(n);
    for (size_t ind=0; ind<n; ++ind) {
        ret[ind] = dims[ind];
    }
    return ret;
}

template<>
float to_float<af::array>(const af::array& arr)
{
    return arr.scalar<float>();
}

template<>
af::array add<af::array>(const af::array& arr1, const af::array& arr2)
{
    return arr1 + arr2;
}

template<>
af::array mul<af::array>(const af::array& arr1, const af::array& arr2)
{
    return arr1 * arr2;
}
    
template<>
af::array mul<af::array>(const af::array& arr1, float val)
{
    return val*arr1;
}

template<>
float norm<af::array>(const af::array& arr)
{
    return to_float(af::max(arr));
}


template<>
af::array real<af::array>(const af::array& arr)
{
    return af::real(arr);
}
template<>
af::array imag<af::array>(const af::array& arr)
{
    return af::imag(arr);
}

template<>
af::array fft2<af::array>(const af::array& arr)
{
//    return af::fft2(arr);
    auto got = af::fftR2C<2>(arr);
    
    auto ret = af::constant(0, arr.dims(), c32); // c32 defined in global namespace!!!
    ret(af::seq(got.dims(0)), af::seq(got.dims(1))) = got;

    return ret;
}
template<>
af::array ifft2<af::array>(const af::array& arr)
{
    return af::ifft2(arr);
}
template<>
af::array median<af::array>(const af::array& arr, int dim)
{
    return af::median(arr, dim);
}
    

template<> af::array& eval<af::array>(af::array& arr)
{
    return af::eval(arr);
}
template<> void sync<af::array>()
{
    af::sync();
}

template<> af::array sort<af::array>(const af::array& arr, int dim)
{
    return af::sort(arr, dim);
}


} // namespace vs

#endif  // USE_ARRAYFIRE

#ifdef USE_LIBTORCH
#include "torch/torch.h"

struct PerfLT : public PerfT<torch::Tensor> {

    std::string device{"cpu"};

    virtual json prerun(json tst) {
        device = getdef<std::string>(tst, "device", "cpu");
        tst["tech"] = "lt";
        return tst;
    }

    virtual void add(const std::string& name, const shape_type& shape) {
        auto data = randu(shape_size(shape));
        std::vector<int64_t> tshape(shape.begin(), shape.end());
        auto borrowed = torch::from_blob((void*)(data.get()), tshape);
        borrowed = borrowed.clone();
        if (device == "gpu") {
            borrowed = borrowed.to(torch::kCUDA);
        }
        arrays[name] = borrowed.clone();
    }

    virtual void add(const std::string& name, const shape_type& shape, const float* data) {
        std::vector<int64_t> tshape(shape.begin(), shape.end());
        auto borrowed = torch::from_blob((void*)(data), tshape);
        borrowed = borrowed.clone();
        if (device == "gpu") {
            borrowed = borrowed.to(torch::kCUDA);
        }
        arrays[name] = borrowed.clone();
    }

};


namespace vs {

template<> shape_type shape(const torch::Tensor& arr)
{
    auto dims = arr.sizes();
    shape_type ret(dims.begin(), dims.end());
    return ret;
}

template<>
float to_float<torch::Tensor>(const torch::Tensor& arr)
{
    return arr.item<float>();
}

template<>
torch::Tensor add<torch::Tensor>(const torch::Tensor& arr1, const torch::Tensor& arr2)
{
    return arr1 + arr2;
}

template<>
torch::Tensor mul<torch::Tensor>(const torch::Tensor& arr1, const torch::Tensor& arr2)
{
    return arr1 * arr2;
}
    
template<>
torch::Tensor mul<torch::Tensor>(const torch::Tensor& arr1, float val)
{
    return val*arr1;
}

template<>
float norm<torch::Tensor>(const torch::Tensor& arr)
{
    return to_float(arr.max());
}


template<>
torch::Tensor real<torch::Tensor>(const torch::Tensor& arr)
{
    return torch::real(arr);
}
template<>
torch::Tensor imag<torch::Tensor>(const torch::Tensor& arr)
{
    return torch::imag(arr);
}

template<>
torch::Tensor fft2<torch::Tensor>(const torch::Tensor& ten)
{
    return torch::fft::fft2(ten);
}
template<>
torch::Tensor ifft2<torch::Tensor>(const torch::Tensor& ten)
{
    return torch::fft::ifft2(ten);
}
template<>
torch::Tensor median<torch::Tensor>(const torch::Tensor& ten, int dim)
{
    auto tup = ten.median(dim);
    return std::get<0>(tup);
}

template<> torch::Tensor sort<torch::Tensor>(const torch::Tensor& arr, int dim)
{
    auto tup = torch::sort(arr, dim);
    return std::get<0>(tup);
}

// Note, torch IS lazy also, at least for GPU.  But, more is needed to implement
// eval/sync than I want to do right now.
// https://pytorch.org/docs/stable/notes/cuda.html#asynchronous-execution
template<> torch::Tensor& eval<torch::Tensor>(torch::Tensor& arr)
{
    return arr;
}

template<> void sync<torch::Tensor>()
{
    return;
}

} // namespace vs

#endif  // USE_LIBTORCH


#ifdef USE_EIGEN
#include <Eigen/Core>
#include <variant>
#include <algorithm>
#include <fftw3.h>

namespace ei {
    // Type mimicry of the af/lt array handle
    
    using scalar = float;
    using oned = Eigen::ArrayXf;
    // Make 2D arrays RowMajor to cheat a little in sort() and median()
    using twod = Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using twodc = Eigen::Array<std::complex<float>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    using array = std::variant<float, oned, twod, twodc>;
}

struct PerfEI : public PerfT<ei::array> {

    virtual json prerun(json tst) {
        auto device = getdef<std::string>(tst, "device", "cpu");
        if (device != "cpu") {
            throw std::runtime_error("unsupported device: " + device);
        }
        tst["tech"] = "ei";
        return tst;
    }

    virtual void add(const std::string& name, const shape_type& shape) {
        auto data = randu(shape_size(shape));
        return add(name, shape, data.get());
    }

    virtual void add(const std::string& name, const shape_type& shape, const float* data) {
        const size_t ndims = shape.size();
        if (ndims == 2) {
            ei::twod arr = Eigen::Map<ei::twod>((float*)data, shape[0], shape[1]);
            arrays[name] = arr;
        }
        else if (ndims == 1) {
            ei::oned arr = Eigen::Map<ei::oned>((float*)data, shape[0]);
            arrays[name] = arr;
        }
        else {
            throw std::runtime_error("ei::add() does not support shape " + shape_str(shape));
        }
    }

};

namespace vs {

template<> shape_type shape(const ei::array& arr)
{
    if (std::holds_alternative<ei::twodc>(arr)) {
        const auto& tarr = std::get<ei::twodc>(arr);
        return shape_type{tarr.rows(), tarr.cols()};
    }
    if (std::holds_alternative<ei::twod>(arr)) {
        const auto& tarr = std::get<ei::twod>(arr);
        return shape_type{tarr.rows(), tarr.cols()};
    }
    if (std::holds_alternative<ei::oned>(arr)) {
        const auto& tarr = std::get<ei::oned>(arr);
        return shape_type{tarr.size()};
    }
    return shape_type{};
}


template<> float to_float<ei::array>(const ei::array& arr)
{
    if (std::holds_alternative<ei::twodc>(arr)) {
        return std::real(std::get<ei::twodc>(arr)(0,0));
    }
    if (std::holds_alternative<ei::twod>(arr)) {
        return std::get<ei::twod>(arr)(0,0);
    }
    if (std::holds_alternative<ei::oned>(arr)) {
        return std::get<ei::oned>(arr)(0);
    }
    return std::get<ei::scalar>(arr);
}

template<>
ei::array add<ei::array>(const ei::array& arr1, const ei::array& arr2)
{
    if (std::holds_alternative<ei::twodc>(arr1)) {
        ei::twodc ret = std::get<ei::twodc>(arr1) + std::get<ei::twodc>(arr2);
        return ret;
    }
    if (std::holds_alternative<ei::twod>(arr1)) {
        ei::twod ret = std::get<ei::twod>(arr1) + std::get<ei::twod>(arr2);
        return ret;
    }
    if (std::holds_alternative<ei::oned>(arr1)) {
        ei::oned ret = std::get<ei::oned>(arr1) + std::get<ei::oned>(arr2);
        return ret;
    }
    return std::get<ei::scalar>(arr1) + std::get<ei::scalar>(arr2);
}

template<>
ei::array mul<ei::array>(const ei::array& arr1, const ei::array& arr2)
{
    if (std::holds_alternative<ei::twodc>(arr1)) {
        ei::twodc ret = std::get<ei::twodc>(arr1) * std::get<ei::twodc>(arr2);
        return ret;
    }
    if (std::holds_alternative<ei::twod>(arr1)) {
        ei::twod ret = std::get<ei::twod>(arr1) * std::get<ei::twod>(arr2);
        return ret;
    }
    if (std::holds_alternative<ei::oned>(arr1)) {
        ei::oned ret = std::get<ei::oned>(arr1) * std::get<ei::oned>(arr2);
        return ret;
    }
    return std::get<ei::scalar>(arr1) * std::get<ei::scalar>(arr2);
}
    
template<>
ei::array mul<ei::array>(const ei::array& arr1, float val)
{
    if (std::holds_alternative<ei::twodc>(arr1)) {
        ei::twodc ret = std::get<ei::twodc>(arr1) * val;
        return ret;
    }
    if (std::holds_alternative<ei::twod>(arr1)) {
        ei::twod ret = std::get<ei::twod>(arr1) * val;
        return ret;
    }
    if (std::holds_alternative<ei::oned>(arr1)) {
        ei::oned ret = std::get<ei::oned>(arr1) * val;
        return ret;
    }
    return std::get<ei::scalar>(arr1) * val;
}

template<> float norm<ei::array>(const ei::array& arr)
{
    if (std::holds_alternative<ei::twodc>(arr)) {
        return std::get<ei::twodc>(arr).abs().maxCoeff();
    }
    if (std::holds_alternative<ei::twod>(arr)) {
        return std::abs(std::get<ei::twod>(arr).maxCoeff());
    }
    if (std::holds_alternative<ei::oned>(arr)) {
        return std::abs(std::get<ei::oned>(arr).maxCoeff());
    }
    return std::abs(std::get<ei::scalar>(arr));
}


template<> 
ei::array real<ei::array>(const ei::array& arr)
{
    if (std::holds_alternative<ei::twodc>(arr)) {
        ei::twod ret = std::get<ei::twodc>(arr).real();
        return ret;
    }
    return arr;
}

template<> 
ei::array imag<ei::array>(const ei::array& arr)
{
    if (std::holds_alternative<ei::twodc>(arr)) {
        ei::twod ret = std::get<ei::twodc>(arr).imag();
        return ret;
    }
    return arr;
}
    
template<> 
ei::array fft2<ei::array>(const ei::array& arr)
{
    if (std::holds_alternative<ei::twodc>(arr)) {
        ei::twodc tarr = std::get<ei::twodc>(arr);

        ei::twodc tout = ei::twodc::Zero(tarr.rows(), tarr.cols());

        fftwf_complex* src = const_cast<fftwf_complex*>( reinterpret_cast<const fftwf_complex*>(tarr.data()) );
        fftwf_complex* dst = reinterpret_cast<fftwf_complex*>(tout.data());

        fftwf_plan plan = fftwf_plan_dft_2d(tarr.cols(), tarr.rows(), src, dst, 
                                            FFTW_FORWARD, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
        fftwf_execute_dft(plan, src, dst);
        return tout;
    }
    if (std::holds_alternative<ei::twod>(arr)) { // R2C
        auto tarr = std::get<ei::twod>(arr);

        ei::twodc tout = ei::twod::Zero(tarr.rows(), tarr.cols());

        float* src = const_cast<float*>( reinterpret_cast<const float*>(tarr.data()) );
        fftwf_complex* dst = reinterpret_cast<fftwf_complex*>(tout.data());

        fftwf_plan plan = fftwf_plan_dft_r2c_2d(tarr.cols(), tarr.rows(), src, dst, 
                                                FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);

        fftwf_execute(plan); // use fftw_execute_dft_r2c() if plan was predefined on another array
        return tout;

    }
    throw std::runtime_error("ei::fft2 unsupported array type");
}
    
template<> 
ei::array ifft2<ei::array>(const ei::array& arr)
{
    ei::twodc tarr;
    if (std::holds_alternative<ei::twodc>(arr)) {
        tarr = std::get<ei::twodc>(arr);
    }
    else if (std::holds_alternative<ei::twod>(arr)) {
        auto twod = std::get<ei::twod>(arr);
        tarr = twod.cast<std::complex<float>>();
    }
    else {
        throw std::runtime_error("ei::fft2 unsupported array type");
    }

    ei::twodc tout = ei::twodc::Zero(tarr.rows(), tarr.cols());

    fftwf_complex* src = const_cast<fftwf_complex*>( reinterpret_cast<const fftwf_complex*>(tarr.data()) );
    fftwf_complex* dst = reinterpret_cast<fftwf_complex*>(tout.data());

    fftwf_plan plan = fftwf_plan_dft_2d(tarr.cols(), tarr.rows(), src, dst, 
                                        FFTW_BACKWARD, FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);

    fftwf_execute_dft(plan, src, dst);
    return tout;
}
    
template<> 
ei::array median<ei::array>(const ei::array& arr, int /*dim*/)
{
    if (std::holds_alternative<ei::scalar>(arr)) {
        return arr;
    }
    if (std::holds_alternative<ei::oned>(arr)) {    
        ei::oned tarr = std::get<ei::oned>(arr); // copy as nth_element reorders
        const size_t siz = tarr.size();
        float* data = tarr.data();
        const double p = 0.5;
        size_t mid = p * siz;
        mid = std::min(mid, siz - 1);
        std::nth_element(data, data+mid, data+siz);
        return *(data+mid);
    }
    if (std::holds_alternative<ei::twod>(arr)) {    
        const ei::twod& tarr = std::get<ei::twod>(arr);
        size_t nrows = tarr.rows();
        ei::oned ret = ei::oned::Zero(nrows);
        for (size_t irow=0; irow < nrows; ++irow) {
            const ei::oned& row = tarr.row(irow);
            const ei::array& rarr = row;
            ei::array marr = median(rarr);
            ret[irow] = std::get<ei::scalar>(marr);
        }
        return ret;
    }

    return arr;
}
    
template<> 
ei::array sort<ei::array>(const ei::array& arr, int /*dim*/)
{
    if (std::holds_alternative<ei::oned>(arr)) {
        ei::oned tarr = std::get<ei::oned>(arr);
        std::sort(tarr.data(), tarr.data() + tarr.size());
        return tarr;
    }
    if (std::holds_alternative<ei::twod>(arr)) {    
        ei::twod tarr = std::get<ei::twod>(arr);
        const size_t nrows = tarr.rows();
        const size_t ncols = tarr.cols();
        for (size_t irow=0; irow<nrows; ++irow) {
            float* data = tarr.row(irow).data();
            std::sort(data, data + ncols);
        }
        return tarr;
    }
    throw std::runtime_error("ei::sort unsupported array type");
}

template<> ei::array& eval<ei::array>(ei::array& arr) { return arr; } // no-op
template<> void sync<ei::array>() { } // no-op
    
} 

#endif // USE_EIGEN


json perf(json tsts)
{
    json results;

    for (auto& tst : tsts) {
        auto tech = getdef<std::string>(tst, "tech", "af");
        std::unique_ptr<Perf> p;

        try {
            if (tech == "af") {
                p = std::make_unique<PerfAF>();
            }
            else if (tech == "lt") {
                p = std::make_unique<PerfLT>();
            }
            else if (tech == "ei") {
                p = std::make_unique<PerfEI>();
            }
        }
        catch (std::runtime_error& err) {
            std::cerr << "skipping failed construction: " << tech << "\n" << tst << "\n";
            continue;
        }
        if (!p) {
            std::cerr << "skipping unknown tech: " << tech << "\n" << tst << "\n";
            continue;
        }

        tst = p->prerun(tst);

        std::cerr << tst << "\n";

        json result;
        std::string kind = tst["kind"];
        if (kind == "arith") {
            result = p->run_arith(tst);
        }
        else if (kind == "convo") {
            result = p->run_convo(tst);
        }
        else if (kind == "median") {
            result = p->run_median(tst);
        }
        else if (kind == "sort") {
            result = p->run_sort(tst);
        }
        else {
            std::cerr << "unknown test: " << kind << " config: " << tst << "\n";
        }
        json jtest;
        jtest["input"] = tst;
        jtest["output"] = result;
        results.push_back(jtest);
    }
    return results;
}

int main(int argc, char* argv[])
{
    std::string iname = "/dev/stdin";
    std::string oname = "/dev/stdout";
    if (argc > 1) {
        iname = argv[1];
        if (iname == "-") iname = "/dev/stdin";
    }
    if (argc > 2) {
        oname = argv[2];
        if (oname == "-") oname = "/dev/stdout";
    }
        
    std::ifstream ifile(iname);
    auto tsts = json::parse(ifile);

    auto res = perf(tsts);

    std::ofstream ofile(oname);
    ofile << std::setw(4) << res << std::endl;
    return 0;    
}
