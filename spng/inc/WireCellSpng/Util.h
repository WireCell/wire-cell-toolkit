/**
   Utility functions for Torch tensors.

   This is meant to collect small utility functions.

   The intention is to hold GENERIC and GENERAL PURPOSE functions that augment
   what Torch provides.

   This should not have WCT specific functions not any functions that are very
   application specific.

   This header should not depend on ITorchTensor!  See TensorTools instead.

 */

#ifndef WIRECELL_SPNG_UTIL
#define WIRECELL_SPNG_UTIL

#include "WireCellSpng/Torch.h"
#include "WireCellUtil/Exceptions.h"

#include <type_traits> // For std::enable_if and std::is_same

namespace WireCell::SPNG {

    /// Return a torch device as a string.
    std::string to_string(const torch::Device& device);

    // Map C++ types to torch::Dtype
    template <typename T, typename Enable=void>
    struct CppTypeToTorchDtype { }; // error about a .value missing?  You are using an usupported dtype.

    // Specializations for types expected in SPNG with some examples where they may be found.
    // "traces"
    template <> struct CppTypeToTorchDtype<float> { static constexpr torch::Dtype value = torch::kFloat32; };
    // "summaries"
    template <> struct CppTypeToTorchDtype<double> { static constexpr torch::Dtype value = torch::kFloat64; };
    // "chids", "chmasks"
    template <> struct CppTypeToTorchDtype<int> { static constexpr torch::Dtype value = torch::kInt32; };
    // indices 
    template <> struct CppTypeToTorchDtype<long> { static constexpr torch::Dtype value = torch::kInt64; };
    // pixel masks
    template <> struct CppTypeToTorchDtype<bool> { static constexpr torch::Dtype value = torch::kBool; };
    // Be thoughtful when extending this list.

    // Cast torch tensor to dtype with C++ type.
    template <typename T>
    torch::Tensor to_dtype(const torch::Tensor& tensor) {
        return tensor.to(CppTypeToTorchDtype<T>::value);
    }


    /// Convert a 1D tensor to a C++ vector of a given type.
    template<typename T>
    std::vector<T> to_vector(const torch::Tensor& tensor) {
        torch::Tensor row_tensor = to_dtype<T>(tensor).to(torch::kCPU).contiguous();
        return std::vector<T>(row_tensor.data_ptr<T>(),
                              row_tensor.data_ptr<T>() + row_tensor.numel());
    }

    // Return a sampled, normalized 1D Gausian pdf.
    torch::Tensor gaussian1d(double mean, double sigma,
                             int64_t npoints, double xmin, double xmax,
                             torch::TensorOptions options = torch::TensorOptions());

    // Return a tensor sizes() as a vector.
    inline
    std::vector<int64_t> vshape(const torch::IntArrayRef& sizes) {
        return std::vector<int64_t>(sizes.begin(), sizes.end());
    }
    inline
    std::vector<int64_t> vshape(const torch::Tensor& ten) {
        return vshape(ten.sizes());
    }

    
    // Note, removed pad() as it was too limiting and was never used.  Use
    // torch::nn::pad().  See spng/test/doctest_simple_convo.cxx for one
    // example.

}

#endif
