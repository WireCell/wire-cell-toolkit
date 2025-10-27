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
#include <complex>

namespace WireCell::SPNG {

    /// Return a torch device as a string.
    std::string to_string(const torch::Device& device);

    /// Return Tensor::sizes() as a string
    std::string to_string(const torch::IntArrayRef& sizes);

    /// Return torch::Tensor as string
    std::string to_string(const torch::Tensor& tensor);

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
    // Caution: torch has layout requirements for complex tensors.
    template <> struct CppTypeToTorchDtype<std::complex<float>> { static constexpr torch::Dtype value = torch::kComplexFloat; };
    template <> struct CppTypeToTorchDtype<std::complex<double>> { static constexpr torch::Dtype value = torch::kComplexDouble; };
    // Be thoughtful when extending this list.

    // Cast torch tensor to dtype with C++ type.
    template <typename T>
    torch::Tensor to_dtype(const torch::Tensor& tensor) {
        return tensor.to(CppTypeToTorchDtype<T>::value);
    }


    /// Convert a 1D tensor to a C++ vector of a given type.
    ///
    /// If you REALLY want to flatten an N-dimensional tensor to vector, pass
    /// nd_okay=true.
    template<typename T>
    std::vector<T> to_vector(const torch::Tensor& tensor, bool nd_okay=false) {
        if (!nd_okay && tensor.dim() != 1) {
            raise<ValueError>("to_vector given %d-dimensional tensor", tensor.dim());
        }

        torch::Tensor row_tensor = to_dtype<T>(tensor).to(torch::kCPU).contiguous();
        return std::vector<T>(row_tensor.data_ptr<T>(),
                              row_tensor.data_ptr<T>() + row_tensor.numel());
    }

    /// Helper to make a 1D tensor from a vector.  This does a copy.
    template<typename T>
    torch::Tensor to_tensor(const std::vector<T>& vec) {
        return torch::from_blob((void*)vec.data(), {(long)vec.size()},
                                CppTypeToTorchDtype<T>::value).clone();
    }


    // Return a sampled, normalized 1D Gausian pdf.
    torch::Tensor gaussian1d(double mean, double sigma,
                             int64_t npoints, double xmin, double xmax,
                             torch::TensorOptions options = torch::TensorOptions());

    // Return a tensor sizes() as a vector.
    // Actually, can just as easily to tensor.sizes().vec();
    inline
    std::vector<int64_t> vshape(const torch::IntArrayRef& sizes) {
        return std::vector<int64_t>(sizes.begin(), sizes.end());
    }
    inline
    std::vector<int64_t> vshape(const torch::Tensor& ten) {
        return vshape(ten.sizes());
    }


    /// Return true if any element is a NaN.
    inline
    bool has_nan(const torch::Tensor& tensor) {
        return tensor.isnan().any().item<bool>();
    }


    /**
    * @brief Calculates the number of unique, non-redundant positive frequency bins in a discrete spectrum of size N.
    *
    * @details Determines the "half size" count of spectral components, excluding both the
    *          Direct Current (DC) component (zero frequency) and the Nyquist frequency bin,
    *          if one exists (i.e., when N is even).
    *          This calculation yields the number of independent positive-frequency bins
    *          and is typically equivalent to floor(N / 2) - 1.
    *
    * @param N The full size or length of the discrete spectrum (e.g., the FFT size).
    * @return size_t The number of unique positive frequency samples (excluding DC and Nyquist).
    */
    inline size_t nhalf(size_t N)
    {
        if (N%2) {
            return (N-1)/2;     // odd
        }
        return (N-2)/2;         // even
    }

    /**
     * @brief Return the "middle" index of a tensor dimension.
     *
     * This is the index of the first element past the "half".
     *
     * Eg, if N is an even size of a spectrum, the middle index is that of the
     * Nyquist bin, odd size it's the sample at the Nyquist frequency.  Or,
     * if it is an odd size interval waveform, the middle index is the upper
     * half, not counting the sample zero.
     */
    inline size_t middle_index(size_t N) { return 1 + nhalf(N); }

    /// Return (positive) "a mod b".
    inline int64_t modulo(int64_t a, int64_t b) {
        return ((a % b) + b) % b;
    }

    /**
     * @brief Resize a tensor along one tensor dimension.
     *
     * @param in A tensor to resize.
     * @param axis The dimension to resize.
     * @param index The index along the dimension at which to apply resizing.
     * @param size The desired size of the resized tensor dimension.
     *
     * The index gives the location at which resizing occurs.  Padding will be
     * inserted just before this index.  Likewise, truncation will start just
     * before this index (element at index will be the first to be removed).
     *
     * When the resize is asked to enlarge at an index that is equal to the size
     * of the input tensor dimension, padding will be appended to the dimension.
     * Truncating in this case is equivalent passing an index of 0.
     *
     * When truncation requires removing more elements than are at or above the
     * given index, the removal will wrap-around and remove elements from the
     * front of the tensor dimension.
     *
     * A negative index is mapped to the positive range to allow Numpy array
     * style indexing.
     *
     * See resize_tensor_head() and resize_tensor_tail() for special case
     * resizing at the front or end of the tensor dimension.
     */
    torch::Tensor resize_tensor(const torch::Tensor& in,
                                int64_t axis, int64_t index, int64_t size);

    /**
     * @brief Resize a tensor along one dimension by operating at its end.
     *
     * This either truncates or pads the end of the tensor dimension.
     *
     * Note, truncation here is different than passing Ns as index to
     * resize_tensor() as that would "wrap around" and actually do a head
     * truncation.  Instead, the index below Ns is found.  Tail padding is not
     * specialized.
     */
    torch::Tensor resize_tensor_tail(const torch::Tensor& in,
                                     int64_t axis, int64_t size);
    /**
     * @brief Resize a tensor along one dimension by operating at its beginning.
     *
     * This will either pad or truncate the beginning of the tensor dimension.
     * It merely passes index=0 to the general resize_tensor.  It is not
     * specialized and is provided for symmetry with resize_tensor_tail().
     */
    torch::Tensor resize_tensor_head(const torch::Tensor& in,
                                     int64_t axis, int64_t size);



    /**
     * @brief Resize the "middle" of a tensor along.
     *
     * This will find the "middle index" and call the general resize_tensor().
     */
    torch::Tensor resize_tensor_middle(const torch::Tensor& in,
                                     int64_t axis, int64_t size);




    // Note, removed pad() as it was too limiting and was never used.  Use
    // torch::nn::pad().  See spng/test/doctest_simple_convo.cxx for one
    // example.

}

#endif
