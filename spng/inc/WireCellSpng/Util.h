/**
   Utility functions for Torch tensors.

   This is meant to collect small utilities that "should" be in libtorch, proper.


   In future, this may be moved into the "pytorch/" sub-package.
 */

#ifndef WIRECELLTORCHUTIL
#define WIRECELLTORCHUTIL

#include "WireCellSpng/Torch.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

#include <type_traits> // For std::enable_if and std::is_same

// Capitalized "Torch" namespace should not collide with any in torch, proper.
namespace WireCell::Torch {

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
        return std::vector<T>(row_tensor.data_ptr<float>(),
                              row_tensor.data_ptr<float>() + row_tensor.numel());
    }

    


    // Return 2D tensor of shape with ten filling lower corner and remaining
    // elements holding the given value.
    torch::Tensor pad(torch::Tensor ten, double value, torch::IntArrayRef shape);

    // Return a sampled, normalized 1D Gausian pdf.
    torch::Tensor gaussian1d(double mean, double sigma,
                             int64_t npoints, double xmin, double xmax,
                             torch::TensorOptions options = torch::TensorOptions());

    // Return a shape large enough to assure linear convolution if the given
    // tensors and another tensor of shape extra_shape were cyclically
    // convolved.
    //
    // Each element of shape is the sum of corresponding dimension size reduced
    // by the number of tensors.  If there is no intention to include an
    // additional tensor of shape extra_shape, this inflates from the strict
    // minimum required size to avoid cyclic artifacts by 1.
    std::vector<int64_t> linear_shape(const std::vector<torch::Tensor>& tens, 
                                      torch::IntArrayRef extra_shape = {0,0});


    // Perform 2D, multi-array cyclic convolution in the given shape.
    //
    // Result is a complex-valued, Fourier-domain 2D array with given shape.
    //
    // In order for the convolution to be linear the shape must be at least as
    // given by linear_shape().
    //
    // The shape of all tens must be inside the shape.
    torch::Tensor convo_spec(const std::vector<torch::Tensor>& tens, 
                             torch::IntArrayRef shape);


    // Perform 2D, multi-array linear convolution / deconvolution.
    //
    // Result is a complex-valued, Fourier-domain 2D array with given shape.
    //
    // Shape should be as from linear_shape() given both sets of filter and
    // response tensors.
    // 
    // All input arrays must be 2D.  To provide a row-array from a 1D tensor,
    // use row.reshape({1,-1}).  To provide a column array from a 1D tensor use
    // col.reshape({-1,1}).
    torch::Tensor filtered_decon_2d(const std::vector<torch::Tensor>& filters,
                                    const std::vector<torch::Tensor>& responses,
                                    torch::IntArrayRef shape);

    // Call filtered_decon_2d with shape calculated from filters, responses and
    // extra_shape.
    torch::Tensor filtered_decon_2d_auto(const std::vector<torch::Tensor>& filters,
                                         const std::vector<torch::Tensor>& responses,
                                         torch::IntArrayRef extra_shape = {0,0});



}

namespace WireCell::SPNG {
    // FIXME: this probably should be moved into WireCell::Torch.
    std::string tensor_shape_string(const torch::Tensor& t);

    // FIXME: Remove.  This already exists Configuration.h as update()
  void metadata_passthrough(
    const WireCell::Configuration & metadata_in,
    WireCell::Configuration & metadata_out,
    const Json::Value & passing_values);

    // FIXME: Remove.  This is not util'ish but io'ish.  It also does not seem to be used anywerhe.
    void save_torchtensor_data(const torch::Tensor& tensor, const std::string& filename);
    void save_simpletensor_data(const ITorchTensorSet::pointer& in, const std::string& filename);

    // FIXME:  Remove, eventually.  Not generic, obsoleted by TDM
    std::vector<torch::IValue> from_itensor(const ITorchTensorSet::pointer& in, bool is_gpu = false);
    ITorchTensorSet::pointer to_itensor(const std::vector<torch::IValue>& inputs);
}
#endif
