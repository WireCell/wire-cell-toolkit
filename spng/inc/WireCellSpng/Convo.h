/** Functions related to convolutions.

 */

#ifndef WIRECELL_SPNG_CONVO
#define WIRECELL_SPNG_CONVO

#include "WireCellSpng/Util.h"
#include "WireCellUtil/Exceptions.h"

namespace WireCell::SPNG {


    /// Return a shape large enough to assure linear convolution if the given
    /// tensors and another tensor of shape extra_shape were cyclically
    /// convolved.
    ///
    /// Scalar form
    template<typename T>
    T linear_shape(const T& a, const T& b) { return a + b - 1; }

    /// Vector form
    template<typename T>
    std::vector<T> linear_shape(const std::vector<T>& a, const std::vector<T>& b)
    {
        const size_t ndima = a.size();
        const size_t ndimb = b.size();
        if (ndima != ndimb) {
            raise<ValueError>("linear_shape ndim mismatch %d != %d", ndima, ndimb);
        }
        std::vector<T> res(ndima);
        for (size_t ind=0; ind<ndima; ++ind) {
            res[ind] = linear_shape(a[ind], b[ind]);
        }
        return res;
    }        

    inline
    std::vector<int64_t> linear_shape(const torch::IntArrayRef& a, const torch::IntArrayRef& b)
    {
        return linear_shape(vshape(a), vshape(b));
    }
    // A cumulative version.
    std::vector<int64_t> linear_shape(const std::vector<torch::Tensor>& tens, 
                                      torch::IntArrayRef extra_shape);
    std::vector<int64_t> linear_shape(const std::vector<torch::Tensor>& tens, 
                                      std::vector<int64_t> extra_shape);

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
#endif
