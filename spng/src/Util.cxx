#include "WireCellSpng/Util.h"
#include "WireCellSpng/Stream.h"
#include <cmath>
#include <sstream>


namespace WireCell::SPNG {

    std::string to_string(const torch::Device& device)
    {
        std::stringstream ss;
        ss << device;
        return ss.str();        
    }

    std::string to_string(const torch::IntArrayRef& sizes)
    {
        std::stringstream ss;
        ss << "(";
        std::string comma="";
        for (const auto& size : sizes) {
            ss << comma << size;
            comma = ", ";
        }
        ss << ")";
        return ss.str();        
    }

    std::string to_string(const torch::Tensor& tensor, bool content)
    {
        std::stringstream ss;
        /// Using this to show min/max is too restrictive, eg bool on cpu lacks abs().
        // auto tabs = torch::abs(tensor);
        // ss << to_string(tensor.sizes()) << " dtype:" << tensor.dtype() << " vmm:["
        //    << torch::min(tabs).item() << ", " << torch::max(tabs).item() << "]";
        ss << to_string(tensor.sizes()) << " dtype:" << tensor.dtype();
        if (content) {
            ss << " " << tensor;
        }
        return ss.str();
    }

    torch::Tensor gaussian1d(double mean, double sigma,
                             int64_t npoints, double xmin, double xmax,
                             torch::TensorOptions options)
    {
        auto x = torch::linspace(xmin, xmax, npoints, options);
        auto rel = (x - mean)/sigma;
        const double norm = sqrt(2*M_PI*sigma*sigma);
        return norm * torch::exp(-0.5*rel*rel);
    }


    torch::Tensor resize_tensor(const torch::Tensor& in,
                                int64_t axis, int64_t index, int64_t size)
    {

        if (in.dim() == 0) return in.clone(); // touches temple

        axis = modulo(axis, in.dim());
            
        std::vector<int64_t> out_size = in.sizes().vec();
        const int64_t Ns = out_size[axis];
        const int64_t Nr = size; // just for symmetry in variable names
        if (in.dim() > 0) {
            out_size[axis] = Nr;
        }

        // Support for Numpy array style negative indices.
        if (index < 0) {
            index = modulo(index, Ns);
        }
        TORCH_CHECK(index <= Ns, "resize_tensor expects axis in [0,Ns] inclusive.");

        torch::Tensor rs = torch::zeros(out_size, in.options());

        if (Ns > Nr) {          // truncation

            const int64_t nloss = Ns - Nr;

            /// Indices of range of original dimension subject to truncation
            const int beg = modulo(index, Ns);
            const int end = modulo(index+nloss, Ns);

            if (beg < end) {    // slice internal, clip middle, keep ends.
                if (beg) { // low-side end, if not empty
                    rs.narrow(axis, 0, beg).copy_(in.narrow(axis, 0, beg));
                }
                if (end < Ns) { // high-side end, if nto empty
                    rs.narrow(axis, beg, Ns-end).copy_(in.narrow(axis, end, Ns-end));
                }
            }
            else if (end < beg) { // slice wraps around, clip the ends, keep middle
                // end here is start, beg is, er, the ending.
                rs.copy_(in.narrow(axis, end, Nr));
            }
        }
        else if (Ns < Nr) {     // padding

            if (index) {        // we are NOT pre-padding, copy first part of input
                // Note, this also captures the post-padding case when index == Ns
                rs.narrow(axis, 0, index).copy_(in.narrow(axis, 0, index));
            }
            if (index < Ns) {   // we are NOT post-padding, copy last part of input
                const int64_t nkeep = Ns - index;
                const int64_t npad = Nr - Ns;
                const int64_t rindex = index + npad;
                const int64_t sindex = index;
                rs.narrow(axis, rindex, nkeep).copy_(in.narrow(axis, sindex, nkeep));
            }
        }
        else {
            return in.clone();
        }

        return rs;
    }

    torch::Tensor resize_tensor_tail(const torch::Tensor& in,
                                     int64_t axis, int64_t size)
    {
        const int64_t Ns = in.size(axis);
        if (size < Ns) { // truncate end
            return resize_tensor(in, axis, size-Ns, size);
        }
        // pad end
        return resize_tensor(in, axis, Ns, size);
    }

    torch::Tensor resize_tensor_head(const torch::Tensor& in,
                                     int64_t axis, int64_t size)
    {
        return resize_tensor(in, axis, 0, size);

    }


    torch::Tensor resize_tensor_middle(const torch::Tensor& in,
                                     int64_t axis, int64_t size)
    {
        return resize_tensor(in, axis, middle_index(in.size(axis)), size);
    }

    std::string tensor_shape_string(const at::Tensor& t) {
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < t.sizes().size(); ++i) {
            oss << t.sizes()[i];
            if (i != t.sizes().size() - 1)
                oss << ", ";
        }
        oss << "]";
        return oss.str();
    }

    void write_torch_to_npy(const torch::Tensor &ten, const std::string &filename){
        torch::save(ten,filename);
    }


}


