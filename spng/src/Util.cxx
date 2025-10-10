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

    std::string to_string(const torch::Tensor& tensor)
    {
        std::stringstream ss;
        ss << to_string(tensor.sizes()) << " <" << tensor.dtype() << ">";
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


}


