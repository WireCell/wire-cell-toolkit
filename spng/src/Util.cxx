#include "WireCellSpng/Util.h"
#include <cmath>
#include <sstream>


namespace WireCell::SPNG {

    std::string to_string(const torch::Device& device)
    {
        std::stringstream ss;
        ss << device;
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



    torch::Tensor pad(torch::Tensor ten, double value, torch::IntArrayRef shape)
    {
        using torch::indexing::Slice;

        torch::Tensor padded = torch::zeros(shape, ten.options()) + value;
        auto s = ten.sizes();
        padded.index_put_({
                Slice(0,std::min(s[0], shape[0])),
                Slice(0,std::min(s[1], shape[1]))
            }, ten);
        return padded;    
    }

}


