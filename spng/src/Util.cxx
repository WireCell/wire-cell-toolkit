#include "WireCellSpng/Util.h"
#include <cmath>
using namespace WireCell;

torch::Tensor Torch::gaussian1d(double mean, double sigma,
                                int64_t npoints, double xmin, double xmax,
                                torch::TensorOptions options)
{
    auto x = torch::linspace(xmin, xmax, npoints, options);
    auto rel = (x - mean)/sigma;
    const double norm = sqrt(2*M_PI*sigma*sigma);
    return norm * torch::exp(-0.5*rel*rel);
}


std::vector<int64_t> Torch::linear_shape(const std::vector<torch::Tensor>& tens, 
                                         torch::IntArrayRef extra_shape)
{
    // Find shape that assures linear convolution.
    std::vector<int64_t> shape = {extra_shape[0], extra_shape[1]};
    for (const auto& ten : tens) {
        auto sizes = ten.sizes();
        shape[0] += sizes[0] - 1;
        shape[1] += sizes[1] - 1;
    }
    return shape;
}


torch::Tensor Torch::pad(torch::Tensor ten, double value, torch::IntArrayRef shape)
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

torch::Tensor Torch::convo_spec(const std::vector<torch::Tensor>& tens, 
                                torch::IntArrayRef shape)
{
    using torch::indexing::Slice;

    // Return value will be complex spectrum in Fourier domain.
    torch::Tensor fourier;
    // Allocate working array in interval domain.
    torch::Tensor interval = torch::zeros(shape, tens[0].options());

    const size_t ntens = tens.size();

    // First, accumulate denominator.
    for (size_t ind=0; ind<ntens; ++ind) {
        // Caveat: zero-padding is not always appropriate for every ten in tens.
        interval.zero_();

        const auto& ten = tens[ind];
        auto s = ten.sizes();
        interval.index_put_({
                Slice(0,std::min(s[0], shape[0])),
                Slice(0,std::min(s[1], shape[1]))
            }, ten);

        if (ind == 0) {
            fourier = torch::fft::fft2(interval);
        }
        else {
            fourier *= torch::fft::fft2(interval);
        }
    }

    return fourier;
}


torch::Tensor Torch::filtered_decon_2d(const std::vector<torch::Tensor>& numerator,
                                       const std::vector<torch::Tensor>& denominator,
                                       torch::IntArrayRef shape)
{
    if (denominator.empty()) {  // graceful degradation 
        return convo_spec(numerator, shape);
    }

    // Note: this suffers holding an extra array which could be avoided by
    // accumulating denominator convolutions, inverting that result and
    // continuing accumulating numerator convolutions.
    auto num = convo_spec(numerator, shape);
    auto den = convo_spec(denominator, shape);
    return torch::divide(num, den); // fixme, divide-by-zero?
}


torch::Tensor Torch::filtered_decon_2d_auto(const std::vector<torch::Tensor>& numerator,
                                          const std::vector<torch::Tensor>& denominator,
                                          torch::IntArrayRef extra_shape)
{
    std::vector<torch::Tensor> all_in(numerator.begin(), numerator.end());
    all_in.insert(all_in.end(), denominator.begin(), denominator.end());

    if (all_in.empty()) {
        return torch::Tensor();
    }

    auto shape = linear_shape(all_in, extra_shape);

    return filtered_decon_2d(numerator, denominator, shape);
}

void WireCell::SPNG::metadata_passthrough(
    const WireCell::Configuration & metadata_in,
    WireCell::Configuration & metadata_out,
    const Json::Value & passing_values) {
        // Throw
        // if (!passing_values.isArray()) {
        // }

        // std::cout << "Passing values:" << passing_values << std::endl;
        for (Json::ArrayIndex i = 0; i < passing_values.size(); ++i) {
            const auto & name = passing_values[i].asString();
            // std::cout << "Passing " << name << std::endl;
            // std::cout << "In input? " << metadata_in.isMember(name) << std::endl;
            metadata_out[name] = metadata_in[name];
        }
}
