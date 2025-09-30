#include "WireCellSpng/Convo.h"

namespace WireCell::SPNG {

    std::vector<int64_t> linear_shape(const std::vector<torch::Tensor>& tens, 
                                      torch::IntArrayRef extra_shape)
    {
        return linear_shape(tens, vshape(extra_shape));
    }
    std::vector<int64_t> linear_shape(const std::vector<torch::Tensor>& tens, 
                                      std::vector<int64_t> shape)
    {
        // start with whatever user may give in input "extra" shape
        for (const auto& ten : tens) {
            auto sizes = ten.sizes();
            const int64_t ndim = sizes.size();
            for (int64_t ind=0; ind<ndim; ++ind) {
                shape[ind] += sizes[0] - 1;
            }
        }
        return std::move(shape);
    }


    torch::Tensor convo_spec(const std::vector<torch::Tensor>& tens, 
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


    torch::Tensor filtered_decon_2d(const std::vector<torch::Tensor>& numerator,
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


    torch::Tensor filtered_decon_2d_auto(const std::vector<torch::Tensor>& numerator,
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
}
