/**
   Rebinning is a sibling to Resampling that operates in the interval domain.

   It only supports integer-factor upsampling or downsampling.  When upsampling,
   all final samples corresponding to an original sample take the same value.
   Both upsamping and downsampling are aligned by the integer factor.  When
   downsampling, any residual inputs samples due to original array size not
   being an integer multiple of the rebin factor are ignored.  When upsampling,
   the resulting array is exactly larger than the original by the rebin factor.
 */

#pragma once

#include "WireCellSpng/TensorFilter.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {
    struct RebinnerConfig {
        // The dimension to rebin
        int dim = -1;

        // The integer rebin factor.  Required (must be non-zero).
        // If positive, a downsampling is done and the new sample period is factor-times larger.
        // If negative, an upsampling is done, and the new sample period is factor-times smaller.
        int factor=0;

        // The rebinning normalization interpretation.
        //
        // This determines how the rebinned values are calculated.
        //
        // - integral :: samples represent an integral (or sum) over the sample period.
        //
        // With "integral", a downsampled sample will be a sum over its original
        // samples.  Upsampled samples take a common value that of the original
        // sample's value divided by the rebin factor.
        //
        // - interpolation :: samples represent instantaneous values.
        //
        // A downsampled value will be the average of its corresponding original
        // sample values.  Upsampled samples all take the common value that of
        // the original sample's value.
        //
        // - maximum :: samples represent a maximum over the sample period
        //
        // A Downsampled value will be the maximum of its corresponding original
        // sample values.  Upsampled samples all take the common valuethat of
        // the original sample's value.  This interpretation is useful when
        // rebinning Boolean or "near boolean" (but real) valued tensors.
        std::string norm = "integral";
    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::RebinnerConfig, dim, factor, norm);

namespace WireCell::SPNG {

    struct Rebinner : public TensorFilter,
                      virtual public IConfigurable {

        Rebinner();
        virtual ~Rebinner() = default ;

        /// Tensorfilter
        virtual ITorchTensor::pointer filter_tensor(const ITorchTensor::pointer& in);

        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        
    private:

        RebinnerConfig m_config;

        using rebin_func = std::function<torch::Tensor(const torch::Tensor&)>;
        rebin_func m_rebin_func;

    };
}
