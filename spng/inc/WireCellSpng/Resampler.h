// Apply LMN resampling to input tensor.
//
// This is a simpler version of WCT's node of the same name in aux.  Besides
// being Torch based, it operates simply on a single input tensor with no
// considering of physical units (eg time).

#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/TorchLMN.h"
#include "WireCellSpng/ITorchTensorFilter.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {
    struct ResamplerConfig {
        // The dimension to resample
        int dim = -1;

        // The ratio of input sample period to output sample period.  This is
        // equivalently the ratio of the output number of samples to input
        // number of samples.
        //
        // A ratio larger than 1.0 is an upsampling, smaller than 1.0 is
        // downsampling.  A ratio of 1.0 will make the node simply pass through
        // the input untouched.
        double ratio=1.0;

        // The LMN normalization interpretation.  May be "interpolation"
        // (default, array holds instantaneous samples of underlying function),
        // "integral" (array holds samples that integrate the function over one
        // period) or "energy" (array holds probability amplitude).  Note, ADC
        // waveforms should likely be interpolated.  Signals are integrated.
        std::string norm = "interpolation";

        // FIXME: additional options to consider: padding, filtering/tapering
        // schemes and handling extra padding due to LMN condition.
    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::ResamplerConfig,
                        dim, ratio, norm);

namespace WireCell::SPNG {

    struct Resampler : public ContextBase,
                       public Logger,
                       public ITorchTensorFilter,
                       virtual public IConfigurable {

        Resampler();
        virtual ~Resampler() = default ;

        /// ITorchTensorFilter
        virtual bool operator()(const input_pointer& in, output_pointer& out);

        // IConfigurable - see KernelConvolveConfig for configuration documentation.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        
    private:

        ResamplerConfig m_config;
        LMN::Normalization m_norm;
    };
}
