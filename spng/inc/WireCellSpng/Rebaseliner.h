#pragma once

#include "WireCellSpng/TensorFilter.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {
    struct RebaselinerConfig {
        // The dimension to rebin
        int dim = -1;

        // The threshold to detect that a sample is inside a region of interest.
        float threshold=0;

    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::RebaselinerConfig, dim, threshold);

namespace WireCell::SPNG {

    struct Rebaseliner : public TensorFilter,
                         virtual public IConfigurable {

        Rebaseliner();
        virtual ~Rebaseliner() = default ;

        /// TensorFilter
        virtual ITorchTensor::pointer filter_tensor(const ITorchTensor::pointer& in);

        // IConfigurable - see KernelConvolveConfig for configuration documentation.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        
    private:

        RebaselinerConfig m_config;

    };
}
