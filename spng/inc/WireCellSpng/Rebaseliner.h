#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/ITorchTensorFilter.h"

#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellIface/IConfigurable.h"

namespace WireCell::SPNG {
    struct RebaselinerConfig {
        // The dimension to rebin
        int dim = -1;

        // The threshold to detect that a sample is inside a region of interest.
        float threshold=0;

        // Set the tensor metadata attribute named 'tag' if not empty.
        std::string tag = "";
    };
}

BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::RebaselinerConfig, dim, threshold, tag);

namespace WireCell::SPNG {

    struct Rebaseliner : public ContextBase,
                         public Logger,
                         public ITorchTensorFilter,
                         virtual public IConfigurable {

        Rebaseliner();
        virtual ~Rebaseliner() = default ;

        /// ITorchTensorFilter
        virtual bool operator()(const input_pointer& in, output_pointer& out);

        // IConfigurable - see KernelConvolveConfig for configuration documentation.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;
        
    private:

        RebaselinerConfig m_config;

    };
}
