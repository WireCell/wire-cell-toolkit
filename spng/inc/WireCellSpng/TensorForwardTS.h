#pragma once

#include "WireCellSpng/ITensorForward.h"
#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/HanaJsonCPP.h"

namespace WireCell::SPNG {

    struct TensorForwardTSConfig {
        /// Required.  File name for a torchscript file.
        std::string ts_filename="";
    };

}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TensorForwardTSConfig, ts_filename);

namespace WireCell::SPNG {
    /// A "service" type component that sends a tensor though torch script.
    ///
    struct TensorForwardTS : public Logger, public ContextBase,
                             public virtual IConfigurable,
                             public virtual ITensorForward
    {
        TensorForwardTS();
        virtual ~TensorForwardTS() = default;
        
        // IConfigurable - see TensorForwardTSConfig for configuration 
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

        // ITensorForward
        virtual ITorchTensor::pointer forward(const ITorchTensor::pointer& input) const;

    private:
        TensorForwardTSConfig m_config;
        mutable torch::jit::script::Module m_module;
    };
                             

}
