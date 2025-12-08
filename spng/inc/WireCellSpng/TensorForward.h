#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/ITorchTensorFilter.h"
#include "WireCellSpng/ITensorForward.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/HanaJsonCPP.h"


namespace WireCell::SPNG {

    struct TensorForwardConfig {
        // type/name of an ITensorForward such as TensorForwardTS
        std::string forward=""; 
    };

}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TensorForwardConfig, forward);
namespace WireCell::SPNG {

    // A tensor filter that sends a tensor through an ITensorForward
    struct TensorForward : public ContextBase,
                           public Logger,
                           public ITorchTensorFilter,
                           virtual public IConfigurable
    {
        TensorForward();
        virtual ~TensorForward() = default;

        /// ITorchTensorFilter
        virtual bool operator()(const input_pointer& in, output_pointer& out);

        // IConfigurable - see KernelConvolveConfig for configuration documentation.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

    private:
        TensorForwardConfig m_config;
        ITensorForward::pointer m_forward;
    };

}
