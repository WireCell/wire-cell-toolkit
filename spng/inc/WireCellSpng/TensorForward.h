#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/ITorchTensorFilter.h"
#include "WireCellSpng/ITensorForward.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/HanaJsonCPP.h"


namespace WireCell::SPNG {

    /// Configuration for TensorForwardConfig
    struct TensorForwardConfig {
        /// The type/name of an ITensorForward such as TensorForwardTS
        std::string forward=""; 

        /// Limit the maximum amount of batches fed forward.  If 0, the full
        /// input tensor is sent forward.  Otherwise, it is chunked to be this size.
        bool nbatch = 0;

        /// The batch dimension.
        int batch_dimension = 0;


        // FIXME: add support for chunking, orthogonal to sequencing
        // bool chunk = false;
        // std::vector<int> shape.  // need to cope with batch dimension

    };

}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TensorForwardConfig, forward);
namespace WireCell::SPNG {

    /// A BARE torch tensor function that sends a BARE tensor through an ITensorForward.
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
