#pragma once

#include "WireCellSpng/TensorFilter.h"
#include "WireCellSpng/ITensorForward.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/HanaJsonCPP.h"


namespace WireCell::SPNG {

    /// Configuration for TensorForwardConfig
    struct TensorForwardConfig {
        /// The type/name of an ITensorForward such as TensorForwardTS
        std::string forward=""; 

        /// Limit the maximum amount of batches to be fed forward.  If nbatch is
        /// positive nonzero, each (possibly sub) tensor that is fed forward
        /// will have a batch dimension given by batch_dimension that has a size
        /// no larger than nbatch.  If nbatch is set to 0, only one batch is
        /// forwarded at a time but with no batch dimension.  Set to nbatch to 1
        /// to forward one batch at a time and keep the batch dimension.
        int nbatch = 0;

        /// The expected number of NON-BATCH dimensions.  
        int ndims = 3;

        /// The batch dimension.  This should be the dimension that is batched
        /// if the input tensor is batched.  If the input tensor is not batched
        /// and if nbatch is nonzero, this dimension will be created in the
        /// tensor that is set forward.
        int batch_dimension = 0;

        // FIXME: add support for chunking?

    };

}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TensorForwardConfig, forward, nbatch, ndims, batch_dimension);
namespace WireCell::SPNG {

    /// A BARE torch tensor function that sends a BARE tensor through an ITensorForward.
    struct TensorForward : public TensorFilter,
                           virtual public IConfigurable
    {
        TensorForward();
        virtual ~TensorForward() = default;

        /// Tensorfilter
        virtual ITorchTensor::pointer filter_tensor(const ITorchTensor::pointer& in);

        // IConfigurable - see KernelConvolveConfig for configuration documentation.
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

    private:
        TensorForwardConfig m_config;
        ITensorForward::pointer m_forward;
    };

}
