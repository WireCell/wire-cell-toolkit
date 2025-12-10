#pragma once

#include "WireCellSpng/Logger.h"
#include "WireCellSpng/ContextBase.h"
#include "WireCellSpng/HanaConfigurable.h"
#include "WireCellSpng/ITorchTensorFilter.h"

namespace WireCell::SPNG {

    struct TransformOpConfig {

        /** Name the transform operation.

            - permute :: Rearrange all dimensions.  Requires "dims" list spanning all dimensions. 

            - transpose :: Rearrange two dimensions.  Requires "dims" of size 2.

            - scale :: Multiply the given "scalar" value to tensor values.

            - offset :: Add the given "scalar" value to tensor values.

            - normalize :: Scale and offset data so that it is in the range [0,1].

            - medsub :: Subtract the median value calculated along dimension given by dims[0].

            - lowmedsub :: Same as medsub but will exclude any values from the
              median that are larger than given by "scalar".

            - threshold :: Convert to boolean with true for all pixels strictly
              above the threshold given by "scalar".

        */

        /// 
        std::string operation = "";

        /// A scalar value used in some operations
        float scalar = 0;

        /// Specify a the indices for a number of dimensions aka axes.
        std::vector<int> dims = {};

    };

    struct TransformConfig {
        std::vector<TransformOpConfig> operations;
    };
}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TransformOpConfig, operation, dims);
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TransformConfig, operations);

namespace WireCell::SPNG {

    using TransformOperation = std::function<torch::Tensor(const torch::Tensor&)>;
 
    // Apply transform operations to input tensor.
    struct Transform : public ITorchTensorFilter, 
                       public Logger,
                       public ContextBase,
                       public virtual IConfigurable
    {
        Transform();
        virtual ~Transform() = default;

        virtual bool operator()(const input_pointer& in, output_pointer& out);

        virtual void configure(const WireCell::Configuration& jconfig);
        virtual WireCell::Configuration default_configuration() const;


    private:
        std::vector<TransformOperation> m_ops;
        TransformConfig m_config;
    };

}

