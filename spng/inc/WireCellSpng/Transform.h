#pragma once

#include "WireCellSpng/TensorFilter.h"
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

            - slice :: Apply a slice() along the dimension given by dims[0] from
              start index in dims[1] to end index in dims[2].
        */

        /// 
        std::string operation = "";

        /// A scalar value used in some operations
        float scalar = 0;

        /// Specify which dimension(s) to operate on.
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
    //
    // This node is basically a bag of torch functions to allow a large variety
    // of transforms to be programmed in the configuration without having to
    // make dedicated DFP nodes.
    struct Transform : public TensorFilter, 
                       public virtual IConfigurable
    {
        Transform();
        virtual ~Transform() = default;

        /// TensorFilter
        virtual ITorchTensor::pointer filter_tensor(const ITorchTensor::pointer& in);

        virtual void configure(const WireCell::Configuration& jconfig);
        virtual WireCell::Configuration default_configuration() const;


    private:
        std::vector<TransformOperation> m_ops;
        TransformConfig m_config;
    };

}

