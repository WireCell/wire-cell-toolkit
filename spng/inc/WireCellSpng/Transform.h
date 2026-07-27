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

            - squeeze :: Apply squeeze() to the tensor on each dimension in dims.

            - unsqueeze :: Apply unsqueeze() to the tensor on each dimension in dims.

            - noop :: Explicit no-operation on the tensor.  If you want to rely
              on Transform to change metadata, set this operation to avoid a warning.

            - to :: casts as the given dtype

            - view :: returns input_tensr.view(dims), where "dims" is a vector of dimension lengths
        */

        /// 
        std::string operation = "";

        /// A scalar value used in some operations
        double scalar = 0;

        /// Specify which dimension(s) to operate on.
        std::vector<int> dims = {};

        std::string dtype = "";
    };

    struct TransformConfig {
        std::vector<TransformOpConfig> operations;
    };

    torch::Dtype resolve_dtype(std::string dtype) {
        if (dtype == "int64") {
            return torch::kInt64;
        } else if (dtype == "int32") {
            return torch::kInt32;
        } else if (dtype == "int16") {
            return torch::kInt16;
        } else if (dtype == "int8") {
            return torch::kInt8;
        } else if (dtype == "float64") {
            return torch::kFloat64;
        } else if (dtype == "float32") {
            return torch::kFloat32;
        } else if (dtype == "float16") {
            return torch::kFloat16;
        } else if (dtype == "bool") {
            return torch::kBool;
        } else {
            raise<ValueError>("unknown dtype: %s", dtype);
        }
        //Default value that will never be reached 
        //In order to quiet compiler warnings
        return torch::kUInt8; 
    };
}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::TransformOpConfig, operation, scalar, dims, dtype);
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

