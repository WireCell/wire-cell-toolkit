#pragma once

// opposite of Reduce.h

#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/Torch.h"



namespace WireCell::SPNG {

    struct ExpandConfig {

        /** The Reduction operation names a function to apply element-wise
            across the input tensors.  Function that operate along a specific,
            given dimension "dim" are:

            - unbind :: Remove the dim and output element of the dimension as
              its own tensor with one less dimension.  Opposite of "cat".

            Operations that do not accept a dimension are:

        */

        /// unbind
        std::string operation = "";

        /// Specify a particular tensor dimension for operations that accept it.
        int dim = 0;
    };
}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::ExpandConfig, operation, dim);

namespace WireCell::SPNG {

    using ExpandOperation = std::function<std::vector<torch::Tensor>(const torch::Tensor&)>;
 
    // Apply a reduction operation to all tensors in the fan-in to produce one tensor out.
    struct Expand : FanoutBase<ITorchTensor>, virtual public IConfigurable
    {
        Expand();
        virtual ~Expand() = default;

        // FanoutBase API.
        virtual void fanout_separate(const input_pointer& in, output_vector& outv);

        virtual void configure(const WireCell::Configuration& jconfig);
        virtual WireCell::Configuration default_configuration() const;


    private:
        ExpandOperation m_op;
        ExpandConfig m_config;
    };

}

