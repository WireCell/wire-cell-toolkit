#pragma once

// Opposite of Expand.h

#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/Torch.h"



namespace WireCell::SPNG {

    struct ReduceConfig {

        /** The Reduction operation names a function to apply element-wise
            across the input tensors.  Function that operate along a specific,
            given dimension "dim" are:

            - cat :: Concatenates tensors along an existing dimension dim.
              Other dimensions must have same size.

            - stack :: Stacks tensors along a NEW dimension (dim). All inputs
              must have the exact same size.

            Operations that do not accept a dimension are:

            - hstack :: Horizontal Stack. Equivalent to cat along dimension 0
              for 1D tensors, and dimension 1 for higher-dimensional tensors.

            - vstack :: Vertical Stack. Equivalent to stack along dimension 0
              for 1D tensors, and cat along dimension 0 for 2D+ tensors.

            - dstack :: Depth Stack. Stacks tensors along the third dimension
              (index 2) after adding a dimension of size 1 if necessary.

            - min :: Minimum elements across the tensors.

            - max :: Maximum elements across the tensors.

            - mean :: Mean of elements across the tensors.  The alias "avg" is accepted.

            - sum :: Summation of elements across the tensors.  The alias "add" is accepted.

            - mul :: Multiplication of elements across the tensors.  The alias "prod" is accepted.

        */

        /// vstack, hstack, max, min, mean, sum, prod.
        std::string operation = "";

        /// Specify a particular tensor dimension for operations that accept it.
        int dim = 0;
    };
}
BOOST_HANA_ADAPT_STRUCT(WireCell::SPNG::ReduceConfig, operation, dim);

namespace WireCell::SPNG {

    using ReduceOperation = std::function<torch::Tensor(const std::vector<torch::Tensor>&)>;
 
    // Apply a reduction operation to all tensors in the fan-in to produce one tensor out.
    struct Reduce : FaninBase<ITorchTensor>, virtual public IConfigurable
    {
        Reduce();
        virtual ~Reduce() = default;

        // FaninBase API.
        virtual void fanin_combine(const input_vector& inv, output_pointer& out);

        virtual void configure(const WireCell::Configuration& jconfig);
        virtual WireCell::Configuration default_configuration() const;


    private:
        ReduceOperation m_op;
        ReduceConfig m_config;
    };

}

