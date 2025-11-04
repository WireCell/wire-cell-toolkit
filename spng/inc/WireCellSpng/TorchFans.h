/** This provides simple fan in/out nodes for ITorchTensor and ITorchTensorSet.

    See the datamodel.org document for more information on the SPNG Tensor Data
    Model and the C++ class support.

 */
#ifndef WIRECELL_SPNG_FANOUTNODE
#define WIRECELL_SPNG_FANOUTNODE

#include "WireCellSpng/FanBase.h"
#include "WireCellSpng/Logger.h"

namespace WireCell::SPNG {

    // This performs a "union".
    struct FaninTensorSets : public FaninBase<ITorchTensorSet>
    {
        FaninTensorSets(const std::string& type_name="SPNGFaninTensorSets",
                        const std::string& group_name="spng");
        virtual ~FaninTensorSets() = default;

        // FaninBase interface.
        virtual void fanin_combine(const input_vector& inv, output_pointer& out);
        virtual void fanin_eos() { };
    };

    // This replicates input to all outputs.
    struct FanoutTensorSets : public FanoutBase<ITorchTensorSet>
    {

        FanoutTensorSets(const std::string& type_name="SPNGFoutTensorSets",
                         const std::string& group_name="spng");
        virtual ~FanoutTensorSets() = default;

        // FaninBase interface.
        virtual void fanout_separate(const input_pointer& in, output_vector& outv);
        virtual void fanout_eos() { };
    };

    // This batches inputs of common size.
    struct FaninTensors : public FaninBase<ITorchTensor>
    {
        FaninTensors(const std::string& type_name="SPNGFaninTensors",
                     const std::string& group_name="spng");
        virtual ~FaninTensors() = default;

        // FaninBase interface.
        virtual void fanin_combine(const input_vector& inv, output_pointer& out);
        virtual void fanin_eos() { };

    };

    // This replicates input to all outputs.
    struct FanoutTensors : public FanoutBase<ITorchTensor>
    {

        FanoutTensors(const std::string& type_name="SPNGFoutTensors",
                      const std::string& group_name="spng");
        virtual ~FanoutTensors() = default;

        // FaninBase interface.
        virtual void fanout_separate(const input_pointer& in, output_vector& outv);
        virtual void fanout_eos() { };

    };

}

#endif
