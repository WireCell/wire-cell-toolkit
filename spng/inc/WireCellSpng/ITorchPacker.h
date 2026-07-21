#ifndef WIRECELL_SPNG_ITORCHPACKER
#define WIRECELL_SPNG_ITORCHPACKER

#include "WireCellIface/IFaninNode.h"
#include "WireCellSpng/ITorchTensor.h"
#include "WireCellSpng/ITorchTensorSet.h"

#include <string>

namespace WireCell {
    /*
     * Pack ITorchTensors into ITorchTensorSet
     */
    class ITorchPacker : public IFaninNode<ITorchTensor, ITorchTensorSet, 0> {
       public:
        virtual ~ITorchPacker() {};

        virtual std::string signature() { return typeid(ITorchPacker).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> input_types() = 0;
        // and the already abstract:
        // virtual bool operator()(const input_pointer& in, output_vector& outv);
    };
}  // namespace WireCell

#endif