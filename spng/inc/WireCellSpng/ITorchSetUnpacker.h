#ifndef WIRECELL_SPNG_ITORCHSETUNPACKER
#define WIRECELL_SPNG_ITORCHSETUNPACKER

#include "WireCellIface/IFanoutNode.h"
#include "WireCellSpng/ITorchTensor.h"
#include "WireCellSpng/ITorchTensorSet.h"

#include <string>

namespace WireCell {
    /*
     * Unpack ITorchTensorSet into tensors
     */
    class ITorchSetUnpacker : public IFanoutNode<ITorchTensorSet, ITorchTensor, 0> {
       public:
        virtual ~ITorchSetUnpacker();

        virtual std::string signature() { return typeid(ITorchSetUnpacker).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> output_types() = 0;
        // and the already abstract:
        // virtual bool operator()(const input_pointer& in, output_vector& outv);
    };
}  // namespace WireCell

#endif