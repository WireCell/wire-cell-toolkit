#ifndef WIRECELL_SPNG_ITORCHTENSORSETFANOUT
#define WIRECELL_SPNG_ITORCHTENSORSETFANOUT

#include "WireCellIface/IFanoutNode.h"
#include "WireCellSpng/ITorchTensorSet.h"

namespace WireCell {
namespace SPNG {
    /** A fan-out component takes 1 input TorchTensorSet and produces N
     * TorchTensorSets.

     */
    class ITorchTensorSetFanout : public IFanoutNode<ITorchTensorSet, ITorchTensorSet, 0> {
       public:
        virtual ~ITorchTensorSetFanout() {};

        virtual std::string signature() { return typeid(ITorchTensorSetFanout).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> output_types() = 0;
        // and the already abstract:
        // virtual bool operator()(const input_pointer& in, output_vector& outv);
    };
}  // namespace SPNG
}  // namespace WireCell

#endif
