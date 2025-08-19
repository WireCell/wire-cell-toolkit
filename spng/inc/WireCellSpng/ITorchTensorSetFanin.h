#ifndef WIRECELL_SPNG_ITORCHTENSORSETFANIN
#define WIRECELL_SPNG_ITORCHTENSORSETFANIN

#include "WireCellIface/IFaninNode.h"
#include "WireCellSpng/ITorchTensorSet.h"

namespace WireCell {
namespace SPNG {
    /** A frame fan-out component takes N input TorchTensorSets and produces one
     * TorchTensorSet.

     */
    class ITorchTensorSetFanin : public IFaninNode<ITorchTensorSet, ITorchTensorSet, 0> {
       public:
        virtual ~ITorchTensorSetFanin() {};

        virtual std::string signature() { return typeid(ITorchTensorSetFanin).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> input_types() = 0;
    };
}  // namespace SPNG
}  // namespace WireCell

#endif