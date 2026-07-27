#ifndef WIRECELL_SPNG_ITORCHTENSORSETTOFRAMEFANIN
#define WIRECELL_SPNG_ITORCHTENSORSETTOFRAMEFANIN

#include "WireCellIface/IFaninNode.h"
#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellIface/IFrame.h"

namespace WireCell {
namespace SPNG {
    /** A frame fan-out component takes N input TorchTensorSets and produces one
     * IFrame.

     */
    class ITorchTensorSetToFrameFanin : public IFaninNode<ITorchTensorSet, IFrame, 0> {
       public:
        virtual ~ITorchTensorSetToFrameFanin() {};

        virtual std::string signature() { return typeid(ITorchTensorSetToFrameFanin).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> input_types() = 0;
    };
}  // namespace SPNG
}  // namespace WireCell

#endif