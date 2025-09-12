#ifndef WIRECELL_SPNG_ITORCHSETTOFRAME
#define WIRECELL_SPNG_ITORCHSETToFRAME

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/IFrame.h"
#include "WireCellSpng/ITorchTensorSet.h"

namespace WireCell::SPNG {

    /** A frame to torch tensor set.
     */
    class ITorchSetToFrame : public IFunctionNode<ITorchTensorSet, IFrame> {
    public:
        virtual ~ITorchSetToFrame() {};

        virtual std::string signature() { return typeid(ITorchSetToFrame).name(); }

        // Subclass must implement:
        // virtual bool operator()(const input_pointer& in, output_pointer& out);
    };
}  // namespace WireCell::SPNG

#endif
