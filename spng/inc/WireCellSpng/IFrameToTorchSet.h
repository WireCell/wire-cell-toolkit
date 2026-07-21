#ifndef WIRECELL_SPNG_IFRAMETOTORCHSET
#define WIRECELL_SPNG_IFRAMETOTORCHSET

#include "WireCellIface/IFunctionNode.h"
#include "WireCellIface/IFrame.h"
#include "WireCellSpng/ITorchTensorSet.h"

namespace WireCell::SPNG {

    /** A frame to torch tensor set.
     */
    class IFrameToTorchSet : public IFunctionNode<IFrame, ITorchTensorSet> {
    public:
        virtual ~IFrameToTorchSet() {};

        virtual std::string signature() { return typeid(IFrameToTorchSet).name(); }

        // Subclass must implement:
        // virtual bool operator()(const input_pointer& in, output_pointer& out);
    };
}  // namespace WireCell::SPNG

#endif
