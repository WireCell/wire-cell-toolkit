#ifndef WIRECELL_SPNG_IFRAMETOTORCHSETFANOUT
#define WIRECELL_SPNG_IFRAMETOTORCHSETFANOUT

#include "WireCellIface/IFanoutNode.h"
#include "WireCellIface/IFrame.h"
#include "WireCellSpng/ITorchTensorSet.h"

namespace WireCell {
namespace SPNG {
    /** A frame fan-out component takes 1 input frame and produces one
     * TorchTensorSet.

     */
    class IFrameToTorchSetFanout : public IFanoutNode<IFrame, ITorchTensorSet, 0> {
       public:
        virtual ~IFrameToTorchSetFanout() {};

        virtual std::string signature() { return typeid(IFrameToTorchSetFanout).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> output_types() = 0;
        // and the already abstract:
        // virtual bool operator()(const input_pointer& in, output_vector& outv);
    };
}  // namespace SPNG
}  // namespace WireCell

#endif