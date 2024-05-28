#ifndef WIRECELL_IFACE_IFRAMEFANIN
#define WIRECELL_IFACE_IFRAMEFANIN

#include "WireCellIface/IFaninNode.h"
#include "WireCellIface/IFrame.h"

namespace WireCell {
    /** A frame fan-in component takes N frames on input ports and
     * produces a single output frame.  The merge policy and number of
     * inputs is left to the implementation.
     */
    class IFrameFanin : public IFaninNode<IFrame, IFrame, 0> {
      public:
        virtual ~IFrameFanin();

        typedef std::shared_ptr<IFrameFanin> pointer;

        virtual std::string signature() { return typeid(IFrameFanin).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> input_types() = 0;
        // and the already abstract:
        // virtual bool operator()(const input_vector& invec, output_pointer& out) = 0;
    };
}  // namespace WireCell

#endif
