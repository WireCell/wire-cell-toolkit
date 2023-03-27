#ifndef WIRECELL_IFACE_IBLOBSETFANOUT
#define WIRECELL_IFACE_IBLOBSETFANOUT

#include "WireCellIface/IFanoutNode.h"
#include "WireCellIface/IBlobSet.h"

namespace WireCell {

    /** A blob set fan-out component takes a vector of blob sets (see
     * IBlobFanin, IBlobPipeline) and produces a fan of blob sets.
     */
    class IBlobSetFanout : public IFanoutNode<IBlobSet, IBlobSet, 0> {
       public:
        virtual ~IBlobSetFanout();

        virtual std::string signature() { return typeid(IBlobSetFanout).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> output_types() = 0;
        // as well as the already abstract:
        // virtual bool operator()(const input_pointer& in, output_vector& outv);
    };
}  // namespace WireCell

#endif
