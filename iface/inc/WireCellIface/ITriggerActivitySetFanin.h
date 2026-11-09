#ifndef WIRECELL_ITRIGGERACTIVITYSETFANIN
#define WIRECELL_ITRIGGERACTIVITYSETFANIN

#include "WireCellIface/IFaninNode.h"
#include "WireCellIface/ITriggerActivitySet.h"

namespace WireCell {
    /** Merge TA sets from N input ports into one TA set.

        For example, the TAs of several anodes, to feed one TC finder.
     */
    class ITriggerActivitySetFanin : public IFaninNode<ITriggerActivitySet, ITriggerActivitySet, 0> {
       public:
        virtual ~ITriggerActivitySetFanin();

        virtual std::string signature() { return typeid(ITriggerActivitySetFanin).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> input_types() = 0;
        // as well as the already abstract:
        // virtual bool operator()(const input_vector& invec, output_pointer& out) = 0;
    };
}  // namespace WireCell

#endif
