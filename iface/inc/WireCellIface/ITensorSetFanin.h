#ifndef WIRECELL_IFACE_ITENSORSETFANIN
#define WIRECELL_IFACE_ITENSORSETFANIN

#include "WireCellIface/IFaninNode.h"
#include "WireCellIface/ITensorSet.h"

namespace WireCell {
    /** 
     * Fanin ITensorSets, output ITensorSet
     */
    class ITensorSetFanin : public IFaninNode<ITensorSet, ITensorSet, 0> {
      public:
        virtual ~ITensorSetFanin();

        virtual std::string signature() { return typeid(ITensorSetFanin).name(); }

        // Subclass must implement:
        virtual std::vector<std::string> input_types() = 0;
        // and the already abstract:
        // virtual bool operator()(const input_vector& invec, output_pointer& out) = 0;
    };
}  // namespace WireCell

#endif
