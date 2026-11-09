#ifndef WIRECELL_ITRIGGERACTIVITYTRUTH
#define WIRECELL_ITRIGGERACTIVITYTRUTH

#include "WireCellIface/IDepoSet.h"
#include "WireCellIface/IJoinNode.h"
#include "WireCellIface/ITensorSet.h"
#include "WireCellIface/ITriggerActivitySet.h"

namespace WireCell {

    /** Join the depos of a frame with the TAs found in it to produce truth
        tables about the TAs, as a tensor set.
     */
    class ITriggerActivityTruth : public IJoinNode<std::tuple<IDepoSet, ITriggerActivitySet>, ITensorSet> {
       public:
        virtual ~ITriggerActivityTruth();

        virtual std::string signature() { return typeid(ITriggerActivityTruth).name(); }

        // supply:
        // virtual bool operator()(const input_tuple_type& intup, output_pointer& out);
    };

}  // namespace WireCell

#endif
