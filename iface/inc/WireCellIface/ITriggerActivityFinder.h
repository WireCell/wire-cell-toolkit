#ifndef WIRECELL_ITRIGGERACTIVITYFINDER
#define WIRECELL_ITRIGGERACTIVITYFINDER

#include "WireCellIface/IQueuedoutNode.h"
#include "WireCellIface/ITriggerPrimitiveSet.h"
#include "WireCellIface/ITriggerActivitySet.h"

namespace WireCell {

    /** Consume a set of TPs and produce a set of the TAs found.

        A finder emits one output set per input set, possibly empty, with the
        input's coverage.  It may hold state across sets.  On end of stream
        it may emit a last set with what it held, then end of stream.  See
        ITriggerPrimitiveSet.h for the stream contract.
     */
    class ITriggerActivityFinder : public IQueuedoutNode<ITriggerPrimitiveSet, ITriggerActivitySet> {
       public:
        typedef std::shared_ptr<ITriggerActivityFinder> pointer;

        virtual ~ITriggerActivityFinder();

        virtual std::string signature() { return typeid(ITriggerActivityFinder).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_queue& outq);
    };

}  // namespace WireCell

#endif
