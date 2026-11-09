#ifndef WIRECELL_ITRIGGERCANDIDATEFINDER
#define WIRECELL_ITRIGGERCANDIDATEFINDER

#include "WireCellIface/IQueuedoutNode.h"
#include "WireCellIface/ITriggerActivitySet.h"
#include "WireCellIface/ITriggerCandidateSet.h"

namespace WireCell {

    /** Consume a set of TAs and produce a set of the TCs found.

        A finder emits one output set per input set, possibly empty, with the
        input's coverage.  It may hold state across sets.  On end of stream
        it may emit a last set with what it held, then end of stream.  See
        ITriggerPrimitiveSet.h for the stream contract.
     */
    class ITriggerCandidateFinder : public IQueuedoutNode<ITriggerActivitySet, ITriggerCandidateSet> {
       public:
        typedef std::shared_ptr<ITriggerCandidateFinder> pointer;

        virtual ~ITriggerCandidateFinder();

        virtual std::string signature() { return typeid(ITriggerCandidateFinder).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_queue& outq);
    };

}  // namespace WireCell

#endif
