#ifndef WIRECELL_ITRIGGERPRIMITIVEFINDER
#define WIRECELL_ITRIGGERPRIMITIVEFINDER

#include "WireCellIface/IQueuedoutNode.h"
#include "WireCellIface/IFrame.h"
#include "WireCellIface/ITriggerPrimitiveSet.h"

namespace WireCell {

    /** Consume a frame of ADC waveforms and produce a set of the TPs found.

        A finder emits one output set per input set, possibly empty, with the
        input's coverage.  It may hold state across sets.  On end of stream
        it may emit a last set with what it held, then end of stream.  See
        ITriggerPrimitiveSet.h for the stream contract.
     */
    class ITriggerPrimitiveFinder : public IQueuedoutNode<IFrame, ITriggerPrimitiveSet> {
       public:
        typedef std::shared_ptr<ITriggerPrimitiveFinder> pointer;

        virtual ~ITriggerPrimitiveFinder();

        virtual std::string signature() { return typeid(ITriggerPrimitiveFinder).name(); }

        // supply:
        // virtual bool operator()(const input_pointer& in, output_queue& outq);
    };

}  // namespace WireCell

#endif
