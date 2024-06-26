// This component produces an output frame which consists of traces
// from the two input frames that have been merged together based on
// tags and a rule.

#ifndef WIRECELL_SIGPROC_FRAMEMERGER
#define WIRECELL_SIGPROC_FRAMEMERGER

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameJoiner.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellAux/Logger.h"

namespace WireCell::SigProc {
    class FrameMerger : public Aux::Logger, public IFrameJoiner, public IConfigurable {
    public:
        FrameMerger();
        virtual ~FrameMerger();

        // IJoinNode
        virtual bool operator()(const input_tuple_type& intup, output_pointer& out);

        // IConfigurable
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

    private:
        Configuration m_cfg;
        size_t m_count{0};
    };

}  // namespace WireCell::SigProc

#endif
