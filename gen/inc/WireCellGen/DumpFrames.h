#ifndef WIRECELLGEN_DUMPFRAMES
#define WIRECELLGEN_DUMPFRAMES

#include "WireCellIface/IFrameSink.h"
#include "WireCellAux/Logger.h"

namespace WireCell::Gen {

    class DumpFrames : public Aux::Logger, public IFrameSink {
    public:
        DumpFrames();
        virtual ~DumpFrames();

        // IFrameSink 
        virtual bool operator()(const IFrame::pointer& frame);

    private:
        size_t m_count{0};

    };
}

#endif
