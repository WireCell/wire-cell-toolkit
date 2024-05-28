#ifndef WIRECELLSIGPROC_FRAMESPLITTER
#define WIRECELLSIGPROC_FRAMESPLITTER

#include "WireCellIface/IFrameSplitter.h"
#include "WireCellAux/Logger.h"

namespace WireCell::SigProc {

    // fixme/note: this class is generic and nothing to do with
    // sigproc per se.  There are a few such frame related nodes
    // also in gen.  They should be moved into some neutral package.
    class FrameSplitter : public Aux::Logger, public WireCell::IFrameSplitter {
    public:
        FrameSplitter();
        virtual ~FrameSplitter();

        virtual bool operator()(const input_pointer& in, output_tuple_type& out);
    private:
        size_t m_count{0};
    };

}  // namespace WireCell::SigProc

#endif
