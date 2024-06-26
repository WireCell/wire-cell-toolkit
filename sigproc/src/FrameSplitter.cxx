#include "WireCellSigProc/FrameSplitter.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/FrameTools.h"


WIRECELL_FACTORY(FrameSplitter,
                 WireCell::SigProc::FrameSplitter,
                 WireCell::INamed, WireCell::IFrameSplitter)

using namespace WireCell::SigProc;

FrameSplitter::FrameSplitter()
  : Aux::Logger("FrameSplitter", "sigproc")
{
}

FrameSplitter::~FrameSplitter() {}

bool FrameSplitter::operator()(const input_pointer& in, output_tuple_type& out)
{
    get<0>(out) = in;
    get<1>(out) = in;

    if (!in) {
        log->debug("EOS at call={}", m_count++);
        return true;
    }

    log->debug("call={} {}", m_count++, WireCell::Aux::taginfo(in));
    return true;
}
