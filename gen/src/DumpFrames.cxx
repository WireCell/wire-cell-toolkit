#include "WireCellGen/DumpFrames.h"

#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/NamedFactory.h"


WIRECELL_FACTORY(DumpFrames, WireCell::Gen::DumpFrames, WireCell::IFrameSink, WireCell::INamed)

using namespace WireCell;

Gen::DumpFrames::DumpFrames()
    : Aux::Logger("DumpFrames", "gen")
{
}

Gen::DumpFrames::~DumpFrames() {}

bool Gen::DumpFrames::operator()(const IFrame::pointer& frame)
{
    if (!frame) {
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    log->debug("frame taginfo at call={}: {}", m_count, Aux::taginfo(frame));
    Aux::dump_frame(frame, log);
    ++m_count;
    return true;
}
