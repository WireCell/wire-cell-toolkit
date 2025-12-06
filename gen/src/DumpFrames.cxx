#include "WireCellGen/DumpFrames.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(DumpFrames, WireCell::Gen::DumpFrames, WireCell::IFrameSink)

using namespace WireCell;

Gen::DumpFrames::DumpFrames()
  : log(Log::logger("glue"))
{
}

Gen::DumpFrames::~DumpFrames() {}

bool Gen::DumpFrames::operator()(const IFrame::pointer& frame)
{
    if (!frame) {
        log->debug("frame sink sees EOS");
        return true;
    }
    auto traces = frame->traces();
    int ntraces = 0;
    if (traces) ntraces = traces->size();

    std::stringstream ss;
    ss << "sink frame: #" << frame->ident() << " @" << frame->time() / units::ms << " with " << ntraces << " traces";
    {
        std::string comma = "";
        ss << ", frame tags:[";
        for (auto ftag : frame->frame_tags()) {
            ss << comma << ftag;
            comma = ", ";
        }
        ss << "]";
    }
    {
        std::string comma = "";
        ss << ", trace tags:[";
        for (auto ftag : frame->trace_tags()) {
            ss << comma << ftag;
            comma = ", ";
        }
        ss << "]";
    }
    log->debug(ss.str());
    return true;
}
