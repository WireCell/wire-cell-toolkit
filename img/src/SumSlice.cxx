#include "WireCellImg/SumSlice.h"
#include "WireCellImg/ImgData.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"

#include "WireCellAux/FrameTools.h"

WIRECELL_FACTORY(SumSlicer, WireCell::Img::SumSlicer,
                 WireCell::INamed,
                 WireCell::IFrameSlicer, WireCell::IConfigurable)

WIRECELL_FACTORY(SumSlices, WireCell::Img::SumSlices,
                 WireCell::INamed,
                 WireCell::IFrameSlices, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

Img::SumSliceBase::SumSliceBase()
    : Aux::Logger("SumSlice", "img")
    , m_tick_span(4)
    , m_tag("")
{
}
Img::SumSliceBase::~SumSliceBase() {}
Img::SumSlicer::~SumSlicer() {}
Img::SumSlices::~SumSlices() {}

WireCell::Configuration Img::SumSliceBase::default_configuration() const
{
    Configuration cfg;

    // Name of an IAnodePlane from which we can resolve channel ident to IChannel
    cfg["anode"] = "";

    // Number of ticks over which each output slice should sum
    cfg["tick_span"] = m_tick_span;

    // If given, use tagged traces.  Otherwise use all traces.
    cfg["tag"] = m_tag;

    // If true, emit an EOS at end of one input frame's slices
    cfg["slice_eos"] = m_slice_eos;

    return cfg;
}

void Img::SumSliceBase::configure(const WireCell::Configuration& cfg)
{
    m_anode = Factory::find_tn<IAnodePlane>(cfg["anode"].asString());  // throws
    m_tick_span = get(cfg, "tick_span", m_tick_span);
    m_tag = get<std::string>(cfg, "tag", m_tag);
    m_slice_eos = get<bool>(cfg, "slice_eos", m_slice_eos);
    m_pad_empty = get<bool>(cfg, "pad_empty", m_pad_empty);
}

ISlice::pointer Img::SumSliceBase::make_empty(const IFrame::pointer& inframe)
{
    const double tick = inframe->tick();
    const double span = tick * m_tick_span;
    return std::make_shared<Img::Data::Slice>(inframe, 0, 0, span);
}

void Img::SumSliceBase::slice(const IFrame::pointer& inframe, slice_map_t& svcmap)
{
    const double tick = inframe->tick();
    const double span = tick * m_tick_span;

    for (auto trace : Aux::tagged_traces(inframe, m_tag)) {
        const int tbin = trace->tbin();
        const int chid = trace->channel();
        IChannel::pointer ich = m_anode->channel(chid);
        const auto& charge = trace->charge();
        const size_t nq = charge.size();
        for (size_t qind = 0; qind != nq; ++qind) {
            const auto q = charge[qind];
            if (q == 0.0) {
                continue;
            }
            size_t slicebin = (tbin + qind) / m_tick_span;
            auto s = svcmap[slicebin];
            if (!s) {
                // ISlice::start() is relative to the frame's time.
                const double start = slicebin * span;
                svcmap[slicebin] = s = new Img::Data::Slice(inframe, slicebin, start, span);
            }
            s->sum(ich, q);
        }
    }
}

bool Img::SumSlicer::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        return true;  // eos
    }

    // Slices will be sparse in general.  Index by a "slice bin" number
    slice_map_t svcmap;
    slice(in, svcmap);

    // intern
    ISlice::vector islices;
    for (auto sit : svcmap) {
        auto s = sit.second;
        islices.push_back(ISlice::pointer(s));
    }
    if (m_slice_eos) {
        islices.push_back(nullptr);
    }
    out = make_shared<Img::Data::SliceFrame>(islices, in->ident(), in->time());

    return true;
}

bool Img::SumSlices::operator()(const input_pointer& inframe, output_queue& slices)
{
    if (!inframe) {
        log->debug("EOS");
        slices.push_back(nullptr);
        return true;  // eos
    }

    // Slices will be sparse in general.  Index by a "slice bin" number
    slice_map_t svcmap;
    slice(inframe, svcmap);

    // intern
    for (auto sit : svcmap) {
        auto s = sit.second;

        /// debug
        // double qtot = 0;
        // for (const auto& a : s->activity()) {
        //     qtot += a.second;
        // }

        slices.push_back(ISlice::pointer(s));
    }

    if (slices.empty()) {
        log->warn("frame={} time={} ms, made no slices",
                  inframe->ident(), inframe->time()/units::ms);
        if (m_pad_empty) {
            slices.push_back(make_empty(inframe));
        }
    }
    else {
        log->debug("frame={}, made {} slices in [{},{}] from {}",
                   inframe->ident(), slices.size(),
                   slices.front()->ident(), slices.back()->ident(),
                   svcmap.size());
    }
    if (m_slice_eos) {
        slices.push_back(nullptr);
    }
    return true;
}
