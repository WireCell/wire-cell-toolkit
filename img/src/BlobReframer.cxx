/*
  This starts by finding the set of all frames accessible via slices
  in a cluster that have at least one blob.



 */

#include "WireCellImg/BlobReframer.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellUtil/Exceptions.h"

#include <algorithm>
#include <cmath>

WIRECELL_FACTORY(BlobReframer, WireCell::Img::BlobReframer,
                 WireCell::INamed,
                 WireCell::IClusterFramer, WireCell::IConfigurable)

using namespace WireCell;

using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

Img::BlobReframer::BlobReframer(const std::string& frame_tag)
    : Aux::Logger("BlobReframer", "img")
    , m_nticks(0)
    , m_period(0)
    , m_frame_tag(frame_tag)
{
}

Img::BlobReframer::~BlobReframer() {}

void Img::BlobReframer::configure(const WireCell::Configuration& cfg)
{
    // m_nticks = get(cfg, "nticks", m_nticks);
    m_frame_tag = get(cfg, "frame_tag", m_frame_tag);
}

WireCell::Configuration Img::BlobReframer::default_configuration() const
{
    Configuration cfg;
    // cfg["nticks"] = m_nticks;
    cfg["frame_tag"] = m_frame_tag;
    return cfg;
}

// cluster goes in, frame goes out.
bool Img::BlobReframer::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        return true;  // EOS
    }

    cluster_indexed_graph_t grind(in->graph());

    // Find the frame and the slices.  The frame may have no traces: the
    // slices from MaskSlicer hold a TracelessFrame whose traces() is
    // nullptr, and frames rebuilt from cluster files have no traces.
    WireCell::IFrame::pointer iframe;
    std::vector<ISlice::pointer> islices;
    for (auto iblob : oftype<IBlob::pointer>(grind)) {
        auto islice = iblob->slice();
        if (!iframe) {
            iframe = islice->frame();
            if (!iframe) {
                raise<ValueError>("cluster %d has a blob whose slice has no frame", in->ident());
            }
            // typical values for protoDUNE
            // period: 0.5us, nticks: 6000
            m_period = iframe->tick();
        }
        islices.push_back(islice);
    }

    if (!iframe) {
        log->debug("cluster {} has no blobs, output empty frame", in->ident());
        auto sframe = std::make_shared<SimpleFrame>(in->ident(), 0.0, ITrace::vector{});
        sframe->tag_frame(m_frame_tag);
        out = sframe;
        return true;
    }

    // Output traces start at tbin 0 (the frame time) and run to the end of
    // the frame's traces, if it has any, and at least to the end of the
    // last slice.
    m_nticks = 0;
    auto ftraces = iframe->traces();
    if (ftraces) {
        for (const auto& trace : *ftraces) {
            m_nticks = std::max(m_nticks, (int) (trace->tbin() + trace->charge().size()));
        }
    }
    for (const auto& islice : islices) {
        m_nticks = std::max(m_nticks, (int) std::round((islice->start() + islice->span()) / m_period));
    }
    log->debug("nticks={} from {} frame traces and {} slices",
               m_nticks, ftraces ? ftraces->size() : 0, islices.size());

    std::unordered_map<int, ITrace::ChargeSequence> map_charges;
    for (const auto& islice : islices) {
        int itick = islice->start() / m_period;
        int nspan = islice->span() / m_period;

        for (auto const& act : islice->activity()) {
            auto ichan = act.first;
            if (map_charges.find(ichan->ident()) == map_charges.end()) {
                map_charges.insert(std::make_pair(ichan->ident(), ITrace::ChargeSequence(m_nticks, 0.0)));
            }
            float q = act.second;
            for (int ns = 0; ns < nspan; ns++) {
                map_charges.at(ichan->ident()).at(itick + ns) += q / (float) nspan;
            }
        }
    }

    ITrace::vector* itraces = new ITrace::vector;  // will become shared_ptr.
    for (auto const& it : map_charges) {
        int ident = it.first;
        auto charge = it.second;
        // Save the waveform densely, including zeros.
        SimpleTrace* trace = new SimpleTrace(ident, 0, charge);
        itraces->push_back(ITrace::pointer(trace));
    }

    SimpleFrame* sframe =
        new SimpleFrame(iframe->ident(), iframe->time(), ITrace::shared_vector(itraces), iframe->tick());
    sframe->tag_frame(m_frame_tag);
    out = IFrame::pointer(sframe);

    return true;
}
