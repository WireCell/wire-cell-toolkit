#include "WireCellImg/BlobSetReframer.h"

#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

#include "WireCellIface/IAnodePlane.h"

#include "WireCellUtil/String.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(BlobSetReframer, WireCell::Img::BlobSetReframer,
                 WireCell::INamed,
                 WireCell::IBlobSetFramer, WireCell::IConfigurable)

using namespace WireCell;

using WireCell::String::startswith;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

Img::BlobSetReframer::BlobSetReframer()
    : Aux::Logger("BlobSetReframer", "img")
{
}

Img::BlobSetReframer::~BlobSetReframer() {}

struct FromBlob {
    double tick;
    bool use_error{false};

    ITrace::vector operator()(const IBlob::pointer& iblob, double time0) {
        ITrace::vector traces;  // our return value

        // Get value to use as frame "charge"
        const double blob_charge = use_error ? iblob->uncertainty() : iblob->value();

        // Unpack info about the blob.
        auto islice = iblob->slice();
        const double start = islice->start();
        const double span = islice->span();
        const auto& blob = iblob->shape();
        const auto iface = iblob->face();
        const int ticks_in_slice = round(span/tick);
        const int tbin = round((start - time0)/tick);

        // check each layer of the blob
        for (const auto& strip : blob.strips()) {

            // Skip the overall active anode boundary layers
            if (strip.layer < 2) {
                continue;
            }

            // construct our constant "charge density" waveform over
            // the slice span
            const auto [wip1, wip2] = strip.bounds;
            const int channels_in_view = wip2-wip1;
            const double q = blob_charge / (ticks_in_slice * channels_in_view);
            const std::vector<float> charge(ticks_in_slice, q);

            // Make one trace for each wire spanning the blob bounds.
            auto iwplane = iface->plane(strip.layer-2);
            const auto& wires = iwplane->wires();
            for (auto iwip = wip1; iwip < wip2; ++iwip) {
                auto iwire = wires[iwip];
                int chid = iwire->channel();

                auto strace = std::make_shared<SimpleTrace>(chid, tbin, charge);
                traces.push_back(strace);
            }
        }
        return traces;
    }
};

struct FromActivity {
    double tick;
    bool use_error{false};

    const ISlice::value_t bogus_value{0, 1e-10};

    // A cache of chid->IChannel lookup.  Will replenish with iwplane if cache miss.
    std::unordered_map<int, IChannel::pointer> chid2ich;
    IChannel::pointer ichan(int chid, const IAnodeFace::pointer& iface)
    {
        auto it = chid2ich.find(chid);
        if (it == chid2ich.end()) { // cache miss
            for (auto iwplane : iface->planes()) {
                for (const auto& ich : iwplane->channels()) {
                    chid2ich[ich->ident()] = ich;
                }
            }
            it = chid2ich.find(chid); // 2nd chance
        }
        if (it == chid2ich.end()) {
            raise<ValueError>("no chid %d in face %d", chid, iface->ident());
        }
        return it->second;        
    }

    ITrace::vector operator()(const IBlob::pointer& iblob, double time0)
    {
        ITrace::vector traces;  // our return value

        // Unpack info about the blob.
        auto islice = iblob->slice();
        const double start = islice->start();
        const double span = islice->span();
        const auto& blob = iblob->shape();
        const auto iface = iblob->face();
        const int ticks_in_slice = round(span/tick);
        const int tbin = round((start - time0)/tick);

        // we do not have easy way to go from wire->chid->ich so go
        // from ich->chid and wire->chid.
        std::unordered_map<int, ISlice::value_t> activity;
        for (const auto& [ich, act] : islice->activity()) {
            activity[ich->ident()] = act;
        }

        // double qtot=0;

        // check each layer of the blob
        for (const auto& strip : blob.strips()) {

            // Skip the overall active anode boundary layers
            if (strip.layer < 2) {
                continue;
            }

            // Walk wires, find channel, get activity, spread if over
            // the slice span for the channel.
            const auto [wip1, wip2] = strip.bounds;
            auto iwplane = iface->plane(strip.layer-2);
            const auto& wires = iwplane->wires();
            for (auto iwip = wip1; iwip < wip2; ++iwip) {
                auto iwire = wires[iwip];
                const int chid = iwire->channel();
                const auto val = activity.at(chid);
                const double act_charge = use_error ? val.uncertainty() : val.value();
                const double q = act_charge / ticks_in_slice;
                const std::vector<float> charge(ticks_in_slice, q);
                auto strace = std::make_shared<SimpleTrace>(chid, tbin, charge);
                traces.push_back(strace);

                // qtot += act_charge;
            }
        }
        return traces;
    }
};

void Img::BlobSetReframer::configure(const WireCell::Configuration& cfg)
{
    m_frame_tag = get(cfg, "frame_tag", m_frame_tag);
    m_tick = get(cfg, "tick", m_tick);

    bool use_error = false;
    {
        const std::string meas = get<std::string>(cfg, "measure", "value");

        if (meas.empty() or startswith(meas, "val")) {
            use_error = false;
        }
        else if (startswith(meas, "err") or startswith(meas, "unc")) {
            use_error = true;
        }
        else {
            raise<ValueError>("unknown type of measure: %s", meas);
        }
    }    

    const std::string source = get<std::string>(cfg, "source", "blob");
    if (source == "blob") {
        m_trace_maker = FromBlob{m_tick, use_error};
        return;
    }
    if (startswith(source, "act")) {
        m_trace_maker = FromActivity{m_tick, use_error};
        return;
    }
    raise<ValueError>("unknown type of source: %s", source);
}

WireCell::Configuration Img::BlobSetReframer::default_configuration() const
{
    Configuration cfg;
    cfg["frame_tag"] = m_frame_tag;
    cfg["tick"] = m_tick;
    return cfg;
}


bool Img::BlobSetReframer::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count++);
        return true;  // EOS
    }

    ITrace::vector traces;

    std::unordered_map<int, IFrame::trace_list_t> anode_traces;

    const int ident = in->ident();
    const auto blobs = in->blobs();

    // Find earliest time of all blobs
    double time0 = 0;
    for (size_t ind=0; ind<blobs.size(); ++ind) {
        if (ind) time0 = std::min(time0, blobs[ind]->slice()->start());
        else time0 = blobs[ind]->slice()->start();
    }

    for (auto iblob : blobs) {
        
        auto face = iblob->face();
        const int aid = face->anode();

        auto blob_traces = m_trace_maker(iblob, time0);
        // log->debug("call={} aid={} bset={} iblob={} ntraces={}",
        //            m_count, aid, ident, iblob->ident(), blob_traces.size());
        for (auto trace : blob_traces) {
            anode_traces[aid].push_back(traces.size());
            traces.push_back(trace);
        }
    }

    auto sframe = std::make_shared<SimpleFrame>(ident, time0, traces, m_tick);
    if (! m_frame_tag.empty()) {
        sframe->tag_frame(m_frame_tag);
    }
    for (const auto& [aid, indices] : anode_traces) {
        std::string tag = "anode" + std::to_string(aid);
        sframe->tag_traces(tag, indices);
    }

    out = sframe;
    log->debug("call={} time0={} output: {}", m_count++, time0, Aux::taginfo(out));

    return true;
}
