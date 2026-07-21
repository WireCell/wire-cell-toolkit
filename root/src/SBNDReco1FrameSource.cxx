#include "WireCellRoot/SBNDReco1FrameSource.h"
#include "SBNDReco1Reader.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include "TTimeStamp.h"

WIRECELL_FACTORY(SBNDReco1FrameSource, WireCell::Root::SBNDReco1FrameSource,
                 WireCell::INamed,
                 WireCell::IFrameSource, WireCell::IConfigurable)

using namespace WireCell;

Root::SBNDReco1FrameSource::SBNDReco1FrameSource()
    : Aux::Logger("SBNDReco1FrameSource", "root")
{
}

Root::SBNDReco1FrameSource::~SBNDReco1FrameSource() {}

WireCell::Configuration Root::SBNDReco1FrameSource::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_filename;
    cfg["entry"] = m_entry;
    cfg["run"] = m_run;
    cfg["subrun"] = m_subrun;
    cfg["event"] = m_event;
    cfg["wire_product"] = m_wire_product;
    cfg["badmask_product"] = m_badmask_product;
    cfg["summary_product"] = m_summary_product;
    cfg["frame_tag"] = m_frame_tag;
    cfg["frame_scale"] = m_frame_scale;
    cfg["summary_scale"] = m_summary_scale;
    cfg["nticks"] = m_nticks;
    cfg["tick"] = m_tick;
    return cfg;
}

void Root::SBNDReco1FrameSource::configure(const WireCell::Configuration& cfg)
{
    m_filename = get(cfg, "filename", m_filename);
    m_entry = get(cfg, "entry", m_entry);
    m_run = get(cfg, "run", m_run);
    m_subrun = get(cfg, "subrun", m_subrun);
    m_event = get(cfg, "event", m_event);
    m_wire_product = get(cfg, "wire_product", m_wire_product);
    m_badmask_product = get(cfg, "badmask_product", m_badmask_product);
    m_summary_product = get(cfg, "summary_product", m_summary_product);
    m_frame_tag = get(cfg, "frame_tag", m_frame_tag);
    m_frame_scale = get(cfg, "frame_scale", m_frame_scale);
    m_summary_scale = get(cfg, "summary_scale", m_summary_scale);
    m_nticks = get(cfg, "nticks", m_nticks);
    m_tick = get(cfg, "tick", m_tick);
    if (m_filename.empty()) {
        THROW(ValueError() << errmsg{"SBNDReco1FrameSource: must supply input \"filename\""});
    }
    // Selection: entry >= 0, or run >= 0 (+event), or neither = stream all.
    if (m_entry < 0 and m_run >= 0 and m_event < 0) {
        THROW(ValueError() << errmsg{"SBNDReco1FrameSource: \"run\" selection also needs \"event\""});
    }
}

bool Root::SBNDReco1FrameSource::operator()(IFrame::pointer& out)
{
    out = nullptr;
    ++m_calls;
    const bool single = (m_entry >= 0 or m_run >= 0);
    if (single) {
        // One configured event, then EOS.
        if (m_calls > 2) {
            log->debug("past EOS at call={}", m_calls - 1);
            return false;
        }
        if (m_calls > 1) {
            log->debug("EOS at call={}", m_calls - 1);
            return true;  // EOS
        }
    }

    SBNDReco1::FilePtr fp;
    SBNDReco1::open_events(fp, m_filename, "SBNDReco1FrameSource");
    const Long64_t nentries = fp.events->GetEntries();

    Long64_t entry = -1;
    if (single) {
        entry = m_entry;
        if (entry < 0) {
            entry = SBNDReco1::find_entry(fp.events, m_run, m_subrun, m_event);
            if (entry < 0) {
                THROW(IOError() << errmsg{"SBNDReco1FrameSource: run/event not found in " + m_filename});
            }
        }
        if (entry >= nentries) {
            THROW(IOError() << errmsg{"SBNDReco1FrameSource: entry out of range in " + m_filename});
        }
    }
    else {
        // Stream all entries, then EOS.
        entry = m_calls - 1;
        if (entry > nentries) {
            log->debug("past EOS at call={}", m_calls - 1);
            return false;
        }
        if (entry == nentries) {
            log->debug("EOS after {} entries", nentries);
            return true;  // EOS
        }
    }

    art::EventAuxiliary aux;
    if (!SBNDReco1::read_event_aux(fp.events, entry, aux)) {
        THROW(IOError() << errmsg{"SBNDReco1FrameSource: no EventAuxiliary branch"});
    }
    const unsigned int run = aux.id_.subRun_.run_.run_;
    const unsigned int event = aux.id_.event_;

    art::Wrapper<std::vector<recob::Wire> > wires;
    if (!SBNDReco1::read_wrapper(fp.events, m_wire_product, entry, wires) or !wires.present) {
        THROW(IOError() << errmsg{"SBNDReco1FrameSource: no product " + m_wire_product});
    }

    // Frame time convention of larwirecell CookedFrameSource: the
    // event-vs-run-start wall-clock difference in SECONDS passed as a
    // bare number (existing dumps show tickinfo[0] = O(100)).  Replicated
    // verbatim so our frames are drop-in compatible; nothing downstream
    // consumes it as a physical time.
    double frame_time = 0.0;
    {
        const auto begin = SBNDReco1::run_begin_time(fp.file, run);
        if (begin.timeHigh_ or begin.timeLow_) {
            TTimeStamp t0(begin.timeHigh_, begin.timeLow_);
            TTimeStamp t1(aux.time_.timeHigh_, aux.time_.timeLow_);
            frame_time = t1.AsDouble() - t0.AsDouble();
        }
        else {
            log->warn("run {} not found in Runs tree; frame time set to 0", run);
        }
    }

    ITrace::vector all_traces;
    IFrame::trace_list_t indices;
    const size_t nwires = wires.obj.size();
    all_traces.reserve(nwires);
    for (size_t iw = 0; iw < nwires; ++iw) {
        const auto& wire = wires.obj[iw];
        std::vector<float> charge(m_nticks, 0.0f);
        for (const auto& range : wire.fSignalROI.ranges) {
            size_t tick = range.offset;
            for (const float val : range.values) {
                if (tick >= (size_t) m_nticks) break;
                charge[tick] = m_frame_scale * val;
                ++tick;
            }
        }
        indices.push_back(all_traces.size());
        all_traces.push_back(std::make_shared<Aux::SimpleTrace>(wire.fChannel, 0, charge));
    }

    // Per-channel summary, parallel to the wire product order.
    IFrame::trace_summary_t summary;
    {
        art::Wrapper<std::vector<double> > sums;
        if (SBNDReco1::read_wrapper(fp.events, m_summary_product, entry, sums) and sums.present) {
            if (sums.obj.size() == nwires) {
                summary.reserve(nwires);
                for (const double val : sums.obj) {
                    summary.push_back(val * m_summary_scale);
                }
            }
            else if (!sums.obj.empty()) {
                log->warn("summary product size {} != nwires {}; dropping summary",
                          sums.obj.size(), nwires);
            }
        }
    }

    // Channel mask map from [channel, first, last] triplets.
    WireCell::Waveform::ChannelMaskMap cmm;
    {
        art::Wrapper<std::vector<int> > masks;
        if (SBNDReco1::read_wrapper(fp.events, m_badmask_product, entry, masks) and masks.present) {
            const size_t ntriplet = masks.obj.size() / 3;
            for (size_t i = 0; i < ntriplet; ++i) {
                const int ch = masks.obj[3 * i];
                const int lo = masks.obj[3 * i + 1];
                const int hi = masks.obj[3 * i + 2];
                cmm["bad"][ch].push_back(std::make_pair(lo, hi));
            }
            log->debug("run {} event {}: {} bad-channel mask triplets", run, event, ntriplet);
        }
        else {
            log->warn("no badmask product {}; empty channel mask", m_badmask_product);
        }
    }

    auto sframe = new Aux::SimpleFrame(event, frame_time, all_traces, m_tick * units::ns, cmm);
    if (summary.size() == indices.size()) {
        sframe->tag_traces(m_frame_tag, indices, summary);
    }
    else {
        sframe->tag_traces(m_frame_tag, indices);
    }
    out = IFrame::pointer(sframe);

    log->debug("run {} event {} (entry {}): emit {} traces x {} ticks tag \"{}\"",
               run, event, entry, nwires, m_nticks, m_frame_tag);
    return true;
}
