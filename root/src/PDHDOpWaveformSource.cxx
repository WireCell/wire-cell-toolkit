#include "WireCellRoot/PDHDOpWaveformSource.h"
#include "PDHDLightTrig.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TH1.h"

#include <algorithm>
#include <map>
#include <vector>

WIRECELL_FACTORY(PDHDOpWaveformSource, WireCell::Root::PDHDOpWaveformSource,
                 WireCell::INamed,
                 WireCell::IFrameSource, WireCell::IConfigurable)

using namespace WireCell;

static const double dts_tick = 16 * units::ns;

namespace {
    // One decoana snippet.
    struct Snippet {
        int channel;
        int64_t tbin;  // ticks relative to the event's earliest snippet
        std::vector<float> values;
    };

    // Read all "waveform_*" TH1Ds under dir/<kind>/ for one channel.
    void read_snippets(TDirectory* chdir, const char* kind, int channel,
                       std::vector<Snippet>& out)
    {
        TDirectory* sub = dynamic_cast<TDirectory*>(chdir->Get(kind));
        if (!sub) return;
        TIter next(sub->GetListOfKeys());
        while (TKey* key = dynamic_cast<TKey*>(next())) {
            TH1* h = dynamic_cast<TH1*>(key->ReadObj());
            if (!h) continue;
            Snippet snip;
            snip.channel = channel;
            // The x-axis low edge encodes ticks numerically (the bin
            // width is 16 ns written in us -- only the offset is used).
            snip.tbin = std::llround(h->GetXaxis()->GetBinLowEdge(1));
            const int nbins = h->GetNbinsX();
            snip.values.resize(nbins);
            for (int i = 0; i < nbins; ++i) {
                snip.values[i] = h->GetBinContent(i + 1);
            }
            out.push_back(std::move(snip));
            delete h;
        }
    }
}

Root::PDHDOpWaveformSource::PDHDOpWaveformSource()
    : Aux::Logger("PDHDOpWaveformSource", "root")
{
}

Root::PDHDOpWaveformSource::~PDHDOpWaveformSource() {}

WireCell::Configuration Root::PDHDOpWaveformSource::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_filename;
    cfg["run"] = m_run;
    cfg["subrun"] = m_subrun;
    cfg["event"] = m_event;
    cfg["nominal_offset_us"] = m_nominal_offset_us;
    cfg["frame_tag"] = m_frame_tag;
    cfg["min_peak"] = m_min_peak;
    cfg["t0_search_window_us"] = m_t0_search_window_us;
    return cfg;
}

void Root::PDHDOpWaveformSource::configure(const WireCell::Configuration& cfg)
{
    m_filename = get(cfg, "filename", m_filename);
    m_run = get(cfg, "run", m_run);
    m_subrun = get(cfg, "subrun", m_subrun);
    m_event = get(cfg, "event", m_event);
    m_nominal_offset_us = get(cfg, "nominal_offset_us", m_nominal_offset_us);
    m_frame_tag = get(cfg, "frame_tag", m_frame_tag);
    m_min_peak = get(cfg, "min_peak", m_min_peak);
    m_t0_search_window_us = get(cfg, "t0_search_window_us", m_t0_search_window_us);
    if (m_filename.empty()) {
        THROW(ValueError() << errmsg{"PDHDOpWaveformSource: must supply input \"filename\""});
    }
    if (m_run < 0 or m_event < 0) {
        THROW(ValueError() << errmsg{"PDHDOpWaveformSource: must supply \"run\" and \"event\""});
    }
}

bool Root::PDHDOpWaveformSource::operator()(IFrame::pointer& out)
{
    out = nullptr;
    ++m_calls;
    if (m_calls > 2) {
        log->debug("past EOS at call={}", m_calls - 1);
        return false;
    }
    if (m_calls > 1) {
        log->debug("EOS at call={}", m_calls - 1);
        return true;  // EOS
    }

    TFile* tfile = TFile::Open(m_filename.c_str());
    if (!tfile or tfile->IsZombie()) {
        THROW(IOError() << errmsg{"PDHDOpWaveformSource: failed to open " + m_filename});
    }

    const auto trig = PDHDLight::read_trigger(tfile, m_run, m_subrun, m_event, m_nominal_offset_us);

    const auto evdir_name = String::format("decoana/run_%d_evt_%d", m_run, m_event);
    TDirectory* evdir = dynamic_cast<TDirectory*>(tfile->Get(evdir_name.c_str()));
    if (!evdir) {
        THROW(IOError() << errmsg{"PDHDOpWaveformSource: no " + evdir_name + " in input file"});
    }

    std::vector<Snippet> raw, deconv;
    TIter next(evdir->GetListOfKeys());
    while (TKey* key = dynamic_cast<TKey*>(next())) {
        const std::string name = key->GetName();
        if (name.rfind("ch", 0) != 0) continue;
        TDirectory* chdir = dynamic_cast<TDirectory*>(key->ReadObj());
        if (!chdir) continue;
        const int channel = std::atoi(name.c_str() + 2);
        read_snippets(chdir, "raw", channel, raw);
        read_snippets(chdir, "deconv", channel, deconv);
    }

    // Recover the absolute DTS time of the event's earliest snippet
    // (t_first, the tick origin of all snippet tbins): for each clear
    // deconv peak, each OpHit on the same channel votes
    //   t_first = PeakTimeAbs - snippet tbin - peak bin
    // and votes near rd_timestamp cluster sharply at the true value.
    TTree* oht = dynamic_cast<TTree*>(tfile->Get("opflashana/PerOpHitTree"));
    if (!oht) {
        THROW(IOError() << errmsg{"PDHDOpWaveformSource: no opflashana/PerOpHitTree"});
    }
    Int_t h_event, h_chan;
    Double_t h_peakabs;
    oht->SetBranchAddress("EventID", &h_event);
    oht->SetBranchAddress("OpChannel", &h_chan);
    oht->SetBranchAddress("PeakTimeAbs", &h_peakabs);
    std::map<int, std::vector<int64_t>> hit_ticks;  // channel -> peak times (abs ticks)
    for (Long64_t ent = 0; ent < oht->GetEntries(); ++ent) {
        oht->GetEntry(ent);
        if (h_event != m_event) continue;
        hit_ticks[h_chan].push_back(std::llround(h_peakabs));
    }
    oht->ResetBranchAddresses();

    const int64_t rd = trig.rd_timestamp;
    const int64_t window = std::llround(m_t0_search_window_us * 62.5);  // ticks
    std::map<int64_t, int> votes;
    for (const auto& snip : deconv) {
        auto it = hit_ticks.find(snip.channel);
        if (it == hit_ticks.end()) continue;
        const auto pk = std::max_element(snip.values.begin(), snip.values.end());
        if (*pk < m_min_peak) continue;
        const int64_t offset = snip.tbin + (pk - snip.values.begin());
        for (const int64_t hit : it->second) {
            const int64_t tf = hit - offset;
            if (std::llabs(tf - rd) < window) {
                ++votes[tf];
            }
        }
    }
    // The true value can split over a few neighboring ticks (the hit
    // peak time is not exactly the deconv argmax), so score each
    // candidate by the votes within +-2 ticks and take the raw-count
    // maximum inside the best neighborhood.
    int64_t t_first = rd;
    int best_votes = 0;
    int best_score = 0;
    for (const auto& [tf, n] : votes) {
        int score = 0;
        for (int64_t dt = -2; dt <= 2; ++dt) {
            auto it = votes.find(tf + dt);
            if (it != votes.end()) score += it->second;
        }
        if (score > best_score or (score == best_score and n > best_votes)) {
            best_score = score;
            best_votes = n;
            t_first = tf;
        }
    }
    if (best_votes < 3) {
        log->warn("run {} event {}: weak t_first recovery ({} votes), falling back near rd_timestamp",
                  m_run, m_event, best_votes);
    }
    log->debug("run {} event {}: t_first - rd = {} ticks ({} votes), {} raw + {} deconv snippets",
               m_run, m_event, t_first - rd, best_votes, raw.size(), deconv.size());

    tfile->Close();
    delete tfile;

    ITrace::vector all_traces;
    IFrame::trace_list_t raw_idx, deconv_idx;
    for (const auto& snip : raw) {
        raw_idx.push_back(all_traces.size());
        all_traces.push_back(std::make_shared<Aux::SimpleTrace>(snip.channel, (int)snip.tbin, snip.values));
    }
    for (const auto& snip : deconv) {
        deconv_idx.push_back(all_traces.size());
        all_traces.push_back(std::make_shared<Aux::SimpleTrace>(snip.channel, (int)snip.tbin, snip.values));
    }

    const double frame_time = static_cast<double>(t_first - static_cast<int64_t>(trig.tc_time)) * dts_tick;
    auto sframe = new Aux::SimpleFrame(m_event, frame_time, all_traces, dts_tick);
    sframe->tag_frame(m_frame_tag);
    sframe->tag_traces("raw", raw_idx);
    sframe->tag_traces("deconv", deconv_idx);
    out = IFrame::pointer(sframe);
    log->debug("output: {}", Aux::taginfo(out));
    return true;
}
