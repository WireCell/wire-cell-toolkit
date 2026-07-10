#include "WireCellRoot/PDVDOpWaveformSource.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"

#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <vector>

WIRECELL_FACTORY(PDVDOpWaveformSource, WireCell::Root::PDVDOpWaveformSource,
                 WireCell::INamed,
                 WireCell::IFrameSource, WireCell::IConfigurable)

using namespace WireCell;

static const double dts_tick = 16 * units::ns;

Root::PDVDOpWaveformSource::PDVDOpWaveformSource()
    : Aux::Logger("PDVDOpWaveformSource", "root")
{
}

Root::PDVDOpWaveformSource::~PDVDOpWaveformSource() {}

WireCell::Configuration Root::PDVDOpWaveformSource::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_filename;
    cfg["run"] = m_run;
    cfg["subrun"] = m_subrun;
    cfg["event"] = m_event;
    cfg["frame_tag"] = m_frame_tag;
    cfg["opch_lo"] = m_opch_lo;
    cfg["opch_hi"] = m_opch_hi;
    cfg["trig_sample"] = m_trig_sample;
    cfg["snippet_nsamp"] = m_snippet_nsamp;
    return cfg;
}

void Root::PDVDOpWaveformSource::configure(const WireCell::Configuration& cfg)
{
    m_filename = get(cfg, "filename", m_filename);
    m_run = get(cfg, "run", m_run);
    m_subrun = get(cfg, "subrun", m_subrun);
    m_event = get(cfg, "event", m_event);
    m_frame_tag = get(cfg, "frame_tag", m_frame_tag);
    m_opch_lo = get(cfg, "opch_lo", m_opch_lo);
    m_opch_hi = get(cfg, "opch_hi", m_opch_hi);
    m_trig_sample = get(cfg, "trig_sample", m_trig_sample);
    m_snippet_nsamp = get(cfg, "snippet_nsamp", m_snippet_nsamp);
    if (m_filename.empty()) {
        THROW(ValueError() << errmsg{"PDVDOpWaveformSource: must supply input \"filename\""});
    }
    if (m_run < 0 or m_event < 0) {
        THROW(ValueError() << errmsg{"PDVDOpWaveformSource: must supply \"run\" and \"event\""});
    }
}

bool Root::PDVDOpWaveformSource::operator()(IFrame::pointer& out)
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
        THROW(IOError() << errmsg{"PDVDOpWaveformSource: failed to open " + m_filename});
    }
    // Newer extractions nest the tree under a TDirectory.
    TTree* tree = dynamic_cast<TTree*>(tfile->Get("raw_waveform"));
    if (!tree) {
        tree = dynamic_cast<TTree*>(tfile->Get("rawdump/raw_waveform"));
    }
    if (!tree) {
        THROW(IOError() << errmsg{"PDVDOpWaveformSource: no raw_waveform tree in " + m_filename});
    }

    Int_t b_run, b_subrun, b_event, b_opch, b_nsamp;
    Double_t b_timestamp;
    std::vector<short>* b_adc = nullptr;
    tree->SetBranchAddress("run", &b_run);
    tree->SetBranchAddress("subrun", &b_subrun);
    tree->SetBranchAddress("event", &b_event);
    tree->SetBranchAddress("opchannel", &b_opch);
    tree->SetBranchAddress("nsamp", &b_nsamp);
    tree->SetBranchAddress("timestamp", &b_timestamp);
    tree->SetBranchAddress("adc", &b_adc);

    // Pass 1 (only the selection/timing branches enabled, so the adc and
    // x/y/z baskets are not deserialized): find the event's records, their
    // start ticks and the common tick origin t0 = min start over ALL
    // channels -- deliberately ignoring the opch selection so that
    // per-population source instances share one origin.
    tree->SetBranchStatus("*", false);
    for (const char* b : {"run", "subrun", "event", "opchannel", "nsamp", "timestamp"}) {
        tree->SetBranchStatus(b, true);
    }
    const Long64_t nent = tree->GetEntries();
    std::vector<std::pair<Long64_t, int64_t>> selected;  // entry, start tick
    int64_t t0 = 0;
    bool have_t0 = false;
    int nrecords = 0;
    for (Long64_t ent = 0; ent < nent; ++ent) {
        tree->GetEntry(ent);
        if (b_run != m_run) continue;
        if (m_subrun >= 0 and b_subrun != m_subrun) continue;
        if (b_event != m_event) continue;
        // timestamp is us on the 16 ns clock; *62.5 is tick-exact in double.
        const int64_t start = std::llround(b_timestamp * 62.5)
            - (b_nsamp <= m_snippet_nsamp ? m_trig_sample : 0);
        if (!have_t0 or start < t0) {
            t0 = start;
            have_t0 = true;
        }
        ++nrecords;
        if (m_opch_lo >= 0 and b_opch < m_opch_lo) continue;
        if (m_opch_hi >= 0 and b_opch > m_opch_hi) continue;
        selected.push_back({ent, start});
    }
    if (!have_t0) {
        THROW(IOError() << errmsg{String::format(
            "PDVDOpWaveformSource: no records for run %d event %d in %s",
            m_run, m_event, m_filename.c_str())});
    }

    // Pass 2: read the selected waveforms.
    tree->SetBranchStatus("adc", true);
    ITrace::vector all_traces;
    IFrame::trace_list_t raw_idx;
    for (const auto& [ent, start] : selected) {
        tree->GetEntry(ent);
        std::vector<float> values(b_adc->begin(), b_adc->end());
        raw_idx.push_back(all_traces.size());
        all_traces.push_back(std::make_shared<Aux::SimpleTrace>(
            b_opch, (int)(start - t0), std::move(values)));
    }
    tree->ResetBranchAddresses();
    tfile->Close();
    delete tfile;

    auto sframe = new Aux::SimpleFrame(m_event, 0.0, all_traces, dts_tick);
    sframe->tag_frame(m_frame_tag);
    sframe->tag_traces("raw", raw_idx);
    out = IFrame::pointer(sframe);
    log->debug("run {} event {}: {} of {} records selected (opch {}..{}), t0 = {} ticks",
               m_run, m_event, selected.size(), nrecords, m_opch_lo, m_opch_hi, t0);
    log->debug("output: {}", Aux::taginfo(out));
    return true;
}
