#include "WireCellRoot/PDHDOpFlashSource.h"
#include "PDHDLightTrig.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include "TFile.h"
#include "TTree.h"

#include <map>
#include <vector>

WIRECELL_FACTORY(PDHDOpFlashSource, WireCell::Root::PDHDOpFlashSource,
                 WireCell::INamed,
                 WireCell::ITensorSetSource, WireCell::IConfigurable)

using namespace WireCell;

// DTS timestamps count 16 ns ticks (62.5 MHz optical clock).
static const double dts_tick = 16 * units::ns;

Root::PDHDOpFlashSource::PDHDOpFlashSource()
    : Aux::Logger("PDHDOpFlashSource", "root")
{
}

Root::PDHDOpFlashSource::~PDHDOpFlashSource() {}

WireCell::Configuration Root::PDHDOpFlashSource::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_filename;
    cfg["run"] = m_run;
    cfg["subrun"] = m_subrun;
    cfg["event"] = m_event;
    cfg["nchan"] = m_nchan;
    cfg["nominal_offset_us"] = m_nominal_offset_us;
    cfg["include_summary"] = m_include_summary;
    cfg["include_ophits"] = m_include_ophits;
    return cfg;
}

void Root::PDHDOpFlashSource::configure(const WireCell::Configuration& cfg)
{
    m_filename = get(cfg, "filename", m_filename);
    m_run = get(cfg, "run", m_run);
    m_subrun = get(cfg, "subrun", m_subrun);
    m_event = get(cfg, "event", m_event);
    m_nchan = get(cfg, "nchan", m_nchan);
    m_nominal_offset_us = get(cfg, "nominal_offset_us", m_nominal_offset_us);
    m_include_summary = get(cfg, "include_summary", m_include_summary);
    m_include_ophits = get(cfg, "include_ophits", m_include_ophits);
    if (m_filename.empty()) {
        THROW(ValueError() << errmsg{"PDHDOpFlashSource: must supply input \"filename\""});
    }
    if (m_run < 0 or m_event < 0) {
        THROW(ValueError() << errmsg{"PDHDOpFlashSource: must supply \"run\" and \"event\""});
    }
}

bool Root::PDHDOpFlashSource::operator()(ITensorSet::pointer& out)
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
        THROW(IOError() << errmsg{"PDHDOpFlashSource: failed to open " + m_filename});
    }

    const auto trig = PDHDLight::read_trigger(tfile, m_run, m_subrun, m_event, m_nominal_offset_us);
    log->debug("run {} event {}: trigger tc_type={} offset_us={} ({} candidates)",
               m_run, m_event, trig.tc_type, trig.offset_us, trig.ncandidates);
    const double tc = static_cast<double>(trig.tc_time);

    // Flash enumeration and summary from PerFlashTree.  The example
    // files hold a single run, so EventID alone identifies the event.
    TTree* pft = dynamic_cast<TTree*>(tfile->Get("opflashana/PerFlashTree"));
    if (!pft) {
        THROW(IOError() << errmsg{"PDHDOpFlashSource: no opflashana/PerFlashTree"});
    }
    Int_t f_event, f_id;
    Double_t f_time;
    Float_t f_totpe, f_yc, f_zc, f_yw, f_zw;
    pft->SetBranchAddress("EventID", &f_event);
    pft->SetBranchAddress("FlashID", &f_id);
    pft->SetBranchAddress("FlashTime", &f_time);  // absolute DTS ticks (as double)
    pft->SetBranchAddress("TotalPE", &f_totpe);
    pft->SetBranchAddress("YCenter", &f_yc);
    pft->SetBranchAddress("ZCenter", &f_zc);
    pft->SetBranchAddress("YWidth", &f_yw);
    pft->SetBranchAddress("ZWidth", &f_zw);

    std::map<int, size_t> fid2row;
    std::vector<double> ftime_ns, ftotpe, fyc, fzc, fyw, fzw, fabs_dts;
    std::vector<int> fids;
    for (Long64_t ent = 0; ent < pft->GetEntries(); ++ent) {
        pft->GetEntry(ent);
        if (f_event != m_event) continue;
        fid2row[f_id] = fids.size();
        fids.push_back(f_id);
        ftime_ns.push_back((f_time - tc) * dts_tick);
        fabs_dts.push_back(f_time);
        ftotpe.push_back(f_totpe);
        fyc.push_back(f_yc * units::cm);
        fzc.push_back(f_zc * units::cm);
        fyw.push_back(f_yw * units::cm);
        fzw.push_back(f_zw * units::cm);
    }
    pft->ResetBranchAddresses();
    const size_t nflash = fids.size();

    // Per-OpDet PE from flashopdet/flash_opdet.  For PDHD OpChannel ==
    // OpDet (ChannelsPerOpDet = 1), 0..159.
    TTree* fot = dynamic_cast<TTree*>(tfile->Get("flashopdet/flash_opdet"));
    if (!fot) {
        THROW(IOError() << errmsg{"PDHDOpFlashSource: no flashopdet/flash_opdet"});
    }
    Int_t o_run, o_subrun, o_event, o_fid, o_opdet;
    Double_t o_pe;
    fot->SetBranchAddress("run", &o_run);
    fot->SetBranchAddress("subrun", &o_subrun);
    fot->SetBranchAddress("event", &o_event);
    fot->SetBranchAddress("flash_id", &o_fid);
    fot->SetBranchAddress("opdet", &o_opdet);
    fot->SetBranchAddress("pe", &o_pe);

    const size_t ncol = 1 + m_nchan;
    std::vector<double> matrix(nflash * ncol, 0.0);
    for (size_t r = 0; r < nflash; ++r) {
        matrix[r * ncol] = ftime_ns[r];
    }
    int nskipped = 0;
    for (Long64_t ent = 0; ent < fot->GetEntries(); ++ent) {
        fot->GetEntry(ent);
        if (o_run != m_run or o_event != m_event) continue;
        if (m_subrun >= 0 and o_subrun != m_subrun) continue;
        auto it = fid2row.find(o_fid);
        if (it == fid2row.end() or o_opdet < 0 or o_opdet >= m_nchan) {
            ++nskipped;
            continue;
        }
        matrix[it->second * ncol + 1 + o_opdet] += o_pe;
    }
    fot->ResetBranchAddresses();
    if (nskipped) {
        log->warn("run {} event {}: skipped {} flash_opdet rows (unknown flash or opdet)",
                  m_run, m_event, nskipped);
    }

    ITensor::vector* tensors = new ITensor::vector;

    {
        Configuration md;
        md["name"] = "opflash";
        tensors->push_back(std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{nflash, ncol}, matrix.data(), md));
    }

    if (m_include_summary) {
        const size_t nsum = 8;
        std::vector<double> summary(nflash * nsum, 0.0);
        for (size_t r = 0; r < nflash; ++r) {
            double* row = &summary[r * nsum];
            row[0] = fids[r];
            row[1] = ftotpe[r];
            row[2] = fyc[r];
            row[3] = fzc[r];
            row[4] = fyw[r];
            row[5] = fzw[r];
            row[6] = fabs_dts[r];
            row[7] = -1;  // nhits unknown for converted flashes
        }
        Configuration md;
        md["name"] = "flash_summary";
        tensors->push_back(std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{nflash, nsum}, summary.data(), md));
    }

    if (m_include_ophits) {
        TTree* oht = dynamic_cast<TTree*>(tfile->Get("opflashana/PerOpHitTree"));
        if (!oht) {
            THROW(IOError() << errmsg{"PDHDOpFlashSource: no opflashana/PerOpHitTree"});
        }
        Int_t h_event, h_chan;
        Double_t h_peakabs;
        Float_t h_width, h_area, h_amp, h_pe, h_ftt;
        oht->SetBranchAddress("EventID", &h_event);
        oht->SetBranchAddress("OpChannel", &h_chan);
        oht->SetBranchAddress("PeakTimeAbs", &h_peakabs);  // absolute DTS ticks
        oht->SetBranchAddress("Width", &h_width);          // us
        oht->SetBranchAddress("Area", &h_area);
        oht->SetBranchAddress("Amplitude", &h_amp);
        oht->SetBranchAddress("PE", &h_pe);
        oht->SetBranchAddress("FastToTotal", &h_ftt);

        const size_t nhcol = 9;
        std::vector<double> hits;
        for (Long64_t ent = 0; ent < oht->GetEntries(); ++ent) {
            oht->GetEntry(ent);
            if (h_event != m_event) continue;
            const double peak_ns = (h_peakabs - tc) * dts_tick;
            const double width_ns = h_width * units::us;
            hits.push_back(h_chan);
            hits.push_back(peak_ns);
            hits.push_back(width_ns);
            hits.push_back(h_area);
            hits.push_back(h_amp);
            hits.push_back(h_pe);
            hits.push_back(peak_ns - 0.5 * width_ns);  // approximate start time
            hits.push_back(-1);                        // flash association not in the dump
            hits.push_back(h_ftt);
        }
        oht->ResetBranchAddresses();
        Configuration md;
        md["name"] = "ophits";
        tensors->push_back(std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{hits.size() / nhcol, nhcol}, hits.data(), md));
        log->debug("run {} event {}: {} ophits", m_run, m_event, hits.size() / nhcol);
    }

    tfile->Close();
    delete tfile;

    Configuration md;
    md["run"] = m_run;
    md["subrun"] = trig.subrun;
    md["event"] = m_event;
    md["tc_type"] = trig.tc_type;
    md["tc_time_dts"] = std::to_string(trig.tc_time);
    md["rd_timestamp_dts"] = std::to_string(trig.rd_timestamp);
    md["offset_us"] = trig.offset_us;
    md["nchan"] = m_nchan;
    md["producer"] = "opflashana";

    out = std::make_shared<Aux::SimpleTensorSet>(m_event, md,
                                                 ITensor::shared_vector(tensors));
    log->debug("run {} event {}: emit {} flashes x {} channels", m_run, m_event, nflash, m_nchan);
    return true;
}
