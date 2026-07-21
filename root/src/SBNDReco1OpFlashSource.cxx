#include "WireCellRoot/SBNDReco1OpFlashSource.h"
#include "SBNDReco1Reader.h"
#include "SBNDReco1FrameShift.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

WIRECELL_FACTORY(SBNDReco1OpFlashSource, WireCell::Root::SBNDReco1OpFlashSource,
                 WireCell::INamed,
                 WireCell::ITensorSetSource, WireCell::IConfigurable)

using namespace WireCell;

Root::SBNDReco1OpFlashSource::SBNDReco1OpFlashSource()
    : Aux::Logger("SBNDReco1OpFlashSource", "root")
{
}

Root::SBNDReco1OpFlashSource::~SBNDReco1OpFlashSource() {}

WireCell::Configuration Root::SBNDReco1OpFlashSource::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_filename;
    cfg["entry"] = m_entry;
    cfg["run"] = m_run;
    cfg["subrun"] = m_subrun;
    cfg["event"] = m_event;
    cfg["flash_product"] = m_flash_product;
    cfg["ptb_product"] = m_ptb_product;
    cfg["tdc_product"] = m_tdc_product;
    cfg["rawheader_prefix"] = m_rawheader_prefix;
    cfg["frameshift_product"] = m_frameshift_product;
    cfg["npmts"] = m_npmts;
    cfg["caf_offset_mode"] = m_caf_offset_mode;
    cfg["caf_offset_override"] = m_caf_offset_override;
    return cfg;
}

void Root::SBNDReco1OpFlashSource::configure(const WireCell::Configuration& cfg)
{
    m_filename = get(cfg, "filename", m_filename);
    m_entry = get(cfg, "entry", m_entry);
    m_run = get(cfg, "run", m_run);
    m_subrun = get(cfg, "subrun", m_subrun);
    m_event = get(cfg, "event", m_event);
    m_flash_product = get(cfg, "flash_product", m_flash_product);
    m_ptb_product = get(cfg, "ptb_product", m_ptb_product);
    m_tdc_product = get(cfg, "tdc_product", m_tdc_product);
    m_rawheader_prefix = get(cfg, "rawheader_prefix", m_rawheader_prefix);
    m_frameshift_product = get(cfg, "frameshift_product", m_frameshift_product);
    m_npmts = get(cfg, "npmts", m_npmts);
    m_caf_offset_mode = get(cfg, "caf_offset_mode", m_caf_offset_mode);
    m_caf_offset_override = get(cfg, "caf_offset_override", m_caf_offset_override);
    if (m_filename.empty()) {
        THROW(ValueError() << errmsg{"SBNDReco1OpFlashSource: must supply input \"filename\""});
    }
    // Selection: entry >= 0, or run >= 0 (+event), or neither = stream all.
    if (m_entry < 0 and m_run >= 0 and m_event < 0) {
        THROW(ValueError() << errmsg{"SBNDReco1OpFlashSource: \"run\" selection also needs \"event\""});
    }
    if (m_caf_offset_mode != "none" and m_caf_offset_mode != "product" and
        m_caf_offset_mode != "auto" and m_caf_offset_mode != "override") {
        THROW(ValueError() << errmsg{"SBNDReco1OpFlashSource: bad caf_offset_mode \"" +
                                     m_caf_offset_mode + "\""});
    }
}

bool Root::SBNDReco1OpFlashSource::operator()(ITensorSet::pointer& out)
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
    SBNDReco1::open_events(fp, m_filename, "SBNDReco1OpFlashSource");
    const Long64_t nentries = fp.events->GetEntries();

    Long64_t entry = -1;
    if (single) {
        entry = m_entry;
        if (entry < 0) {
            entry = SBNDReco1::find_entry(fp.events, m_run, m_subrun, m_event);
            if (entry < 0) {
                THROW(IOError() << errmsg{"SBNDReco1OpFlashSource: run/event not found in " + m_filename});
            }
        }
        if (entry >= nentries) {
            THROW(IOError() << errmsg{"SBNDReco1OpFlashSource: entry out of range in " + m_filename});
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
        THROW(IOError() << errmsg{"SBNDReco1OpFlashSource: no EventAuxiliary branch"});
    }
    const unsigned int run = aux.id_.subRun_.run_.run_;
    const unsigned int subrun = aux.id_.subRun_.subRun_;
    const unsigned int event = aux.id_.event_;

    art::Wrapper<std::vector<recob::OpFlash> > flashes;
    if (!SBNDReco1::read_wrapper(fp.events, m_flash_product, entry, flashes) or !flashes.present) {
        THROW(IOError() << errmsg{"SBNDReco1OpFlashSource: no product " + m_flash_product});
    }

    // Matrix layout of larwirecell OpFlashSource: col 0 = Time() [us]
    // scaled to WCT units (ns), cols 1.. = PE per OpDet.
    const size_t nflash = flashes.obj.size();
    const size_t ncol = 1 + m_npmts;
    std::vector<double> matrix(nflash * ncol, 0.0);
    for (size_t r = 0; r < nflash; ++r) {
        const auto& flash = flashes.obj[r];
        matrix[r * ncol] = flash.fTime * units::microsecond;
        const size_t npe = flash.fPEperOpDet.size();
        if (npe > (size_t) m_npmts) {
            THROW(ValueError() << errmsg{"SBNDReco1OpFlashSource: flash has more PE entries than npmts"});
        }
        for (size_t ipmt = 0; ipmt < npe; ++ipmt) {
            matrix[r * ncol + 1 + ipmt] = flash.fPEperOpDet[ipmt];
        }
    }

    ITensor::vector* tensors = new ITensor::vector;
    {
        Configuration md;
        md["name"] = "opflash";
        tensors->push_back(std::make_shared<Aux::SimpleTensor>(
            ITensor::shape_t{nflash, ncol}, matrix.data(), md));
    }

    Configuration setmd;
    setmd["run"] = run;
    setmd["subrun"] = subrun;
    setmd["event"] = event;
    setmd["nchan"] = m_npmts;

    if (m_caf_offset_mode == "override") {
        setmd["frame_apply_at_caf"] = m_caf_offset_override;
        log->debug("run {} event {}: frame_apply_at_caf override {} ns",
                   run, event, m_caf_offset_override);
    }
    else if (m_caf_offset_mode == "product") {
        // The authoritative per-event value, written by the sbndcode
        // FrameShift producer.  Fail loudly rather than emit un-offset
        // flashes: a file without the product needs mode "auto"/"none".
        art::Wrapper<sbnd::timing::FrameShiftInfo> fsi;
        if (!SBNDReco1::read_wrapper(fp.events, m_frameshift_product, entry, fsi) or !fsi.present) {
            THROW(IOError() << errmsg{"SBNDReco1OpFlashSource: caf_offset_mode=product but no product " +
                                      m_frameshift_product + " in " + m_filename});
        }
        setmd["frame_apply_at_caf"] = fsi.obj.fFrameApplyAtCaf;
        setmd["caf_timing_type"] = fsi.obj.fTimingType;
        log->debug("run {} event {}: timing_type {} frame_apply_at_caf {} ns (product)",
                   run, event, fsi.obj.fTimingType, fsi.obj.fFrameApplyAtCaf);
    }
    else if (m_caf_offset_mode == "auto") {
        const uint64_t raw_ts = SBNDReco1::raw_header_timestamp(fp.events, entry, m_rawheader_prefix);
        art::Wrapper<std::vector<raw::ptb::sbndptb> > ptbs;
        art::Wrapper<std::vector<sbnd::timing::DAQTimestamp> > tdcs;
        const bool have_ptb =
            SBNDReco1::read_wrapper(fp.events, m_ptb_product, entry, ptbs) and ptbs.present;
        const bool have_tdc =
            SBNDReco1::read_wrapper(fp.events, m_tdc_product, entry, tdcs) and tdcs.present;
        if (!raw_ts or !have_ptb or !have_tdc) {
            log->warn("run {} event {}: timing products missing (raw={} ptb={} tdc={}); "
                      "no frame_apply_at_caf",
                      run, event, raw_ts != 0, have_ptb, have_tdc);
        }
        else {
            const auto fs = SBNDReco1::compute_frame_shift(raw_ts, ptbs.obj, tdcs.obj,
                                                           SBNDReco1::FrameShiftConfig{});
            if (!fs.valid) {
                log->warn("run {} event {}: frame shift failed ({}); no frame_apply_at_caf",
                          run, event, fs.status);
            }
            else {
                // Approximation for pre-FrameShift files.  The true
                // FrameApplyAtCaf (caf_offset_mode=product) equals the
                // decoded-frame-to-SPEC-TDC-RWM shift; on the 48-event
                // dev sample frame_etrig - frame_default sits 43..482 ns
                // (mean 262) below it -- close enough to land the in-time
                // beam flash in the +0.3..1.9 us window, not exact.
                // Prefer "product" whenever the producer has run.
                const double offset_ns = static_cast<double>(fs.frame_etrig) -
                                         static_cast<double>(fs.frame_default);
                setmd["frame_apply_at_caf"] = offset_ns;
                // Diagnostics for the offset-validation study.
                setmd["caf_stream"] = fs.status;
                setmd["caf_hlt_etrig"] = fs.hlt_etrig;
                log->debug("run {} event {}: stream {} hlt_etrig {} frame_apply_at_caf {} ns",
                           run, event, fs.status, fs.hlt_etrig, offset_ns);
            }
        }
    }

    out = std::make_shared<Aux::SimpleTensorSet>(event, setmd,
                                                 ITensor::shared_vector(tensors));
    log->debug("run {} event {} (entry {}): emit {} flashes x {} channels from {}",
               run, event, entry, nflash, m_npmts, m_flash_product);
    return true;
}
