#include "WireCellImg/MaskSlice.h"
#include "WireCellImg/ImgData.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/PlaneTools.h"

// Test to see if we can make slice time absolute
#undef SLICE_START_TIME_IS_RELATIVE

WIRECELL_FACTORY(MaskSlicer, WireCell::Img::MaskSlicer,
                 WireCell::INamed,
                 WireCell::IFrameSlicer, WireCell::IConfigurable)

WIRECELL_FACTORY(MaskSlices, WireCell::Img::MaskSlices,
                 WireCell::INamed,
                 WireCell::IFrameSlices, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

namespace WireCell {
    /// fixme: why invent TracelessFrame instead of using SimpleFrame with no traces?
    class TracelessFrame : public IFrame {
       public:
        TracelessFrame(int ident, double time = 0, double tick = 0.5 * units::microsecond)
          : m_ident(ident)
          , m_time(time)
          , m_tick(tick)
        {
        }
        virtual ~TracelessFrame(){};
        virtual const tag_list_t& frame_tags() const { return m_frame_tags; };
        virtual const tag_list_t& trace_tags() const { return m_trace_tags; };
        virtual const trace_list_t& tagged_traces(const tag_t& tag) const { return m_empty_trace_list; };
        virtual const trace_summary_t& trace_summary(const tag_t& tag) const { return m_empty_trace_summary; };
        virtual ITrace::shared_vector traces() const { return nullptr; };
        virtual Waveform::ChannelMaskMap masks() const { return Waveform::ChannelMaskMap(); }
        virtual int ident() const { return m_ident; }
        virtual double time() const { return m_time; };
        virtual double tick() const { return m_tick; };

       private:
        int m_ident;
        double m_time, m_tick;
        IFrame::tag_list_t m_frame_tags, m_trace_tags;
        trace_list_t m_empty_trace_list;
        trace_summary_t m_empty_trace_summary;
    };
}  // namespace WireCell

Img::MaskSliceBase::MaskSliceBase()
    : Aux::Logger("MaskSlice", "img")
{
}
Img::MaskSliceBase::~MaskSliceBase() {}
Img::MaskSlicer::~MaskSlicer() {}
Img::MaskSlices::~MaskSlices() {}

WireCell::Configuration Img::MaskSliceBase::default_configuration() const
{
    Configuration cfg;

    // Name of an IAnodePlane from which we can resolve channel ident to IChannel
    cfg["anode"] = "";

    // Number of ticks over which each output slice should sum
    cfg["tick_span"] = m_tick_span;

    // Used to judge active or not
    cfg["wiener_tag"] = m_wiener_tag;

    // Used to judge active or not
    cfg["summary_tag"] = m_summary_tag;

    // Assign charge for activities
    cfg["charge_tag"] = m_charge_tag;

    // use these traces for activity error
    cfg["error_tag"] = m_error_tag;

    // 0, 1, 2 for U, V, W
    // WirePlaneId::index()
    cfg["active_planes"][0] = 0;
    cfg["active_planes"][1] = 1;
    cfg["active_planes"][2] = 2;

    //
    cfg["dummy_planes"] = Json::arrayValue;

    //
    cfg["masked_planes"] = Json::arrayValue;

    //
    cfg["dummy_charge"] = m_dummy_charge;
    cfg["dummy_error"] = m_dummy_error;
    cfg["masked_charge"] = m_masked_charge;
    cfg["masked_error"] = m_masked_error;

    cfg["nthreshold"][0] = m_nthreshold[0];
    cfg["nthreshold"][1] = m_nthreshold[1];
    cfg["nthreshold"][2] = m_nthreshold[2];
    cfg["default_threshold"][0] = m_default_threshold[0];
    cfg["default_threshold"][1] = m_default_threshold[1];
    cfg["default_threshold"][2] = m_default_threshold[2];

    // if both are zero, (default) determine tbin range from input
    // frame.
    cfg["min_tbin"] = m_min_tbin;
    cfg["max_tbin"] = m_max_tbin;

    return cfg;
}

void Img::MaskSliceBase::configure(const WireCell::Configuration& cfg)
{
    m_anode = Factory::find_tn<IAnodePlane>(cfg["anode"].asString());  // throws
    m_tick_span = get(cfg, "tick_span", m_tick_span);
    m_wiener_tag = get<std::string>(cfg, "wiener_tag", m_wiener_tag);
    m_summary_tag = get<std::string>(cfg, "summary_tag", m_wiener_tag); // if not given, use wiener_tag
    m_charge_tag = get<std::string>(cfg, "charge_tag", m_charge_tag);
    m_error_tag = get<std::string>(cfg, "error_tag", m_error_tag);
    m_dummy_charge = get<double>(cfg, "dummy_charge", m_dummy_charge);
    m_dummy_error = get<double>(cfg, "dummy_error", m_dummy_error);
    log->debug("double dummy_error = {}", m_dummy_error);
    m_masked_charge = get<double>(cfg, "masked_charge", m_masked_charge);
    m_masked_error = get<double>(cfg, "masked_error", m_masked_error);
    m_min_tbin = get<int>(cfg, "min_tbin", m_min_tbin);
    m_max_tbin = get<int>(cfg, "max_tbin", m_max_tbin);
    if (cfg.isMember("active_planes")) {
        m_active_planes.clear();
        for (auto id : cfg["active_planes"]) {
            m_active_planes.push_back(id.asInt());
        }
    }
    if (cfg.isMember("dummy_planes")) {
        m_dummy_planes.clear();
        for (auto id : cfg["dummy_planes"]) {
            m_dummy_planes.push_back(id.asInt());
        }
    }
    if (cfg.isMember("masked_planes")) {
        m_masked_planes.clear();
        for (auto id : cfg["masked_planes"]) {
            m_masked_planes.push_back(id.asInt());
        }
    }
    for (auto id : m_active_planes) {
        log->debug("active planes: {}", id);
    }
    for (auto id : m_dummy_planes) {
        log->debug("dummy planes: {}", id);
    }
    for (auto id : m_masked_planes) {
        log->debug("masked planes: {}", id);
    }
    if (cfg.isMember("nthreshold")) {
        m_nthreshold.clear();
        for (auto var : cfg["nthreshold"]) {
            m_nthreshold.push_back(var.asDouble());
        }
    }
    if (cfg.isMember("default_threshold")) {
        m_default_threshold.clear();
        for (auto var : cfg["default_threshold"]) {
            m_default_threshold.push_back(var.asDouble());
        }
    }
}

bool Img::MaskSliceBase::thresholding(const WireCell::ITrace::ChargeSequence& wiener_charge,
                                   const WireCell::ITrace::ChargeSequence& gauss_charge, const size_t qind,
                                   const double threshold, const int tick_span, const bool verbose)
{
    const auto q_wiener = wiener_charge[qind];
    const auto q_gauss = gauss_charge[qind];
    if (q_wiener > threshold) {
        if (verbose == true) {
            log->debug("wiener: {} gauss: {} threshold: {} active: {}", q_wiener, q_gauss, threshold, true);
        }
        return true;
    }
    const int nq = gauss_charge.size();
    const int sbin = qind / tick_span;
    const int sbin_next = sbin + 1;
    const int sbin_prev = sbin - 1;
    double q_next = 0;
    double q_prev = 0;
    if (sbin_next * tick_span > 0 && sbin_next * tick_span < nq) {
        int count = 0;
        for (int i = sbin_next * tick_span; i < (sbin_next + 1) * tick_span && i < nq; ++i) {
            q_next += wiener_charge[i];
            ++count;
        }
        q_next /= count;
    }
    if (sbin_prev * tick_span > 0 && sbin_prev * tick_span < nq) {
        int count = 0;
        for (int i = sbin_prev * tick_span; i < (sbin_prev + 1) * tick_span && i < nq; ++i) {
            q_prev += wiener_charge[i];
            ++count;
        }
        q_prev /= count;
    }
    bool is_active = false;
    if ((q_gauss > q_next / 3. && q_next > threshold) || (q_gauss > q_prev / 3. && q_prev > threshold)) {
        is_active =  true;
    }
    if (verbose == true) {
        log->debug("gauss: {} wiener: {} prev: {} next: {} threshold: {} active: {}", q_gauss, q_wiener, q_prev, q_next, threshold,
                   is_active);
    }
    return is_active;
}

void Img::MaskSliceBase::slice(const IFrame::pointer& in, slice_map_t& svcmap)
{
    log->debug("call={} input frame: {}", m_count++, Aux::taginfo(in));

    const double tick = in->tick();
    const double span = tick * m_tick_span;

    // get charge traces
    auto charge_traces = Aux::tagged_traces(in, m_charge_tag);
    const size_t ntraces = charge_traces.size();

    if (!ntraces) {
        log->debug("no traces {} from frame {}", m_charge_tag, in->ident());
        return;
    }

    /// fixme: this should be created once in the calling function and given to each slice.
    // needed by ISlice
    auto tlframe = new TracelessFrame(in->ident(), in->time(), in->tick());
    auto tlframe_ptr = IFrame::pointer(tlframe);

    // min_max tbin from charge trace
    int min_tbin = m_min_tbin;
    int max_tbin = m_max_tbin;
    if (min_tbin == 0 and max_tbin == 0) {
        std::vector<int> tbins(ntraces, 0), tends(ntraces, 0);
        for (size_t tind = 0; tind<ntraces; ++tind) {
            const auto& tr = charge_traces[tind];
            tends[tind] = tbins[tind] = tr->tbin();
            tends[tind] += tr->charge().size();
        }
        min_tbin = *std::min_element(tbins.begin(), tbins.end());
        max_tbin = *std::max_element(tends.begin(), tends.end());
    }

    // [min, max) slice bins
    const size_t min_slicebin = floor((double)min_tbin/m_tick_span);
    const size_t max_slicebin = ceil((double)max_tbin/m_tick_span);

    // Need to create slices for all possible tbins to make sure syncing in BlobSetMerge
    // TODO: make this optional?
    // FIXME: how to handle same slice ident? Currently slicebin is used as slice ident
    for (size_t slicebin=min_slicebin; slicebin < max_slicebin; ++slicebin) {
        const double start = slicebin * m_tick_span * tick;
        const double span = m_tick_span * tick;
        auto s = new Img::Data::Slice(tlframe_ptr, slicebin, start, span);
        svcmap[slicebin] = s;
    }

    // get wiener traces
    auto wiener_traces = Aux::tagged_traces(in, m_wiener_tag);
    if(charge_traces.size()!=wiener_traces.size()) {
        log->error("trace size mismatch: {}: {} and {}: {}",
                   m_charge_tag, charge_traces.size(), m_wiener_tag, wiener_traces.size());
        THROW(RuntimeError() << errmsg{"charge_traces.size()!=wiener_traces.size()"});
    }

    // get RMS for traces
    auto const& summary = in->trace_summary(m_summary_tag);
    if (summary.size()!=wiener_traces.size()) {
        log->error("size unmatched for tag \"{}\", trace: {}, summary: {}. needed for threshold calc.", m_summary_tag, wiener_traces.size(), summary.size());        
        THROW(RuntimeError() << errmsg{"size unmatched"});
    }

    // get charge error traces
    auto charge_err_traces = Aux::tagged_traces(in, m_error_tag);
    if(charge_traces.size()!=charge_err_traces.size()) {
        log->error("trace size mismatch: {}: {} and {}: {}",
                   m_charge_tag, charge_traces.size(), m_error_tag, charge_err_traces.size());
        THROW(RuntimeError() << errmsg{"charge_traces.size()!=charge_err_traces.size()"});
    }

    // process active planes 
    for (size_t idx=0; idx < charge_traces.size(); ++idx) {
        auto trace = charge_traces[idx];
        auto wiener_trace = wiener_traces[idx];
        auto err_trace = charge_err_traces[idx];
        const int tbin = trace->tbin();
        const int chid = trace->channel();
        if (tbin != wiener_trace->tbin() || chid != wiener_trace->channel()) {
            THROW(RuntimeError() << errmsg{"tbin != wiener_trace->tbin() || chid != wiener_trace->channel()"});
        }
        if (tbin != err_trace->tbin() || chid != err_trace->channel()) {
            THROW(RuntimeError() << errmsg{"tbin != err_trace->tbin() || chid != err_trace->channel()"});
        }
        IChannel::pointer ich = m_anode->channel(chid);
        auto planeid = ich->planeid();
        if (std::find(m_active_planes.begin(),m_active_planes.end(),planeid.index())==m_active_planes.end()) {
            continue;
        }
        const auto& charge = trace->charge();
        const auto& wiener_charge = wiener_trace->charge();
        const auto& error = err_trace->charge();
        if (charge.size() != wiener_charge.size() || charge.size() != error.size()) {
            THROW(RuntimeError() << errmsg{"charge.size() != wiener_charge.size() || charge.size() != error.size()"});
        }
        const size_t nq = charge.size();
        for (size_t qind = 0; qind != nq; ++qind) {
            if ((tbin + (int) qind) < min_tbin || (tbin + (int) qind) >= max_tbin) {
                log->warn("trace {} {} exceeds given range [{},{}), breaking.", trace->channel(), trace->tbin(),
                          min_tbin, max_tbin);
                break;
            }
            const auto q = charge[qind];
            const auto e = error[qind];
            // threshold
            double threshold = m_nthreshold[planeid.index()] * summary[idx];
            if (threshold == 0) {
                threshold = m_default_threshold[planeid.index()];
            }
            bool verbose = false;
            // TODO: remove debug code
            // if ( (tbin + qind == 824 && ich->ident() == 5575) ||
            // (tbin + qind == 1132 && ich->ident() == 5526)
            // ) {
            //     log->debug("tick: {} channel: {}", tbin + qind, ich->ident());
            //     verbose = true;
            // }
            if (thresholding(wiener_charge, charge, qind, threshold, m_tick_span, verbose) != true) continue;
            size_t slicebin = (tbin + qind) / m_tick_span;
            auto s = svcmap[slicebin];
            if (!s) {
#ifdef SLICE_START_TIME_IS_RELATIVE
                const double start = slicebin * span;  // thus relative to slice frame's time.
#else
                // Slice start time is absolute with frame time as origin
                const double start = in->time() + slicebin * span;  
#endif
                s = new Img::Data::Slice(tlframe_ptr, slicebin, start, span);
                svcmap[slicebin] = s;
            }
            // TODO: how to handle error?
            s->sum(ich, {q, e});
            // if (tbin+qind == 0 && chid > 1517 && chid < 1523) {
            //     log->debug("chid: {} slicebin: {} charge: {} error {}", chid, slicebin, s->activity()[ich].value(),
            //                s->activity()[ich].uncertainty());
            // }
        }
    }

    // dummy: to form slices active through all ticks/channels with special charge/error
    for (int plane_index : m_dummy_planes) {
        auto ichans = Aux::plane_channels(m_anode, plane_index);
        log->debug("dummy plane: {} size {}", plane_index, ichans.size());
        for (auto ich : ichans) {
            for (size_t slicebin = min_slicebin; slicebin < max_slicebin; ++slicebin) {
                auto s = svcmap[slicebin];
                if (!s) {
#ifdef SLICE_START_TIME_IS_RELATIVE
                    const double start = slicebin * span;  // thus relative to slice frame's time.
#else
                    // Slice start time is absolute with frame time as origin
                    const double start = in->time() + slicebin * span;  
#endif
                    s = new Img::Data::Slice(tlframe_ptr, slicebin, start, span);
                    svcmap[slicebin] = s;
                }
                s->assign(ich, {(float)m_dummy_charge, (float)m_dummy_error});
            }
        }
    }

    // masked slices
    auto cmm = in->masks()["bad"];
    log->debug("cmm.size(): {}", cmm.size());
    for (auto ch_tbins : cmm) {
        const int chid = ch_tbins.first;
        const auto& tbins = ch_tbins.second;
        // if (chid == 802) {
        //     for (auto tbin : tbins) {
        //         log->debug("chid {} tbin {} {}", chid, tbin.first, tbin.second);
        //     }
        // }
        IChannel::pointer ich = m_anode->channel(chid);
        if (!ich) {
            continue;
        }
        auto planeid = ich->planeid();
        if (std::find(m_masked_planes.begin(),m_masked_planes.end(),planeid.index())==m_masked_planes.end()) {
            continue;
        }
        /// slicebin based, should be faster for fewer slices
        for (size_t slicebin = min_slicebin; slicebin < max_slicebin; ++slicebin) {
            for (auto tbin : tbins) {
                if (slicebin * m_tick_span > (size_t)tbin.second || (slicebin + 1) * m_tick_span <= (size_t)tbin.first) continue;
                auto s = svcmap[slicebin];
                if (!s) {
#ifdef SLICE_START_TIME_IS_RELATIVE
                    const double start = slicebin * span;  // thus relative to slice frame's time.
#else
                    // Slice start time is absolute with frame time as origin
                    const double start = in->time() + slicebin * span;  
#endif
                    s = new Img::Data::Slice(tlframe_ptr, slicebin, start, span);
                    svcmap[slicebin] = s;
                }
                s->assign(ich, {(float) m_masked_charge, (float) m_masked_error});
            }
        }
        /// tbin based
        // for (auto tbin : tbins) {
        //     // log->debug("t: {} {}", tbin.first, tbin.second);
        //     for (auto t = tbin.first; t != tbin.second; ++t) {
        //         if (t < min_tbin || t >= max_tbin) {
        //             log->warn("cmm {} {} exceeds given range [{},{}), breaking.", chid, t, min_tbin, max_tbin);
        //             break;
        //         }
        //         size_t slicebin = t / m_tick_span;
        //         auto s = svcmap[slicebin];
        //         if (!s) {
        //             const double start = slicebin * span;  // thus relative to slice frame's time.
        //             s = new Img::Data::Slice(tlframe_ptr, slicebin, start, span);
        //             svcmap[slicebin] = s;
        //         }
        //         s->assign(ich, {(float)m_masked_charge, (float)m_masked_error});
        //     }
        // }
    }

    log->debug("nslices={} from ntraces={}, tbin=[{}, {}]",
               svcmap.size(), ntraces, min_tbin, max_tbin);
}

bool Img::MaskSlicer::operator()(const input_pointer& in, output_pointer& out)
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
    out = make_shared<Img::Data::SliceFrame>(islices, in->ident(), in->time());

    return true;
}

bool Img::MaskSlices::operator()(const input_pointer& in, output_queue& slices)
{
    if (!in) {
        log->debug("EOS");
        slices.push_back(nullptr);
        return true;  // eos
    }

    // Slices will be sparse in general.  Index by a "slice bin" number
    slice_map_t svcmap;
    slice(in, svcmap);

    // intern
    for (auto sit : svcmap) {
        auto s = sit.second;

        /// debug
        // double qtot = 0;
        // for (const auto& a : s->activity()) {
        //     qtot += a.second;
        // }

        // log->debug("slice: id={} t={} activity={}", s->ident(), s->start(), s->activity().size());
        slices.push_back(ISlice::pointer(s));
    }

    return true;
}
