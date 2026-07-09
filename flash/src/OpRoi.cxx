#include "WireCellFlash/OpRoi.h"

#include "WireCellAux/DftTools.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"

#include <algorithm>
#include <cmath>
#include <complex>

WIRECELL_FACTORY(OpRoi, WireCell::Flash::OpRoi,
                 WireCell::INamed,
                 WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;

Flash::OpRoi::OpRoi()
    : Aux::Logger("OpRoi", "flash")
{
}

Flash::OpRoi::~OpRoi() {}

WireCell::Configuration Flash::OpRoi::default_configuration() const
{
    Configuration cfg;
    cfg["intag"] = m_intag;
    cfg["outtag"] = m_outtag;
    cfg["hpf_tau_mhz"] = m_hpf_tau_mhz;
    cfg["roi_seed_nsigma"] = m_roi_seed_nsigma;
    cfg["roi_ext_nsigma"] = m_roi_ext_nsigma;
    cfg["roi_pad_pre"] = m_roi_pad_pre;
    cfg["roi_post_peak"] = m_roi_post_peak;
    cfg["veto_sigma"] = m_veto_sigma;
    cfg["apply_baseline"] = m_apply_baseline;
    cfg["veto_channels"] = Json::arrayValue;
    cfg["dft"] = "FftwDFT";
    return cfg;
}

void Flash::OpRoi::configure(const WireCell::Configuration& cfg)
{
    m_intag = get(cfg, "intag", m_intag);
    m_outtag = get(cfg, "outtag", m_outtag);
    m_hpf_tau_mhz = get(cfg, "hpf_tau_mhz", m_hpf_tau_mhz);
    m_roi_seed_nsigma = get(cfg, "roi_seed_nsigma", m_roi_seed_nsigma);
    m_roi_ext_nsigma = get(cfg, "roi_ext_nsigma", m_roi_ext_nsigma);
    m_roi_pad_pre = get(cfg, "roi_pad_pre", m_roi_pad_pre);
    m_roi_post_peak = get(cfg, "roi_post_peak", m_roi_post_peak);
    m_veto_sigma = get(cfg, "veto_sigma", m_veto_sigma);
    m_apply_baseline = get(cfg, "apply_baseline", m_apply_baseline);

    if (cfg.isMember("veto_channels")) {
        m_veto_channels.clear();
        for (const auto& jch : cfg["veto_channels"]) m_veto_channels.insert(jch.asInt());
    }

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);
}

void Flash::OpRoi::ensure_hpf(int n)
{
    if (n == m_hpf_n) return;
    // Sample frequency 62.5 MHz (16 ns ticks).  H(f) = 1 - exp(-(f/tau)^2),
    // built on the full-size (mirrored) spectrum: bin k -> f = df*min(k, n-k).
    const double sample_freq = 62.5;  // MHz
    const double df = sample_freq / n;
    m_hpf.resize(n);
    for (int k = 0; k < n; ++k) {
        const double f = df * (k <= n / 2 ? k : n - k);
        const double r = f / m_hpf_tau_mhz;
        m_hpf[k] = 1.0 - std::exp(-r * r);
    }
    m_hpf_n = n;
}

// Median via nth_element on the caller's vector, which is reordered.
static double vmedian_inplace(std::vector<float>& v)
{
    if (v.empty()) return 0.0;
    const size_t mid = v.size() / 2;
    std::nth_element(v.begin(), v.begin() + mid, v.end());
    return v[mid];
}

std::vector<float> Flash::OpRoi::clean(const std::vector<float>& decon)
{
    const int n = decon.size();
    std::vector<float> out(n, 0.0f);
    if (n == 0) return out;

    // 1. High-pass filter -> "ROI-finding" waveform h.
    auto D = Aux::DftTools::fwd_r2c(m_dft, decon);
    for (int k = 0; k < n; ++k) D[k] *= std::complex<float>((float) m_hpf[k], 0.0f);
    auto h = Aux::DftTools::inv_c2r(m_dft, D);  // real, size n

    // 2. Baseline removal: subtract median(h).
    m_scratch.assign(h.begin(), h.end());
    const double hmed = vmedian_inplace(m_scratch);
    for (int i = 0; i < n; ++i) h[i] -= (float) hmed;

    // 3. Noise rms = 1.4826 * MAD(h), and the ringing-channel veto.
    m_scratch.assign(h.begin(), h.end());
    const double hmed2 = vmedian_inplace(m_scratch);
    m_dev.resize(n);
    for (int i = 0; i < n; ++i) m_dev[i] = std::abs(h[i] - (float) hmed2);
    const double rms = 1.4826 * vmedian_inplace(m_dev);
    if (rms > m_veto_sigma) return out;  // zeroed channel

    // 4. ROIs (hysteresis): contiguous runs of h > ext (low) that reach seed
    // (high) somewhere inside.  Pad m_roi_pad_pre before the run and extend to at
    // least m_roi_post_peak ticks past the decon pulse peak (the late-light
    // window); a brighter tail above ext runs further on its own.  Merge overlaps.
    const double ext = m_roi_ext_nsigma * rms;
    const double seed = m_roi_seed_nsigma * rms;
    std::vector<std::pair<int, int>> rois;
    int i = 0;
    while (i < n) {
        if (h[i] > ext) {
            const int s = i;
            float hmax = h[i];
            int peak = s;
            float dpeak = decon[s];
            while (i < n && h[i] > ext) {
                if (h[i] > hmax) hmax = h[i];
                if (decon[i] > dpeak) { dpeak = decon[i]; peak = i; }
                ++i;
            }
            const int e = i - 1;
            if (hmax < seed) continue;   // no real seed inside the run -> drop
            const int ps = std::max(0, s - m_roi_pad_pre);
            const int pe = std::min(n - 1, std::max(e, peak + m_roi_post_peak));
            if (!rois.empty() && ps <= rois.back().second) {
                rois.back().second = std::max(rois.back().second, pe);
            }
            else {
                rois.emplace_back(ps, pe);
            }
        }
        else {
            ++i;
        }
    }

    // 5. Apply to the ORIGINAL decon: zero outside ROIs; per ROI subtract the
    // line through the endpoints so the ROI starts and ends exactly at zero.
    for (const auto& roi : rois) {
        const int s = roi.first, e = roi.second;
        const double v0 = decon[s], v1 = decon[e];
        for (int j = s; j <= e; ++j) {
            const double base = (e > s) ? v0 + (v1 - v0) * double(j - s) / double(e - s) : v0;
            out[j] = m_apply_baseline ? (float) (decon[j] - base) : decon[j];
        }
    }
    return out;
}

bool Flash::OpRoi::operator()(const IFrame::pointer& in, IFrame::pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS");
        return true;
    }

    auto traces = Aux::tagged_traces(in, m_intag);
    ITrace::vector all_traces(in->traces()->begin(), in->traces()->end());
    IFrame::trace_list_t out_idx;
    for (const auto& trace : traces) {
        const auto& d = trace->charge();
        ensure_hpf(d.size());
        // Hard per-channel veto: zero the trace outright (skip ROI cleaning).
        auto cleaned = m_veto_channels.count(trace->channel())
                           ? std::vector<float>(d.size(), 0.0f)
                           : clean(d);
        out_idx.push_back(all_traces.size());
        all_traces.push_back(std::make_shared<Aux::SimpleTrace>(trace->channel(), trace->tbin(), std::move(cleaned)));
    }

    // Forward any incoming masks (e.g. OpDecon's "saturation" flags) so the
    // full-stream branch carries them to OpHitFinder.  Empty by default.
    auto sframe = new Aux::SimpleFrame(in->ident(), in->time(), all_traces, in->tick(), in->masks());
    for (const auto& tag : in->frame_tags()) {
        sframe->tag_frame(tag);
    }
    for (const auto& tag : in->trace_tags()) {
        sframe->tag_traces(tag, in->tagged_traces(tag), in->trace_summary(tag));
    }
    sframe->tag_traces(m_outtag, out_idx);
    out = IFrame::pointer(sframe);
    log->debug("frame {}: ROI-cleaned {} of {} '{}' traces -> '{}'",
               in->ident(), out_idx.size(), traces.size(), m_intag, m_outtag);
    return true;
}
