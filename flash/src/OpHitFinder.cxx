#include "WireCellFlash/OpHitFinder.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <algorithm>
#include <cmath>

WIRECELL_FACTORY(OpHitFinder, WireCell::Flash::OpHitFinder,
                 WireCell::INamed,
                 WireCell::IFrameTensorSet, WireCell::IConfigurable)

using namespace WireCell;

Flash::OpHitFinder::OpHitFinder()
    : Aux::Logger("OpHitFinder", "flash")
{
}

Flash::OpHitFinder::~OpHitFinder() {}

WireCell::Configuration Flash::OpHitFinder::default_configuration() const
{
    Configuration cfg;
    cfg["intag"] = m_intag;
    cfg["scale"] = m_scale;
    cfg["spe_area"] = m_spe_area;
    cfg["hit_threshold"] = m_hit_threshold;
    cfg["ped_nsamples"] = m_ped_nsamples;
    // Robust per-channel baseline for the continuous full stream (default OFF
    // -> head method, bit-identical).  See the header / fullstream doc.
    cfg["robust_baseline"] = m_robust_baseline;
    cfg["robust_nsigma"] = m_robust_nsigma;
    cfg["robust_veto_sigma"] = m_robust_veto_sigma;
    cfg["fixed_ped_sigma"] = m_fixed_ped_sigma;
    // AlgoSlidingWindow parameters, dune_ophit_finder_deco values.
    Configuration algo;
    algo["adc_threshold"] = 3.0;
    algo["nsigma_threshold"] = 1.0;
    algo["tail_adc_threshold"] = 1.0;
    algo["tail_nsigma_threshold"] = 1.0;
    algo["end_adc_threshold"] = 1.0;
    algo["end_nsigma_threshold"] = 1.0;
    algo["min_pulse_width"] = 1;
    algo["num_presample"] = 2;
    algo["num_postsample"] = 2;
    // Overlapping-pulse splitting (post-processing of each found pulse).
    // Default OFF -> bit-identical to the larana AlgoSlidingWindow output.
    algo["split_enable"] = false;            // master gate
    algo["split_min_prominence"] = 0.4;      // valley depth as a fraction of the smaller flanking peak
    algo["split_min_prominence_abs"] = 0.0;  // absolute valley-depth floor (scaled units); 0 = relative only
    algo["split_min_peak"] = 3.0;            // both sub-peaks must exceed this (scaled, ped-subtracted)
    algo["split_min_separation"] = 2;        // min ticks between two accepted maxima (noise guard)
    cfg["algo"] = algo;
    return cfg;
}

void Flash::OpHitFinder::configure(const WireCell::Configuration& cfg)
{
    auto defs = default_configuration();
    m_intag = get(cfg, "intag", m_intag);
    m_scale = get(cfg, "scale", m_scale);
    m_spe_area = get(cfg, "spe_area", m_spe_area);
    m_hit_threshold = get(cfg, "hit_threshold", m_hit_threshold);
    m_ped_nsamples = get(cfg, "ped_nsamples", m_ped_nsamples);
    m_robust_baseline = get(cfg, "robust_baseline", m_robust_baseline);
    m_robust_nsigma = get(cfg, "robust_nsigma", m_robust_nsigma);
    m_robust_veto_sigma = get(cfg, "robust_veto_sigma", m_robust_veto_sigma);
    m_fixed_ped_sigma = get(cfg, "fixed_ped_sigma", m_fixed_ped_sigma);
    m_algo = defs["algo"];
    if (cfg.isMember("algo")) {
        for (const auto& key : cfg["algo"].getMemberNames()) {
            m_algo[key] = cfg["algo"][key];
        }
    }
}

std::vector<Flash::OpHitFinder::Pulse> Flash::OpHitFinder::sliding_window(
    const std::vector<short>& wf, double ped_mean, double ped_sigma, const Configuration& pars)
{
    // Direct port of pmtana::AlgoSlidingWindow::RecoPulse (positive
    // polarity) with the constant pedestal of PedAlgoEdges (kHEAD).
    const double adc_thres = pars["adc_threshold"].asDouble();
    const double nsigma = pars["nsigma_threshold"].asDouble();
    const double tail_adc_thres = pars["tail_adc_threshold"].asDouble();
    const double tail_nsigma = pars["tail_nsigma_threshold"].asDouble();
    const double end_adc_thres = pars["end_adc_threshold"].asDouble();
    const double end_nsigma = pars["end_nsigma_threshold"].asDouble();
    const size_t min_width = pars["min_pulse_width"].asUInt();
    const size_t num_presample = pars["num_presample"].asUInt();
    const size_t num_postsample = pars["num_postsample"].asUInt();

    const double start_threshold = std::max(adc_thres, nsigma * ped_sigma);
    const double tail_threshold = std::max(tail_adc_thres, tail_nsigma * ped_sigma);
    const double end_threshold = std::max(end_adc_thres, end_nsigma * ped_sigma);

    std::vector<Pulse> pulses;
    Pulse pulse;
    bool fire = false, in_tail = false, in_post = false;
    int post_integration = 0;

    auto register_pulse = [&](int t_end) {
        pulse.t_end = t_end;
        if (size_t(pulse.t_end - pulse.t_start) >= min_width) {
            pulses.push_back(pulse);
        }
        pulse = Pulse{};
    };

    for (size_t i = 0; i < wf.size(); ++i) {
        const double value = double(wf[i]) - ped_mean;

        if ((!fire || in_tail || in_post) && value > start_threshold) {
            if (in_tail) {
                register_pulse(int(i) - 1);
            }
            int buffer_num_index = 0;
            if (!pulses.empty()) {
                buffer_num_index = int(i) - pulses.back().t_end - 1;
            }
            else {
                buffer_num_index = std::min(num_presample, i);
            }
            if (buffer_num_index > int(num_presample)) buffer_num_index = num_presample;
            if (in_post) {
                int t_end = int(i) - buffer_num_index;
                if (t_end > 0) --t_end;  // leave a gap, if we can
                register_pulse(t_end);
            }
            pulse.t_start = i - buffer_num_index;
            for (int pre = pulse.t_start; pre < int(i); ++pre) {
                const double pre_adc = double(wf[pre]) - ped_mean;
                if (pre_adc > 0.) pulse.area += pre_adc;
            }
            fire = true;
            in_tail = in_post = false;
        }

        if (fire && value < tail_threshold) {
            fire = false;
            in_tail = true;
            in_post = false;
        }

        if ((fire || in_tail) && value < end_threshold) {
            in_post = true;
            fire = in_tail = false;
            post_integration = num_postsample;
        }

        if (in_post && post_integration < 1) {
            register_pulse(int(i) - 1);
            fire = in_tail = in_post = false;
        }

        if (fire || in_tail || in_post) {
            pulse.area += value;
            if (pulse.peak < value) {
                pulse.peak = value;
                pulse.t_max = i;
            }
            if (in_post) --post_integration;
        }
    }
    if (fire || in_tail || in_post) {
        register_pulse(int(wf.size()) - 1);
    }
    return pulses;
}

std::vector<Flash::OpHitFinder::Pulse> Flash::OpHitFinder::split_pulse(
    const std::vector<short>& wf, double ped_mean, const Pulse& pulse,
    const Configuration& pars)
{
    // Disabled, or a degenerate window: pass the pulse through verbatim so
    // the caller (and the emitted hit row) are bit-identical to the
    // un-split larana behaviour.  Do NOT recompute any field here.
    if (!pars["split_enable"].asBool() || pulse.t_end <= pulse.t_start) {
        return {pulse};
    }
    const double frac = pars["split_min_prominence"].asDouble();
    const double absdepth = pars["split_min_prominence_abs"].asDouble();
    const double minpk = pars["split_min_peak"].asDouble();
    const int minsep = pars["split_min_separation"].asInt();

    const int a = pulse.t_start, b = pulse.t_end;
    auto v = [&](int i) { return double(wf[i]) - ped_mean; };

    // Local maxima over the closed window [a, b].  Samples outside the
    // window are treated as -inf so the window edges can themselves be
    // peaks (matching sliding_window, whose t_max may sit at an edge).  A
    // flat-top run of equal maxima collapses to its centre index.
    std::vector<int> peaks;
    for (int i = a; i <= b; ) {
        const double left = (i > a) ? v(i - 1) : -1e300;
        if (v(i) < left) { ++i; continue; }
        int j = i;
        while (j < b && v(j + 1) == v(i)) ++j;  // plateau run [i, j]
        const double right = (j < b) ? v(j + 1) : -1e300;
        if (v(i) >= left && v(i) >= right) {
            const int c = (i + j) / 2;           // plateau centre
            if (v(c) > minpk) peaks.push_back(c);
        }
        i = j + 1;
    }
    if (peaks.size() <= 1) return {pulse};

    // Valley (argmin) on the open interval (p, q).
    auto valley_between = [&](int p, int q) {
        int vmin = p + 1;
        for (int i = p + 1; i < q; ++i) {
            if (v(i) < v(vmin)) vmin = i;
        }
        return vmin;
    };

    // Adjacent-pair prominence test to a fixed point: a shallow shoulder
    // between two peaks is not a real second pulse, so merge the pair
    // (drop the lower peak) and re-evaluate.
    bool changed = true;
    while (changed && peaks.size() > 1) {
        changed = false;
        for (size_t k = 0; k + 1 < peaks.size(); ++k) {
            const int p = peaks[k], q = peaks[k + 1];
            const int val = valley_between(p, q);
            const double smaller = std::min(v(p), v(q));
            const double depth = smaller - v(val);
            // The valley must be deep both relative to the smaller flanking
            // peak (separates comparable pulses) and in absolute terms (so a
            // shallow ripple on the slow scintillation tail is not split off
            // as a spurious hit -- mirrors the prototype's absolute PE margin).
            const bool prominent = depth >= frac * smaller && depth >= absdepth
                                   && (q - p) >= minsep;
            if (!prominent) {
                peaks.erase(peaks.begin() + (v(p) < v(q) ? k : k + 1));
                changed = true;
                break;
            }
        }
    }
    if (peaks.size() <= 1) return {pulse};

    // Cut at each accepted valley.  The valley sample belongs to the left
    // sub-pulse; the next sub-pulse starts at valley+1, so the union of the
    // sub-windows is exactly [t_start, t_end] with no double counting.
    std::vector<Pulse> subs;
    int start = a;
    for (size_t k = 0; k + 1 < peaks.size(); ++k) {
        const int cut = valley_between(peaks[k], peaks[k + 1]);
        subs.push_back(Pulse{start, cut, 0, 0.0, 0.0});
        start = cut + 1;
    }
    subs.push_back(Pulse{start, b, 0, 0.0, 0.0});

    // Recompute fields the way operator() consumes them.  area is the
    // unclamped baseline-subtracted sum over the sub-window (matching the
    // in-pulse accumulation in sliding_window); it differs from the parent
    // area only by the parent's presample clamping, which is intentional
    // (the parent pulse is never emitted when it is split).  t_max is the
    // first index reaching the max, matching sliding_window's strict-`<`
    // tie-break.
    for (auto& s : subs) {
        s.t_max = s.t_start;
        s.peak = v(s.t_start);
        s.area = 0.0;
        for (int i = s.t_start; i <= s.t_end; ++i) {
            const double val = v(i);
            s.area += val;
            if (s.peak < val) { s.peak = val; s.t_max = i; }
        }
    }
    return subs;
}

bool Flash::OpHitFinder::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count);
        return true;
    }
    ++m_count;

    const double t0 = in->time();
    const double tick = in->tick();

    auto traces = Aux::tagged_traces(in, m_intag);
    const size_t ncol = 9;
    std::vector<double> hits;
    for (const auto& trace : traces) {
        const auto& charge = trace->charge();

        // Scale and cast to short, as RunHitFinder_deco does before
        // pulse finding (thresholds and areas are in these units).
        std::vector<short> wf(charge.size());
        for (size_t i = 0; i < charge.size(); ++i) {
            wf[i] = static_cast<short>(m_scale * charge[i]);
        }

        // Pedestal and noise.  Default (m_robust_baseline false): PedAlgoEdges
        // head method (mean/std of the first samples) -- bit-identical to the
        // self-trigger snippet path and every existing config.
        double ped_mean = 0, ped_sigma = 0;
        Configuration algo = m_algo;
        if (m_fixed_ped_sigma > 0) {
            // ROI-cleaned input: baseline is 0 (endpoint-zeroed) and the noise
            // floor is the known clean value (the in-ROI samples are too
            // signal-dominated to estimate from).  See the header.
            ped_mean = 0;
            ped_sigma = m_fixed_ped_sigma;
            algo["nsigma_threshold"] = m_robust_nsigma;
        }
        else if (m_robust_baseline) {
            // Robust per-channel baseline for the CONTINUOUS full stream, where
            // the head method is meaningless: ped_mean = median, ped_sigma =
            // MAD (both over the whole waveform; signal is sparse so the median
            // is the baseline -- removes a per-channel DC offset).  Channels
            // whose MAD is far above the noise floor are ringing data-quality
            // channels and emit no hits.  The start gate is raised to
            // robust_nsigma * MAD (high for noisy, ~unchanged for clean).  See
            // pdhd/docs/pdhd-fullstream-light-reco.md.
            if (wf.empty()) continue;
            std::vector<short> tmp(wf);
            const size_t mid = tmp.size() / 2;
            std::nth_element(tmp.begin(), tmp.begin() + mid, tmp.end());
            ped_mean = tmp[mid];
            std::vector<double> dev(tmp.size());
            for (size_t i = 0; i < tmp.size(); ++i) dev[i] = std::abs(tmp[i] - ped_mean);
            std::nth_element(dev.begin(), dev.begin() + mid, dev.end());
            ped_sigma = 1.4826 * dev[mid];
            if (ped_sigma >= m_robust_veto_sigma) continue;  // veto ringing channel
            algo["nsigma_threshold"] = m_robust_nsigma;
        }
        else {
            const int nped = std::min<int>(m_ped_nsamples, wf.size());
            for (int i = 0; i < nped; ++i) ped_mean += wf[i];
            ped_mean /= nped;
            for (int i = 0; i < nped; ++i) ped_sigma += (wf[i] - ped_mean) * (wf[i] - ped_mean);
            ped_sigma = std::sqrt(ped_sigma / nped);
        }

        for (const auto& pulse : sliding_window(wf, ped_mean, ped_sigma, algo)) {
          // Split each found pulse at prominent sub-peaks (no-op and
          // bit-identical when split_enable is false: one sub == pulse).
          for (const auto& sub : split_pulse(wf, ped_mean, pulse, algo)) {
            if (sub.peak < m_hit_threshold) continue;
            const double peak_time = t0 + tick * (trace->tbin() + sub.t_max);
            const double start_time = t0 + tick * (trace->tbin() + sub.t_start);
            const double width = (sub.t_end - sub.t_start) * tick;
            hits.push_back(trace->channel());
            hits.push_back(peak_time);
            hits.push_back(width);
            hits.push_back(sub.area);
            hits.push_back(sub.peak);
            hits.push_back(sub.area / m_spe_area);
            hits.push_back(start_time);
            hits.push_back(-1);  // flash_id, assigned by OpFlashFinder
            hits.push_back(0);   // fast_to_total
          }
        }
    }

    ITensor::vector* tensors = new ITensor::vector;
    Configuration tmd;
    tmd["name"] = "ophits";
    tensors->push_back(std::make_shared<Aux::SimpleTensor>(
        ITensor::shape_t{hits.size() / ncol, ncol}, hits.data(), tmd));

    Configuration md;
    md["event"] = in->ident();
    md["producer"] = "wct-flash";
    out = std::make_shared<Aux::SimpleTensorSet>(in->ident(), md,
                                                 ITensor::shared_vector(tensors));
    log->debug("frame {}: {} ophits from {} '{}' traces",
               in->ident(), hits.size() / ncol, traces.size(), m_intag);
    return true;
}
