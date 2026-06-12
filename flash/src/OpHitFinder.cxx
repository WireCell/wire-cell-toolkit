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

        // PedAlgoEdges, head method: mean/std of the first samples.
        double ped_mean = 0, ped_sigma = 0;
        const int nped = std::min<int>(m_ped_nsamples, wf.size());
        for (int i = 0; i < nped; ++i) ped_mean += wf[i];
        ped_mean /= nped;
        for (int i = 0; i < nped; ++i) ped_sigma += (wf[i] - ped_mean) * (wf[i] - ped_mean);
        ped_sigma = std::sqrt(ped_sigma / nped);

        for (const auto& pulse : sliding_window(wf, ped_mean, ped_sigma, m_algo)) {
            if (pulse.peak < m_hit_threshold) continue;
            const double peak_time = t0 + tick * (trace->tbin() + pulse.t_max);
            const double start_time = t0 + tick * (trace->tbin() + pulse.t_start);
            const double width = (pulse.t_end - pulse.t_start) * tick;
            hits.push_back(trace->channel());
            hits.push_back(peak_time);
            hits.push_back(width);
            hits.push_back(pulse.area);
            hits.push_back(pulse.peak);
            hits.push_back(pulse.area / m_spe_area);
            hits.push_back(start_time);
            hits.push_back(-1);  // flash_id, assigned by OpFlashFinder
            hits.push_back(0);   // fast_to_total
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
