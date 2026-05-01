#include "WireCellSigProc/L1SPFilterPD.h"

#include "WireCellAux/DftTools.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"

#include "WireCellIface/IFieldResponse.h"

#include "WireCellUtil/LassoModel.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Eigen.h"
#include "WireCellUtil/Response.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/cnpy.h"

#include <filesystem>
#include <numeric>
#include <cmath>

#ifdef __clang__
#  if defined(__has_warning)
#    define HAS_WARNING(warning) __has_warning(warning)
#  else
#    define HAS_WARNING(warning) 1
#  endif
#else
#  define HAS_WARNING(warning) 1
#endif

WIRECELL_FACTORY(L1SPFilterPD,
                 WireCell::SigProc::L1SPFilterPD,
                 WireCell::INamed, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace Eigen;
using namespace WireCell;
using namespace WireCell::SigProc;

using WireCell::Aux::DftTools::fwd_r2c;
using WireCell::Aux::DftTools::inv_c2r;

// ── anonymous helpers ─────────────────────────────────────────────────────────
namespace {

// Build N×2N response matrix for one segment.
// Column j           ← basis0(dt + overall)
// Column N+j         ← basis1(dt + overall + basis1_offset)
// Response window: dt ∈ (t_lo, t_hi) relative to overall_time_offset.
// Tick size assumed 0.5 µs (2 MHz ADC).
MatrixXd build_G(int nbin,
                 double t_lo, double t_hi,
                 double overall_time_offset,
                 double basis1_offset,
                 double scaling, double resp_scale,
                 linterp<double>* basis0,
                 linterp<double>* basis1)
{
    MatrixXd G = MatrixXd::Zero(nbin, nbin * 2);
    for (int i = 0; i < nbin; i++) {
        double t_meas = i * 0.5 * units::us;
        for (int j = 0; j < nbin; j++) {
            double t_sig = j * 0.5 * units::us;
            double dt = t_meas - t_sig;
            if (dt > t_lo && dt < t_hi) {
                double dt_adj = dt + overall_time_offset;
#if HAS_WARNING("-Wstringop-overread")
#pragma GCC diagnostic push
#pragma GCC diagnostic warning "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
                G(i, j)        = (*basis0)(dt_adj)                 * scaling * resp_scale;
                G(i, nbin + j) = (*basis1)(dt_adj + basis1_offset) * scaling * resp_scale;
#if HAS_WARNING("-Wstringop-overread")
#pragma GCC diagnostic pop
#endif
            }
        }
    }
    return G;
}

// Run LASSO on one segment.  Returns beta of length 2*nbin.
VectorXd lasso_solve(const MatrixXd& G, const VectorXd& W,
                     double lambda, int niter, double eps)
{
    WireCell::LassoModel m(lambda, niter, eps);
    m.SetData(G, W);
    m.Fit();
    return m.Getbeta();
}

// Per-ROI asymmetry statistics used by the trigger and the calibration dump.
struct AsymRecord {
    int    nbin_fit{0};
    double temp_sum{0}, temp1_sum{0}, temp2_sum{0};
    double max_val{-1e30}, min_val{1e30};
};

// Compute ADC asymmetry quantities for ticks [start_tick, end_tick).
// adc/sig charges are indexed with the same tbin offset (both traces assumed
// to share the same tbin, which is the case for all pdhd/pdvd frames).
AsymRecord compute_asym(const WireCell::ITrace::ChargeSequence& adc,
                        const WireCell::ITrace::ChargeSequence& sig,
                        int tbin,
                        int start_tick, int end_tick,
                        double threshold)
{
    AsymRecord r;
    r.nbin_fit = end_tick - start_tick;
    for (int i = 0; i < r.nbin_fit; i++) {
        int idx = i + start_tick - tbin;
        double w = adc.at(idx);
        double b = sig.at(idx);
        if (w > r.max_val) r.max_val = w;
        if (w < r.min_val) r.min_val = w;
        if (std::fabs(w) > threshold) {
            r.temp_sum  += w;
            r.temp1_sum += std::fabs(w);
            r.temp2_sum += std::fabs(b);
        }
    }
    return r;
}

}  // namespace

// ── L1SPFilterPD ─────────────────────────────────────────────────────────────

L1SPFilterPD::L1SPFilterPD(double gain, double shaping, double postgain, double ADC_mV,
                            double fine_time_offset, double coarse_time_offset)
  : Aux::Logger("L1SPFilterPD", "sigproc")
  , m_gain(gain)
  , m_shaping(shaping)
  , m_postgain(postgain)
  , m_ADC_mV(ADC_mV)
  , m_fine_time_offset(fine_time_offset)
  , m_coarse_time_offset(coarse_time_offset)
{
}

WireCell::Configuration L1SPFilterPD::default_configuration() const
{
    Configuration cfg;

    // Name of component providing the bipolar induction-plane field response.
    cfg["fields"] = "FieldResponse";
    cfg["bipolar_plane"] = 1;   // plane index within that response (0=U, 1=V, 2=W)

    // Optional: field responses for the unipolar bases.  Empty = not configured;
    // the component acts as a pass-through until these are provided.
    cfg["fields_pos_unipolar"] = "";   // arriving-only (electrons collected on induction)
    cfg["fields_neg_unipolar"] = "";   // leaving-only  (anode-induction)
    cfg["unipolar_plane"] = 1;

    cfg["filter"] = Json::arrayValue;

    cfg["adctag"] = "raw";
    cfg["sigtag"] = "gauss";
    cfg["outtag"] = "l1sp";

    cfg["raw_ROI_th_nsigma"] = 4;
    cfg["raw_ROI_th_adclimit"] = 10;
    cfg["overall_time_offset"] = 0;
    // Time offset of the unipolar basis relative to the bipolar basis (µs).
    // Analogous to collect_time_offset in L1SPFilter.  Tune once the
    // unipolar field-response files are available.
    cfg["unipolar_time_offset"] = 3.0;

    cfg["roi_pad"] = 3;
    cfg["raw_pad"] = 15;

    cfg["adc_l1_threshold"] = 6;
    cfg["adc_sum_threshold"] = 160;
    cfg["adc_sum_rescaling"] = 90.;
    cfg["adc_sum_rescaling_limit"] = 50.;
    cfg["adc_ratio_threshold"] = 0.2;

    cfg["l1_seg_length"] = 120;
    cfg["l1_scaling_factor"] = 500;
    cfg["l1_lambda"] = 5;
    cfg["l1_epsilon"] = 0.05;
    cfg["l1_niteration"] = 100000;
    cfg["l1_decon_limit"] = 100;
    cfg["l1_resp_scale"] = 0.5;
    cfg["l1_basis0_scale"] = 1.15;  // weight for bipolar component
    cfg["l1_basis1_scale"] = 0.5;   // weight for unipolar component

    cfg["peak_threshold"] = 1000;
    cfg["mean_threshold"] = 500;

    cfg["gain"] = m_gain;
    cfg["shaping"] = m_shaping;
    cfg["postgain"] = m_postgain;
    cfg["ADC_mV"] = m_ADC_mV;
    cfg["fine_time_offset"] = m_fine_time_offset;
    cfg["coarse_time_offset"] = m_coarse_time_offset;

    cfg["dft"] = "FftwDFT";

    // Plane-scope filter: IAnodePlane typename + list of plane indices to process.
    // Leave "anode" empty to disable plane filtering (all channels processed).
    cfg["anode"] = "";
    cfg["process_planes"][0] = 0;   // U
    cfg["process_planes"][1] = 1;   // V

    // Calibration dump mode: write per-ROI asymmetry records to NPZ files.
    cfg["dump_mode"] = false;
    cfg["dump_path"] = "";    // directory; one NPZ per operator() call written here
    cfg["dump_tag"] = "";     // label baked into filename (e.g. "apa1")

    return cfg;
}

void L1SPFilterPD::configure(const WireCell::Configuration& cfg)
{
    m_gain = get(cfg, "gain", m_gain);
    m_shaping = get(cfg, "shaping", m_shaping);
    m_postgain = get(cfg, "postgain", m_postgain);
    m_ADC_mV = get(cfg, "ADC_mV", m_ADC_mV);
    m_fine_time_offset = get(cfg, "fine_time_offset", m_fine_time_offset);
    m_coarse_time_offset = get(cfg, "coarse_time_offset", m_coarse_time_offset);

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    m_adctag = get<std::string>(cfg, "adctag", "raw");
    m_sigtag = get<std::string>(cfg, "sigtag", "gauss");
    m_outtag = get<std::string>(cfg, "outtag", "l1sp");

    m_roi_pad = get(cfg, "roi_pad", m_roi_pad);
    m_raw_pad = get(cfg, "raw_pad", m_raw_pad);
    m_raw_ROI_th_nsigma = get(cfg, "raw_ROI_th_nsigma", m_raw_ROI_th_nsigma);
    m_raw_ROI_th_adclimit = get(cfg, "raw_ROI_th_adclimit", m_raw_ROI_th_adclimit);

    // fixme: the use of units here is broken (same issue as in L1SPFilter)
    m_overall_time_offset = get(cfg, "overall_time_offset", 0.0) * units::us;
    m_unipolar_time_offset = get(cfg, "unipolar_time_offset", 3.0) * units::us;

    m_adc_l1_threshold = get(cfg, "adc_l1_threshold", m_adc_l1_threshold);
    m_adc_sum_threshold = get(cfg, "adc_sum_threshold", m_adc_sum_threshold);
    m_adc_sum_rescaling = get(cfg, "adc_sum_rescaling", m_adc_sum_rescaling);
    m_adc_sum_rescaling_limit = get(cfg, "adc_sum_rescaling_limit", m_adc_sum_rescaling_limit);
    m_adc_ratio_threshold = get(cfg, "adc_ratio_threshold", m_adc_ratio_threshold);

    m_l1_seg_length = get(cfg, "l1_seg_length", m_l1_seg_length);
    m_l1_scaling_factor = get(cfg, "l1_scaling_factor", m_l1_scaling_factor);
    m_l1_lambda = get(cfg, "l1_lambda", m_l1_lambda);
    m_l1_epsilon = get(cfg, "l1_epsilon", m_l1_epsilon);
    m_l1_niteration = get(cfg, "l1_niteration", m_l1_niteration);
    m_l1_decon_limit = get(cfg, "l1_decon_limit", m_l1_decon_limit);
    m_l1_resp_scale = get(cfg, "l1_resp_scale", m_l1_resp_scale);
    m_l1_basis0_scale = get(cfg, "l1_basis0_scale", m_l1_basis0_scale);
    m_l1_basis1_scale = get(cfg, "l1_basis1_scale", m_l1_basis1_scale);
    m_peak_threshold = get(cfg, "peak_threshold", m_peak_threshold);
    m_mean_threshold = get(cfg, "mean_threshold", m_mean_threshold);

    m_smearing_vec = get<std::vector<double>>(cfg, "filter");

    m_bipolar_plane = get(cfg, "bipolar_plane", m_bipolar_plane);
    m_unipolar_plane = get(cfg, "unipolar_plane", m_unipolar_plane);

    m_cfg_fields     = get<std::string>(cfg, "fields", "FieldResponse");
    m_cfg_fields_pos = get<std::string>(cfg, "fields_pos_unipolar", "");
    m_cfg_fields_neg = get<std::string>(cfg, "fields_neg_unipolar", "");

    // Plane-scope filter
    m_cfg_anode = get<std::string>(cfg, "anode", "");
    if (!m_cfg_anode.empty()) {
        m_anode = Factory::find_tn<IAnodePlane>(m_cfg_anode);
    }
    if (cfg.isMember("process_planes") && cfg["process_planes"].isArray()) {
        m_process_planes.clear();
        for (auto const& v : cfg["process_planes"]) {
            m_process_planes.push_back(v.asInt());
        }
    }

    // Calibration dump mode
    m_dump_mode = get(cfg, "dump_mode", m_dump_mode);
    m_dump_path = get<std::string>(cfg, "dump_path", m_dump_path);
    m_dump_tag  = get<std::string>(cfg, "dump_tag", m_dump_tag);

    // Reset interpolators so init_resp() rebuilds them on next operator() call.
    m_lin_bipolar.reset();
    m_lin_pos_unipolar.reset();
    m_lin_neg_unipolar.reset();
}

std::unique_ptr<linterp<double>> L1SPFilterPD::build_response(const std::string& fields_tn,
                                                               int plane_index)
{
    auto ifr = Factory::find_tn<IFieldResponse>(fields_tn);
    Response::Schema::FieldResponse fr = ifr->field_response();
    Response::Schema::FieldResponse fravg = Response::average_1D(fr);

    WireCell::Binning tbins(Response::as_array(fravg.planes[0]).cols(), 0,
                            Response::as_array(fravg.planes[0]).cols() * fravg.period);
    Response::ColdElec ce(m_gain, m_shaping);
    auto ewave = ce.generate(tbins);
    Waveform::scale(ewave, m_postgain * m_ADC_mV * (-1));
    auto elec = fwd_r2c(m_dft, ewave);

    std::complex<float> fine_period(fravg.period, 0);

    WireCell::Waveform::realseq_t resp = fravg.planes[plane_index].paths[0].current;
    auto spectrum = fwd_r2c(m_dft, resp);
    Waveform::scale(spectrum, elec);
    Waveform::scale(spectrum, fine_period);
    resp = inv_c2r(m_dft, spectrum);

    double intrinsic_time_offset = fravg.origin / fravg.speed;
    double x0   = (-intrinsic_time_offset - m_coarse_time_offset + m_fine_time_offset);
    double xstep = fravg.period;

    return std::make_unique<linterp<double>>(resp.begin(), resp.end(), x0, xstep);
}

bool L1SPFilterPD::channel_in_scope(int channel) const
{
    if (m_process_planes.empty() || !m_anode) return true;
    int plane = m_anode->resolve(channel).index();
    for (int p : m_process_planes) {
        if (p == plane) return true;
    }
    return false;
}

void L1SPFilterPD::init_resp()
{
    if (!m_lin_bipolar) {
        m_lin_bipolar = build_response(m_cfg_fields, m_bipolar_plane);
    }
    if (!m_lin_pos_unipolar && !m_cfg_fields_pos.empty()) {
        m_lin_pos_unipolar = build_response(m_cfg_fields_pos, m_unipolar_plane);
    }
    if (!m_lin_neg_unipolar && !m_cfg_fields_neg.empty()) {
        m_lin_neg_unipolar = build_response(m_cfg_fields_neg, m_unipolar_plane);
    }
}

int L1SPFilterPD::l1_fit(std::shared_ptr<Aux::SimpleTrace>& newtrace,
                          const std::shared_ptr<const WireCell::ITrace>& adctrace,
                          int start_tick, int end_tick,
                          int hint_polarity)
{
    const int nbin_fit = end_tick - start_tick;

    VectorXd init_W    = VectorXd::Zero(nbin_fit);
    VectorXd init_beta = VectorXd::Zero(nbin_fit);

    double temp_sum = 0, temp1_sum = 0, temp2_sum = 0;
    double max_val = -1, min_val = 1;
    for (int i = 0; i < nbin_fit; i++) {
        init_W(i)    = adctrace->charge().at(i + start_tick - newtrace->tbin());
        init_beta(i) = newtrace->charge().at(i + start_tick - newtrace->tbin());

        if (init_W(i) > max_val) max_val = init_W(i);
        if (init_W(i) < min_val) min_val = init_W(i);

        if (std::fabs(init_W(i)) > m_adc_l1_threshold) {
            temp_sum  += init_W(i);
            temp1_sum += std::fabs(init_W(i));
            temp2_sum += std::fabs(init_beta(i));
        }
    }

    // These are computed above for the trigger; suppress unused-variable
    // warnings while the stub below is in place.
    (void)temp_sum; (void)temp1_sum; (void)temp2_sum;
    (void)max_val;  (void)min_val;

    // ── STUB: trigger classification ──────────────────────────────────────────
    // Always returns 0 (pass-through) until trigger cuts are determined.
    //
    // When implementing Strategy B (see sigproc/docs/l1sp/README.md):
    //
    //   Polarity asymmetry ratio = temp_sum / (temp1_sum * m_adc_sum_rescaling / nbin_fit)
    //   Guard against temp1_sum == 0 before dividing.
    //
    //   Positive polarity (electrons collected on induction, +unipolar):
    //     if (temp1_sum > m_adc_sum_threshold && temp1_sum > 0 &&
    //         temp_sum / (temp1_sum * m_adc_sum_rescaling / nbin_fit) > m_adc_ratio_threshold)
    //         flag_l1 = 1;
    //
    //   Negative polarity (anode-induction, -unipolar):
    //     if (temp1_sum > m_adc_sum_threshold && temp1_sum > 0 &&
    //         -temp_sum / (temp1_sum * m_adc_sum_rescaling / nbin_fit) > m_adc_ratio_threshold)
    //         flag_l1 = -1;
    //
    //   Artifact (flag=2) — same logic as L1SPFilter:
    //     else if (temp1_sum * m_adc_sum_rescaling / nbin_fit < m_adc_sum_rescaling_limit)
    //         flag_l1 = 2;
    //     else if (temp2_sum > 30*nbin_fit && temp1_sum < 2.0*nbin_fit && max_val-min_val < 22)
    //         flag_l1 = 2;
    //
    //   Pass-through (0): otherwise.
    // ──────────────────────────────────────────────────────────────────────────
    int flag_l1 = 0;

    // When propagation hints a polarity, override the stub.
    if (hint_polarity != 0) flag_l1 = hint_polarity;

    // Select the appropriate unipolar basis according to polarity.
    linterp<double>* lin_unipolar = (flag_l1 > 0) ? m_lin_pos_unipolar.get()
                                                   : m_lin_neg_unipolar.get();

    if ((flag_l1 == 1 || flag_l1 == -1) && lin_unipolar == nullptr) {
        // Unipolar basis not yet configured — fall back to pass-through.
        log->warn("l1_fit: polarity {} triggered but unipolar response not configured; "
                  "falling back to pass-through", flag_l1);
        flag_l1 = 0;
    }

    if (flag_l1 == 1 || flag_l1 == -1) {
        // ── LASSO fit ──────────────────────────────────────────────────────
        int n_section = std::max(1, (int)std::round(nbin_fit / m_l1_seg_length));
        std::vector<int> bounds;
        for (int i = 0; i < n_section; i++) bounds.push_back(int(i * nbin_fit / n_section));
        bounds.push_back(nbin_fit);

        VectorXd final_beta = VectorXd::Zero(nbin_fit * 2);

        for (int s = 0; s < n_section; s++) {
            int sn = bounds[s + 1] - bounds[s];
            VectorXd W_seg = VectorXd::Zero(sn);
            for (int j = 0; j < sn; j++) W_seg(j) = init_W(bounds[s] + j);

            MatrixXd G = build_G(sn,
                                 -15 * units::us - m_overall_time_offset,
                                  10 * units::us - m_overall_time_offset,
                                 m_overall_time_offset,
                                 m_unipolar_time_offset,
                                 m_l1_scaling_factor, m_l1_resp_scale,
                                 m_lin_bipolar.get(), lin_unipolar);

            VectorXd beta = lasso_solve(G, W_seg, m_l1_lambda,
                                        (int)m_l1_niteration, m_l1_epsilon);
            for (int j = 0; j < sn; j++) {
                final_beta(bounds[s] + j)           = beta(j);
                final_beta(nbin_fit + bounds[s] + j) = beta(sn + j);
            }
        }

        // Check that there is a non-trivial total reconstructed charge.
        double sum_beta = 0;
        for (int i = 0; i < nbin_fit * 2; i++) sum_beta += final_beta(i);

        if (sum_beta > m_adc_l1_threshold) {
            // Combine basis components and apply smearing filter.
            Waveform::realseq_t l1_signal(nbin_fit, 0);
            for (int j = 0; j < nbin_fit; j++) {
                l1_signal[j] = final_beta(j)           * m_l1_basis0_scale
                             + final_beta(nbin_fit + j) * m_l1_basis1_scale;
            }

            Waveform::realseq_t l2_signal(nbin_fit, 0);
            int mid_bin = ((int)m_smearing_vec.size() - 1) / 2;
            for (int j = 0; j < nbin_fit; j++) {
                if (l1_signal[j] > 0) {
                    for (int k = 0; k < (int)m_smearing_vec.size(); k++) {
                        int bin = j + k - mid_bin;
                        if (bin >= 0 && bin < nbin_fit)
                            l2_signal[bin] += l1_signal[j] * m_smearing_vec[k];
                    }
                }
            }

            for (int j = 0; j < nbin_fit; j++) {
                l1_signal[j] = (l2_signal[j] < m_l1_decon_limit / m_l1_scaling_factor)
                               ? 0.0
                               : l2_signal[j] * m_l1_scaling_factor;
            }

            // Remove small isolated peaks.
            std::vector<std::pair<int, int>> peak_rois;
            {
                bool in_roi = false;
                int start_bin = -1, end_bin = -1;
                for (int j = 0; j < nbin_fit; j++) {
                    if (l1_signal[j] > 0) {
                        if (!in_roi) { start_bin = end_bin = j; in_roi = true; }
                        else         { end_bin = j; }
                    } else if (in_roi) {
                        peak_rois.push_back({start_bin, end_bin});
                        in_roi = false;
                    }
                }
                if (in_roi) peak_rois.push_back({start_bin, end_bin});
            }
            for (auto& roi : peak_rois) {
                double mx = -1, mean_v = 0;
                for (int k = roi.first; k <= roi.second; k++) {
                    if (l1_signal[k] > mx) mx = l1_signal[k];
                    mean_v += l1_signal[k];
                }
                mean_v /= (roi.second - roi.first + 1);
                if (mx < m_peak_threshold && mean_v < m_mean_threshold) {
                    for (int k = roi.first; k <= roi.second; k++) l1_signal[k] = 0;
                }
            }

            for (int t = start_tick; t < end_tick; t++)
                newtrace->charge().at(t - newtrace->tbin()) = l1_signal[t - start_tick];
        }
    }
    else if (flag_l1 == 2) {
        for (int t = start_tick; t < end_tick; t++)
            newtrace->charge().at(t - newtrace->tbin()) = 0;
    }

    return flag_l1;
}

bool L1SPFilterPD::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS at call={}", m_count++);
        return true;
    }

    init_resp();

    auto adctraces = Aux::tagged_traces(in, m_adctag);
    auto sigtraces = Aux::tagged_traces(in, m_sigtag);

    if (adctraces.empty() || sigtraces.empty() || adctraces.size() != sigtraces.size()) {
        log->error("unexpected input: {} ADC traces, {} signal traces at call={}",
                   adctraces.size(), sigtraces.size(), m_count++);
        raise<RuntimeError>("L1SPFilterPD: unexpected input");
    }

    // Pre-compute total tick extent before iterating over ADC traces.
    int ntot_ticks = 0;
    for (auto trace : adctraces) {
        int n = (int)trace->charge().size();
        if (n > ntot_ticks) ntot_ticks = n;
    }

    // Collect ticks with positive decon signal per channel.
    std::map<int, std::set<int>> init_map;
    for (auto trace : sigtraces) {
        int ch = trace->channel();
        int tbin = trace->tbin();
        auto const& charges = trace->charge();
        std::set<int>& ticks = init_map[ch];
        for (int qi = 0; qi < (int)charges.size(); qi++) {
            if (charges[qi] > 0) ticks.insert(tbin + qi);
        }
    }

    // Augment with ticks above the raw-ADC noise threshold.
    std::map<int, std::shared_ptr<const WireCell::ITrace>> adctrace_ch_map;
    for (auto trace : adctraces) {
        int ch = trace->channel();
        adctrace_ch_map[ch] = trace;
        int tbin = trace->tbin();
        auto const& charges = trace->charge();
        const int ntbins = (int)charges.size();
        std::set<int>& ticks = init_map[ch];

        Waveform::realseq_t tmp(charges);
        double mean       = Waveform::percentile(tmp, 0.5);
        double mean_p1sig = Waveform::percentile(tmp, 0.5 + 0.34);
        double mean_n1sig = Waveform::percentile(tmp, 0.5 - 0.34);
        double cut = m_raw_ROI_th_nsigma *
                     std::sqrt((std::pow(mean_p1sig - mean, 2) + std::pow(mean_n1sig - mean, 2)) / 2.);
        if (cut < m_raw_ROI_th_adclimit) cut = m_raw_ROI_th_adclimit;

        for (int qi = 0; qi < ntbins; qi++) {
            if (std::fabs(charges[qi]) > cut) {
                for (int qii = -m_raw_pad; qii <= m_raw_pad; qii++) {
                    int t = tbin + qi + qii;
                    if (t >= 0 && t < ntot_ticks) ticks.insert(t);
                }
            }
        }
    }

    // Build merged, padded ROIs per channel.
    std::map<int, std::vector<std::pair<int, int>>> map_ch_rois;
    for (auto& [wire_index, tick_set] : init_map) {
        if (tick_set.empty()) continue;
        std::vector<int> ts(tick_set.begin(), tick_set.end());

        std::vector<std::pair<int, int>> rois;
        rois.push_back({ts.front(), ts.front()});
        for (size_t i = 1; i < ts.size(); i++) {
            if (ts[i] - rois.back().second <= m_roi_pad * 2)
                rois.back().second = ts[i];
            else
                rois.push_back({ts[i], ts[i]});
        }

        for (auto& r : rois) {
            r.first  = std::max(0, r.first  - m_roi_pad);
            r.second = std::min(ntot_ticks - 1, r.second + m_roi_pad);
        }

        // Merge overlapping after padding.
        std::vector<std::pair<int, int>> merged;
        for (auto& r : rois) {
            if (merged.empty() || r.first > merged.back().second)
                merged.push_back(r);
            else
                merged.back().second = std::max(merged.back().second, r.second);
        }

        map_ch_rois[wire_index] = merged;
    }

    // Main per-channel processing: classify ROIs and run L1 fits (or dump).
    std::map<int, std::vector<int>> map_ch_flag_rois;
    ITrace::vector out_traces;

    // Calibration dump: per-ROI parallel vectors accumulated across all channels.
    std::vector<int32_t> d_channel, d_roi_start, d_roi_end, d_nbin_fit;
    std::vector<double>  d_temp_sum, d_temp1_sum, d_temp2_sum, d_max_val, d_min_val;
    std::vector<int32_t> d_prev_roi_end, d_next_roi_start, d_prev_gap, d_next_gap;

    for (auto trace : sigtraces) {
        auto newtrace = std::make_shared<Aux::SimpleTrace>(
            trace->channel(), trace->tbin(), trace->charge());

        int ch = trace->channel();

        // Skip channels with no ROIs or outside the configured plane scope.
        if (map_ch_rois.find(ch) == map_ch_rois.end() || !channel_in_scope(ch)) {
            out_traces.push_back(newtrace);
            continue;
        }

        auto& rois_save = map_ch_rois[ch];
        auto& adctrace  = adctrace_ch_map[ch];

        if (m_dump_mode) {
            // Calibration / dump path: compute asymmetry stats per ROI, record
            // them, then pass the trace through unchanged (no LASSO, no zeroing).
            for (size_t i = 0; i < rois_save.size(); i++) {
                auto rec = compute_asym(adctrace->charge(), trace->charge(),
                                        trace->tbin(),
                                        rois_save[i].first, rois_save[i].second + 1,
                                        m_adc_l1_threshold);
                d_channel.push_back(ch);
                d_roi_start.push_back(rois_save[i].first);
                d_roi_end.push_back(rois_save[i].second);
                d_nbin_fit.push_back(rec.nbin_fit);
                d_temp_sum.push_back(rec.temp_sum);
                d_temp1_sum.push_back(rec.temp1_sum);
                d_temp2_sum.push_back(rec.temp2_sum);
                d_max_val.push_back(rec.max_val);
                d_min_val.push_back(rec.min_val);

                int32_t prev_end   = (i > 0) ? rois_save[i - 1].second : -1;
                int32_t next_start = (i + 1 < rois_save.size())
                                     ? rois_save[i + 1].first : -1;
                d_prev_roi_end.push_back(prev_end);
                d_next_roi_start.push_back(next_start);
                d_prev_gap.push_back(prev_end >= 0
                                     ? (int32_t)(rois_save[i].first - prev_end) : -1);
                d_next_gap.push_back(next_start >= 0
                                     ? (int32_t)(next_start - rois_save[i].second) : -1);
            }
        } else {
            // Normal processing path: classify ROIs and run L1 fits.
            std::vector<int> flag_rois;
            flag_rois.reserve(rois_save.size());

            // First pass: classify and fit each ROI.
            for (auto& roi : rois_save) {
                int flag = l1_fit(newtrace, adctrace, roi.first, roi.second + 1);
                flag_rois.push_back(flag);
                // Zero any negative decon values within the ROI (same as L1SPFilter).
                for (int t = roi.first; t <= roi.second; t++) {
                    auto& v = newtrace->charge().at(t - trace->tbin());
                    if (v < 0) v = 0;
                }
            }
            map_ch_flag_rois[ch] = flag_rois;

            // Forward propagation: if a flag-1 or flag-(-1) ROI is within 20 ticks
            // of a flag-2 ROI, re-run the flag-2 ROI with the neighbour's polarity.
            {
                int active_polarity = 0;
                int prev_end = -2000;
                for (size_t i = 0; i < rois_save.size(); i++) {
                    if (rois_save[i].first - prev_end > 20) active_polarity = 0;
                    int f = flag_rois[i];
                    if (f == 1 || f == -1) {
                        active_polarity = f;
                    } else if (f == 2 && active_polarity != 0) {
                        l1_fit(newtrace, adctrace, rois_save[i].first, rois_save[i].second + 1,
                               active_polarity);
                        flag_rois[i] = 0;
                    } else if (f == 0) {
                        active_polarity = 0;
                    }
                    prev_end = rois_save[i].second;
                }
            }

            // Reverse propagation.
            if (!rois_save.empty()) {
                int active_polarity = 0;
                int prev_start = rois_save.back().second + 2000;
                for (size_t i = 0; i < rois_save.size(); i++) {
                    size_t ri = rois_save.size() - 1 - i;
                    if (prev_start - rois_save[ri].second > 20) active_polarity = 0;
                    int f = flag_rois[ri];
                    if (f == 1 || f == -1) {
                        active_polarity = f;
                    } else if (f == 2 && active_polarity != 0) {
                        l1_fit(newtrace, adctrace, rois_save[ri].first, rois_save[ri].second + 1,
                               active_polarity);
                        flag_rois[ri] = 0;
                    } else if (f == 0) {
                        active_polarity = 0;
                    }
                    prev_start = rois_save[ri].first;
                }
            }
        }

        out_traces.push_back(newtrace);
    }

    // Calibration dump: write NPZ for this frame.
    if (m_dump_mode && !m_dump_path.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(m_dump_path, ec);
        const std::string fname =
            fmt::format("{}/{}_{:04d}_{}.npz", m_dump_path, m_dump_tag, m_count, in->ident());
        std::filesystem::remove(fname, ec);

        auto save_i32 = [&](const std::string& key, const std::vector<int32_t>& v) {
            if (v.empty()) { int32_t d = 0; cnpy::npz_save(fname, key, &d, {0}, "a"); }
            else cnpy::npz_save(fname, key, v.data(), {v.size()}, "a");
        };
        auto save_f64 = [&](const std::string& key, const std::vector<double>& v) {
            if (v.empty()) { double d = 0; cnpy::npz_save(fname, key, &d, {0}, "a"); }
            else cnpy::npz_save(fname, key, v.data(), {v.size()}, "a");
        };
        auto save_i32s = [&](const std::string& key, int32_t val) {
            cnpy::npz_save(fname, key, &val, {1}, "a");
        };
        auto save_f64s = [&](const std::string& key, double val) {
            cnpy::npz_save(fname, key, &val, {1}, "a");
        };

        save_i32s("frame_ident",  (int32_t)in->ident());
        save_f64s("frame_time",   in->time());
        save_i32s("call_count",   (int32_t)m_count);
        save_i32s("n_rois",       (int32_t)d_channel.size());
        save_i32("channel",       d_channel);
        save_i32("roi_start",     d_roi_start);
        save_i32("roi_end",       d_roi_end);
        save_i32("nbin_fit",      d_nbin_fit);
        save_f64("temp_sum",      d_temp_sum);
        save_f64("temp1_sum",     d_temp1_sum);
        save_f64("temp2_sum",     d_temp2_sum);
        save_f64("max_val",       d_max_val);
        save_f64("min_val",       d_min_val);
        save_i32("prev_roi_end",  d_prev_roi_end);
        save_i32("next_roi_start", d_next_roi_start);
        save_i32("prev_gap",      d_prev_gap);
        save_i32("next_gap",      d_next_gap);

        log->debug("call={} dump_mode: {} ROIs -> {}", m_count, d_channel.size(), fname);
    }

    // Layer 4 cross-channel cleaning is intentionally omitted.
    // The uBooNE shorted-wire rationale (paired channels) does not apply
    // to PDHD/PDVD unipolar-induction geometry.

    auto sf = new Aux::SimpleFrame(in->ident(), in->time(), out_traces, in->tick());
    if (!m_outtag.empty()) {
        IFrame::trace_list_t tl(out_traces.size());
        std::iota(tl.begin(), tl.end(), 0);
        sf->tag_traces(m_outtag, tl);
    }
    out = IFrame::pointer(sf);

    log->debug("call={} adctag={} sigtag={} outtag={}", m_count, m_adctag, m_sigtag, m_outtag);
    log->debug("call={} in frame: {}", m_count, Aux::taginfo(in));
    log->debug("call={} out frame: {}", m_count, Aux::taginfo(out));
    ++m_count;

    return true;
}
