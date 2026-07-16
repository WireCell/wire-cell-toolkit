#include "WireCellFlash/OpDecon.h"

#include "WireCellAux/DftTools.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Units.h"

#include <cmath>
#include <complex>

WIRECELL_FACTORY(OpDecon, WireCell::Flash::OpDecon,
                 WireCell::INamed,
                 WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;

Flash::OpDecon::OpDecon()
    : Aux::Logger("OpDecon", "flash")
{
}

Flash::OpDecon::~OpDecon() {}

WireCell::Configuration Flash::OpDecon::default_configuration() const
{
    Configuration cfg;
    cfg["intag"] = m_intag;
    cfg["outtag"] = m_outtag;
    cfg["spe_file"] = m_spe_file;
    cfg["noise_file"] = m_noise_file;
    cfg["samples"] = m_samples;
    cfg["pre_trigger"] = m_pre_trigger;
    cfg["pedestal_buffer"] = m_pedestal_buffer;
    cfg["line_noise_rms"] = m_line_noise_rms;
    cfg["input_polarity"] = m_input_polarity;
    cfg["auto_scale"] = m_auto_scale;
    cfg["scale"] = m_scale;
    cfg["fixed_snr"] = m_fixed_snr;
    cfg["wiener_inspired"] = m_wiener_inspired;
    cfg["wi_sigma_mhz"] = m_wi_sigma_mhz;
    cfg["wi_power"] = m_wi_power;
    cfg["wi_eps_rel"] = m_wi_eps_rel;
    cfg["apply_postfilter"] = m_apply_postfilter;
    cfg["postfilter_cutoff"] = m_postfilter_cutoff;
    cfg["apply_post_blcorr"] = m_apply_post_blcorr;
    cfg["use_real_dft"] = m_use_real_dft;
    cfg["fold_postfilter"] = m_fold_postfilter;
    cfg["detect_saturation"] = m_detect_saturation;
    cfg["saturation_adc"] = m_saturation_adc;
    cfg["saturation_min_samples"] = m_saturation_min_samples;
    cfg["saturation_pad"] = m_saturation_pad;
    cfg["saturation_repair"] = m_saturation_repair;
    cfg["repair_fit_samples"] = m_repair_fit_samples;
    cfg["overflow_to_rail"] = m_overflow_to_rail;
    cfg["overflow_adc"] = m_overflow_adc;
    cfg["overflow_min_samples"] = m_overflow_min_samples;
    cfg["overflow_min_neighbor"] = m_overflow_min_neighbor;
    cfg["dft"] = "FftwDFT";
    return cfg;
}

void Flash::OpDecon::configure(const WireCell::Configuration& cfg)
{
    m_intag = get(cfg, "intag", m_intag);
    m_outtag = get(cfg, "outtag", m_outtag);
    m_spe_file = get(cfg, "spe_file", m_spe_file);
    m_noise_file = get(cfg, "noise_file", m_noise_file);
    m_samples = get(cfg, "samples", m_samples);
    m_pre_trigger = get(cfg, "pre_trigger", m_pre_trigger);
    m_pedestal_buffer = get(cfg, "pedestal_buffer", m_pedestal_buffer);
    m_line_noise_rms = get(cfg, "line_noise_rms", m_line_noise_rms);
    m_input_polarity = get(cfg, "input_polarity", m_input_polarity);
    m_auto_scale = get(cfg, "auto_scale", m_auto_scale);
    m_scale = get(cfg, "scale", m_scale);
    m_fixed_snr = get(cfg, "fixed_snr", m_fixed_snr);
    m_wiener_inspired = get(cfg, "wiener_inspired", m_wiener_inspired);
    m_wi_sigma_mhz = get(cfg, "wi_sigma_mhz", m_wi_sigma_mhz);
    m_wi_power = get(cfg, "wi_power", m_wi_power);
    m_wi_eps_rel = get(cfg, "wi_eps_rel", m_wi_eps_rel);
    m_apply_postfilter = get(cfg, "apply_postfilter", m_apply_postfilter);
    m_postfilter_cutoff = get(cfg, "postfilter_cutoff", m_postfilter_cutoff);
    m_apply_post_blcorr = get(cfg, "apply_post_blcorr", m_apply_post_blcorr);
    m_use_real_dft = get(cfg, "use_real_dft", m_use_real_dft);
    m_fold_postfilter = get(cfg, "fold_postfilter", m_fold_postfilter);
    m_detect_saturation = get(cfg, "detect_saturation", m_detect_saturation);
    m_saturation_adc = get(cfg, "saturation_adc", m_saturation_adc);
    m_saturation_min_samples = get(cfg, "saturation_min_samples", m_saturation_min_samples);
    m_saturation_pad = get(cfg, "saturation_pad", m_saturation_pad);
    m_saturation_repair = get(cfg, "saturation_repair", m_saturation_repair);
    m_repair_fit_samples = get(cfg, "repair_fit_samples", m_repair_fit_samples);
    m_overflow_to_rail = get(cfg, "overflow_to_rail", m_overflow_to_rail);
    m_overflow_adc = get(cfg, "overflow_adc", m_overflow_adc);
    m_overflow_min_samples = get(cfg, "overflow_min_samples", m_overflow_min_samples);
    m_overflow_min_neighbor = get(cfg, "overflow_min_neighbor", m_overflow_min_neighbor);
    if (m_overflow_to_rail && !m_detect_saturation) {
        log->warn("overflow_to_rail requires detect_saturation; it is off, so the "
                  "remapped runs would never be flagged -- overflow_to_rail disabled");
        m_overflow_to_rail = false;
    }
    if (m_overflow_to_rail) {
        SPDLOG_LOGGER_DEBUG(log, "overflow_to_rail on: adc<={} len>={} neighbor>={} -> {}",
                            m_overflow_adc, m_overflow_min_samples,
                            m_overflow_min_neighbor, m_saturation_adc);
    }

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    // SPE templates + channel map.
    auto jspe = Persist::load(m_spe_file);
    m_templates.clear();
    for (const auto& jt : jspe["templates"]) {
        SPETemplate spe;
        for (const auto& jv : jt["values"]) {
            spe.wave.push_back(jv.asFloat());
        }
        if ((int)spe.wave.size() > m_samples) {
            spe.wave.resize(m_samples);
        }
        // Zero-padding to m_samples does not change the max (amplitude is
        // clamped to >= 1), so it can be computed on the unpadded wave; the
        // padded spectrum itself is built lazily in ensure_fft() on the
        // template's first use.
        spe.amplitude = std::max(1.0f, *std::max_element(spe.wave.begin(), spe.wave.end()));
        // The wiener-inspired PE normalization is 1/H(0) = 1/template-area:
        // a non-positive area makes the DC gain negative and the channel's
        // PE scale meaningless (seen for a template built from noise
        // self-triggers with an unrepaired AC undershoot).  Warn only --
        // behavior unchanged.
        {
            double dc = 0.0;
            for (float v : spe.wave) dc += v;
            if (dc <= 0.0) {
                log->warn("SPE template {} has non-positive area {:.1f}: "
                          "its deconvolved PE scale is invalid (fix the "
                          "template; see pdvd-lightpattern-sp-investigation.md)",
                          m_templates.size(), dc);
            }
        }
        if (m_saturation_repair) {
            // Exponential time constants for the rail repair, fit from the
            // template itself (same fits as the Python study --
            // pdvd/docs/qlmatch/saturation_recovery_study.py Channel).
            const int n = (int)spe.wave.size();
            const int ipk = (int)std::distance(
                spe.wave.begin(), std::max_element(spe.wave.begin(), spe.wave.end()));
            auto fit_logslope = [&](int t0, int t1) -> double {
                // least-squares slope of ln(wave) vs t over [t0, t1), y > 2% amp
                double st = 0, sl = 0, stt = 0, stl = 0;
                int np = 0;
                for (int t = std::max(0, t0); t < std::min(n, t1); ++t) {
                    const double y = spe.wave[t];
                    if (y <= 0.02 * spe.amplitude) continue;
                    const double l = std::log(y);
                    st += t; sl += l; stt += double(t) * t; stl += double(t) * l;
                    ++np;
                }
                if (np < 2) return 0.0;
                const double den = stt - st * st / np;
                return den != 0.0 ? (stl - st * sl / np) / den : 0.0;
            };
            const double s_fall = fit_logslope(ipk + 25, ipk + 200);
            if (s_fall < 0) spe.tau_fall = -1.0 / s_fall;
            const double s_rise = fit_logslope(ipk - 4, ipk);
            if (s_rise > 0) spe.tau_rise = 1.0 / s_rise;
        }
        m_templates.push_back(std::move(spe));
    }
    m_chan2tmpl.clear();
    const auto& jch = jspe["channels"];
    const auto& jti = jspe["template_index"];
    for (Json::ArrayIndex i = 0; i < jch.size(); ++i) {
        m_chan2tmpl[jch[i].asInt()] = jti[i].asUInt();
    }
    log->debug("loaded {} SPE templates, {} mapped channels from {}",
               m_templates.size(), m_chan2tmpl.size(), m_spe_file);

    // Optional per-channel noise power spectra (half-spectrum bins,
    // zero-padded/truncated to samples/2+1 as in the LArSoft module).
    m_noise_templates.clear();
    m_chan2noise.clear();
    if (!m_noise_file.empty()) {
        auto jnoise = Persist::load(m_noise_file);
        for (const auto& jt : jnoise["templates"]) {
            std::vector<double> power;
            for (const auto& jv : jt["values"]) {
                power.push_back(jv.asDouble());
            }
            power.resize(m_samples / 2 + 1, 0.0);
            m_noise_templates.push_back(std::move(power));
        }
        const auto& jnch = jnoise["channels"];
        const auto& jnti = jnoise["template_index"];
        for (Json::ArrayIndex i = 0; i < jnch.size(); ++i) {
            m_chan2noise[jnch[i].asInt()] = jnti[i].asUInt();
        }
        log->debug("loaded {} noise templates, {} mapped channels from {}",
                   m_noise_templates.size(), m_chan2noise.size(), m_noise_file);
    }

    // Gauss post-filter spectrum, built exactly as
    // Deconvolution::BuildExtraFilter: a normalized time-domain
    // Gaussian at the window center, transformed and phase-shifted
    // back.  Sample frequency is 62.5 MHz (16 ns ticks).
    if (m_apply_postfilter) {
        const double sample_freq = 62.5;  // MHz
        const double df = sample_freq / m_samples;
        const double cutoff_bins = m_postfilter_cutoff / df;
        const double sigma = m_samples * std::sqrt(std::log(2.0)) / (2 * M_PI * cutoff_bins);
        const int mu = m_samples / 2;
        std::vector<float> xf(m_samples);
        const double norm = 1.0 / (sigma * std::sqrt(2 * M_PI));
        for (int i = 0; i < m_samples; ++i) {
            const double u = (i - mu) / sigma;
            xf[i] = norm * std::exp(-0.5 * u * u);
        }
        auto F = dft_fwd(xf);
        m_postfilter.resize(m_samples);
        for (int k = 0; k < m_samples; ++k) {
            const double ph = -2 * M_PI * double(k) * mu / m_samples;
            m_postfilter[k] = F[k] * std::complex<float>(std::cos(ph), std::sin(ph));
        }
    }

    // Spectral weights for the fold_postfilter head-pedestal: the mean of
    // the first nped samples of the (unfiltered) inverse transform is
    // mean = Re sum_k X_k w_k with w_k = (1/(N nped)) sum_{i<nped}
    // exp(+2 pi i k i / N) (geometric sum in closed form).
    if (m_apply_postfilter && m_fold_postfilter && m_apply_post_blcorr) {
        const int N = m_samples;
        const int nped = m_pre_trigger - m_pedestal_buffer;
        m_ped_w.assign(N, {0.0, 0.0});
        for (int k = 0; k < N; ++k) {
            if (k == 0) {
                m_ped_w[k] = 1.0 / N;
                continue;
            }
            const double th = 2 * M_PI * double(k) / N;
            const std::complex<double> num = 1.0 - std::polar(1.0, th * nped);
            const std::complex<double> den = 1.0 - std::polar(1.0, th);
            m_ped_w[k] = num / den / double(N) / double(nped);
        }
    }

    // Wiener-inspired band filter F(f) = exp(-0.5 (f/sigma)^power),
    // F(0) = 1 exactly, on the full-spectrum frequency grid.
    if (m_wiener_inspired) {
        const double sample_freq = 62.5;  // MHz
        const double df = sample_freq / m_samples;
        m_wi_filter.resize(m_samples);
        for (int k = 0; k < m_samples; ++k) {
            const double f = df * std::min(k, m_samples - k);
            m_wi_filter[k] = std::exp(-0.5 * std::pow(f / m_wi_sigma_mhz, m_wi_power));
        }
        log->debug("wiener-inspired filter: sigma={} MHz power={} eps_rel={}",
                   m_wi_sigma_mhz, m_wi_power, m_wi_eps_rel);
    }
}

std::vector<std::complex<float>> Flash::OpDecon::dft_fwd(const std::vector<float>& wave) const
{
    return m_use_real_dft ? Aux::DftTools::fwd_r2c_real(m_dft, wave)
                          : Aux::DftTools::fwd_r2c(m_dft, wave);
}

std::vector<float> Flash::OpDecon::dft_inv(const std::vector<std::complex<float>>& spec) const
{
    return m_use_real_dft ? Aux::DftTools::inv_c2r_real(m_dft, spec)
                          : Aux::DftTools::inv_c2r(m_dft, spec);
}

void Flash::OpDecon::ensure_fft(SPETemplate& spe)
{
    if (!spe.fft.empty()) {
        return;
    }
    std::vector<float> padded(spe.wave);
    padded.resize(m_samples, 0);
    spe.fft = dft_fwd(padded);
    if (m_wiener_inspired) {
        double hmax = 0;
        for (const auto& c : spe.fft) hmax = std::max(hmax, std::abs(std::complex<double>(c)));
        spe.wi_eps = std::pow(m_wi_eps_rel * hmax, 2);
    }
    // With a fixed SNR and flat noise the filter -- hence the AutoScale
    // normalization -- is the same for every record of this template.
    if (m_auto_scale && !m_wiener_inspired && m_fixed_snr > 0.0 && m_noise_file.empty()) {
        const int N = m_samples;
        const double N2_flat = m_line_noise_rms * m_line_noise_rms * N;
        const double S2 = m_fixed_snr * N2_flat;
        std::vector<std::complex<float>> xG(N);
        for (int k = 0; k < N; ++k) {
            const std::complex<double> H = spe.fft[k];
            const double H2 = std::norm(H);
            const std::complex<double> g = std::conj(H) * S2 / (H2 * S2 + N2_flat);
            xG[k] = std::complex<float>(g);
        }
        spe.cached_scale = auto_scale(spe, xG);
        spe.scale_cached = true;
    }
}

// AutoScale normalization: apply the filter to the SPE response, center it,
// and integrate the positive region around the peak.
double Flash::OpDecon::auto_scale(const SPETemplate& spe,
                                  const std::vector<std::complex<float>>& xG) const
{
    const int N = m_samples;
    std::vector<std::complex<float>> xGH(N);
    for (int k = 0; k < N; ++k) {
        // phase shift by half window: exp(+i 2 pi k (N/2) / N) = (-1)^k
        const float sign = (k % 2 == 0) ? 1.0f : -1.0f;
        xGH[k] = xG[k] * spe.fft[k] * sign;
    }
    auto x = dft_inv(xGH);
    const int imax = std::distance(x.begin(), std::max_element(x.begin(), x.end()));
    int ileft = imax, iright = imax;
    while (ileft > 0 && x[ileft] > 0) --ileft;
    while (iright < N && x[iright] > 0) ++iright;
    double norm = 0;
    for (int k = ileft; k <= iright && k < N; ++k) norm += x[k];
    if (norm > 1.0 || norm <= 0.0) norm = 1.0;
    return 1.0 / norm;
}

// Fill each railed run [i,j) with the two-sided exponential bridge: the
// falling-edge back-extrapolation anchored on the first repair_fit_samples
// samples after the run, intersected (min) with the rising-edge
// extrapolation anchored on the last 4 samples before it, clamped >= the
// measured (railed) samples.  Anchors use the robust median of the
// per-sample amplitude estimates.  A run with < 2 positive anchor samples
// (or a template whose tau fits failed) is left clipped.
// Mirrors repair_runs of pdvd/docs/qlmatch/saturation_recovery_study.py.
void Flash::OpDecon::repair_runs(std::vector<float>& w,
                                 const std::vector<std::pair<int, int>>& runs,
                                 double pedestal, const SPETemplate& spe) const
{
    if (spe.tau_fall <= 0) return;
    const int n = (int)w.size();
    auto median = [](std::vector<double>& v) {
        std::sort(v.begin(), v.end());
        const size_t m = v.size() / 2;
        return v.size() % 2 ? v[m] : 0.5 * (v[m - 1] + v[m]);
    };
    for (const auto& [i, j] : runs) {
        std::vector<double> ya;
        for (int k = 0; k < m_repair_fit_samples && j + k < n; ++k) {
            const double y = w[j + k] - pedestal;
            if (y > 0) ya.push_back(y * std::exp(k / spe.tau_fall));
        }
        if (ya.size() < 2) continue;
        const double A = median(ya);
        double B = -1;
        if (spe.tau_rise > 0) {
            std::vector<double> yb;
            for (int t = std::max(0, i - 4); t < i; ++t) {
                const double y = w[t] - pedestal;
                if (y > 0) yb.push_back(y * std::exp((i - t) / spe.tau_rise));
            }
            if (!yb.empty()) B = median(yb);
        }
        for (int t = i; t < j; ++t) {
            double fill = A * std::exp(-(t - j) / spe.tau_fall);
            if (B > 0) fill = std::min(fill, B * std::exp((t - i) / spe.tau_rise));
            w[t] = std::max(w[t], (float)(pedestal + fill));
        }
    }
}

// Rewrite floor-pinned OVERFLOW runs to the rail so the ordinary rail scan
// sees them.  A run of >= m_overflow_min_samples samples at <= m_overflow_adc
// is an overflow ONLY if the true signal was above the ceiling across it; the
// evidence for that is the pair of IMMEDIATE neighbours, which must BOTH reach
// m_overflow_min_neighbor (the pulse enters and leaves the run near the
// ceiling).  The same floor pin also occurs in the deep post-pulse undershoot,
// where the signal is below 0 and the neighbours are low -- those must be left
// alone, since raising them would fabricate a rail-height pulse.  Immediate
// neighbours only: widening the window lets a one-sample spike next to an
// undershoot run masquerade as the exit of an overflow.
// A run touching either trace edge has no neighbour pair and is left alone, as
// is a run reaching into the pedestal window (rewriting there would wreck the
// head pedestal that deconvolve() derives).  See the doc cited in OpDecon.h.
int Flash::OpDecon::unclip_overflow(std::vector<float>& w) const
{
    const int n = (int)w.size();
    const int nped = m_pre_trigger - m_pedestal_buffer;
    int nfix = 0;
    int i = 0;
    while (i < n) {
        if (w[i] > m_overflow_adc) { ++i; continue; }
        int j = i;
        while (j < n && w[j] <= m_overflow_adc) ++j;
        if (j - i >= m_overflow_min_samples && i > 0 && j < n && i >= nped) {
            if (w[i - 1] >= m_overflow_min_neighbor && w[j] >= m_overflow_min_neighbor) {
                for (int t = i; t < j; ++t) w[t] = (float) m_saturation_adc;
                ++nfix;
            }
        }
        i = j;
    }
    return nfix;
}

std::vector<float> Flash::OpDecon::deconvolve(const std::vector<float>& adc, const SPETemplate& spe,
                                              const std::vector<double>* noise) const
{
    const int N = m_samples;

    // Pedestal from the first (PreTrigger - PedestalBuffer) samples.
    const int nped = m_pre_trigger - m_pedestal_buffer;
    double pedestal = 0;
    for (int i = 0; i < nped; ++i) pedestal += adc[i];
    pedestal /= nped;

    // Baseline-subtracted, positive-polarity waveform.
    std::vector<float> xv(N, 0);
    const int nin = std::min<int>(N, adc.size());
    for (int i = 0; i < nin; ++i) {
        xv[i] = m_input_polarity * (adc[i] - pedestal);
    }
    // (LArSoft pads a short waveform with random noise; snippets here
    // are always full-length, shorter input is zero-padded.)

    // Expected input: delta of SPE_Max p.e. -> flat |S|^2.
    const double N2_flat = m_line_noise_rms * m_line_noise_rms * N;
    // Signal level S2.  Adaptive (default): from this waveform's own peak,
    // so the filter follows the brightest pulse in the window.  Fixed
    // (m_fixed_snr > 0): S2 = R * N^2 with R the chosen S2/N^2 ratio, giving
    // a filter independent of signal amplitude and record length.
    double S2 = 0;
    if (m_wiener_inspired) {
        // unused: the wiener-inspired filter has no signal/noise model.
    }
    else if (m_fixed_snr > 0.0) {
        S2 = m_fixed_snr * N2_flat;
    }
    else {
        const double spe_max = *std::max_element(xv.begin(), xv.end()) / spe.amplitude;
        S2 = spe_max * spe_max;
    }

    // Wiener filter G = conj(H) S2 / (|H|^2 S2 + N2), full spectrum.
    // N2 per bin from the channel's noise power spectrum (half-spectrum
    // indexed, mirrored onto the upper bins) when available.
    // Wiener-inspired mode: G = conj(H) F / (|H|^2 + eps) instead.
    auto xV = dft_fwd(xv);
    const bool fold = m_apply_postfilter && m_fold_postfilter;
    std::vector<std::complex<float>> xG(N), xY(N);
    std::complex<double> ped_acc(0.0, 0.0);
    for (int k = 0; k < N; ++k) {
        const std::complex<double> H = spe.fft[k];
        const double H2 = std::norm(H);
        std::complex<double> g;
        if (m_wiener_inspired) {
            g = std::conj(H) * m_wi_filter[k] / (H2 + spe.wi_eps);
        }
        else {
            const double N2 = noise ? (*noise)[k <= N / 2 ? k : N - k] : N2_flat;
            g = std::conj(H) * S2 / (H2 * S2 + N2);
        }
        xG[k] = std::complex<float>(g);
        xY[k] = xG[k] * xV[k];
        if (fold) {
            // Head pedestal of the UNfiltered decon, evaluated spectrally.
            if (m_apply_post_blcorr) ped_acc += std::complex<double>(xY[k]) * m_ped_w[k];
            xY[k] *= m_postfilter[k];
        }
    }
    auto xvdec = dft_inv(xY);  // normalized inverse

    // AutoScale normalization (cached per template when the filter is
    // record-independent; see ensure_fft).
    // Wiener-inspired: F(0)=1 already gives unit 1-PE area; keep m_scale.
    double scale = m_scale;
    if (m_auto_scale && !m_wiener_inspired) {
        scale = spe.scale_cached ? spe.cached_scale : auto_scale(spe, xG);
    }

    // Baseline of the deconvolved waveform, computed before the
    // post-filter, subtracted after (as in the LArSoft module).  With
    // fold_postfilter the same pre-filter baseline comes from the
    // spectral accumulation above.
    double dec_pedestal = 0;
    if (m_apply_post_blcorr) {
        if (fold) {
            dec_pedestal = ped_acc.real();
        }
        else {
            for (int i = 0; i < nped; ++i) dec_pedestal += xvdec[i];
            dec_pedestal /= nped;
        }
    }

    if (m_apply_postfilter && !fold) {
        auto Y = dft_fwd(xvdec);
        for (int k = 0; k < N; ++k) Y[k] *= m_postfilter[k];
        xvdec = dft_inv(Y);
    }

    std::vector<float> out(N);
    for (int i = 0; i < N; ++i) {
        out[i] = (xvdec[i] - dec_pedestal) * scale;
    }
    return out;
}

bool Flash::OpDecon::operator()(const IFrame::pointer& in, IFrame::pointer& out)
{
    out = nullptr;
    if (!in) {
        log->debug("EOS");
        return true;
    }

    auto traces = Aux::tagged_traces(in, m_intag);
    ITrace::vector all_traces(in->traces()->begin(), in->traces()->end());
    IFrame::trace_list_t out_idx;
    int nskipped = 0;
    // Saturation flags collected over the snippets (empty unless enabled).
    Waveform::ChannelMaskMap cmm;
    int nsaturated = 0;
    int noverflow = 0;   // floor-pinned overflow runs rewritten to the rail
    for (const auto& trace : traces) {
        const int chan = trace->channel();
        auto it = m_chan2tmpl.find(chan);
        if (it == m_chan2tmpl.end() or it->second >= m_templates.size()) {
            ++nskipped;
            continue;
        }
        // overflow_to_rail: rewrite floor-pinned overflow runs to the rail on a
        // COPY, before the rail scan below, so detect/flag/repair treat them as
        // ordinary saturation.  Off (or nothing to rewrite) => `src` stays the
        // input trace and every path below is bit-identical.
        const std::vector<float>* src = &trace->charge();
        std::vector<float> unclipped;
        if (m_detect_saturation && m_overflow_to_rail) {
            unclipped = trace->charge();
            const int nfix = unclip_overflow(unclipped);
            if (nfix) {
                src = &unclipped;
                noverflow += nfix;
            }
        }
        std::vector<std::pair<int, int>> sat_runs;  // unpadded, trace-local
        if (m_detect_saturation) {
            // Flag each contiguous run of >= saturation_min_samples raw samples
            // at/above the rail as a saturated tick sub-range.  Marking the run
            // (not the whole trace) keeps a long full-stream channel from being
            // vetoed wholesale on one stray sample: the broad over-integrated
            // hit a clipped flat-top produces overlaps the run and is dropped,
            // while real light elsewhere on the trace survives.
            const auto& q = *src;
            const int tb = trace->tbin();
            const int n = (int)q.size();
            int i = 0;
            while (i < n) {
                if (q[i] >= m_saturation_adc) {
                    int j = i;
                    while (j < n && q[j] >= m_saturation_adc) ++j;
                    if (j - i >= m_saturation_min_samples) {
                        const int lo = std::max(0, i - m_saturation_pad);
                        const int hi = std::min(n, j + m_saturation_pad);
                        cmm["saturation"][chan].push_back({tb + lo, tb + hi});
                        sat_runs.emplace_back(i, j);
                        ++nsaturated;
                    }
                    i = j;
                }
                else {
                    ++i;
                }
            }
        }
        const std::vector<double>* noise = nullptr;
        auto nit = m_chan2noise.find(chan);
        if (nit != m_chan2noise.end() and nit->second < m_noise_templates.size()) {
            noise = &m_noise_templates[nit->second];
        }
        ensure_fft(m_templates[it->second]);
        // saturation_repair: deconvolve a repaired COPY; the flagged mask
        // ranges above are emitted unchanged (repair AND flag, not
        // repair-instead-of-flag).  Default off -> original waveform.
        const std::vector<float>* wf = src;
        std::vector<float> repaired;
        if (m_saturation_repair && !sat_runs.empty()) {
            repaired = *src;
            const int nped = m_pre_trigger - m_pedestal_buffer;
            double pedestal = 0;
            for (int k = 0; k < nped; ++k) pedestal += repaired[k];
            pedestal /= nped;
            repair_runs(repaired, sat_runs, pedestal, m_templates[it->second]);
            wf = &repaired;
        }
        auto dec = deconvolve(*wf, m_templates[it->second], noise);
        out_idx.push_back(all_traces.size());
        all_traces.push_back(std::make_shared<Aux::SimpleTrace>(chan, trace->tbin(), std::move(dec)));
    }
    if (nskipped) {
        log->warn("frame {}: skipped {} traces with no SPE template", in->ident(), nskipped);
    }

    // When saturation detection is off, build the frame exactly as before
    // (bit-identical).  When on, forward any incoming masks and add the
    // collected "saturation" ranges.
    Aux::SimpleFrame* sframe;
    if (m_detect_saturation) {
        Waveform::ChannelMaskMap outcmm = in->masks();
        for (const auto& [label, chmasks] : cmm) {
            for (const auto& [chan, ranges] : chmasks) {
                auto& dst = outcmm[label][chan];
                dst.insert(dst.end(), ranges.begin(), ranges.end());
            }
        }
        sframe = new Aux::SimpleFrame(in->ident(), in->time(), all_traces, in->tick(), outcmm);
        log->debug("frame {}: {} saturated tick-runs flagged ({} floor-pinned overflow runs "
                   "remapped to the rail)", in->ident(), nsaturated, noverflow);
    }
    else {
        sframe = new Aux::SimpleFrame(in->ident(), in->time(), all_traces, in->tick());
    }
    for (const auto& tag : in->frame_tags()) {
        sframe->tag_frame(tag);
    }
    for (const auto& tag : in->trace_tags()) {
        sframe->tag_traces(tag, in->tagged_traces(tag), in->trace_summary(tag));
    }
    sframe->tag_traces(m_outtag, out_idx);
    out = IFrame::pointer(sframe);
    log->debug("frame {}: deconvolved {} of {} '{}' traces -> '{}'",
               in->ident(), out_idx.size(), traces.size(), m_intag, m_outtag);
    return true;
}
