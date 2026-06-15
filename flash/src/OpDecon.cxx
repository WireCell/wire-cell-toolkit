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
    cfg["apply_postfilter"] = m_apply_postfilter;
    cfg["postfilter_cutoff"] = m_postfilter_cutoff;
    cfg["apply_post_blcorr"] = m_apply_post_blcorr;
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
    m_apply_postfilter = get(cfg, "apply_postfilter", m_apply_postfilter);
    m_postfilter_cutoff = get(cfg, "postfilter_cutoff", m_postfilter_cutoff);
    m_apply_post_blcorr = get(cfg, "apply_post_blcorr", m_apply_post_blcorr);

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
        spe.wave.resize(m_samples, 0);
        spe.fft = Aux::DftTools::fwd_r2c(m_dft, spe.wave);
        spe.amplitude = std::max(1.0f, *std::max_element(spe.wave.begin(), spe.wave.end()));
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
        auto F = Aux::DftTools::fwd_r2c(m_dft, xf);
        m_postfilter.resize(m_samples);
        for (int k = 0; k < m_samples; ++k) {
            const double ph = -2 * M_PI * double(k) * mu / m_samples;
            m_postfilter[k] = F[k] * std::complex<float>(std::cos(ph), std::sin(ph));
        }
    }
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
    double S2;
    if (m_fixed_snr > 0.0) {
        S2 = m_fixed_snr * N2_flat;
    }
    else {
        const double spe_max = *std::max_element(xv.begin(), xv.end()) / spe.amplitude;
        S2 = spe_max * spe_max;
    }

    // Wiener filter G = conj(H) S2 / (|H|^2 S2 + N2), full spectrum.
    // N2 per bin from the channel's noise power spectrum (half-spectrum
    // indexed, mirrored onto the upper bins) when available.
    auto xV = Aux::DftTools::fwd_r2c(m_dft, xv);
    std::vector<std::complex<float>> xG(N), xY(N);
    for (int k = 0; k < N; ++k) {
        const std::complex<double> H = spe.fft[k];
        const double H2 = std::norm(H);
        const double N2 = noise ? (*noise)[k <= N / 2 ? k : N - k] : N2_flat;
        const std::complex<double> g = std::conj(H) * S2 / (H2 * S2 + N2);
        xG[k] = std::complex<float>(g);
        xY[k] = xG[k] * xV[k];
    }
    auto xvdec = Aux::DftTools::inv_c2r(m_dft, xY);  // normalized inverse

    // AutoScale normalization: apply the filter to the SPE response,
    // center it, and integrate the positive region around the peak.
    double scale = m_scale;
    if (m_auto_scale) {
        std::vector<std::complex<float>> xGH(N);
        for (int k = 0; k < N; ++k) {
            // phase shift by half window: exp(+i 2 pi k (N/2) / N) = (-1)^k
            const float sign = (k % 2 == 0) ? 1.0f : -1.0f;
            xGH[k] = xG[k] * spe.fft[k] * sign;
        }
        auto x = Aux::DftTools::inv_c2r(m_dft, xGH);
        const int imax = std::distance(x.begin(), std::max_element(x.begin(), x.end()));
        int ileft = imax, iright = imax;
        while (ileft > 0 && x[ileft] > 0) --ileft;
        while (iright < N && x[iright] > 0) ++iright;
        double norm = 0;
        for (int k = ileft; k <= iright && k < N; ++k) norm += x[k];
        if (norm > 1.0 || norm <= 0.0) norm = 1.0;
        scale = 1.0 / norm;
    }

    // Baseline of the deconvolved waveform, computed before the
    // post-filter, subtracted after (as in the LArSoft module).
    double dec_pedestal = 0;
    if (m_apply_post_blcorr) {
        for (int i = 0; i < nped; ++i) dec_pedestal += xvdec[i];
        dec_pedestal /= nped;
    }

    if (m_apply_postfilter) {
        auto Y = Aux::DftTools::fwd_r2c(m_dft, xvdec);
        for (int k = 0; k < N; ++k) Y[k] *= m_postfilter[k];
        xvdec = Aux::DftTools::inv_c2r(m_dft, Y);
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
    for (const auto& trace : traces) {
        const int chan = trace->channel();
        auto it = m_chan2tmpl.find(chan);
        if (it == m_chan2tmpl.end() or it->second >= m_templates.size()) {
            ++nskipped;
            continue;
        }
        const std::vector<double>* noise = nullptr;
        auto nit = m_chan2noise.find(chan);
        if (nit != m_chan2noise.end() and nit->second < m_noise_templates.size()) {
            noise = &m_noise_templates[nit->second];
        }
        auto dec = deconvolve(trace->charge(), m_templates[it->second], noise);
        out_idx.push_back(all_traces.size());
        all_traces.push_back(std::make_shared<Aux::SimpleTrace>(chan, trace->tbin(), dec));
    }
    if (nskipped) {
        log->warn("frame {}: skipped {} traces with no SPE template", in->ident(), nskipped);
    }

    auto sframe = new Aux::SimpleFrame(in->ident(), in->time(), all_traces, in->tick());
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
