// ==============================
// UncorrelatedAddNoise.cxx
// ==============================
/**
 * @file UncorrelatedAddNoise.cxx
 * @brief Add per-channel uncorrelated electronics noise to an IFrame.
 *
 * This is the PSD-only counterpart to CorrelatedAddNoise and is intentionally kept
 * parallel in structure, scaling, and FFT conventions to ease comparison.
 *
 * Key points
 * ----------
 * - Noise is generated in frequency space with independent CSCG coefficients per wire/bin.
 * - Coefficients are scaled using Rayleigh-mean scaling to match a stored mean-magnitude
 *   target spectrum avg_mag(w,f).
 * - DC is forced to zero.
 * - If N is even, Nyquist bin is real-only (required for real-valued time-domain signal).
 * - One-sided spectrum is converted to time domain using a NumPy-compatible irfft
 *   convention implemented via FFTW: bin-doubling + 1/N + optional ifft_scale.
 *
 * @author Avik Ghosh
 * @version 1.0.0
 * @date 2026-02-04
 */

#include "WireCellGen/UncorrelatedAddNoise.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"

#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

#include <json/json.h>
#include <fftw3.h>
#include <bzlib.h>

#include <algorithm>
#include <complex>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <memory>

namespace {

// -----------------------------------------------------------------------------
// Versioning (local to this translation unit)
// -----------------------------------------------------------------------------
static constexpr const char* kComponentVersion = "1.0.0";

// -----------------------------------------------------------------------------
// Statistical constants used in Rayleigh-mean scaling
//
// Convention used by this code:
//   Z = (N(0,1) + i N(0,1)) / sqrt(2)  =>  Z ~ CN(0,1)
//   => E|Z| = sqrt(pi)/2
//
// For a real Gaussian:
//   n ~ N(0,1) => E|n| = sqrt(2/pi)
//
// NOTE on naming: SQRT_PI_OVER_2 == sqrt(pi)/2, *not* sqrt(pi/2).
// -----------------------------------------------------------------------------
static const double SQRT_PI_OVER_2 = 0.5 * std::sqrt(M_PI);     // E|CN(0,1)|
static const double SQRT_2_OVER_PI = std::sqrt(2.0 / M_PI);     // E|N(0,1)|

// -----------------------------------------------------------------------------
// Utility helpers (kept local; mirrored in CorrelatedAddNoise for symmetry)
// -----------------------------------------------------------------------------
inline bool ends_with(const std::string& s, const std::string& suf)
{
    return s.size() >= suf.size() &&
           s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

inline std::string slurp_file_text(const std::string& fname)
{
    std::ifstream in(fname, std::ios::in | std::ios::binary);
    if (!in.good()) {
        THROW(WireCell::IOError()
              << WireCell::errmsg{"UncorrelatedAddNoise: cannot open model_file=" + fname});
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

inline std::string slurp_file_bz2_text(const std::string& fname)
{
    FILE* fp = std::fopen(fname.c_str(), "rb");
    if (!fp) {
        THROW(WireCell::IOError()
              << WireCell::errmsg{"UncorrelatedAddNoise: cannot open model_file=" + fname});
    }

    int bzerr = BZ_OK;
    BZFILE* bz = BZ2_bzReadOpen(&bzerr, fp, 0, 0, nullptr, 0);
    if (bzerr != BZ_OK || !bz) {
        std::fclose(fp);
        THROW(WireCell::RuntimeError()
              << WireCell::errmsg{"UncorrelatedAddNoise: BZ2_bzReadOpen failed for " + fname});
    }

    std::string out;
    out.reserve(1024 * 1024);

    const int CHUNK = 1 << 15;
    std::vector<char> buf((size_t)CHUNK);

    while (true) {
        int nread = BZ2_bzRead(&bzerr, bz, buf.data(), CHUNK);

        if (bzerr == BZ_OK || bzerr == BZ_STREAM_END) {
            if (nread > 0) out.append(buf.data(), (size_t)nread);
            if (bzerr == BZ_STREAM_END) break;
            continue;
        }

        BZ2_bzReadClose(&bzerr, bz);
        std::fclose(fp);
        THROW(WireCell::RuntimeError()
              << WireCell::errmsg{"UncorrelatedAddNoise: BZ2_bzRead failed for " + fname});
    }

    BZ2_bzReadClose(&bzerr, bz);
    std::fclose(fp);
    return out;
}

inline Json::Value parse_json_string(const std::string& text,
                                     const std::string& label_for_errors)
{
    Json::CharReaderBuilder rb;
    rb["collectComments"] = false;

    Json::Value root;
    std::string errs;

    const std::unique_ptr<Json::CharReader> reader(rb.newCharReader());
    const bool ok = reader->parse(text.data(), text.data() + text.size(), &root, &errs);
    if (!ok) {
        THROW(WireCell::ValueError()
              << WireCell::errmsg{"UncorrelatedAddNoise: JSON parse failed for "
                                  + label_for_errors + " : " + errs});
    }
    return root;
}

// -----------------------------------------------------------------------------
// FFTW helper: NumPy-like irfft with one-sided spectrum
//
// Implements a NumPy-compatible inverse real FFT convention:
//
// Input:
//   Xpos : one-sided spectrum of length K = N/2 + 1
// Behavior:
//   - If N even, enforce Nyquist bin imaginary part = 0
//   - Double bins that have conjugate partners:
//       * N even : k = 1 .. K-2
//       * N odd  : k = 1 .. K-1
//   - FFTW c2r is unnormalized => apply 1/N explicitly
//   - Apply final_scale after 1/N (this corresponds to ifft_scale)
// -----------------------------------------------------------------------------
inline void irfft_numpy_like_fftwf(const std::vector<std::complex<float>>& Xpos,
                                  std::vector<float>& x_time,
                                  int N,
                                  float final_scale)
{
    const int Kt = N / 2 + 1;
    if ((int)Xpos.size() != Kt) {
        THROW(WireCell::ValueError()
              << WireCell::errmsg{"irfft_numpy_like_fftwf: Xpos size != N/2+1"});
    }

    x_time.assign((size_t)N, 0.0f);

    fftwf_complex* in  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)Kt);
    float*         out = (float*)fftwf_malloc(sizeof(float) * (size_t)N);

    if (!in || !out) {
        if (in)  fftwf_free(in);
        if (out) fftwf_free(out);
        THROW(WireCell::RuntimeError()
              << WireCell::errmsg{"irfft_numpy_like_fftwf: fftwf_malloc failed"});
    }

    // Copy one-sided spectrum into FFTW input layout.
    for (int k = 0; k < Kt; ++k) {
        in[k][0] = Xpos[(size_t)k].real();
        in[k][1] = Xpos[(size_t)k].imag();
    }

    const bool is_even = ((N % 2) == 0);
    if (is_even) {
        // For real output, the Nyquist bin must be purely real when N is even.
        in[Kt - 1][1] = 0.0f;
    }

    // Bin-doubling to account for omitted negative-frequency half.
    const int k_max = is_even ? (Kt - 2) : (Kt - 1);
    for (int k = 1; k <= k_max; ++k) {
        in[k][0] *= 2.0f;
        in[k][1] *= 2.0f;
    }

    fftwf_plan plan = fftwf_plan_dft_c2r_1d(N, in, out, FFTW_ESTIMATE);
    if (!plan) {
        fftwf_free(in);
        fftwf_free(out);
        THROW(WireCell::RuntimeError()
              << WireCell::errmsg{"irfft_numpy_like_fftwf: fftwf_plan_dft_c2r_1d failed"});
    }

    fftwf_execute(plan);

    // Explicit 1/N to match NumPy irfft conventions.
    const float invN = 1.0f / (float)N;
    for (int t = 0; t < N; ++t) {
        x_time[(size_t)t] = out[t] * invN * final_scale;
    }

    fftwf_destroy_plan(plan);
    fftwf_free(in);
    fftwf_free(out);
}

} // namespace

using namespace WireCell;
using namespace WireCell::Gen;

// Factory registration
WIRECELL_FACTORY(UncorrelatedAddNoise,
                 WireCell::Gen::UncorrelatedAddNoise,
                 WireCell::INamed,
                 WireCell::IFrameFilter,
                 WireCell::IConfigurable)

// -----------------------------------------------------------------------------
// Ctors / dtor
// -----------------------------------------------------------------------------
UncorrelatedAddNoise::UncorrelatedAddNoise()
    : Aux::Logger("UncorrelatedAddNoise")
{
}
UncorrelatedAddNoise::~UncorrelatedAddNoise() = default;

// -----------------------------------------------------------------------------
// INamed
// -----------------------------------------------------------------------------
void UncorrelatedAddNoise::set_name(const std::string& name) { m_name = name; }
std::string UncorrelatedAddNoise::get_name() const { return m_name; }

// -----------------------------------------------------------------------------
// Default configuration
// -----------------------------------------------------------------------------
Configuration UncorrelatedAddNoise::default_configuration() const
{
    Configuration cfg;
    cfg["rng"]        = "Random:default";
    cfg["model_file"] = "uncorrelated_noise_model.json.bz2";
    cfg["nsamples"]   = 2128;
    cfg["dt"]         = 0.5 * units::us;
    cfg["ifft_scale"] = 1.0;
    return cfg;
}

// -----------------------------------------------------------------------------
// Configure: resolve RNG, read config, load model, basic sanity checks
// -----------------------------------------------------------------------------
void UncorrelatedAddNoise::configure(const Configuration& cfg)
{
    // RNG tool
    const auto rng_tn = get<std::string>(cfg, "rng", "Random:default");
    m_rng = Factory::find_tn<IRandom>(rng_tn);
    if (!m_rng) {
        THROW(ValueError() << errmsg{"UncorrelatedAddNoise: failed to get RNG tool: " + rng_tn});
    }

    // Primary config
    m_nsamples   = static_cast<size_t>(get<int>(cfg, "nsamples", 2128));
    m_dt         = get<double>(cfg, "dt", 0.5 * units::us);
    m_ifft_scale = get<double>(cfg, "ifft_scale", 1.0);
    m_model_file = get<std::string>(cfg, "model_file", "uncorrelated_noise_model.json.bz2");

    // Convenience: if dt was provided as a small raw seconds value, convert to WCT units.
    if (m_dt > 0.0 && m_dt < 1.0) {
        const double dt_seconds = m_dt;
        m_dt = dt_seconds * units::second;
        log->info("dt provided as seconds ({}) -> converted to WCT base time units ({} ns)", dt_seconds, m_dt);
    }

    // Load model file
    load_model(m_model_file);

    // Minimal sanity checks
    if (m_avg_mag.cols() == 0) {
        THROW(ValueError() << errmsg{"UncorrelatedAddNoise: avg_mag has zero columns"});
    }
    if (m_avg_mag.cols() != m_freq.size()) {
        THROW(ValueError() << errmsg{"UncorrelatedAddNoise: avg_mag N_freq != freq.size()"});
    }
    if (m_avg_mag.rows() == 0) {
        THROW(ValueError() << errmsg{"UncorrelatedAddNoise: avg_mag has zero rows"});
    }

    log->info("Configured UncorrelatedAddNoise v{}: nsamples={} dt={}ns ifft_scale={} model_file={}",
              kComponentVersion, (int)m_nsamples, m_dt, m_ifft_scale, m_model_file);
}

// -----------------------------------------------------------------------------
// Model loading
//
// Model keys supported:
//   - freq_hz   (preferred) or freq_ghz
//   - avg_mag   [Nw x Nf]
//   - live_mask (optional) [Nw]
//   - meta.stored_units.avg_mag = "MV" or "mV" (optional; otherwise assume WCT base)
// -----------------------------------------------------------------------------
void UncorrelatedAddNoise::load_model(const std::string& fname)
{
    const std::string text =
        ends_with(fname, ".bz2") ? slurp_file_bz2_text(fname) : slurp_file_text(fname);

    const Json::Value root = parse_json_string(text, fname);

    // ---- Frequency axis ----
    Json::Value jfreq;
    double freq_scale = 1.0;

    if (root.isMember("freq_hz")) {
        jfreq = root["freq_hz"];
        freq_scale = units::hertz;
    }
    else if (root.isMember("freq_ghz")) {
        jfreq = root["freq_ghz"];
        freq_scale = (1.0e9 * units::hertz);
    }
    else {
        THROW(ValueError() << errmsg{"UncorrelatedAddNoise: missing freq_hz or freq_ghz"});
    }

    const int Nf = (int)jfreq.size();
    if (Nf <= 0) {
        THROW(ValueError() << errmsg{"UncorrelatedAddNoise: frequency array has zero length"});
    }

    m_freq.resize(Nf);
    for (int i = 0; i < Nf; ++i) {
        m_freq(i) = jfreq[(Json::ArrayIndex)i].asDouble() * freq_scale;
    }

    // ---- avg_mag units ----
    // We store m_avg_mag in WCT voltage units; interpret meta if present.
    double mag_scale = 1.0;
    if (root.isMember("meta")) {
        const auto& meta = root["meta"];

        // Optional cross-checks (kept as in your working code)
        if (meta.isMember("dt_ns")) {
            const double dt_model = meta["dt_ns"].asDouble();
            if (std::fabs(dt_model - m_dt) > 1e-6) {
                log->warn("Model meta.dt_ns={} differs from configured dt={} (both expected WCT base ns)",
                          dt_model, m_dt);
            }
        }
        if (meta.isMember("n_ticks")) {
            const size_t nt_model = (size_t)meta["n_ticks"].asLargestUInt();
            if (nt_model != m_nsamples) {
                log->warn("Model meta.n_ticks={} differs from configured nsamples={}", nt_model, m_nsamples);
            }
        }

        if (meta.isMember("stored_units") && meta["stored_units"].isMember("avg_mag")) {
            const std::string u = meta["stored_units"]["avg_mag"].asString();
            if (u == "MV")      mag_scale = units::megavolt;
            else if (u == "mV") mag_scale = units::millivolt;
            else {
                log->warn("Unknown meta.stored_units.avg_mag='{}' -> assuming WCT base units", u);
                mag_scale = 1.0;
            }
        }
    }

    // ---- avg_mag ----
    if (!root.isMember("avg_mag")) {
        THROW(ValueError() << errmsg{"UncorrelatedAddNoise: missing avg_mag"});
    }

    const auto& jam = root["avg_mag"];
    const int Nw = (int)jam.size();
    if (Nw <= 0) {
        THROW(ValueError() << errmsg{"UncorrelatedAddNoise: avg_mag has zero rows"});
    }

    m_avg_mag.resize(Nw, Nf);
    for (int w = 0; w < Nw; ++w) {
        const auto& row = jam[(Json::ArrayIndex)w];
        if ((int)row.size() != Nf) {
            THROW(ValueError() << errmsg{
                "UncorrelatedAddNoise: avg_mag row size mismatch (row="
                + std::to_string((int)row.size()) + " != N_freq=" + std::to_string(Nf) + ")"
            });
        }
        for (int k = 0; k < Nf; ++k) {
            m_avg_mag(w, k) = row[(Json::ArrayIndex)k].asDouble() * mag_scale;
        }
    }

    // ---- live_mask (optional) ----
    m_live_mask.resize(Nw);
    m_live_mask.setOnes();

    if (root.isMember("live_mask")) {
        const auto& jlive = root["live_mask"];
        if ((int)jlive.size() == Nw) {
            for (int w = 0; w < Nw; ++w) {
                m_live_mask(w) = jlive[(Json::ArrayIndex)w].asInt();
            }
        }
        else {
            log->warn("live_mask size {} != N_wires {}, ignoring", (int)jlive.size(), Nw);
        }
    }

    log->info("Loaded uncorrelated noise model from {} with N_wires={}, N_freq={}",
              fname, (int)m_avg_mag.rows(), (int)m_freq.size());
}

// -----------------------------------------------------------------------------
// Frame utility: sort traces by channel id
// -----------------------------------------------------------------------------
ITrace::vector UncorrelatedAddNoise::sorted_traces_by_channel(const IFrame::pointer& frame) const
{
    ITrace::vector out;
    auto tvp = frame->traces();
    if (tvp && !tvp->empty()) {
        out = *tvp;
        std::sort(out.begin(), out.end(),
                  [](const ITrace::pointer& a, const ITrace::pointer& b) {
                      return a->channel() < b->channel();
                  });
    }
    return out;
}

// -----------------------------------------------------------------------------
// Helper: interpolate target mean magnitude avg_mag(w,f) along model frequency axis
// -----------------------------------------------------------------------------
double UncorrelatedAddNoise::interp_avg_mag_w(int w, double f) const
{
    const int Nf = (int)m_freq.size();
    if (Nf <= 1) return m_avg_mag(w, 0);

    // Clamp to endpoints
    if (f <= m_freq(0))    return m_avg_mag(w, 0);
    if (f >= m_freq(Nf-1)) return m_avg_mag(w, Nf-1);

    // Find first bin >= f
    const double* begin = m_freq.data();
    const double* end   = m_freq.data() + Nf;
    const double* it = std::lower_bound(begin, end, f);

    int i1 = (int)(it - begin);
    if (i1 <= 0)  return m_avg_mag(w, 0);
    if (i1 >= Nf) return m_avg_mag(w, Nf-1);

    const int i0 = i1 - 1;

    // Linear interpolation
    const double f0 = m_freq(i0);
    const double f1 = m_freq(i1);
    const double a  = (f1 > f0) ? ((f - f0) / (f1 - f0)) : 0.0;

    return (1.0 - a) * m_avg_mag(w, i0) + a * m_avg_mag(w, i1);
}

// -----------------------------------------------------------------------------
// Core synthesis routine: build per-wire time-domain noise
//
// Output: noise_time[w, t] in WCT voltage units (consistent with avg_mag units).
//
// Notes:
// - Uses the configured dt (WCT ns) to convert bin index -> physical frequency.
// - DC bin is forced to zero (consistent with pedestal-subtracted inputs).
// - Nyquist bin is handled only if N is even (real-only).
// -----------------------------------------------------------------------------
Eigen::MatrixXd UncorrelatedAddNoise::make_uncorrelated_noise(int nwires_frame,
                                                              size_t nsamp) const
{
    const int N_wires_model = (int)m_avg_mag.rows();

    // If the frame has a different wire count than the model, keep the behavior:
    // warn and truncate to the smaller of the two.
    if (nwires_frame != N_wires_model) {
        log->warn("UncorrelatedAddNoise: nwires_frame={} != model N_wires={}, truncating/minor mismatch",
                  nwires_frame, N_wires_model);
        nwires_frame = std::min(nwires_frame, N_wires_model);
    }

    const int    N   = (int)nsamp;      // time samples
    const int    Kt  = N/2 + 1;         // one-sided bins
    const double dfT = 1.0 / (double(N) * m_dt); // bin width in Hz (WCT units)

    const bool is_even = ((N % 2) == 0);
    const int  nyq_k   = is_even ? (Kt - 1) : -1;

    // Build one-sided spectrum per wire
    Eigen::MatrixXcd Xpos = Eigen::MatrixXcd::Zero(nwires_frame, Kt);

    // DC forced to 0
    if (Kt > 0) Xpos.col(0).setZero();

    // ---- Complex bins (exclude DC and Nyquist) ----
    for (int kt = 1; kt < Kt; ++kt) {
        if (is_even && kt == nyq_k) continue;  // Nyquist handled separately below

        const double f = kt * dfT;

        for (int w = 0; w < nwires_frame; ++w) {
            if (m_live_mask(w) == 0) continue;

            const double target = interp_avg_mag_w(w, f);

            // Draw Z ~ CN(0,1)
            const double xr = m_rng->normal(0.0, 1.0);
            const double xi = m_rng->normal(0.0, 1.0);
            const std::complex<double> Z =
                std::complex<double>(xr, xi) / std::sqrt(2.0);

            // Rayleigh-mean scaling: E|Z| = sqrt(pi)/2
            const double scale = (SQRT_PI_OVER_2 > 0.0) ? (target / SQRT_PI_OVER_2) : 0.0;
            Xpos(w, kt) = Z * scale;
        }
    }

    // ---- Nyquist bin (real-only) for even N ----
    if (is_even && nyq_k >= 0) {
        const double fnyq = nyq_k * dfT;

        for (int w = 0; w < nwires_frame; ++w) {
            if (m_live_mask(w) == 0) continue;

            const double target = interp_avg_mag_w(w, fnyq);

            // Draw n ~ N(0,1) (real)
            const double n = m_rng->normal(0.0, 1.0);

            // Scale using E|N(0,1)| = sqrt(2/pi)
            const double scale = (SQRT_2_OVER_PI > 0.0) ? (target / SQRT_2_OVER_PI) : 0.0;
            Xpos(w, nyq_k) = std::complex<double>(n * scale, 0.0);
        }
    }

    // ---- IFFT per wire using NumPy-like irfft convention ----
    Eigen::MatrixXd noise_time(nwires_frame, N);

    std::vector<std::complex<float>> one_sided((size_t)Kt);
    std::vector<float> x_time;

    for (int w = 0; w < nwires_frame; ++w) {
        // Copy Eigen one-sided spectrum into std::vector<complex<float>>
        for (int k = 0; k < Kt; ++k) {
            const auto xd = Xpos(w, k);
            one_sided[(size_t)k] = std::complex<float>((float)xd.real(), (float)xd.imag());
        }

        // Enforce Nyquist imag=0 (safety; helper also enforces for even N)
        if (is_even) {
            one_sided[(size_t)(Kt - 1)] =
                std::complex<float>(one_sided[(size_t)(Kt - 1)].real(), 0.0f);
        }

        irfft_numpy_like_fftwf(one_sided, x_time, N, (float)m_ifft_scale);

        // Store into Eigen matrix
        for (int t = 0; t < N; ++t) {
            noise_time(w, t) = (double)x_time[(size_t)t];
        }
    }

    return noise_time;
}

// -----------------------------------------------------------------------------
// Main filter operator(): read frame, generate noise, add to traces, write new frame
// -----------------------------------------------------------------------------
bool UncorrelatedAddNoise::operator()(const input_pointer& inframe,
                                      output_pointer& outframe)
{
    // End-of-stream handling
    if (!inframe) {
        outframe = nullptr;
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    // IMPORTANT: trace ordering assumption
    // We sort traces by channel id, and then assume:
    //   model row i corresponds to sorted_traces[i].
    auto traces = sorted_traces_by_channel(inframe);
    const size_t nwires_frame = traces.size();

    if (nwires_frame == 0) {
        outframe = inframe;
        ++m_count;
        return true;
    }

    // Use the trace length of the first trace as the per-frame sample count.
    const size_t ncharge = traces[0]->charge().size();
    if (ncharge != m_nsamples) {
        log->warn("UncorrelatedAddNoise: frame nsamples={} differs from configured nsamples={}, adapting per-frame",
                  ncharge, m_nsamples);
    }

    // Generate a noise matrix (time-domain) matching the frame's length.
    Eigen::MatrixXd noise = make_uncorrelated_noise((int)nwires_frame, ncharge);

    // Build new traces with noise added.
    ITrace::vector outtraces;
    outtraces.reserve(nwires_frame);

    for (size_t i = 0; i < nwires_frame; ++i) {
        const auto& intrace = traces[i];
        const int chid = intrace->channel();
        auto charge = intrace->charge();

        // Add noise row i to this trace (if available)
        if ((int)i < noise.rows()) {
            for (size_t t = 0; t < charge.size(); ++t) {
                charge[t] += noise((int)i, (int)t);
            }
        }

        outtraces.push_back(std::make_shared<Aux::SimpleTrace>(chid, intrace->tbin(), charge));
    }

    // Output a new frame with modified traces. Keep ident/time/tick consistent.
    outframe = std::make_shared<Aux::SimpleFrame>(inframe->ident(),
                                                  inframe->time(),
                                                  outtraces,
                                                  inframe->tick());

    log->debug("UncorrelatedAddNoise v{}: call={} frame={} {} traces",
               kComponentVersion, m_count, inframe->ident(), outtraces.size());

    ++m_count;
    return true;
}