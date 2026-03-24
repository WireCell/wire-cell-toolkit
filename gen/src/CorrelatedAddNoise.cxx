// ==============================
// CorrelatedAddNoise.cxx
// ==============================
/**
 * @file CorrelatedAddNoise.cxx
 * @brief Add frequency-banded, inter-channel correlated electronics noise to an IFrame.
 *
 * Key points
 * ----------
 * - Per frequency bin, generate a correlated complex Gaussian vector across wires:
 *     U = A_b Z, where Z ~ CN(0,I) and A_b is the band "colorer".
 * - Per wire, scale U(w) so that the mean magnitude matches avg_mag(w,f).
 * - DC is forced to 0.
 * - If N is even, Nyquist bin is handled as real-only and correlated via A_b on a
 *   real normal vector, with scaling based on E|N(0,1)|.
 * - One-sided spectrum -> time domain via NumPy-like irfft convention:
 *     FFTW c2r + 1/N + optional ifft_scale.
 *
 * @author Avik Ghosh
 * @version 1.0.0
 * @date 2026-02-04
 */

#include "WireCellGen/CorrelatedAddNoise.h"

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
// Statistical constants used in magnitude scaling
//
// For CSCG Z~CN(0,1): E|Z| = sqrt(pi)/2
// For real n~N(0,1): E|n| = sqrt(2/pi)
//
// NOTE on naming: SQRT_PI_OVER_2 == sqrt(pi)/2, *not* sqrt(pi/2).
// -----------------------------------------------------------------------------
static constexpr double SQRT_PI_OVER_2 = 0.5 * std::sqrt(M_PI);     // E|CN(0,1)|
static constexpr double SQRT_2_OVER_PI = std::sqrt(2.0 / M_PI);     // E|N(0,1)|

// -----------------------------------------------------------------------------
// Utility helpers (kept local; mirrored in UncorrelatedAddNoise for symmetry)
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
              << WireCell::errmsg{"CorrelatedAddNoise: cannot open model_file=" + fname});
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
              << WireCell::errmsg{"CorrelatedAddNoise: cannot open model_file=" + fname});
    }

    int bzerr = BZ_OK;
    BZFILE* bz = BZ2_bzReadOpen(&bzerr, fp, 0, 0, nullptr, 0);
    if (bzerr != BZ_OK || !bz) {
        std::fclose(fp);
        THROW(WireCell::RuntimeError()
              << WireCell::errmsg{"CorrelatedAddNoise: BZ2_bzReadOpen failed for " + fname});
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
              << WireCell::errmsg{"CorrelatedAddNoise: BZ2_bzRead failed for " + fname});
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
              << WireCell::errmsg{"CorrelatedAddNoise: JSON parse failed for "
                                  + label_for_errors + " : " + errs});
    }
    return root;
}

// -----------------------------------------------------------------------------
// FFTW helper: NumPy-like irfft with one-sided spectrum
//
// Same convention as in UncorrelatedAddNoise.
// -----------------------------------------------------------------------------
inline void irfft_fftwf(const std::vector<std::complex<float>>& Xpos,
                        std::vector<float>& x_time,
                        int N,
                        float ifft_scale)
{
    const int Kt = N / 2 + 1;
    if ((int)Xpos.size() != Kt) {
        THROW(WireCell::ValueError()
              << WireCell::errmsg{"irfft_fftwf: Xpos size != N/2+1"});
    }

    x_time.resize((size_t)N);

    fftwf_complex* in  = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * (size_t)Kt);
    float*         out = (float*)fftwf_malloc(sizeof(float) * (size_t)N);

    if (!in || !out) {
        if (in)  fftwf_free(in);
        if (out) fftwf_free(out);
        THROW(WireCell::RuntimeError()
              << WireCell::errmsg{"irfft_fftwf: fftwf_malloc failed"});
    }

    for (int k = 0; k < Kt; ++k) {
        in[k][0] = Xpos[(size_t)k].real();
        in[k][1] = Xpos[(size_t)k].imag();
    }

    // Even-N Nyquist bin must be purely real
    if ((N % 2) == 0) {
        in[Kt - 1][1] = 0.0f;
    }

    fftwf_plan plan = fftwf_plan_dft_c2r_1d(N, in, out, FFTW_ESTIMATE);
    if (!plan) {
        fftwf_free(in);
        fftwf_free(out);
        THROW(WireCell::RuntimeError()
              << WireCell::errmsg{"irfft_fftwf: fftwf_plan_dft_c2r_1d failed"});
    }

    fftwf_execute(plan);

    // FFTW backward transform is unnormalized: apply 1/N to match NumPy irfft
    const float norm = ifft_scale / (float)N;
    for (int t = 0; t < N; ++t) {
        x_time[(size_t)t] = out[t] * norm;
    }

    fftwf_destroy_plan(plan);
    fftwf_free(in);
    fftwf_free(out);
}

} // namespace

using namespace WireCell;
using namespace WireCell::Gen;

// Factory registration
WIRECELL_FACTORY(CorrelatedAddNoise,
                 WireCell::Gen::CorrelatedAddNoise,
                 WireCell::INamed,
                 WireCell::IFrameFilter,
                 WireCell::IConfigurable)

// -----------------------------------------------------------------------------
// Ctors / dtor
// -----------------------------------------------------------------------------
CorrelatedAddNoise::CorrelatedAddNoise()
    : Aux::Logger("CorrelatedAddNoise")
{
}
CorrelatedAddNoise::~CorrelatedAddNoise() = default;

// -----------------------------------------------------------------------------
// INamed
// -----------------------------------------------------------------------------
void CorrelatedAddNoise::set_name(const std::string& name) { m_name = name; }
std::string CorrelatedAddNoise::get_name() const { return m_name; }

// -----------------------------------------------------------------------------
// Default configuration
// -----------------------------------------------------------------------------
Configuration CorrelatedAddNoise::default_configuration() const
{
    Configuration cfg;
    cfg["rng"]        = "Random:default";
    cfg["model_file"] = "correlated_noise_model.json.bz2";
    cfg["nsamples"]   = 2128;
    cfg["dt"]         = 0.5 * units::us;
    cfg["ifft_scale"] = 1.0;
    return cfg;
}

// -----------------------------------------------------------------------------
// Configure: resolve RNG, read config, load model, basic sanity checks
// -----------------------------------------------------------------------------
void CorrelatedAddNoise::configure(const Configuration& cfg)
{
    // RNG tool
    const auto rng_tn = get<std::string>(cfg, "rng", "Random:default");
    m_rng = Factory::find_tn<IRandom>(rng_tn);
    if (!m_rng) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: failed to get RNG tool: " + rng_tn});
    }

    // Primary config
    m_nsamples   = static_cast<size_t>(get<int>(cfg, "nsamples", 2128));
    m_dt         = get<double>(cfg, "dt", 0.5 * units::us);
    m_ifft_scale = get<double>(cfg, "ifft_scale", 1.0);
    m_model_file = get<std::string>(cfg, "model_file", "correlated_noise_model.json.bz2");

    // Convenience: allow small raw seconds inputs (kept as-is from your working code)
    if (m_dt > 0.0 && m_dt < 1.0) {
        const double dt_seconds = m_dt;
        m_dt = dt_seconds * units::second;
        log->info("dt provided as seconds ({}) -> converted to WCT base time units ({} ns)", dt_seconds, m_dt);
    }

    // Load model file
    load_model(m_model_file);

    // Minimal sanity checks
    if (m_avg_mag.cols() == 0) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: avg_mag has zero columns"});
    }
    if (m_avg_mag.cols() != m_freq.size()) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: avg_mag N_freq != freq.size()"});
    }
    if (m_avg_mag.rows() == 0) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: avg_mag has zero rows"});
    }
    if (m_A_band.empty()) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: model contains zero correlation bands (A_band empty)"});
    }

    log->info("Configured CorrelatedAddNoise v{}: nsamples={} dt={}ns ifft_scale={} model_file={}",
              kComponentVersion, (int)m_nsamples, m_dt, m_ifft_scale, m_model_file);
}

// -----------------------------------------------------------------------------
// Model loading
//
// Model keys required for correlated mode:
//   - freq_hz   or freq_ghz
//   - avg_mag   [Nw x Nf]
//   - bands_khz [Nb x 2]
//   - A_band    [Nb x Nw x Nw]
//
// Optional:
//   - live_mask [Nw]
//   - meta.stored_units.avg_mag = "MV" or "mV"
// -----------------------------------------------------------------------------
void CorrelatedAddNoise::load_model(const std::string& fname)
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
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: missing freq_hz or freq_ghz"});
    }

    const int Nf = (int)jfreq.size();
    if (Nf <= 0) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: frequency array has zero length"});
    }

    m_freq.resize(Nf);
    for (int i = 0; i < Nf; ++i) {
        m_freq(i) = jfreq[(Json::ArrayIndex)i].asDouble() * freq_scale;
    }

    // ---- Bands ----
    if (!root.isMember("bands_khz")) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: missing bands_khz"});
    }
    const auto& jbands = root["bands_khz"];
    const int nbands = (int)jbands.size();
    if (nbands <= 0) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: bands_khz is empty"});
    }

    // ---- avg_mag units ----
    double mag_scale = 1.0;
    if (root.isMember("meta")) {
        const auto& meta = root["meta"];

        // Optional cross-checks 
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
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: missing avg_mag"});
    }
    const auto& jam = root["avg_mag"];
    const int Nw = (int)jam.size();
    if (Nw <= 0) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: avg_mag has zero rows"});
    }

    m_avg_mag.resize(Nw, Nf);
    for (int w = 0; w < Nw; ++w) {
        const auto& row = jam[(Json::ArrayIndex)w];
        if ((int)row.size() != Nf) {
            THROW(ValueError() << errmsg{
                "CorrelatedAddNoise: avg_mag row size mismatch (row="
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

    // ---- A_band (required) ----
    if (!root.isMember("A_band")) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: missing A_band"});
    }
    const auto& jAb = root["A_band"];
    if ((int)jAb.size() != nbands) {
        THROW(ValueError() << errmsg{"CorrelatedAddNoise: A_band length != bands_khz length"});
    }

    m_A_band.clear();
    m_A_band.reserve((size_t)nbands);

    for (int b = 0; b < nbands; ++b) {
        const auto& jmat = jAb[(Json::ArrayIndex)b];
        if ((int)jmat.size() != Nw) {
            THROW(ValueError() << errmsg{"CorrelatedAddNoise: A_band[b] wrong N_wires"});
        }

        Eigen::MatrixXd A(Nw, Nw);
        for (int i = 0; i < Nw; ++i) {
            const auto& row = jmat[(Json::ArrayIndex)i];
            if ((int)row.size() != Nw) {
                THROW(ValueError() << errmsg{"CorrelatedAddNoise: A_band[b] row wrong N_wires"});
            }
            for (int j = 0; j < Nw; ++j) {
                A(i, j) = row[(Json::ArrayIndex)j].asDouble();
            }
        }
        m_A_band.push_back(std::move(A));
    }

    // ---- Derived band info: convert band edges to WCT Hz and compute diagC ----
    //
    // diagC(w) = sum_j A(w,j)^2  = (A A^T)_{ww}  (row-wise energy)
    // This is used to compute the expected mean magnitude of U(w) = (A Z)(w):
    //   Var[U(w)] = diagC(w)  (for Z ~ CN(0,I))
    //   E|U(w)| = (sqrt(pi)/2) * sqrt(diagC(w))
    //
    // For Nyquist (real draw), we use E|N(0,1)| = sqrt(2/pi).
    //
    m_bands.clear();
    m_bands.reserve((size_t)nbands);

    for (int b = 0; b < nbands; ++b) {
        const double lo_khz = jbands[(Json::ArrayIndex)b][0].asDouble();
        const double hi_khz = jbands[(Json::ArrayIndex)b][1].asDouble();

        BandInfo bi;
        bi.lo_f = lo_khz * units::kilohertz;
        bi.hi_f = hi_khz * units::kilohertz;
        bi.band_index = (size_t)b;

        const Eigen::MatrixXd& A = m_A_band[(size_t)b];
        bi.diagC = A.array().square().rowwise().sum().matrix();

        // Prevent zeros that could lead to divide-by-zero in scaling.
        for (int i = 0; i < bi.diagC.size(); ++i) {
            if (bi.diagC(i) < 1e-30) bi.diagC(i) = 1e-30;
        }

        m_bands.push_back(std::move(bi));
    }

    log->info("Loaded correlated noise model from {} with N_wires={}, N_freq={}, nbands={}",
              fname, (int)m_avg_mag.rows(), (int)m_freq.size(), (int)m_A_band.size());
}

// -----------------------------------------------------------------------------
// Frame utility: sort traces by channel id
// -----------------------------------------------------------------------------
ITrace::vector CorrelatedAddNoise::sorted_traces_by_channel(const IFrame::pointer& frame) const
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
// Helper: find band index for frequency f
// -----------------------------------------------------------------------------
inline int CorrelatedAddNoise::band_for_freq(double f) const
{
    for (size_t bi = 0; bi < m_bands.size(); ++bi) {
        if (f >= m_bands[bi].lo_f && f < m_bands[bi].hi_f) return (int)bi;
    }
    return -1;
}

// -----------------------------------------------------------------------------
// Helper: interpolate avg_mag(w,f) along model frequency axis
// -----------------------------------------------------------------------------
inline double CorrelatedAddNoise::interp_avg_mag_w(int w, double f) const
{
    const int Nf = (int)m_freq.size();
    if (Nf <= 1) return m_avg_mag(w, 0);

    if (f <= m_freq(0))    return m_avg_mag(w, 0);
    if (f >= m_freq(Nf-1)) return m_avg_mag(w, Nf-1);

    const double* begin = m_freq.data();
    const double* end   = m_freq.data() + Nf;
    const double* it = std::lower_bound(begin, end, f);

    int i1 = (int)(it - begin);
    if (i1 <= 0)  return m_avg_mag(w, 0);
    if (i1 >= Nf) return m_avg_mag(w, Nf-1);

    const int i0 = i1 - 1;

    const double f0 = m_freq(i0);
    const double f1 = m_freq(i1);
    const double a  = (f1 > f0) ? ((f - f0) / (f1 - f0)) : 0.0;

    return (1.0 - a) * m_avg_mag(w, i0) + a * m_avg_mag(w, i1);
}

// -----------------------------------------------------------------------------
// Core synthesis routine: build correlated per-wire time-domain noise
//
// Output: noise_time[w, t] in WCT voltage units (consistent with avg_mag units).
//
// Notes:
// - Uses configured dt to convert bin index -> physical frequency.
// - DC forced to 0.
// - Nyquist handled only if N even, and is real-only.
// -----------------------------------------------------------------------------
Eigen::MatrixXd CorrelatedAddNoise::make_correlated_noise(int nwires_frame,
                                                         size_t nsamp) const
{
    const int N_wires_model = (int)m_avg_mag.rows();

    // Keep existing behavior: warn and truncate if frame != model wire counts.
    if (nwires_frame != N_wires_model) {
        log->warn("CorrelatedAddNoise: nwires_frame={} != model N_wires={}, truncating/minor mismatch",
                  nwires_frame, N_wires_model);
        nwires_frame = std::min(nwires_frame, N_wires_model);
    }

    const int    N   = (int)nsamp;
    const int    Kt  = N/2 + 1;
    const double dfT = 1.0 / (double(N) * m_dt);

    // One-sided spectrum per wire
    Eigen::MatrixXcd Xpos = Eigen::MatrixXcd::Zero(nwires_frame, Kt);
    if (Kt > 0) Xpos.col(0).setZero(); // DC forced 0

    const bool is_even = ((N % 2) == 0);
    const int  nyq_k   = is_even ? (Kt - 1) : -1;

    // -------------------------------------------------------------------------
    // Complex bins (exclude DC and Nyquist)
    // -------------------------------------------------------------------------
    for (int kt = 1; kt < Kt; ++kt) {
        if (is_even && kt == nyq_k) continue;

        const double f = kt * dfT;
        const int b = band_for_freq(f);
        if (b < 0) continue; // out-of-band => leave as zero

        const auto& band  = m_bands[(size_t)b];
        const auto& A     = m_A_band[band.band_index];
        const auto& diagC = band.diagC;

        // Draw Z ~ CN(0, I) across all model wires (not just frame wires)
        Eigen::VectorXcd Z(N_wires_model);
        for (int w = 0; w < N_wires_model; ++w) {
            const double xr = m_rng->normal(0.0, 1.0);
            const double xi = m_rng->normal(0.0, 1.0);
            Z(w) = std::complex<double>(xr, xi) / std::sqrt(2.0);
        }

        // Color it: U = A Z  (imposes covariance ~ A A^T)
        const Eigen::VectorXcd U = A.cast<std::complex<double>>() * Z;

        // Per-wire scaling to match target mean magnitude
        for (int w = 0; w < nwires_frame; ++w) {
            if (m_live_mask(w) == 0) continue;

            const double target   = interp_avg_mag_w(w, f);

            // For U(w): Var = diagC(w), so E|U(w)| = (sqrt(pi)/2)*sqrt(diagC(w))
            const double baseline = SQRT_PI_OVER_2 * std::sqrt(diagC(w));
            const double scale    = (baseline > 0.0) ? (target / baseline) : 0.0;

            Xpos(w, kt) = U(w) * scale;
        }
    }

    // -------------------------------------------------------------------------
    // Nyquist bin (real) for even N
    //
    // At Nyquist, the one-sided coefficient must be purely real to ensure a
    // real-valued time-domain waveform. To keep correlations consistent, we:
    //   - draw a real normal vector nvec ~ N(0,I)
    //   - correlate it via y = A nvec
    //   - scale per wire using E|N(0,1)| = sqrt(2/pi)
    // -------------------------------------------------------------------------
    if (is_even && nyq_k >= 0) {
        const double fnyq = nyq_k * dfT;
        const int b = band_for_freq(fnyq);
        if (b >= 0) {
            const auto& band  = m_bands[(size_t)b];
            const auto& A     = m_A_band[band.band_index];
            const auto& diagC = band.diagC;

            Eigen::VectorXd nvec(N_wires_model);
            for (int w = 0; w < N_wires_model; ++w) {
                nvec(w) = m_rng->normal(0.0, 1.0);
            }

            // Correlate real vector
            const Eigen::VectorXd y = A * nvec;

            for (int w = 0; w < nwires_frame; ++w) {
                if (m_live_mask(w) == 0) continue;

                const double target   = interp_avg_mag_w(w, fnyq);

                // For real y(w): stddev is sqrt(diagC(w)) in same A*A^T sense,
                // and E|N(0,1)| = sqrt(2/pi).
                const double baseline = SQRT_2_OVER_PI * std::sqrt(diagC(w));
                const double scale    = (baseline > 0.0) ? (target / baseline) : 0.0;

                Xpos(w, nyq_k) = std::complex<double>(y(w) * scale, 0.0);
            }
        }
    }

    // -------------------------------------------------------------------------
    // IFFT per wire using NumPy-like irfft convention
    // -------------------------------------------------------------------------
    Eigen::MatrixXd noise_time(nwires_frame, N);

    std::vector<std::complex<float>> one_sided((size_t)Kt);
    std::vector<float> x_time;

    for (int w = 0; w < nwires_frame; ++w) {
        for (int k = 0; k < Kt; ++k) {
            const auto xd = Xpos(w, k);
            one_sided[(size_t)k] = std::complex<float>((float)xd.real(), (float)xd.imag());
        }

        // Safety: enforce Nyquist imag=0
        if (is_even) {
            one_sided[(size_t)(Kt - 1)] =
                std::complex<float>(one_sided[(size_t)(Kt - 1)].real(), 0.0f);
        }

        irfft_fftwf(one_sided, x_time, N, (float)m_ifft_scale);

        for (int t = 0; t < N; ++t) {
            noise_time(w, t) = (double)x_time[(size_t)t];
        }
    }

    return noise_time;
}

// -----------------------------------------------------------------------------
// Main filter operator(): read frame, generate noise, add to traces, write new frame
// -----------------------------------------------------------------------------
bool CorrelatedAddNoise::operator()(const input_pointer& inframe,
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
        log->warn("CorrelatedAddNoise: frame nsamples={} differs from configured nsamples={}, adapting per-frame",
                  ncharge, m_nsamples);
    }

    // Generate correlated noise matrix
    Eigen::MatrixXd noise = make_correlated_noise((int)nwires_frame, ncharge);

    // Build new traces with noise added.
    ITrace::vector outtraces;
    outtraces.reserve(nwires_frame);

    for (size_t i = 0; i < nwires_frame; ++i) {
        const auto& intrace = traces[i];
        const int chid = intrace->channel();
        auto charge = intrace->charge();

        if ((int)i < noise.rows()) {
            for (size_t t = 0; t < charge.size(); ++t) {
                charge[t] += noise((int)i, (int)t);
            }
        }

        outtraces.push_back(std::make_shared<Aux::SimpleTrace>(chid, intrace->tbin(), charge));
    }

    // Output a new frame with modified traces
    outframe = std::make_shared<Aux::SimpleFrame>(inframe->ident(),
                                                  inframe->time(),
                                                  outtraces,
                                                  inframe->tick());

    log->debug("CorrelatedAddNoise v{}: call={} frame={} {} traces",
               kComponentVersion, m_count, inframe->ident(), outtraces.size());

    ++m_count;
    return true;
}