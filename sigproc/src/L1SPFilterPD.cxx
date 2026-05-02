#include "WireCellSigProc/L1SPFilterPD.h"

#include "WireCellAux/DftTools.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"

#include "WireCellIface/IFilterWaveform.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/LassoModel.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Eigen.h"
#include "WireCellUtil/Persist.h"
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
    // Split-sign accumulators (same threshold gate as temp_sum / temp1_sum)
    double temp_sum_pos{0}, temp_sum_neg{0};
    int    n_above_pos{0}, n_above_neg{0};
    // Absolute tick of max_val / min_val within the ROI
    int    argmax_tick{-1}, argmin_tick{-1};
    // Decon (gauss) peak and integral over the ROI, ungated
    double sig_peak{-1e30}, sig_integral{0};
    // Gauss-side shape features (in-ROI)
    double gmax{0};                 // max(|gauss[t]|) for t in [start, end)
    double gauss_abs_sum_roi{0};    // Σ|gauss[t]| over [start, end)
    double gauss_fill{0};           // gauss_abs_sum_roi / (gmax · nbin_fit)
    double gauss_fwhm_frac{0};      // fraction of ticks with |gauss| > gmax/2
    // Wide-window features (require pad outside the ROI)
    double roi_energy_frac{0};      // gauss_abs_sum_roi / Σ|gauss| over padded
    double raw_asym_wide{0};        // (pos+neg)/(pos-neg) of raw ADC over padded
    // Core sub-window (contiguous run of |gauss|>core_g_thr around argmax|gauss|).
    // Aligns the trigger features with iter-7's `|gauss|>g_thr=50` ROIs even when
    // the C++ gauss>0 + raw-noise ROI extraction reports a wider window.
    int    core_lo{-1}, core_hi{-1};   // tick bounds (absolute), -1 if no core
    int    core_length{0};
    double core_fill{0};               // Σ|gauss[core]|/(gmax·core_length)
    double core_fwhm_frac{0};
    double core_raw_asym_wide{0};      // (pos+neg)/(pos-neg) on raw over core ± pad
};

// Compute per-ROI shape and asymmetry quantities for ticks [start_tick, end_tick).
// adc/sig charges are indexed with the same tbin offset (both traces assumed
// to share the same tbin, which is the case for all pdhd/pdvd frames).
//
// The wide-window features (roi_energy_frac, raw_asym_wide) are computed over
// padded windows clamped to the trace boundary.  Wide-padded ranges that run
// past the trace edges are silently truncated; downstream divides by zero are
// guarded.
//
//   threshold     — per-tick |ADC| gate for temp_sum / temp1_sum / temp2_sum
//   energy_pad    — ticks added on each side for roi_energy_frac denominator
//   raw_asym_pad  — ticks added on each side for raw_asym_wide / core_raw_asym_wide
//   raw_eps       — per-tick raw ADC threshold for the asymmetry sum (sign-gated)
//   core_g_thr    — per-tick |gauss| gate defining the core sub-window
//                   (matches iter-7 g_thr=50 ADC; pass 0 to disable the core
//                   computation, in which case core fields stay at defaults)
AsymRecord compute_asym(const WireCell::ITrace::ChargeSequence& adc,
                        const WireCell::ITrace::ChargeSequence& sig,
                        int tbin,
                        int start_tick, int end_tick,
                        double threshold,
                        int    energy_pad,
                        int    raw_asym_pad,
                        double raw_eps,
                        double core_g_thr)
{
    AsymRecord r;
    r.nbin_fit = end_tick - start_tick;
    if (r.nbin_fit <= 0) return r;

    const int ntot_adc = (int)adc.size();
    const int ntot_sig = (int)sig.size();

    // Pass 1: in-ROI accumulators (existing + gmax + gauss_abs_sum_roi).
    for (int i = 0; i < r.nbin_fit; i++) {
        const int idx = i + start_tick - tbin;
        const double w = adc.at(idx);
        const double b = sig.at(idx);
        if (w > r.max_val) { r.max_val = w; r.argmax_tick = start_tick + i; }
        if (w < r.min_val) { r.min_val = w; r.argmin_tick = start_tick + i; }
        if (std::fabs(w) > threshold) {
            r.temp_sum  += w;
            r.temp1_sum += std::fabs(w);
            r.temp2_sum += std::fabs(b);
            if (w > 0) { r.temp_sum_pos += w; ++r.n_above_pos; }
            else       { r.temp_sum_neg += w; ++r.n_above_neg; }
        }
        if (b > r.sig_peak) r.sig_peak = b;
        r.sig_integral += b;
        const double absb = std::fabs(b);
        if (absb > r.gmax) r.gmax = absb;
        r.gauss_abs_sum_roi += absb;
    }

    // Pass 2: gauss_fill + gauss_fwhm_frac (need gmax from pass 1).
    if (r.gmax > 0) {
        const double half = 0.5 * r.gmax;
        int n_above_half = 0;
        for (int i = 0; i < r.nbin_fit; i++) {
            const int idx = i + start_tick - tbin;
            if (std::fabs(sig.at(idx)) > half) ++n_above_half;
        }
        r.gauss_fwhm_frac = (double)n_above_half / (double)r.nbin_fit;
        r.gauss_fill      = r.gauss_abs_sum_roi
                          / (r.gmax * (double)r.nbin_fit);
    }

    // Pass 3: roi_energy_frac over [start - energy_pad, end + energy_pad).
    {
        const int wide_lo = std::max(0, (start_tick - tbin) - energy_pad);
        const int wide_hi = std::min(ntot_sig, (end_tick - tbin) + energy_pad);
        double wide_sum = 0;
        for (int idx = wide_lo; idx < wide_hi; idx++) {
            wide_sum += std::fabs(sig.at(idx));
        }
        if (wide_sum > 0) {
            r.roi_energy_frac = r.gauss_abs_sum_roi / wide_sum;
        }
    }

    // Pass 4: raw_asym_wide over [start - raw_asym_pad, end + raw_asym_pad).
    {
        const int wide_lo = std::max(0, (start_tick - tbin) - raw_asym_pad);
        const int wide_hi = std::min(ntot_adc, (end_tick - tbin) + raw_asym_pad);
        double pos = 0, neg = 0;
        for (int idx = wide_lo; idx < wide_hi; idx++) {
            const double w = adc.at(idx);
            if      (w >  raw_eps) pos += w;
            else if (w < -raw_eps) neg += w;   // neg is ≤ 0 by construction
        }
        const double denom = pos - neg;        // pos - neg = pos + |neg|
        if (denom > 0) {
            r.raw_asym_wide = (pos + neg) / denom;
        }
    }

    // Pass 5: scan ALL contiguous runs of |gauss|>core_g_thr inside the C++
    // ROI and pick the one most likely to be the artifact body.
    //
    // Why iterate all runs: the C++ ROI extraction (gauss>0 + raw-noise merge)
    // can contain multiple separate iter-7 ROIs (each is `|gauss|>g_thr`).
    // Anchoring only on argmax(|gauss|) misses the others.
    //
    // Selection score: `len · |aw|` — rewards runs that are simultaneously
    // long *and* asymmetric.  A short +1.0-asym noise sliver loses to a
    // medium-length 0.6-asym artifact body.  Ties broken by length.
    if (core_g_thr > 0 && r.gmax > core_g_thr) {
        double best_score = -1.0;
        int   best_lo = -1, best_hi = -1, best_len = 0;
        double best_aw = 0, best_fill = 0, best_fwhm = 0;

        int run_lo = -1;
        for (int i = 0; i <= r.nbin_fit; i++) {
            const bool above = (i < r.nbin_fit) &&
                std::fabs(sig.at(i + start_tick - tbin)) > core_g_thr;
            if (above) {
                if (run_lo < 0) run_lo = i;
                continue;
            }
            if (run_lo < 0) continue;
            const int lo_i = run_lo, hi_i = i - 1;
            run_lo = -1;
            const int run_start = lo_i + start_tick;
            const int run_end   = hi_i + start_tick + 1;
            const int run_len   = run_end - run_start;

            double abs_sum = 0;
            int n_above_half = 0;
            const double half_core = 0.5 * r.gmax;
            for (int j = lo_i; j <= hi_i; j++) {
                const double absb = std::fabs(sig.at(j + start_tick - tbin));
                abs_sum += absb;
                if (absb > half_core) ++n_above_half;
            }
            const double fill = (r.gmax > 0)
                ? abs_sum / (r.gmax * (double)run_len) : 0.0;
            const double fwhm = (double)n_above_half / (double)run_len;

            const int wlo = std::max(0, (run_start - tbin) - raw_asym_pad);
            const int whi = std::min(ntot_adc, (run_end - tbin) + raw_asym_pad);
            double pos = 0, neg = 0;
            for (int idx = wlo; idx < whi; idx++) {
                const double w = adc.at(idx);
                if      (w >  raw_eps) pos += w;
                else if (w < -raw_eps) neg += w;
            }
            const double denom = pos - neg;
            const double aw    = (denom > 0) ? (pos + neg) / denom : 0.0;
            const double score = (double)run_len * std::fabs(aw);

            if (score > best_score ||
                (score == best_score && run_len > best_len)) {
                best_score = score;
                best_aw    = aw;
                best_fill  = fill;
                best_fwhm  = fwhm;
                best_lo    = run_start;
                best_hi    = run_end - 1;
                best_len   = run_len;
            }
        }
        if (best_len > 0) {
            r.core_lo            = best_lo;
            r.core_hi            = best_hi;
            r.core_length        = best_len;
            r.core_fill          = best_fill;
            r.core_fwhm_frac     = best_fwhm;
            r.core_raw_asym_wide = best_aw;
        }
    }

    return r;
}

// ── Per-ROI trigger gate (Strategy B retuned) ─────────────────────────────────
// Knobs bundled into a small struct so the same gate can be applied in
// l1_fit() (drives the LASSO) and in the dump path (records the decision
// alongside per-ROI features for offline cross-checking).
struct TriggerCfg {
    int    min_length;
    double gmax_min;
    double energy_frac_thr;
    double asym_strong, asym_mod, asym_loose;
    int    len_long_mod, len_long_loose, len_fill_shape;
    double fill_shape_fill_thr, fill_shape_fwhm_thr;
};

// Per-sub-window trigger walk.  Iterates every contiguous run of
// |gauss|>core_g_thr inside [start_tick, end_tick) and tests the multi-arm
// gate against each run's own features.  Fires (returns +1/-1) on the first
// run that passes — matches iter-7's per-candidate gating where each
// |gauss|>g_thr ROI is its own trigger candidate.
//
// Why per-sub-window rather than per-aggregate-or-per-best:
//   - The C++ ROI extraction (gauss>0+raw-noise-merge) often spans multiple
//     iter-7 candidates; aggregating mixes their features ("longest length"
//     from one sub-window with "extreme asym" from another) which lets the
//     L_long arm fire on weak-asym backgrounds.
//   - Picking one "best" sub-window misses real artifacts whose features are
//     split across several |gauss|>50 runs (e.g. evt 12 apa 1 ch=203-212).
//
// Returns: -1, 0, or +1.  Polarity = sign(raw_asym_wide) of the firing run.
int decide_trigger(const WireCell::ITrace::ChargeSequence& adc,
                   const WireCell::ITrace::ChargeSequence& sig,
                   int tbin,
                   int start_tick, int end_tick,
                   double gmax,
                   const TriggerCfg& cfg,
                   double core_g_thr,
                   int    energy_pad,
                   int    raw_asym_pad,
                   double raw_eps)
{
    if (gmax < cfg.gmax_min) return 0;
    if (core_g_thr <= 0)     return 0;

    const int nbin = end_tick - start_tick;
    if (nbin <= 0) return 0;

    const int ntot_adc = (int)adc.size();
    const int ntot_sig = (int)sig.size();
    const double half_g = 0.5 * gmax;

    // First pass: collect every |gauss|>core_g_thr sub-window plus its
    // per-sub-window (fill, fwhm, ef, aw).  Each sub-window is one
    // iter-7-style trigger candidate.
    struct SubInfo {
        int    lo_i, hi_i, run_len;
        double fill, fwhm, ef, aw;
    };
    std::vector<SubInfo> subs;
    subs.reserve(8);
    int run_lo = -1;
    for (int i = 0; i <= nbin; i++) {
        const bool above = (i < nbin) &&
            std::fabs(sig.at(i + start_tick - tbin)) > core_g_thr;
        if (above) { if (run_lo < 0) run_lo = i; continue; }
        if (run_lo < 0) continue;
        const int lo_i = run_lo, hi_i = i - 1;
        run_lo = -1;
        const int run_start = lo_i + start_tick;
        const int run_end   = hi_i + start_tick + 1;
        const int run_len   = run_end - run_start;
        if (run_len < cfg.min_length) continue;

        double abs_sum = 0;
        int n_above_half = 0;
        for (int j = lo_i; j <= hi_i; j++) {
            const double absb = std::fabs(sig.at(j + start_tick - tbin));
            abs_sum += absb;
            if (absb > half_g) ++n_above_half;
        }
        const double fill = abs_sum / (gmax * (double)run_len);
        const double fwhm = (double)n_above_half / (double)run_len;

        double wide_sum = 0;
        const int elo = std::max(0, (run_start - tbin) - energy_pad);
        const int ehi = std::min(ntot_sig, (run_end - tbin) + energy_pad);
        for (int idx = elo; idx < ehi; idx++) {
            wide_sum += std::fabs(sig.at(idx));
        }
        const double ef = (wide_sum > 0) ? abs_sum / wide_sum : 0.0;

        double pos = 0, neg = 0;
        const int wlo = std::max(0, (run_start - tbin) - raw_asym_pad);
        const int whi = std::min(ntot_adc, (run_end - tbin) + raw_asym_pad);
        for (int idx = wlo; idx < whi; idx++) {
            const double w = adc.at(idx);
            if      (w >  raw_eps) pos += w;
            else if (w < -raw_eps) neg += w;
        }
        const double denom = pos - neg;
        const double aw    = (denom > 0) ? (pos + neg) / denom : 0.0;

        subs.push_back({lo_i, hi_i, run_len, fill, fwhm, ef, aw});
    }
    if (subs.empty()) return 0;

    // Per-sub-window gate: each candidate sub-window must individually pass
    // the energy-fraction (isolated-lobe) precondition before the arm tests.
    for (const auto& s : subs) {
        if (s.ef < cfg.energy_frac_thr) continue;
        const double aabs = std::fabs(s.aw);
        const bool fire =
            (aabs >= cfg.asym_strong) ||
            (s.run_len >= cfg.len_long_mod   && aabs >= cfg.asym_mod)   ||
            (s.run_len >= cfg.len_long_loose && aabs >= cfg.asym_loose) ||
            (s.run_len >= cfg.len_fill_shape &&
             s.fill <= cfg.fill_shape_fill_thr &&
             s.fwhm <= cfg.fill_shape_fwhm_thr &&
             aabs >= cfg.asym_mod);
        if (fire) return (s.aw > 0.0) ? +1 : -1;
    }
    return 0;
}

}  // namespace

// ── L1SPFilterPD ─────────────────────────────────────────────────────────────

L1SPFilterPD::L1SPFilterPD()
  : Aux::Logger("L1SPFilterPD", "sigproc")
{
}

WireCell::Configuration L1SPFilterPD::default_configuration() const
{
    Configuration cfg;

    // Path (resolved via WIRECELL_PATH) to the JSON+bz2 file holding the
    // pre-built L1SP response kernels.  Generated offline with:
    //   wirecell-sigproc gen-l1sp-kernels --gain '14*mV/fC' --shaping '2.2*us'
    //     --postgain 1.2 --adc-per-mv 2.048 --coarse-time-offset '-8*us'
    //     <field-response.json.bz2>  <out>_l1sp_kernels.json.bz2
    // See wire-cell-python/wirecell/sigproc/l1sp.py for the schema.
    cfg["kernels_file"] = "";
    // Multiplier applied to every loaded kernel amplitude at init_resp() time.
    // Kernels in ``kernels_file`` are in ADC/electron at the reference 14 mV/fC
    // FE gain; set this to params.elec.gain / (14*mV/fC) when the detector runs
    // at a different gain (same gain_scale as for ADC-domain thresholds).
    cfg["kernels_scale"] = 1.0;

    cfg["filter"] = Json::arrayValue;
    // Auto-derived smearing kernel: if "filter" is empty, look up this
    // IFilterWaveform (by "TypeName:InstanceName"), IFFT it, take a centered
    // window of taps above kernel_threshold*peak, and sum-normalize.
    // Set to "" to leave m_smearing_vec empty (disables L1SP smearing).
    cfg["gauss_filter"]      = "HfFilter:Gaus_wide";
    cfg["kernel_threshold"]  = 1.0e-3;   // relative amplitude cutoff
    cfg["kernel_max_half"]   = 64;        // max half-width in ticks (safety cap)
    cfg["kernel_nticks"]     = 4096;      // IFFT length (>> expected kernel width)

    cfg["adctag"] = "raw";
    cfg["sigtag"] = "gauss";
    cfg["outtag"] = "l1sp";

    cfg["raw_ROI_th_nsigma"] = 4;
    cfg["raw_ROI_th_adclimit"] = 10;
    cfg["overall_time_offset"] = 0;

    cfg["roi_pad"] = 3;
    cfg["raw_pad"] = 15;

    cfg["adc_l1_threshold"] = 6;
    cfg["adc_sum_threshold"] = 160;
    cfg["adc_sum_rescaling"] = 90.;
    cfg["adc_ratio_threshold"] = 0.2;

    // Per-ROI trigger gate (PDHD/PDVD Strategy B retuned).  See header
    // for the exact gate definition.  Defaults seeded from the iter-7
    // offline detector (find_long_decon_artifacts.py).
    cfg["l1_min_length"]            = m_l1_min_length;
    cfg["l1_gmax_min"]              = m_l1_gmax_min;
    cfg["l1_energy_frac_thr"]       = m_l1_energy_frac_thr;
    cfg["l1_energy_pad_ticks"]      = m_l1_energy_pad_ticks;
    cfg["l1_raw_asym_pad_ticks"]    = m_l1_raw_asym_pad_ticks;
    cfg["l1_raw_asym_eps"]          = m_l1_raw_asym_eps;
    cfg["l1_core_g_thr"]            = m_l1_core_g_thr;
    cfg["l1_asym_strong"]           = m_l1_asym_strong;
    cfg["l1_asym_mod"]              = m_l1_asym_mod;
    cfg["l1_asym_loose"]            = m_l1_asym_loose;
    cfg["l1_len_long_mod"]          = m_l1_len_long_mod;
    cfg["l1_len_long_loose"]        = m_l1_len_long_loose;
    cfg["l1_len_fill_shape"]        = m_l1_len_fill_shape;
    cfg["l1_fill_shape_fill_thr"]   = m_l1_fill_shape_fill_thr;
    cfg["l1_fill_shape_fwhm_thr"]   = m_l1_fill_shape_fwhm_thr;

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

    cfg["dft"] = "FftwDFT";

    // Plane-scope filter: IAnodePlane typename + list of plane indices to process.
    // Leave "anode" empty to disable plane filtering (all channels processed).
    cfg["anode"] = "";
    cfg["process_planes"][0] = 0;   // U
    cfg["process_planes"][1] = 1;   // V

    // Optional per-channel eligibility whitelist (channel ID ints).
    // Empty = all channels in process_planes are eligible (default).
    cfg["eligible_channels"] = Json::arrayValue;

    // Calibration dump mode: write per-ROI asymmetry records to NPZ files.
    cfg["dump_mode"] = false;
    cfg["dump_path"] = "";    // directory; one NPZ per operator() call written here
    cfg["dump_tag"] = "";     // label baked into filename (e.g. "apa1")
    // Waveform dump: write per-triggered-ROI NPZ (raw/decon/lasso/smeared).
    // Non-empty path enables it; files go under <waveform_dump_path>/<dump_tag>_<frame_ident>/.
    cfg["waveform_dump_path"] = "";

    return cfg;
}

void L1SPFilterPD::configure(const WireCell::Configuration& cfg)
{
    m_kernels_file  = get<std::string>(cfg, "kernels_file", "");
    m_kernels_scale = get(cfg, "kernels_scale", m_kernels_scale);

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

    m_adc_l1_threshold = get(cfg, "adc_l1_threshold", m_adc_l1_threshold);
    m_adc_sum_threshold = get(cfg, "adc_sum_threshold", m_adc_sum_threshold);
    m_adc_sum_rescaling = get(cfg, "adc_sum_rescaling", m_adc_sum_rescaling);
    m_adc_ratio_threshold = get(cfg, "adc_ratio_threshold", m_adc_ratio_threshold);

    m_l1_min_length          = get(cfg, "l1_min_length",          m_l1_min_length);
    m_l1_gmax_min            = get(cfg, "l1_gmax_min",            m_l1_gmax_min);
    m_l1_energy_frac_thr     = get(cfg, "l1_energy_frac_thr",     m_l1_energy_frac_thr);
    m_l1_energy_pad_ticks    = get(cfg, "l1_energy_pad_ticks",    m_l1_energy_pad_ticks);
    m_l1_raw_asym_pad_ticks  = get(cfg, "l1_raw_asym_pad_ticks",  m_l1_raw_asym_pad_ticks);
    m_l1_raw_asym_eps        = get(cfg, "l1_raw_asym_eps",        m_l1_raw_asym_eps);
    m_l1_core_g_thr          = get(cfg, "l1_core_g_thr",          m_l1_core_g_thr);
    m_l1_asym_strong         = get(cfg, "l1_asym_strong",         m_l1_asym_strong);
    m_l1_asym_mod            = get(cfg, "l1_asym_mod",            m_l1_asym_mod);
    m_l1_asym_loose          = get(cfg, "l1_asym_loose",          m_l1_asym_loose);
    m_l1_len_long_mod        = get(cfg, "l1_len_long_mod",        m_l1_len_long_mod);
    m_l1_len_long_loose      = get(cfg, "l1_len_long_loose",      m_l1_len_long_loose);
    m_l1_len_fill_shape      = get(cfg, "l1_len_fill_shape",      m_l1_len_fill_shape);
    m_l1_fill_shape_fill_thr = get(cfg, "l1_fill_shape_fill_thr", m_l1_fill_shape_fill_thr);
    m_l1_fill_shape_fwhm_thr = get(cfg, "l1_fill_shape_fwhm_thr", m_l1_fill_shape_fwhm_thr);

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
    if (m_smearing_vec.empty()) {
        m_gauss_filter_tn  = get<std::string>(cfg, "gauss_filter",     m_gauss_filter_tn);
        m_kernel_threshold = get(cfg, "kernel_threshold", m_kernel_threshold);
        m_kernel_max_half  = get(cfg, "kernel_max_half",  m_kernel_max_half);
        m_kernel_nticks    = get(cfg, "kernel_nticks",    m_kernel_nticks);

        if (!m_gauss_filter_tn.empty()) {
            // The IFFT bin spacing is 1/(2·max_freq); for the standard
            // HfFilter max_freq=1 MHz this implies a 500 ns tick — the SP
            // tick on uBooNE and on PDHD post-resampler.  Same assumption
            // already baked into build_G (t_meas = i * 0.5 * units::us).
            auto hf = Factory::find_tn<IFilterWaveform>(m_gauss_filter_tn);
            const int N = m_kernel_nticks;
            auto freq_wf = hf->filter_waveform(N);   // vector<float>, length N, freq-domain

            Aux::DftTools::complex_vector_t spec(N);
            for (int i = 0; i < N; ++i) spec[i] = {freq_wf[i], 0.0f};
            auto time_wf = inv_c2r(m_dft, spec);     // vector<float>, peak at [0], wraps for t<0

            // Find half-width by scanning outward from the peak until both
            // the positive-time and negative-time (wrapped) sides drop below threshold.
            const double peak = std::fabs(time_wf[0]);
            const double thr  = m_kernel_threshold * peak;
            int n_half = 0;
            for (int k = 1; k <= m_kernel_max_half; ++k) {
                if (std::fabs(time_wf[k]) < thr && std::fabs(time_wf[N - k]) < thr) break;
                n_half = k;
            }

            // Build centered kernel: time index i ∈ [-n_half, n_half]
            // maps to array index (i + N) % N in time_wf (circular).
            m_smearing_vec.assign(2 * n_half + 1, 0.0);
            for (int i = -n_half; i <= n_half; ++i)
                m_smearing_vec[i + n_half] = time_wf[(i + N) % N];

            // Sum-normalize → kernel sums to 1 (absorbs small DC=0 offset).
            double s = 0.0;
            for (double v : m_smearing_vec) s += v;
            if (s > 0.0)
                for (double& v : m_smearing_vec) v /= s;

            log->debug("smearing kernel from {} n_half={} ntaps={} peak={:.5f}",
                       m_gauss_filter_tn, n_half, (int)m_smearing_vec.size(),
                       m_smearing_vec[n_half]);
        }
        // If gauss_filter_tn is empty, m_smearing_vec stays empty → L1SP smearing disabled.
    }

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
    m_eligible_channels.clear();
    if (cfg.isMember("eligible_channels") && cfg["eligible_channels"].isArray()) {
        for (auto const& v : cfg["eligible_channels"]) {
            m_eligible_channels.insert(v.asInt());
        }
    }

    // Calibration dump mode
    m_dump_mode   = get(cfg, "dump_mode", m_dump_mode);
    m_dump_path   = get<std::string>(cfg, "dump_path", m_dump_path);
    m_dump_tag    = get<std::string>(cfg, "dump_tag", m_dump_tag);
    m_wf_dump_path = get<std::string>(cfg, "waveform_dump_path", m_wf_dump_path);

    // Reset interpolators so init_resp() reloads them on next operator() call.
    m_lin_bipolar.clear();
    m_lin_pos_unipolar.clear();
    m_lin_neg_unipolar.clear();
    m_unipolar_toff_pos.clear();
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

bool L1SPFilterPD::channel_eligible(int ch) const
{
    if (m_eligible_channels.empty()) return true;
    return m_eligible_channels.count(ch) > 0;
}

void L1SPFilterPD::init_resp()
{
    if (!m_lin_bipolar.empty()) return;   // already loaded

    if (m_kernels_file.empty()) {
        THROW(ValueError() << errmsg{"L1SPFilterPD: 'kernels_file' is required. "
              "Generate one with: wirecell-sigproc gen-l1sp-kernels "
              "<field-response> <output>.json.bz2"});
    }

    auto top = Persist::load(m_kernels_file);   // resolves WIRECELL_PATH; .json.bz2 OK

    if (!top.isMember("meta") || !top.isMember("planes")) {
        THROW(ValueError() << errmsg{"L1SPFilterPD: malformed kernels_file '"
              + m_kernels_file + "' (missing 'meta' or 'planes')"});
    }

    const auto& meta = top["meta"];
    const double period_ns = meta["period_ns"].asDouble();
    const double t0_us     = meta["t0_us"].asDouble();
    const int    n_samples = meta["n_samples"].asInt();
    const double xstep     = period_ns * units::ns;
    const double x0        = t0_us * units::us;

    auto load_array = [&](const Json::Value& jv) {
        Waveform::realseq_t v;
        v.reserve(jv.size());
        for (const auto& x : jv) v.push_back(x.asDouble());
        return v;
    };

    auto make_lin = [&](const Waveform::realseq_t& v) {
        return std::make_unique<linterp<double>>(v.begin(), v.end(), x0, xstep);
    };

    for (const auto& jpl : top["planes"]) {
        const int plane = jpl["plane_index"].asInt();

        auto k_bip_pos = load_array(jpl["positive"]["bipolar"]);
        auto k_uni_pos = load_array(jpl["positive"]["unipolar"]);
        auto k_uni_neg = load_array(jpl["negative"]["unipolar"]);
        const double toff_pos_us = jpl["positive"]["unipolar_time_offset_us"].asDouble();

        if ((int)k_bip_pos.size() != n_samples ||
            (int)k_uni_pos.size() != n_samples ||
            (int)k_uni_neg.size() != n_samples) {
            THROW(ValueError() << errmsg{"L1SPFilterPD: kernel length mismatch for plane "
                  + std::to_string(plane) + " in '" + m_kernels_file + "'"});
        }

        if (m_kernels_scale != 1.0) {
            for (auto& v : k_bip_pos) v *= m_kernels_scale;
            for (auto& v : k_uni_pos) v *= m_kernels_scale;
            for (auto& v : k_uni_neg) v *= m_kernels_scale;
        }

        m_lin_bipolar[plane]      = make_lin(k_bip_pos);
        m_lin_pos_unipolar[plane] = make_lin(k_uni_pos);
        m_lin_neg_unipolar[plane] = make_lin(k_uni_neg);
        m_unipolar_toff_pos[plane] = toff_pos_us * units::us;

        log->debug("loaded plane {} kernels from {}: n={} period={} ns t0={:.3f} us "
                   "W-shift(pos)={:+.3f} us kernels_scale={:.4f}",
                   plane, m_kernels_file, n_samples, period_ns, t0_us, toff_pos_us,
                   m_kernels_scale);
    }
}

int L1SPFilterPD::l1_fit(std::shared_ptr<Aux::SimpleTrace>& newtrace,
                          const std::shared_ptr<const WireCell::ITrace>& adctrace,
                          const std::shared_ptr<const WireCell::ITrace>& sigtrace,
                          int start_tick, int end_tick, int plane,
                          std::vector<double>* lasso_unsmeared)
{
    const int nbin_fit = end_tick - start_tick;

    // Compute all per-ROI features once.  sigtrace carries the unmodified
    // gauss data for the whole frame; reading from it (rather than newtrace)
    // ensures the wide-window energy fraction is not corrupted by L1 fits
    // already written into prior ROIs on the same channel.
    const AsymRecord rec = compute_asym(adctrace->charge(),
                                        sigtrace->charge(),
                                        newtrace->tbin(),
                                        start_tick, end_tick,
                                        m_adc_l1_threshold,
                                        m_l1_energy_pad_ticks,
                                        m_l1_raw_asym_pad_ticks,
                                        m_l1_raw_asym_eps,
                                        m_l1_core_g_thr);

    // Per-ROI multi-arm gate, walked per |gauss|>core_g_thr sub-window.
    const TriggerCfg tcfg{
        m_l1_min_length, m_l1_gmax_min, m_l1_energy_frac_thr,
        m_l1_asym_strong, m_l1_asym_mod, m_l1_asym_loose,
        m_l1_len_long_mod, m_l1_len_long_loose, m_l1_len_fill_shape,
        m_l1_fill_shape_fill_thr, m_l1_fill_shape_fwhm_thr,
    };
    int flag_l1 = decide_trigger(adctrace->charge(), sigtrace->charge(),
                                 newtrace->tbin(),
                                 start_tick, end_tick,
                                 rec.gmax, tcfg,
                                 m_l1_core_g_thr,
                                 m_l1_energy_pad_ticks,
                                 m_l1_raw_asym_pad_ticks,
                                 m_l1_raw_asym_eps);

    // Build the per-tick LASSO input from raw ADC.  Kept as a separate pass
    // (rather than wiring into compute_asym) to keep the helper data-only.
    VectorXd init_W = VectorXd::Zero(nbin_fit);
    for (int i = 0; i < nbin_fit; i++) {
        init_W(i) = adctrace->charge().at(i + start_tick - newtrace->tbin());
    }

    // Select the appropriate per-plane bases.  Bipolar is the same for
    // positive and negative cases; unipolar differs (W kernel vs neg-half).
    auto bip_it = m_lin_bipolar.find(plane);
    if (bip_it == m_lin_bipolar.end()) {
        log->warn("l1_fit: no kernels loaded for plane {}; passing through", plane);
        return 0;
    }
    linterp<double>* lin_bipolar  = bip_it->second.get();
    linterp<double>* lin_unipolar = nullptr;
    double basis1_toff = 0.0;     // negative case: no shift
    if (flag_l1 > 0) {
        auto it = m_lin_pos_unipolar.find(plane);
        if (it != m_lin_pos_unipolar.end()) lin_unipolar = it->second.get();
        auto tit = m_unipolar_toff_pos.find(plane);
        if (tit != m_unipolar_toff_pos.end()) basis1_toff = tit->second;
    }
    else if (flag_l1 < 0) {
        auto it = m_lin_neg_unipolar.find(plane);
        if (it != m_lin_neg_unipolar.end()) lin_unipolar = it->second.get();
    }

    if ((flag_l1 == 1 || flag_l1 == -1) && lin_unipolar == nullptr) {
        log->warn("l1_fit: polarity {} triggered on plane {} but unipolar kernel "
                  "not loaded; falling back to pass-through", flag_l1, plane);
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
                                 basis1_toff,
                                 m_l1_scaling_factor, m_l1_resp_scale,
                                 lin_bipolar, lin_unipolar);

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
            // Snapshot the unsmeared LASSO output before smearing overwrites l1_signal.
            if (lasso_unsmeared) {
                lasso_unsmeared->assign(l1_signal.begin(), l1_signal.end());
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

    return flag_l1;
}

void L1SPFilterPD::dump_roi_waveforms(int frame_ident, int channel, int plane,
                                       int start_tick, int end_tick, int polarity,
                                       const std::shared_ptr<const WireCell::ITrace>& adctrace,
                                       const std::shared_ptr<const WireCell::ITrace>& sigtrace,
                                       const std::shared_ptr<Aux::SimpleTrace>& newtrace,
                                       const std::vector<double>& lasso_unsmeared)
{
    const int nbin = end_tick - start_tick;
    const int tbin = sigtrace->tbin();

    // Build output directory and filename.
    const std::string subdir = fmt::format("{}/{}_{:04d}_{}", m_wf_dump_path,
                                           m_dump_tag, m_count, frame_ident);
    std::error_code ec;
    std::filesystem::create_directories(subdir, ec);
    const std::string polsign = (polarity > 0) ? "pos" : "neg";
    const std::string fname = fmt::format("{}/wf_p{}_c{}_t{}_{}.npz",
                                          subdir, plane, channel, start_tick, polsign);

    // Slice the four waveforms over [start_tick, end_tick).
    std::vector<float> raw_arr(nbin), decon_arr(nbin), smeared_arr(nbin);
    for (int i = 0; i < nbin; ++i) {
        const int t = start_tick + i;
        raw_arr[i]    = adctrace->charge().at(t - adctrace->tbin());
        decon_arr[i]  = sigtrace->charge().at(t - tbin);
        smeared_arr[i] = newtrace->charge().at(t - newtrace->tbin());
    }

    bool first = true;
    auto save_f32 = [&](const std::string& key, const std::vector<float>& v) {
        cnpy::npz_save(fname, key, v.data(), {v.size()}, first ? "w" : "a");
        first = false;
    };
    auto save_f64v = [&](const std::string& key, const std::vector<double>& v) {
        cnpy::npz_save(fname, key, v.data(), {v.size()}, first ? "w" : "a");
        first = false;
    };
    auto save_i32s = [&](const std::string& key, int32_t val) {
        cnpy::npz_save(fname, key, &val, {1}, first ? "w" : "a");
        first = false;
    };

    save_f32("raw",    raw_arr);
    save_f32("decon",  decon_arr);
    save_f64v("lasso", lasso_unsmeared);
    save_f32("smeared", smeared_arr);
    save_i32s("channel",     (int32_t)channel);
    save_i32s("plane",       (int32_t)plane);
    save_i32s("start_tick",  (int32_t)start_tick);
    save_i32s("end_tick",    (int32_t)end_tick);
    save_i32s("polarity",    (int32_t)polarity);
    save_i32s("frame_ident", (int32_t)frame_ident);
    save_i32s("call_count",  (int32_t)m_count);
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
    // Only in-scope, eligible channels need ROI build work.
    std::map<int, std::set<int>> init_map;
    for (auto trace : sigtraces) {
        int ch = trace->channel();
        if (!channel_in_scope(ch) || !channel_eligible(ch)) continue;
        int tbin = trace->tbin();
        auto const& charges = trace->charge();
        std::set<int>& ticks = init_map[ch];
        for (int qi = 0; qi < (int)charges.size(); qi++) {
            if (charges[qi] > 0) ticks.insert(tbin + qi);
        }
    }

    // Augment with ticks above the raw-ADC noise threshold.
    // adctrace_ch_map is populated for all channels; the expensive
    // percentile sort and tick addition are skipped for out-of-scope ones.
    std::map<int, std::shared_ptr<const WireCell::ITrace>> adctrace_ch_map;
    for (auto trace : adctraces) {
        int ch = trace->channel();
        adctrace_ch_map[ch] = trace;
        if (!channel_in_scope(ch) || !channel_eligible(ch)) continue;
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
    ITrace::vector out_traces;

    // Calibration dump: per-ROI parallel vectors accumulated across all channels.
    std::vector<int32_t> d_channel, d_roi_start, d_roi_end, d_nbin_fit;
    std::vector<double>  d_temp_sum, d_temp1_sum, d_temp2_sum, d_max_val, d_min_val;
    std::vector<int32_t> d_prev_roi_end, d_next_roi_start, d_prev_gap, d_next_gap;
    // Tier 1: pre-computed flag/ratio and split-sign shape discriminants
    std::vector<int32_t> d_flag;
    std::vector<double>  d_ratio, d_temp_sum_pos, d_temp_sum_neg;
    std::vector<int32_t> d_n_above_pos, d_n_above_neg;
    // Tier 2: peak locations and decon-side scalars
    std::vector<int32_t> d_argmax_tick, d_argmin_tick;
    std::vector<double>  d_sig_peak, d_sig_integral;
    // Tier 3: per-ROI shape features driving the new trigger gate
    std::vector<double>  d_gmax, d_gauss_fill, d_gauss_fwhm_frac;
    std::vector<double>  d_roi_energy_frac, d_raw_asym_wide;
    // Tier 3 (cont.): core sub-window features actually used by the trigger.
    std::vector<int32_t> d_core_lo, d_core_hi, d_core_length;
    std::vector<double>  d_core_fill, d_core_fwhm_frac, d_core_raw_asym_wide;
    // Tier 3 (cont.): the new-gate trigger decision actually applied.  Kept
    // alongside the legacy 'flag' (above, derived from the old ratio test) so
    // offline analyses can compare both decisions on the same ROI.
    std::vector<int32_t> d_flag_l1;

    // Trigger gate config snapshot used by both dump-side decide_trigger calls
    // and the live l1_fit path.  Kept identical so the dump's flag_l1 matches
    // what l1_fit() would have produced if dump_mode were off.
    const TriggerCfg tcfg_dump{
        m_l1_min_length, m_l1_gmax_min, m_l1_energy_frac_thr,
        m_l1_asym_strong, m_l1_asym_mod, m_l1_asym_loose,
        m_l1_len_long_mod, m_l1_len_long_loose, m_l1_len_fill_shape,
        m_l1_fill_shape_fill_thr, m_l1_fill_shape_fwhm_thr,
    };

    for (auto trace : sigtraces) {
        auto newtrace = std::make_shared<Aux::SimpleTrace>(
            trace->channel(), trace->tbin(), trace->charge());

        int ch = trace->channel();

        // Skip channels with no ROIs (includes out-of-scope/ineligible channels
        // that were not added to init_map above).
        if (map_ch_rois.find(ch) == map_ch_rois.end()) {
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
                                        m_adc_l1_threshold,
                                        m_l1_energy_pad_ticks,
                                        m_l1_raw_asym_pad_ticks,
                                        m_l1_raw_asym_eps,
                                        m_l1_core_g_thr);
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

                // Tier 1: flag + ratio (same formula as l1_fit trigger)
                double ratio = (rec.temp1_sum > 0)
                             ? rec.temp_sum / (rec.temp1_sum * m_adc_sum_rescaling / rec.nbin_fit)
                             : 0.0;
                int flag = 0;
                if (rec.temp1_sum > m_adc_sum_threshold) {
                    if      (ratio >  m_adc_ratio_threshold) flag = +1;
                    else if (ratio < -m_adc_ratio_threshold) flag = -1;
                }
                d_flag.push_back(flag);
                d_ratio.push_back(ratio);
                d_temp_sum_pos.push_back(rec.temp_sum_pos);
                d_temp_sum_neg.push_back(rec.temp_sum_neg);
                d_n_above_pos.push_back(rec.n_above_pos);
                d_n_above_neg.push_back(rec.n_above_neg);
                // Tier 2: peak locations and decon scalars
                d_argmax_tick.push_back(rec.argmax_tick);
                d_argmin_tick.push_back(rec.argmin_tick);
                d_sig_peak.push_back(rec.sig_peak);
                d_sig_integral.push_back(rec.sig_integral);
                // Tier 3: shape features driving the trigger gate
                d_gmax.push_back(rec.gmax);
                d_gauss_fill.push_back(rec.gauss_fill);
                d_gauss_fwhm_frac.push_back(rec.gauss_fwhm_frac);
                d_roi_energy_frac.push_back(rec.roi_energy_frac);
                d_raw_asym_wide.push_back(rec.raw_asym_wide);
                d_core_lo.push_back(rec.core_lo);
                d_core_hi.push_back(rec.core_hi);
                d_core_length.push_back(rec.core_length);
                d_core_fill.push_back(rec.core_fill);
                d_core_fwhm_frac.push_back(rec.core_fwhm_frac);
                d_core_raw_asym_wide.push_back(rec.core_raw_asym_wide);
                d_flag_l1.push_back((int32_t)decide_trigger(
                    adctrace->charge(), trace->charge(), trace->tbin(),
                    rois_save[i].first, rois_save[i].second + 1,
                    rec.gmax, tcfg_dump,
                    m_l1_core_g_thr,
                    m_l1_energy_pad_ticks,
                    m_l1_raw_asym_pad_ticks,
                    m_l1_raw_asym_eps));
            }
        } else {
            // Normal processing path: classify ROIs and run L1 fits.
            // Resolve channel → plane once (kernels are indexed by plane).
            const int plane = m_anode ? m_anode->resolve(ch).index() : 0;
            for (auto& roi : rois_save) {
                std::vector<double> lasso_unsmeared_buf;
                const int polarity = l1_fit(newtrace, adctrace, trace,
                                            roi.first, roi.second + 1, plane,
                                            m_wf_dump_path.empty() ? nullptr : &lasso_unsmeared_buf);
                if (!m_wf_dump_path.empty() && polarity != 0) {
                    dump_roi_waveforms(in->ident(), ch, plane,
                                       roi.first, roi.second + 1, polarity,
                                       adctrace, trace, newtrace, lasso_unsmeared_buf);
                }
                // Zero any negative decon values within the ROI.
                for (int t = roi.first; t <= roi.second; t++) {
                    auto& v = newtrace->charge().at(t - trace->tbin());
                    if (v < 0) v = 0;
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
        // Tier 1
        save_i32("flag",          d_flag);
        save_f64("ratio",         d_ratio);
        save_f64("temp_sum_pos",  d_temp_sum_pos);
        save_f64("temp_sum_neg",  d_temp_sum_neg);
        save_i32("n_above_pos",   d_n_above_pos);
        save_i32("n_above_neg",   d_n_above_neg);
        // Tier 2
        save_i32("argmax_tick",   d_argmax_tick);
        save_i32("argmin_tick",   d_argmin_tick);
        save_f64("sig_peak",      d_sig_peak);
        save_f64("sig_integral",  d_sig_integral);
        // Tier 3 — features driving the per-ROI trigger gate
        save_f64("gmax",            d_gmax);
        save_f64("gauss_fill",      d_gauss_fill);
        save_f64("gauss_fwhm_frac", d_gauss_fwhm_frac);
        save_f64("roi_energy_frac", d_roi_energy_frac);
        save_f64("raw_asym_wide",   d_raw_asym_wide);
        save_i32("core_lo",            d_core_lo);
        save_i32("core_hi",            d_core_hi);
        save_i32("core_length",        d_core_length);
        save_f64("core_fill",          d_core_fill);
        save_f64("core_fwhm_frac",     d_core_fwhm_frac);
        save_f64("core_raw_asym_wide", d_core_raw_asym_wide);
        save_i32("flag_l1",         d_flag_l1);

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
