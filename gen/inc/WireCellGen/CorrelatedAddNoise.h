// ==============================
// CorrelatedAddNoise.h
// ==============================
/**
 * @file CorrelatedAddNoise.h
 * @brief Wire-Cell Toolkit: model-driven, banded, inter-channel correlated noise injector.
 *
 * Overview
 * --------
 * This component reads a correlated-noise model from JSON (optionally .bz2) and
 * injects noise in the frequency domain while imposing cross-wire correlations
 * in user-defined frequency bands.
 *
 * The model provides:
 *   - freq axis (freq_ghz)
 *   - avg_mag[w,f] : target mean one-sided FFT magnitude per wire/channel and frequency bin
 *   - bands_khz[b] : list of [lo, hi] band edges in kHz
 *   - A_band[b]    : per-band "colorer" matrix such that A*A^T ≈ correlation matrix
 *   - optional live_mask[w]
 *
 * Algorithm (high level)
 * ----------------------
 * For each positive-frequency bin f_k:
 *   1) Determine the band b that contains f_k
 *   2) Draw Z ~ CN(0,I) across wires
 *   3) Color: U = A_b Z   (imposes cross-wire covariance)
 *   4) Per wire, scale U(w) so E|scaled U(w)| matches avg_mag(w,f_k)
 *      using Rayleigh-mean scaling and the row-wise variance of A_b
 *   5) Handle DC (force to 0) and Nyquist (real-only if N even, correlated draw)
 *   6) Convert one-sided spectrum to time domain using NumPy-like irfft convention:
 *        - bin-doubling + 1/N + optional ifft_scale
 *
 * IMPORTANT ASSUMPTION ABOUT CHANNEL ORDERING
 * ------------------------------------------
 * The model rows are assumed to correspond to the *sorted-by-channel* trace order:
 * after sorting frame traces by channel ID:
 *   - model row 0 corresponds to traces[0]
 *   - model row 1 corresponds to traces[1]
 *   - ...
 *
 *
 * @author Avik Ghosh
 * @version 1.0.0
 * @date 2026-02-04
 */

#ifndef WIRECELLGEN_CORRELATEDADDNOISE_H
#define WIRECELLGEN_CORRELATEDADDNOISE_H

#include "WireCellAux/Logger.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IRandom.h"

#include <Eigen/Dense>

#include <string>
#include <vector>
#include <cstddef>

namespace WireCell {
namespace Gen {

    class CorrelatedAddNoise : public Aux::Logger,
                              public IFrameFilter,
                              public IConfigurable
    {
      public:
        CorrelatedAddNoise();
        virtual ~CorrelatedAddNoise();

        // IConfigurable interface
        virtual WireCell::Configuration default_configuration() const override;
        virtual void configure(const WireCell::Configuration& cfg) override;

        // IFrameFilter interface
        virtual bool operator()(const input_pointer& inframe,
                                output_pointer& outframe) override;

        // INamed (provided via INode in the IFrameFilter chain)
        virtual void set_name(const std::string& name) override;
        virtual std::string get_name() const override;

      private:
        // ---- Model loading/parsing ----
        void load_model(const std::string& fname);

        // ---- Frame utilities ----
        WireCell::ITrace::vector
        sorted_traces_by_channel(const WireCell::IFrame::pointer& frame) const;

        // ---- Small helpers ----
        // Return the index of the correlation band that contains frequency f, or -1 if none.
        int band_for_freq(double f) const;

        // Linear interpolation of avg_mag for wire w at frequency f (in WCT base units).
        double interp_avg_mag_w(int w, double f) const;

        // ---- Main synthesis routine ----
        // Build a [nwires_frame x nsamp] time-domain correlated noise matrix.
        Eigen::MatrixXd make_correlated_noise(int nwires_frame, size_t nsamp) const;

      private:
        // ---- Configuration state ----
        std::string m_name{"CorrelatedAddNoise"};
        std::string m_model_file{"correlated_noise_model.json.bz2"};
        WireCell::IRandom::pointer m_rng{nullptr};

        // Expected number of samples and sample spacing.
        // m_dt is in WCT base time units (ns) after configure().
        size_t m_nsamples{2128};
        double m_dt{0.5 * WireCell::units::us};

        // Optional final amplitude tweak applied after explicit 1/N normalization.
        double m_ifft_scale{1.0};

        // ---- Model data ----
        Eigen::VectorXd m_freq;                 // WCT base freq units
        Eigen::MatrixXd m_avg_mag;              // [Nw x Nf] in WCT voltage units
        Eigen::VectorXi m_live_mask;            // [Nw] 0/1
        std::vector<Eigen::MatrixXd> m_A_band;  // per band [Nw x Nw]

        // ---- Derived band info ----
        struct BandInfo {
            double lo_f{0.0};       // lower band edge in WCT Hz
            double hi_f{0.0};       // upper band edge in WCT Hz
            size_t band_index{0};   // index into m_A_band
            Eigen::VectorXd diagC;  // row-wise sum of squares of A (diag of A*A^T)
        };
        std::vector<BandInfo> m_bands;

        // ---- Bookkeeping ----
        mutable size_t m_count{0};
    };

} // namespace Gen
} // namespace WireCell

#endif // WIRECELLGEN_CORRELATEDADDNOISE_H