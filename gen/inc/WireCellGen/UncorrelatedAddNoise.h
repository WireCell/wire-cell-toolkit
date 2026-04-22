// ==============================
// UncorrelatedAddNoise.h
// ==============================
/**
 * @file UncorrelatedAddNoise.h
 * @brief Wire-Cell Toolkit: model-driven, per-channel uncorrelated noise injector.
 *
 * Overview
 * --------
 * This component reads an uncorrelated noise model from JSON (optionally .bz2) and
 * injects noise independently per channel in the frequency domain.
 *
 * The model provides:
 *   - freq axis (freq_ghz)
 *   - avg_mag[w,f] : target mean one-sided FFT magnitude per wire/channel and frequency bin
 *   - optional live_mask[w] : 0/1 mask to disable dead/missing channels
 *
 * Algorithm (high level)
 * ----------------------
 * For each channel (wire) and each positive-frequency bin:
 *   1) Draw a CSCG random coefficient Z ~ CN(0,1)
 *   2) Scale Z so that E|scaled Z| matches the target mean magnitude avg_mag(w,f)
 *      using Rayleigh-mean scaling: E|CN(0,1)| = sqrt(pi)/2
 *   3) Handle DC (force to 0) and Nyquist (real-only if N is even)
 *   4) Convert one-sided spectrum to time-domain noise via FFTW c2r:
 *        - provide the one-sided spectrum (DC..Nyquist)
 *        - FFTW c2r (unnormalized)
 *        - apply 1/N normalization to match NumPy irfft
 *        - apply optional ifft_scale
 *
 * IMPORTANT ASSUMPTION ABOUT CHANNEL ORDERING
 * ------------------------------------------
 * The model rows are assumed to correspond to the *sorted-by-channel* trace order:
 * after sorting frame traces by channel ID:
 *   - model row 0 corresponds to traces[0]
 *   - model row 1 corresponds to traces[1]
 *   - ...
 *
 * @author Avik Ghosh
 * @version 1.0.0
 * @date 2026-02-04
 */

#ifndef WIRECELLGEN_UNCORRELATEDADDNOISE_H
#define WIRECELLGEN_UNCORRELATEDADDNOISE_H

#include "WireCellAux/Logger.h"

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IRandom.h"

#include <Eigen/Dense>

#include <string>
#include <vector>
#include <cstddef>

namespace WireCell {
namespace Gen {

    class UncorrelatedAddNoise : public Aux::Logger,
                                public IFrameFilter,
                                public IConfigurable
    {
      public:
        UncorrelatedAddNoise();
        virtual ~UncorrelatedAddNoise();

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
        // Return a local copy of traces sorted by channel id.
        WireCell::ITrace::vector
        sorted_traces_by_channel(const WireCell::IFrame::pointer& frame) const;

        // ---- Small helpers ----
        // Linear interpolation of avg_mag for wire w at frequency f (in WCT base units).
        double interp_avg_mag_w(int w, double f) const;

        // ---- Main synthesis routine ----
        // Build a [nwires_frame x nsamp] time-domain noise matrix.
        // NOTE: This is a pure helper; it does NOT touch frames/traces.
        Eigen::MatrixXd make_uncorrelated_noise(int nwires_frame, size_t nsamp) const;

      private:
        // ---- Configuration state ----
        std::string m_name{"UncorrelatedAddNoise"};
        std::string m_model_file{"uncorrelated_noise_model.json.bz2"};
        WireCell::IRandom::pointer m_rng{nullptr};

        // Expected number of samples and sample spacing.
        // m_dt is in WCT base time units (ns) after configure().
        size_t m_nsamples{2128};
        double m_dt{0.5 * WireCell::units::us};

        // Optional final amplitude tweak applied after the explicit 1/N normalization
        // in the NumPy-like irfft helper.
        double m_ifft_scale{1.0};

        // ---- Model data ----
        // Frequency axis in WCT base freq units (Hz in WireCell units).
        Eigen::VectorXd m_freq;

        // Target mean one-sided FFT magnitudes: [Nw x Nf] (stored in WCT voltage units)
        Eigen::MatrixXd m_avg_mag;

        // Optional live mask: 0/1 per wire, defaults to all-ones if absent.
        Eigen::VectorXi m_live_mask;

        // ---- Bookkeeping ----
        mutable size_t m_count{0};
    };

} // namespace Gen
} // namespace WireCell

#endif // WIRECELLGEN_UNCORRELATEDADDNOISE_H