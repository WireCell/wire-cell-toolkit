/** Optical waveform deconvolution for PDHD light reconstruction.
 *
 * A dependency-free port of the duneopdet `Deconvolution` module
 * (Wiener filter built from a single-p.e. template, optional Gauss
 * post-filter and post baseline correction), operating on tagged
 * optical-snippet traces of an IFrame.  See
 * flash/docs/stage2-reconstruction.md for the fcl correspondence.
 *
 * Input: traces tagged `intag` (raw ADC snippets, channel = OpChannel).
 * Output: the input frame plus SPE-normalized traces tagged `outtag`.
 */
#ifndef WIRECELLFLASH_OPDECON
#define WIRECELLFLASH_OPDECON

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDFT.h"
#include "WireCellAux/Logger.h"

#include <map>
#include <vector>

namespace WireCell {
    namespace Flash {
        class OpDecon : public Aux::Logger, public IFrameFilter, public IConfigurable {
          public:
            OpDecon();
            virtual ~OpDecon();

            virtual bool operator()(const IFrame::pointer& in, IFrame::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            std::string m_intag{"raw"};
            std::string m_outtag{"decon"};
            // JSON file with SPE templates and the channel -> template map.
            std::string m_spe_file{"pgrapher/experiment/pdhd/pdhd-spe-templates.json"};
            // Optional JSON file with per-channel noise power spectra
            // (half-spectrum |FFT|^2, LArSoft NoiseTemplateFiles).  Empty
            // means the flat LineNoiseRMS^2 * Samples default.
            std::string m_noise_file{""};
            // protodunehd_deconvolution fcl values.
            int m_samples{1024};
            int m_pre_trigger{50};
            int m_pedestal_buffer{30};
            double m_line_noise_rms{4.5};
            int m_input_polarity{-1};
            bool m_auto_scale{true};       // Wiener: scale from filtered SPE response
            double m_scale{1.0};           // used when auto_scale is false
            // Fixed Wiener signal-to-noise ratio R = S2/N^2.  By default
            // (<= 0) S2 is taken from each waveform's own peak (adaptive,
            // signal- and length-dependent).  When > 0, S2 = R * N^2 so the
            // filter G = conj(H) R / (|H|^2 R + 1) is independent of signal
            // amplitude AND record length -- the same filter for the 1024-tick
            // snippets and the 343808-tick full stream.  Assumes flat noise.
            double m_fixed_snr{-1.0};
            // Wiener-INSPIRED mode.  Default OFF -> bit-identical to every
            // existing config.  When ON the deconvolution spectrum becomes
            // the pure inverse times a band filter pinned to 1 at zero
            // frequency:
            //   G = conj(H) F(f) / (|H|^2 + eps),  F(f) = exp(-0.5 (f/sigma)^p)
            // with eps = (wi_eps_rel * max_k|H_k|)^2 guarding template
            // spectral nulls.  F(0) = 1 makes the deconvolved 1-PE area
            // exactly 1 independent of the filter parameters, so AutoScale
            // is skipped (scale = `scale`, default 1) and fixed_snr /
            // line_noise_rms are unused.  The Gauss post-filter should be
            // disabled in this mode (F *is* the band filter).  See
            // pdvd/docs/pdvd-light-filter.md.
            bool m_wiener_inspired{false};
            double m_wi_sigma_mhz{1.0};   // F half form: exp(-0.5 (f/sigma)^power)
            double m_wi_power{2.0};
            double m_wi_eps_rel{1e-3};
            bool m_apply_postfilter{true};
            double m_postfilter_cutoff{1.5};  // MHz, Gauss post-filter
            bool m_apply_post_blcorr{true};
            // Perf knobs, default OFF -> bit-identical to every existing
            // config.  use_real_dft: run the transforms as true real FFTs
            // (FFTW r2c/c2r plans) -- same results up to round-off (~1e-7),
            // about half the FFT flops and transform memory.
            // fold_postfilter: multiply the Gauss post-filter into the
            // deconvolution spectrum instead of a separate fwd+inv pair per
            // record.  The post baseline correction keeps its exact
            // semantics (mean of the UNfiltered decon head): that mean is a
            // linear functional of the spectrum and is evaluated with the
            // precomputed m_ped_w weights, so only round-off differs.
            bool m_use_real_dft{false};
            bool m_fold_postfilter{false};
            // 14-bit ADC saturation detection.  Default OFF -> the output frame
            // carries no masks, bit-identical to every existing config.  When ON:
            // any input snippet with >= saturation_min_samples raw samples at or
            // above saturation_adc (the DAPHNE 14-bit rail, 16383) is flagged by
            // adding its [tbin, tbin+len) range to a "saturation" ChannelMaskMap
            // on the output frame.  Deconvolving a clipped flat-top yields a broad
            // plateau that OpHitFinder over-integrates (one ~16 us hit) and
            // fragments (many spurious wide hits); OpHitFinder's veto_saturation
            // consumes this mask to drop such snippets.  See
            // pdhd/docs/run29107-evt1015-light-anomaly.md.
            bool m_detect_saturation{false};
            int m_saturation_adc{16383};
            int m_saturation_min_samples{1};
            // Pad (ticks) added each side of a railed run.  Deconvolving a
            // clipped flat-top yields a plateau wider than the clip itself, so
            // the over-integrated hits extend beyond the railed samples; the pad
            // widens the vetoed range to cover them.
            int m_saturation_pad{0};

            IDFT::pointer m_dft;

            // Per-template precomputed state.  The spectrum is built lazily
            // on a template's first use: a branch instance configured with
            // the full template file only pays the (m_samples-long) FFT and
            // spectrum storage for the channels it actually deconvolves.
            struct SPETemplate {
                std::vector<float> wave;            // time domain, unpadded (<= m_samples)
                std::vector<std::complex<float>> fft;  // full-size spectrum, lazy
                double amplitude;                   // max(wave), >= 1
                double wi_eps{0.0};                 // (wi_eps_rel * max|fft|)^2, lazy
                // AutoScale normalization cached when the Wiener filter is
                // record-independent (fixed_snr > 0, flat noise): the filtered
                // SPE response then depends only on the template, so the extra
                // inverse FFT per record collapses to one per template.
                double cached_scale{0.0};
                bool scale_cached{false};
            };
            std::vector<SPETemplate> m_templates;
            void ensure_fft(SPETemplate& spe);
            double auto_scale(const SPETemplate& spe,
                              const std::vector<std::complex<float>>& xG) const;
            // Transform via the complex (default) or real (use_real_dft) path.
            std::vector<std::complex<float>> dft_fwd(const std::vector<float>& wave) const;
            std::vector<float> dft_inv(const std::vector<std::complex<float>>& spec) const;
            std::map<int, size_t> m_chan2tmpl;
            std::vector<std::complex<float>> m_postfilter;  // full-size spectrum
            std::vector<double> m_wi_filter;  // full-size F(f), wiener_inspired
            // Spectral weights for the head-pedestal mean, used by
            // fold_postfilter: mean_{i<nped} x_i = Re sum_k X_k m_ped_w[k].
            std::vector<std::complex<double>> m_ped_w;

            // Noise power spectra, half-spectrum (samples/2+1) bins.
            std::vector<std::vector<double>> m_noise_templates;
            std::map<int, size_t> m_chan2noise;

            std::vector<float> deconvolve(const std::vector<float>& adc, const SPETemplate& spe,
                                          const std::vector<double>* noise) const;
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
