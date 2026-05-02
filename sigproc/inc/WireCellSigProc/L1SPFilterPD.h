/** L1-norm sparse deconvolution for PDHD/PDVD induction-plane channels.
 *
 * Handles unipolar induction-plane signals caused by:
 *   - anode-induction: electrons only leaving  → negative unipolar raw signal
 *   - collection-on-induction: electrons only arriving → positive unipolar raw signal
 *
 * Strategy B (per-ROI polarity detection): each ROI is classified by the
 * signed/unsigned ADC ratio.  Triggered ROIs are processed with a 2-column
 * LASSO using {bipolar, +unipolar} (positive case) or {bipolar, -unipolar}
 * (negative case) response bases.
 *
 * ┌──────────────────────────────────────────────────────────────────────────┐
 * │ Response kernels are NOT computed in C++.  They are pre-built offline by │
 * │ the wire-cell-python tool and shipped as a JSON+bz2 file (one per        │
 * │ detector, e.g. wire-cell-data/pdhd_l1sp_kernels.json.bz2):               │
 * │                                                                          │
 * │   wirecell-sigproc gen-l1sp-kernels \                                    │
 * │     --gain '14*mV/fC' --shaping '2.2*us' \                               │
 * │     --postgain 1.2 --adc-per-mv 2.048 \                                  │
 * │     --coarse-time-offset '-8*us' \                                       │
 * │     <field-response.json.bz2>  <out>_l1sp_kernels.json.bz2               │
 * │                                                                          │
 * │ See wirecell/sigproc/l1sp.py for the schema and the reference            │
 * │ algorithm.  pdhd/nf_plot/track_response_l1sp_kernels.py validates a      │
 * │ generated file (--from-file) by re-building from FR and asserting        │
 * │ bitwise agreement.                                                       │
 * └──────────────────────────────────────────────────────────────────────────┘
 *
 * The file provides per induction plane (U=0, V=1):
 *   - positive case: bipolar kernel + W-plane unipolar kernel + per-plane
 *     unipolar_time_offset (W peak ↔ this plane's bipolar zero crossing).
 *   - negative case: bipolar kernel + neg-half(bipolar) unipolar, no shift.
 *
 * See sigproc/docs/l1sp/README.md, Strategy B.
 */
#ifndef WIRECELLSIGPROC_L1SPFILTERPD
#define WIRECELLSIGPROC_L1SPFILTERPD

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDFT.h"
#include "WireCellIface/IAnodePlane.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Interpolate.h"

#include <map>
#include <memory>
#include <set>
#include <vector>

namespace WireCell {
    namespace SigProc {

        class L1SPFilterPD : public Aux::Logger, public IFrameFilter, public IConfigurable {
        public:
            L1SPFilterPD();
            virtual ~L1SPFilterPD() = default;

            virtual bool operator()(const input_pointer& in, output_pointer& out);
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

        private:
            // Lazy-load the per-plane bipolar and unipolar interpolators from
            // the pre-built JSON+bz2 kernel file (see header comment).
            void init_resp();

            // Classify one ROI by per-ROI shape features computed from raw ADC
            // (adctrace) and the unmodified gauss decon signal (sigtrace), then
            // if triggered run the LASSO fit on newtrace samples in
            // [start_tick, end_tick).  ``plane`` is the induction-plane index
            // (0=U, 1=V) used to select the appropriate kernel set.  Returns:
            // 0=pass-through, +1=L1-positive, -1=L1-negative.
            int l1_fit(std::shared_ptr<WireCell::Aux::SimpleTrace>& newtrace,
                       const std::shared_ptr<const WireCell::ITrace>& adctrace,
                       const std::shared_ptr<const WireCell::ITrace>& sigtrace,
                       int start_tick, int end_tick, int plane,
                       std::vector<double>* lasso_unsmeared = nullptr);

            // Write per-ROI waveform NPZ when m_wf_dump_path is non-empty.
            void dump_roi_waveforms(int frame_ident, int channel, int plane,
                                    int start_tick, int end_tick, int polarity,
                                    const std::shared_ptr<const WireCell::ITrace>& adctrace,
                                    const std::shared_ptr<const WireCell::ITrace>& sigtrace,
                                    const std::shared_ptr<WireCell::Aux::SimpleTrace>& newtrace,
                                    const std::vector<double>& lasso_unsmeared);

            // True if channel belongs to a plane in m_process_planes.
            // Always returns true when m_process_planes is empty or m_anode is null.
            bool channel_in_scope(int channel) const;

            // True if channel is in m_eligible_channels (or that set is empty).
            bool channel_eligible(int ch) const;

            // Per-plane interpolators loaded from the kernel file (m_kernels_file).
            // Keys are induction-plane indices (typically 0=U, 1=V).
            //
            //   m_lin_bipolar[plane]      — plane bipolar kernel (same kernel
            //                               used for positive and negative cases).
            //   m_lin_pos_unipolar[plane] — collection (W) kernel; queried at
            //                               (t − m_unipolar_toff_pos[plane]) so
            //                               its peak lands at the bipolar zero
            //                               crossing.  The stored offset is
            //                               (zero_crossing − W_peak); subtracting
            //                               it puts W_peak at LASSO dt =
            //                               zero_crossing.
            //   m_lin_neg_unipolar[plane] — neg-half(bipolar); queried at t (no shift).
            std::map<int, std::unique_ptr<linterp<double>>> m_lin_bipolar;
            std::map<int, std::unique_ptr<linterp<double>>> m_lin_pos_unipolar;
            std::map<int, std::unique_ptr<linterp<double>>> m_lin_neg_unipolar;
            // Per-plane W-shift in WCT time units (positive case).  Negative
            // case has no shift (always 0) — handled inline in l1_fit().
            std::map<int, double> m_unipolar_toff_pos;
            // Amplitude multiplier applied to every loaded kernel sample at
            // init_resp() time.  Kernels in ``kernels_file`` are in ADC/electron
            // at the reference 14 mV/fC FE gain; set this to
            //   params.elec.gain / (14*mV/fC)
            // (i.e. the same gain_scale used for ADC-domain thresholds in
            // sp.jsonnet) so the LASSO bases track the runtime FE gain.
            // Default 1.0 = no scaling (reference gain).
            double m_kernels_scale{1.0};

            IDFT::pointer m_dft;
            size_t m_count{0};

            // Cached config scalars (populated in configure())
            std::string m_adctag;
            std::string m_sigtag;
            std::string m_outtag;

            int m_roi_pad{3};
            int m_raw_pad{15};
            double m_raw_ROI_th_nsigma{4};
            double m_raw_ROI_th_adclimit{10};

            // Global LASSO frame origin loaded from kernel file
            // meta.frame_origin_us (= V-plane bipolar zero-crossing in
            // kernel native time).  In WCT time units (ns).  Used uniformly
            // for all induction planes as the "source signal at t = 0"
            // anchor; build_G calls pass (m_frame_origin + m_overall_time_offset)
            // as the overall_time_offset argument.
            double m_frame_origin{0};

            // Additive cfg override on top of m_frame_origin.  Default 0;
            // typically only used for tuning/diagnostics.
            double m_overall_time_offset{0};

            double m_adc_l1_threshold{6};
            double m_adc_sum_threshold{160};
            double m_adc_sum_rescaling{90.0};
            // Legacy uBooNE asymmetry-ratio knob; kept for the dump-record
            // 'flag' / 'ratio' fields but no longer drives the trigger.  See
            // m_l1_asym_strong / m_l1_asym_mod / m_l1_asym_loose below.
            double m_adc_ratio_threshold{0.2};

            // ── Per-ROI trigger gate (PDHD/PDVD Strategy B, retuned) ─────────
            // Convention: m_l1_raw_asym_eps and the adc_* gates above are in
            // raw ADC counts at the 14 mV/fC reference FE gain.  When the
            // detector runs at a different gain, the jsonnet multiplies these
            // by gain_scale = params.elec.gain / 14 (cfg/.../pdhd/sp.jsonnet),
            // mirroring the chndb-base pattern.  The response kernel amplitudes
            // (m_lin_bipolar/pos_unipolar/neg_unipolar) are scaled by
            // m_kernels_scale (= gain_scale) at init_resp() time.
            // Deconvolved-domain knobs (m_l1_gmax_min, asym ratios, lengths,
            // energy fractions) operate
            // on gain-normalised signals and do NOT scale with FE gain.
            //
            // The trigger requires three pre-filters
            //   nbin_fit         >= m_l1_min_length
            //   gmax             >= m_l1_gmax_min
            //   roi_energy_frac  >= m_l1_energy_frac_thr
            // followed by ANY of four arms:
            //   |raw_asym_wide| >= m_l1_asym_strong, OR
            //   nbin_fit >= m_l1_len_long_mod   AND |raw_asym_wide| >= m_l1_asym_mod, OR
            //   nbin_fit >= m_l1_len_long_loose AND |raw_asym_wide| >= m_l1_asym_loose, OR
            //   nbin_fit >= m_l1_len_fill_shape AND gauss_fill <= m_l1_fill_shape_fill_thr
            //                                   AND gauss_fwhm_frac <= m_l1_fill_shape_fwhm_thr
            //                                   AND |raw_asym_wide| >= m_l1_asym_mod.
            // Polarity is sign(raw_asym_wide).  Defaults seeded from the iter-7
            // offline detector (find_long_decon_artifacts.py).
            int    m_l1_min_length{30};
            double m_l1_gmax_min{1500.0};
            double m_l1_energy_frac_thr{0.66};
            int    m_l1_energy_pad_ticks{500};
            int    m_l1_raw_asym_pad_ticks{20};
            double m_l1_raw_asym_eps{20.0};
            // Per-tick |gauss| gate defining the core sub-window used by the
            // trigger feature computation.  Matches iter-7 g_thr=50.
            double m_l1_core_g_thr{50.0};
            double m_l1_asym_strong{0.65};
            double m_l1_asym_mod{0.40};
            double m_l1_asym_loose{0.30};
            int    m_l1_len_long_mod{100};
            int    m_l1_len_long_loose{200};
            int    m_l1_len_fill_shape{50};
            double m_l1_fill_shape_fill_thr{0.38};
            double m_l1_fill_shape_fwhm_thr{0.30};

            double m_l1_seg_length{120};
            double m_l1_scaling_factor{500};
            double m_l1_lambda{10};        // doubled from 5 to compensate for 2× G (resp_scale 0.5→1)
            double m_l1_epsilon{0.05};
            double m_l1_niteration{100000};
            double m_l1_decon_limit{100};
            double m_l1_resp_scale{1.0};   // must equal 1.0 when kernels are in ADC/electron
            // Output weights for each basis component of the reconstructed signal.
            // basis0 = bipolar (first half of beta), basis1 = unipolar (second half).
            // Default 1.0: both components are in electrons; override to down-weight one.
            double m_l1_basis0_scale{1.0};
            double m_l1_basis1_scale{1.0};
            double m_peak_threshold{1000};
            double m_mean_threshold{500};
            std::vector<double> m_smearing_vec;
            // Auto-derived smearing kernel config (used when "filter" is empty).
            std::string m_gauss_filter_tn{"HfFilter:Gaus_wide"};
            double m_kernel_threshold{1.0e-3};
            int    m_kernel_max_half{64};
            int    m_kernel_nticks{4096};

            // Path to the JSON+bz2 file holding the pre-built L1SP response
            // kernels.  Resolved via WIRECELL_PATH.  See class header comment
            // for the wirecell-sigproc gen-l1sp-kernels invocation.
            std::string m_kernels_file;

            // Anode + plane-scope filter
            std::string m_cfg_anode;
            IAnodePlane::pointer m_anode;
            std::vector<int> m_process_planes{0, 1};  // 0=U, 1=V; skip W

            // Optional eligibility whitelist (channel IDs). Empty = all in scope.
            std::set<int> m_eligible_channels;

            // Calibration dump mode
            bool m_dump_mode{false};
            std::string m_dump_path;
            std::string m_dump_tag;
            // Waveform dump: write per-triggered-ROI NPZ files with the four
            // waveforms (raw, decon, lasso_unsmeared, lasso_smeared) when
            // non-empty.  Independent of m_dump_mode.  Default "" = OFF.
            std::string m_wf_dump_path;
        };

    }  // namespace SigProc
}  // namespace WireCell

#endif
