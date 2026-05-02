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
 * Trigger thresholds are wired to default uBooNE values and need re-tuning
 * from PDHD/PDVD dump-mode hand-scan data.  Unipolar response bases are
 * configured via "fields_pos_unipolar" and "fields_neg_unipolar"; if empty
 * the component acts as a pass-through.
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

#include <memory>
#include <set>
#include <vector>

namespace WireCell {
    namespace SigProc {

        class L1SPFilterPD : public Aux::Logger, public IFrameFilter, public IConfigurable {
        public:
            L1SPFilterPD(double gain = 14.0 * units::mV / units::fC,
                         double shaping = 2.2 * units::microsecond,
                         double postgain = 1.2,
                         double ADC_mV = 4096 / (2000. * units::mV),
                         double fine_time_offset = 0.0 * units::microsecond,
                         double coarse_time_offset = -8.0 * units::microsecond);
            virtual ~L1SPFilterPD() = default;

            virtual bool operator()(const input_pointer& in, output_pointer& out);
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;

        private:
            // Build one linterp<double> from the named IFieldResponse, using the
            // response of plane `plane_index` (0=U, 1=V, 2=W) convolved with
            // the electronics response and current gain/shaping settings.
            std::unique_ptr<linterp<double>> build_response(const std::string& fields_tn,
                                                            int plane_index);

            // Populate m_lin_bipolar (always) and m_lin_pos/neg_unipolar (if configured).
            // Called lazily on first operator() invocation.
            void init_resp();

            // Classify one ROI by per-ROI shape features computed from raw ADC
            // (adctrace) and the unmodified gauss decon signal (sigtrace), then
            // if triggered run the LASSO fit on newtrace samples in
            // [start_tick, end_tick).  sigtrace is the original gauss-side trace
            // for this channel; it is read for the wide-window energy fraction
            // and is never mutated.  newtrace holds the writable output and is
            // initially a copy of sigtrace.  Returns: 0=pass-through,
            // +1=L1-positive, -1=L1-negative.
            int l1_fit(std::shared_ptr<WireCell::Aux::SimpleTrace>& newtrace,
                       const std::shared_ptr<const WireCell::ITrace>& adctrace,
                       const std::shared_ptr<const WireCell::ITrace>& sigtrace,
                       int start_tick, int end_tick);

            // True if channel belongs to a plane in m_process_planes.
            // Always returns true when m_process_planes is empty or m_anode is null.
            bool channel_in_scope(int channel) const;

            // True if channel is in m_eligible_channels (or that set is empty).
            bool channel_eligible(int ch) const;

            // Bipolar induction-plane response (always built from "fields" config).
            // Unipolar responses built from "fields_pos_unipolar" / "fields_neg_unipolar"
            // if those config keys are non-empty; otherwise they remain null.
            std::unique_ptr<linterp<double>> m_lin_bipolar;
            std::unique_ptr<linterp<double>> m_lin_pos_unipolar;
            std::unique_ptr<linterp<double>> m_lin_neg_unipolar;

            IDFT::pointer m_dft;
            size_t m_count{0};

            // Detector response parameters
            double m_gain;
            double m_shaping;
            double m_postgain;
            double m_ADC_mV;
            double m_fine_time_offset;
            double m_coarse_time_offset;

            // Cached config scalars (populated in configure())
            std::string m_adctag;
            std::string m_sigtag;
            std::string m_outtag;

            int m_roi_pad{3};
            int m_raw_pad{15};
            double m_raw_ROI_th_nsigma{4};
            double m_raw_ROI_th_adclimit{10};

            double m_overall_time_offset{0};
            // Time offset of the unipolar basis relative to the bipolar basis.
            // Analogous to collect_time_offset in L1SPFilter.  Tune once the
            // unipolar field-response files are available.
            double m_unipolar_time_offset{3.0 * units::microsecond};

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
            // mirroring the chndb-base pattern.  Deconvolved-domain knobs
            // (m_l1_gmax_min, asym ratios, lengths, energy fractions) operate
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
            double m_l1_lambda{5};
            double m_l1_epsilon{0.05};
            double m_l1_niteration{100000};
            double m_l1_decon_limit{100};
            double m_l1_resp_scale{0.5};
            // Output weights for each basis component of the reconstructed signal.
            // basis0 = bipolar (first half of beta), basis1 = unipolar (second half).
            double m_l1_basis0_scale{1.15};
            double m_l1_basis1_scale{0.5};
            double m_peak_threshold{1000};
            double m_mean_threshold{500};
            std::vector<double> m_smearing_vec;
            // Auto-derived smearing kernel config (used when "filter" is empty).
            std::string m_gauss_filter_tn{"HfFilter:Gaus_wide"};
            double m_kernel_threshold{1.0e-3};
            int    m_kernel_max_half{64};
            int    m_kernel_nticks{4096};

            int m_bipolar_plane{1};    // field-response plane index for bipolar basis
            int m_unipolar_plane{1};   // field-response plane index for unipolar bases

            // Field-response type-names stored at configure() time for lazy init_resp()
            std::string m_cfg_fields;
            std::string m_cfg_fields_pos;
            std::string m_cfg_fields_neg;

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
        };

    }  // namespace SigProc
}  // namespace WireCell

#endif
