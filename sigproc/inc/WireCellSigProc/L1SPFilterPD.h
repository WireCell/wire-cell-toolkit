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
 * STUB STATUS: trigger classification always returns 0 (pass-through).
 * Fill in classify_roi_flag() in L1SPFilterPD.cxx after dedicated event-scan
 * analysis.  Unipolar response bases are configured via "fields_pos_unipolar"
 * and "fields_neg_unipolar"; if empty the component is a no-op.
 * See sigproc/docs/l1sp/README.md, Strategy B.
 */
#ifndef WIRECELLSIGPROC_L1SPFILTERPD
#define WIRECELLSIGPROC_L1SPFILTERPD

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDFT.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Interpolate.h"

#include <memory>
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

            // Classify one ROI and, if triggered, run the LASSO fit on
            // newtrace samples in [start_tick, end_tick).
            //
            // hint_polarity: 0=auto-classify, 1=force positive, -1=force negative.
            // Returns: 0=pass-through, 1=L1-positive, -1=L1-negative, 2=zero-out.
            int l1_fit(std::shared_ptr<WireCell::Aux::SimpleTrace>& newtrace,
                       const std::shared_ptr<const WireCell::ITrace>& adctrace,
                       int start_tick, int end_tick,
                       int hint_polarity = 0);

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
            double m_adc_sum_rescaling_limit{50.0};
            // Asymmetry ratio threshold.  The trigger fires when
            //   |temp_sum| / (temp1_sum * rescaling / nbin) > adc_ratio_threshold
            // with sign(temp_sum) selecting positive vs negative polarity.
            double m_adc_ratio_threshold{0.2};

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

            int m_bipolar_plane{1};    // field-response plane index for bipolar basis
            int m_unipolar_plane{1};   // field-response plane index for unipolar bases

            // Field-response type-names stored at configure() time for lazy init_resp()
            std::string m_cfg_fields;
            std::string m_cfg_fields_pos;
            std::string m_cfg_fields_neg;
        };

    }  // namespace SigProc
}  // namespace WireCell

#endif
