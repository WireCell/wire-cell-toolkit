/** ROI identification + ROI-based cleaning for the PDHD full-stream decon.
 *
 * The -x full-stream decon (opch 120-159, ~343808 ticks @16ns) is one long
 * record dominated by baseline between sparse scintillation pulses, with
 * per-channel DC offsets and a ringing channel, so a threshold OpHit finder
 * over-produces.  This component cleans the deconvolved waveform fed to
 * OpHitFinder, mirroring the WireCell charge signal-processing induction ROI
 * chain (sigproc/src/ROI_formation.cxx):
 *
 *   1. high-pass filter H(f) = 1 - exp(-(f/tau)^2) (the sigproc LfFilter form;
 *      scintillation is < 20 us so the corner tau ~ 1/20us = 0.05 MHz passes
 *      the pulses and removes slower baseline wander) -> a "ROI-finding"
 *      waveform;
 *   2. subtract its median (baseline);
 *   3. rms = 1.4826 * MAD; channels with rms > veto_sigma (ringing) are zeroed;
 *   4. ROIs = contiguous runs above roi_nsigma * rms, padded (pad_pre, pad_post)
 *      and merged;
 *   5. apply the ROIs to the ORIGINAL decon: zero everything outside ROIs and,
 *      per ROI [s,e], subtract the line through (s,d[s]),(e,d[e]) so the ROI
 *      starts and ends exactly at zero (linear baseline correction).
 *
 * The prototype + tuning live in pdhd/pd_plot/fullstream_roi_proto.py; see
 * pdhd/docs/pdhd-fullstream-light-reco.md.
 *
 * Input:  traces tagged `intag` (deconvolved, channel = OpChannel).
 * Output: the input frame plus cleaned traces tagged `outtag`.
 */
#ifndef WIRECELLFLASH_OPROI
#define WIRECELLFLASH_OPROI

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IDFT.h"
#include "WireCellAux/Logger.h"

#include <vector>

namespace WireCell {
    namespace Flash {
        class OpRoi : public Aux::Logger, public IFrameFilter, public IConfigurable {
          public:
            OpRoi();
            virtual ~OpRoi();

            virtual bool operator()(const IFrame::pointer& in, IFrame::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            std::string m_intag{"decon"};
            std::string m_outtag{"decon_roi"};
            // High-pass corner of H(f) = 1 - exp(-(f/tau)^2), MHz.  0.05 ~ 1/20us
            // (scintillation is < 20 us; slower wander is removed).
            double m_hpf_tau_mhz{0.05};
            // ROI start threshold = roi_nsigma * (1.4826 * MAD of the HPF wave).
            // 5 -> ~0.10 decon at the ~0.02 clean noise floor: low but above the
            // electronic noise (and ~ the OpHit 0.11 hit threshold).
            double m_roi_nsigma{5.0};
            int m_roi_pad_pre{50};     // ticks padded before each ROI (~0.8 us)
            int m_roi_pad_post{700};   // ticks padded after  each ROI (~11 us tail)
            // MAD cap (decon units) above which a channel is treated as ringing
            // and zeroed entirely (carries the OpHitFinder full-stream veto).
            double m_veto_sigma{0.1};
            // Per-ROI linear endpoint-zeroing.  When false the in-ROI decon is
            // kept verbatim (only the outside-ROI zeroing applies).
            bool m_apply_baseline{true};

            IDFT::pointer m_dft;

            // High-pass multiplier per full-size FFT bin, built for the trace
            // length actually seen (rebuilt only if the length changes).
            std::vector<double> m_hpf;
            int m_hpf_n{0};
            void ensure_hpf(int n);

            // ROI-clean one deconvolved waveform.
            std::vector<float> clean(const std::vector<float>& decon) const;
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
