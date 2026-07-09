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
 *   4. ROIs (HYSTERESIS) = contiguous runs above roi_ext_nsigma * rms (low) that
 *      reach roi_seed_nsigma * rms somewhere inside (high, ~ the hit threshold so
 *      a real pulse is required), padded roi_pad_pre before the run and extended to
 *      at least roi_post_peak ticks past the decon pulse peak (the late-light
 *      window); a brighter tail above the extend threshold runs further on its own
 *      and overlapping windows merge.  This keeps the slow scintillation tail
 *      inside the ROI without the 60-75 us ballooning of a fixed blanket pad;
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

#include <set>
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
            // Hysteresis thresholds in units of (1.4826 * MAD of the HPF wave).
            // seed (high): a run must reach this to be a real ROI; 5 -> ~0.10
            // decon at the ~0.02 clean noise floor ~ the OpHit 0.11 hit threshold,
            // so dim ~1-PE pulses still seed but pure noise does not.
            // ext (low): the ROI spans the contiguous run down to this (3 -> ~0.06
            // decon), capturing the rise/tail so the linear baseline lands near
            // true baseline.
            double m_roi_seed_nsigma{5.0};
            double m_roi_ext_nsigma{3.0};
            int m_roi_pad_pre{50};     // pad before the run start (~0.8 us); no early extension
            // ROI reaches at least this many ticks past the DECON pulse peak (the
            // late-light window; 300 ~ 4.8 us ~ 3x the 1.6 us LAr late-light tau).
            // A brighter tail that stays above the extend threshold runs further on
            // its own; overlapping windows merge adjacent pulses.
            int m_roi_post_peak{300};
            // MAD cap (decon units) above which a channel is treated as ringing
            // and zeroed entirely (carries the OpHitFinder full-stream veto).
            double m_veto_sigma{0.1};
            // Hard per-channel veto list (OpChannel ids): these channels are zeroed
            // unconditionally, independent of the MAD test above.  For known-bad
            // data-quality channels that a hand scan flags but whose MAD does not
            // always clear m_veto_sigma.  Default empty -> bit-identical to every
            // existing config.
            std::set<int> m_veto_channels{};
            // Per-ROI linear endpoint-zeroing.  When false the in-ROI decon is
            // kept verbatim (only the outside-ROI zeroing applies).
            bool m_apply_baseline{true};
            // Perf knob, default OFF -> bit-identical: run the HPF transforms
            // as true real FFTs (r2c/c2r plans), same results up to round-off.
            bool m_use_real_dft{false};

            IDFT::pointer m_dft;

            // High-pass multiplier per full-size FFT bin, built for the trace
            // length actually seen (rebuilt only if the length changes).
            std::vector<double> m_hpf;
            int m_hpf_n{0};
            void ensure_hpf(int n);

            // ROI-clean one deconvolved waveform.
            std::vector<float> clean(const std::vector<float>& decon);

            // Reused per-trace scratch (median input copies), avoiding two
            // full-length allocations per full-stream trace.
            std::vector<float> m_scratch;
            std::vector<float> m_dev;
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
