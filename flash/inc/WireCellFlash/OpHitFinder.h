/** Optical hit finding for PDHD light reconstruction.
 *
 * A dependency-free port of the essentials of the duneopdet
 * `OpHitFinderDeco` module in its PDHD configuration
 * (dune_ophit_finder_deco): larana `PedAlgoEdges` (head method) +
 * `AlgoSlidingWindow` pulse finding on the (scaled, short-cast)
 * deconvolved snippets, PE = pulse area / SPEArea.
 *
 * Input: IFrame with traces tagged `intag` (SPE-normalized snippets).
 * Output: ITensorSet with one "ophits" tensor, f8 [nhit, 9]:
 *   channel, peak_time_ns, width_ns, area, amplitude, pe,
 *   start_time_ns, flash_id (-1), fast_to_total (0)
 * (see flash/docs/design.md §3.4; times trigger-relative WCT ns).
 */
#ifndef WIRECELLFLASH_OPHITFINDER
#define WIRECELLFLASH_OPHITFINDER

#include "WireCellIface/IFrameTensorSet.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <vector>

namespace WireCell {
    namespace Flash {
        class OpHitFinder : public Aux::Logger, public IFrameTensorSet, public IConfigurable {
          public:
            OpHitFinder();
            virtual ~OpHitFinder();

            virtual bool operator()(const input_pointer& in, output_pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

            // One reconstructed pulse (larana pmtana::pulse_param essentials).
            struct Pulse {
                int t_start{0}, t_end{0}, t_max{0};
                double peak{0}, area{0};
            };
            // Port of AlgoSlidingWindow::RecoPulse with the constant
            // (head-estimated) pedestal of PedAlgoEdges.  Public and
            // static for unit testing.
            static std::vector<Pulse> sliding_window(const std::vector<short>& wf,
                                                     double ped_mean, double ped_sigma,
                                                     const Configuration& pars);

            // Post-process one (possibly merged) pulse, splitting it at
            // prominent valleys between sub-peaks so two optical pulses
            // riding a common elevated scintillation tail become separate
            // hits.  Returns {pulse} unchanged when splitting is disabled
            // (pars["split_enable"] false) or nothing qualifies.  Public
            // and static for unit testing.
            static std::vector<Pulse> split_pulse(const std::vector<short>& wf,
                                                  double ped_mean, const Pulse& pulse,
                                                  const Configuration& pars);

          private:
            std::string m_intag{"decon"};
            // dune_ophit_finder_deco values.
            double m_scale{100.0};       // ScalingFactor applied before short cast
            double m_spe_area{100.0};    // SPEArea: PE = area / spe_area
            double m_hit_threshold{3.0}; // HitThreshold on pulse peak (scaled units)
            int m_ped_nsamples{3};       // PedAlgoEdges NumSampleFront (head method)
            // Robust per-channel baseline for the CONTINUOUS full stream, where
            // the head method (a few leading samples) is meaningless.  Default
            // OFF -> the head method above, bit-identical to the self-trigger
            // (snippet) path and every existing config.  When ON: ped_mean =
            // per-channel median, ped_sigma = per-channel MAD (both over the
            // whole waveform), the sliding-window start gate is raised to
            // robust_nsigma * MAD, and channels whose MAD exceeds
            // robust_veto_sigma (ringing data-quality channels) emit no hits.
            // See pdhd/docs/pdhd-fullstream-light-reco.md.
            bool m_robust_baseline{false};
            double m_robust_nsigma{3.0};       // start-gate nsigma on the robust MAD
            double m_robust_veto_sigma{10.0};  // veto channel if MAD >= this (scaled units)
            // Fixed noise floor (scaled units) for ROI-cleaned input (OpRoi).  The
            // hysteresis ROIs hug each pulse, so the in-ROI samples are
            // signal-dominated and a median/MAD over them would bias ped_sigma
            // high and close the start gate (clean-channel hits crater to ~45%).
            // When > 0 this overrides robust/head: ped_mean = 0 (the ROIs are
            // endpoint-zeroed) and ped_sigma = this KNOWN clean noise floor (the
            // HPF rms ~0.02 decon -> ~2 scaled), with the start gate at
            // robust_nsigma * this.  Ringing channels are already zeroed by OpRoi,
            // so no separate veto is needed.  Default 0 (off) -> bit-identical to
            // every existing config.  See pdhd/docs/pdhd-fullstream-light-reco.md.
            double m_fixed_ped_sigma{0.0};
            // Drop snippets flagged as ADC-saturated by OpDecon (the
            // "saturation" ChannelMaskMap on the input frame).  A clipped
            // flat-top deconvolves into a broad plateau that over-integrates
            // into one giant ~16 us hit and fragments into many spurious wide
            // hits; skipping the whole saturated snippet removes both at once
            // (the real bright flash survives via the unsaturated channels).
            // Default OFF -> reads no masks, bit-identical to every existing
            // config.  See pdhd/docs/run29107-evt1015-light-anomaly.md.
            bool m_veto_saturation{false};
            // KEEP saturated hits but mark them: emit a 10th ophit column =
            // 1 when the hit overlaps a "saturation" tick sub-range (same
            // overlap test as the veto), 0 otherwise.  OpFlashFinder turns
            // the column into a per-flash per-OpDet "flash_sat" tensor so
            // Q/L matching can mask railed channels per flash instead of
            // receiving them as zero (veto) or a clipped underestimate
            // (veto off).  Default OFF -> 9 columns, bit-identical.
            // See pdvd/docs/qlmatch/pdvd-saturation-recovery.md.
            bool m_flag_saturation{false};
            Configuration m_algo;        // SlidingWindow parameters

            int m_count{0};
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
