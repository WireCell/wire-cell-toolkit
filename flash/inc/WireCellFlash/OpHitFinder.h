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

          private:
            std::string m_intag{"decon"};
            // dune_ophit_finder_deco values.
            double m_scale{100.0};       // ScalingFactor applied before short cast
            double m_spe_area{100.0};    // SPEArea: PE = area / spe_area
            double m_hit_threshold{3.0}; // HitThreshold on pulse peak (scaled units)
            int m_ped_nsamples{3};       // PedAlgoEdges NumSampleFront (head method)
            Configuration m_algo;        // SlidingWindow parameters

            int m_count{0};
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
