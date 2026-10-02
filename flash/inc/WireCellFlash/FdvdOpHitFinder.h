/** OpHit finder for FD-VD (DUNE far detector vertical drift) SiPM waveforms.
 *
 * A port of the LArSoft `ophit10ppm` hit finder (larana OpHitAlg with
 * AlgoSiPM and PedAlgoEdges kHEAD, the FD-VD dunesw settings) by way of its
 * python port fdvd_sim/stageB/overlay_light.py `hitfind` (wcp
 * WireCell/wcp-porting-validation), which was gated hit for hit against
 * LArSoft (fdvd_sim doc 08a).  This component reproduces the python port bit
 * for bit (fdvd_sim doc 16 gate H).
 *
 * Input: ITensorSet holding the digitised waveform snippets as four tensors,
 * found by metadata "name":
 *   "ch"   int   [nsnip]     OpChannel of each snippet
 *   "t"    f8    [nsnip]     snippet start time [us] (the digitiser's TimeStamp)
 *   "off"  i8    [nsnip+1]   offset of each snippet in "adc" (last = total)
 *   "adc"  i2    [nsample]   samples, snippets concatenated
 * The snippet start times are kept as f8 microseconds (not frame tbins) so
 * the hit times are bit-identical to the python port.
 *
 * Output: ITensorSet with one "ophits" tensor (f8 [nhit, 9]) in the
 * OpHitFinder / OpFlashFinder column order (times in WCT ns):
 *   0 OpChannel  1 peak_time  2 width  3 area  4 amplitude  5 PE
 *   6 start_time  7 flash_id (-1)  8 fast_to_total (0)
 * which is what fdvd_sim/stageB/ophits_to_tensor.py writes from the python
 * hits.  Hits with |peak time| > max_abs_peak_us are dropped, as there.
 */
#ifndef WIRECELLFLASH_FDVDOPHITFINDER
#define WIRECELLFLASH_FDVDOPHITFINDER

#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <cstdint>
#include <vector>

namespace WireCell {
    namespace Flash {

        // One OpHit in python-port units: times [us], widths [us].
        struct FdvdOpHit {
            int channel;
            double peak_time_us, start_time_us, width_us, area, amplitude, pe;
        };

        // The ophit10ppm settings (overlay_light.py line 33).
        struct FdvdOpHitParams {
            int nped{3};                 // PedAlgoEdges kHEAD samples
            double threshold{15.0};      // AlgoSiPM ADCThreshold: peak must reach it
            double threshold2{1.0};      // AlgoSiPM 2nd threshold: pulse = samples >= it
            int min_width{60};           // AlgoSiPM MinWidth [ticks], end - start
            double hit_threshold{0.2};   // OpHitAlg HitThreshold on the peak
            double spe_area{130.0};      // SPEArea
            double spe_shift{0.43};      // SPEShift
            double tick_us{0.016};       // 62.5 MHz
        };

        // The python hitfind on concatenated snippets; hits in sample order.
        std::vector<FdvdOpHit> fdvd_find_ophits(const std::vector<int>& ch, const std::vector<double>& t,
                                                const std::vector<int64_t>& off, const std::vector<double>& adc,
                                                const FdvdOpHitParams& par);

        // numpy's float64 add.reduce on a contiguous array (pairwise, 8
        // accumulators, block 128), so sums match numpy bit for bit.
        double fdvd_numpy_sum(const double* a, size_t n);

        class FdvdOpHitFinder : public Aux::Logger, public ITensorSetFilter, public IConfigurable {
          public:
            FdvdOpHitFinder();
            virtual ~FdvdOpHitFinder();

            virtual bool operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            FdvdOpHitParams m_par;
            double m_max_abs_peak_us{1e6};   // ophits_to_tensor.py drop rule
            size_t m_count{0};
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
