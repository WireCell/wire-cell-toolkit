/** Optical flash assembly for PDHD light reconstruction.
 *
 * A dependency-free port of the larana `OpFlashAlg` flash finder
 * (RunFlashFinder: double-offset accumulators, hit claiming, width
 * refinement, late-light removal) in the `protodune_opflash`
 * configuration, with flash y/z from the OpDet geometry.
 *
 * Input: ITensorSet with an "ophits" tensor (f8 [nhit, 9], see
 * OpHitFinder).  Output: the opflash tensor set of
 * flash/docs/design.md §3.4 ("opflash" matrix + "flash_summary" +
 * "ophits" with flash_id assigned).
 */
#ifndef WIRECELLFLASH_OPFLASHFINDER
#define WIRECELLFLASH_OPFLASHFINDER

#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <vector>

namespace WireCell {
    namespace Flash {
        class OpFlashFinder : public Aux::Logger, public ITensorSetFilter, public IConfigurable {
          public:
            OpFlashFinder();
            virtual ~OpFlashFinder();

            virtual bool operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

          private:
            int m_nchan{160};
            std::string m_geom_file{"pgrapher/experiment/pdhd/pdhd-opdet-geom.json"};
            // standard_opflash / protodune_opflash values (times in WCT ns).
            double m_bin_width{1000.0};      // BinWidth, 1 us
            double m_flash_threshold{3.5};   // FlashThreshold [PE] (dunefd/protodune)
            double m_width_tolerance{0.5};   // WidthTolerance
            bool m_remove_late_light{true};

            std::vector<double> m_opdet_y, m_opdet_z;  // [nchan], mm

            int m_count{0};
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
