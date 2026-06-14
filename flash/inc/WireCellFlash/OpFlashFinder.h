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
            // Build flashes independently per side of the cathode (one
            // drift volume each) rather than across all OpDets.  The PDHD
            // cathode is opaque to scintillation light, so the two volumes
            // are optically independent; off by default (all-TPC, larana
            // behaviour), on for PDHD.  See stage2-reconstruction.md.
            bool m_group_by_side{false};

            // Optional post-construction "flash refinement": merge a later,
            // dim, few-PD flash into an earlier one whose lit OpDets are
            // physically adjacent -- the same flash over-split by the
            // per-channel OpHit splitter + 1us accumulators into a bright
            // primary plus small satellites riding its scintillation tail.
            // Cascading and per cathode side.  Off by default => bit-identical
            // to larana (and SBND); on for PDHD via flash.jsonnet.  See
            // pdhd-light-raw-data.md §4.
            bool   m_flash_refine{false};
            double m_refine_window_us{8.0};   // sliding merge window [us]
            double m_refine_pe_ratio{0.5};    // later_pe <= ratio * earlier_pe (run-27305 tuned)
            int    m_refine_max_fired{2};     // later flash fired-PD count cap
            double m_refine_fired_pe{0.5};    // pes[od] >= this counts as "fired" [PE]
            // Subset escape: bypass the few-PD cap when every lit OpDet of the
            // later flash is also lit in the earlier one (the tail fires only
            // PMTs the prompt already fired -> unambiguously the same flash).
            // Catches bright, spatially-extended parents whose tail spreads over
            // more than max_fired of the SAME PDs.  Strictly additive to the
            // existing merges.  Off by default; on for PDHD via flash.jsonnet.
            bool   m_refine_subset_merge{false};
            // Per-OpDet grid coordinate within its cathode side, built once in
            // configure(): row = y-rank 0..9, col = z-rank 0..7, side = 0 (x>=0)
            // / 1 (x<0).  Drives the 8-neighbour (Chebyshev<=1) adjacency test.
            std::vector<int> m_opdet_row, m_opdet_col, m_opdet_side;  // [nchan]

            // Per-event readout-vs-trigger offset (us): (tc - rd_timestamp)*16ns
            // from the ROOT trigoff tree (~250).  The reconstruction can't recover
            // this exactly from the self-triggered snippets (their origin sits ~tens
            // of us inside the readout window), so it is supplied by config and
            // stamped verbatim into the opflash tensor-set metadata as "offset_us"
            // for downstream charge clustering / Q-L matching.  0 => key written as
            // 0 (no offset applied downstream).
            double m_offset_us{0.0};

            std::vector<double> m_opdet_x, m_opdet_y, m_opdet_z;  // [nchan], mm

            int m_count{0};
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
