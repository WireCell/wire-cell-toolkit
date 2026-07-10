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

#include <map>
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
            // Optional OpChannel -> OpDet ganging map, JSON
            // {"channels": [{"opch": <id>, "opdet": <0..nchan-1>}, ...]}.
            // When set, incoming hit channel ids are remapped to OpDet
            // columns before flash building, so several readout channels
            // (e.g. the two DAPHNE channels of a PDVD X-ARAPUCA) sum into
            // one opflash PE column; hits on unmapped channels are ignored.
            // The output "ophits" tensor keeps the ORIGINAL channel ids.
            // Empty (default) = identity, bit-identical to before (hit
            // channel ids are already OpDet indices for PDHD/SBND).
            std::string m_channel_map_file{""};
            std::map<int, int> m_chmap;
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

            // Optional post-construction multiplicity / total-PE quality cut:
            // drop any flash with fewer than m_min_fired_pds lit OpDets
            // (pes[od] >= m_refine_fired_pe) or less than m_min_total_pe total
            // PE.  Removes the single-/few-PD low-PE flood (no nPD cut in the
            // larana port; cf. the prototype's pe>=6 && mult>=3).  Both 0 by
            // default => no cut => bit-identical to larana (and SBND); set for
            // PDHD via flash.jsonnet.  See run27980-processing-status.md §6.
            int    m_min_fired_pds{0};   // keep flash iff nPD >= this
            double m_min_total_pe{0.0};  // keep flash iff total_pe >= this
            // Per-PD PE threshold used ONLY for the m_min_fired_pds count (a PD
            // counts as "fired" iff pes[od] >= this).  Default -1 = use
            // m_refine_fired_pe (0.5), so the count is bit-identical to before
            // regardless of refine_fired_pe; the refinement-merge logic keeps
            // using m_refine_fired_pe unchanged.  0.5 PE is sub-single-
            // photoelectron; PDHD sets 1.0 (= one detected p.e.) via flash.jsonnet
            // so a PD must really fire to count.  See
            // run29107-evt1015-light-anomaly.md.
            double m_min_fired_pe{-1.0};

            // Per-event readout-vs-trigger offset (us): (tc - rd_timestamp)*16ns
            // from the ROOT trigoff tree (~250).  The reconstruction can't recover
            // this exactly from the self-triggered snippets (their origin sits ~tens
            // of us inside the readout window), so it is supplied by config and
            // stamped verbatim into the opflash tensor-set metadata as "offset_us"
            // for downstream charge clustering / Q-L matching.  0 => key written as
            // 0 (no offset applied downstream).
            double m_offset_us{0.0};

            // Extra key/values merged verbatim into the output tensor-set
            // metadata (e.g. PDVD's per-crate offsets offset_bot_us /
            // offset_top_us).  Null (default) => nothing stamped, output
            // byte-identical.
            Configuration m_metadata_extra;

            std::vector<double> m_opdet_x, m_opdet_y, m_opdet_z;  // [nchan], mm

            int m_count{0};
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
