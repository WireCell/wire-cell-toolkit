/** Optical flash finder for SBND reco1 PMT hits.
 *
 * SBND's reco1 OpHits are mostly single photons (tens of ns wide), so a flash is thousands
 * of narrow hits spread over the ~1.6 us argon slow light, unlike the few wide hits per
 * channel that OpFlashFinder (the larana OpFlashAlg port used for ProtoDUNE) expects.  This
 * component keeps OpFlashFinder's hit accumulation and claiming and replaces the rest with
 * steps that suit narrow hits, run in this order:
 *
 *  1. accumulate: two sets of bin_width (8 us) time bins, the second shifted by half a bin,
 *     anchored at the earliest hit; a bin with >= flash_threshold PE is a candidate;
 *  2. claim: candidates biggest first; each takes the hits no bigger candidate took and is
 *     kept if they hold >= flash_threshold PE (no width refinement);
 *  3. pulse split: cut a flash where its hits hold a second prompt burst of light (see the
 *     split_* keys), repeatedly on the later part;
 *  4. join (the tail pieces an 8 us bin edge cut off): a flash joins the nearest earlier flash
 *     when it starts <= join_max_gap_us after that flash's last hit, has >= join_lit_frac of its
 *     PE on OpDets the earlier flash lights (>= join_fired_pe; not all of it, since overlapping
 *     tails of other flashes add a few hits), the earlier flash's slow tail explains its PE
 *     (expected = PE_i x the share of an exponential with late_tau_us, starting at i's onset,
 *     that falls in the piece's time span, normalised to the span i covers; (PE - expected) /
 *     sqrt(expected) < late_nsigma.  Lenient on purpose: larana's step-5 formula underestimates
 *     SBND data tails ~2x) -- this glow test is applied only with join_glow "only"; the default
 *     "any" skips it (data: pieces touching a flash often carry extra light beyond the glow, and
 *     SBND's own 8 us windows count it in the flash; 91 % vs 67 % one-to-one with SBND on 63 data
 *     events), "never" joins only pieces with MORE light than the glow (glow pieces then go to
 *     step 5 and are deleted: loses light), the two are not parts of the same split, and the later
 *     flash has no sharp burst of its own (a split_burst_ns window with >= split_min_pe without its
 *     brightest OpDet, after a dip; fixed thresholds -- a real new flash has one, a tail piece does
 *     not).  No PE cap by default (join_pe_ratio <= 0): a flash starting shortly before a bin edge
 *     can have more light in its tail piece than in its early part;
 *  5. remove late light: a later flash whose PE the slow tail of an earlier one explains
 *     (larana: expected = PE_i * width_j / width_i * exp(-(t_j - t_i) / late_tau_us), deleted
 *     if (PE_j - expected) / sqrt(expected) < late_nsigma) is deleted; its hits go to no flash.
 *     Pieces touching their flash were already joined in step 4, so this only deletes the rest.
 *     Never applied between two parts of one split, nor to a flash at least as bright as the
 *     earlier one (the width ratio in the formula can otherwise "explain" a bright second flash
 *     right after a short split-off first part);
 *  6. quality cut: >= min_fired_pds OpDets with >= min_fired_pe PE, and >= min_total_pe;
 *  7. flash time = SBND's rule (SimpleFlashAlgo candidate bin + FlashT0SelectedChannels):
 *     brightest prompt_bin_ns bin, then the mean time of the largest hits within
 *     -prompt_pre_ns/+prompt_post_ns with PE > prompt_min_hit_pe, until they hold
 *     > prompt_pe_fraction of that PE.
 *
 * Input: ITensorSet with an "ophits" tensor (f8 [nhit, 9], OpHitFinder schema: channel,
 * time, width, area, amplitude, PE, start time, flash id, fast/total), e.g. from
 * SBNDReco1OpHitSource, one TPC's PMTs per set.  Output: the OpFlashFinder tensor set
 * ("opflash" [nflash, 1 + nchan] = time + PE per OpDet, "flash_summary", "ophits" with the
 * flash id filled in).  Times in WCT ns.  Defaults are the SBND settings.
 */
#ifndef WIRECELLFLASH_SBNDOPFLASHFINDER
#define WIRECELLFLASH_SBNDOPFLASHFINDER

#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/Logger.h"

#include <string>
#include <vector>

namespace WireCell {
    namespace Flash {
        class SBNDOpFlashFinder : public Aux::Logger, public ITensorSetFilter, public IConfigurable {
          public:
            SBNDOpFlashFinder();
            virtual ~SBNDOpFlashFinder();

            virtual bool operator()(const ITensorSet::pointer& in, ITensorSet::pointer& out);

            virtual WireCell::Configuration default_configuration() const;
            virtual void configure(const WireCell::Configuration& config);

            struct Params {
                double bin_width{8000.0};       // ns
                double flash_threshold{20.0};   // PE
                bool   pulse_split{true};
                double split_bin_ns{10.0};
                double split_burst_ns{100.0};
                double split_min_gap_ns{500.0};
                double split_spike_pe{30.0};
                double split_spike_frac{0.01};
                double split_min_ratio{0.05};
                double split_min_pe{100.0};
                double split_dip_ns{500.0};
                double split_dip_frac{0.5};
                bool   split_drop_brightest{true};
                bool   join{true};
                double join_max_gap_ns{1000.0};
                double join_pe_ratio{0.0};     // > 0: joined piece PE <= ratio * PE; <= 0: no cap
                double join_fired_pe{0.5};
                double join_lit_frac{0.9};
                std::string join_glow{"any"};  // any | only | never (see step 4)
                bool   remove_late_light{true};
                double late_tau_ns{1600.0};
                double late_nsigma{3.0};
                int    min_fired_pds{3};
                double min_fired_pe{1.0};
                double min_total_pe{0.0};
                double prompt_bin_ns{10.0};
                double prompt_pre_ns{20.0};
                double prompt_post_ns{10.0};
                double prompt_min_hit_pe{1.0};
                double prompt_pe_fraction{0.6};
            };

          private:
            int m_nchan{312};
            std::string m_geom_file{""};  // OpDet positions [mm]; empty => y/z centroids 0
            Params m_par;
            double m_offset_us{0.0};
            Configuration m_metadata_extra;
            std::vector<double> m_opdet_y, m_opdet_z;  // [nchan], mm
            int m_count{0};
        };
    }  // namespace Flash
}  // namespace WireCell

#endif
