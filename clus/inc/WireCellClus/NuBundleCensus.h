#ifndef WIRECELL_CLUS_NUBUNDLECENSUS_H
#define WIRECELL_CLUS_NUBUNDLECENSUS_H

/// What TaggerCheckNeutrino's candidate selection saw and decided, per event,
/// per in-beam-window flash bundle and per optical flash
/// (sbnd_xin/docs/109, groups 1 and 2 of doc 108).
///
/// Filled only when TaggerCheckNeutrino's `nu_provenance` knob is on, and
/// carried to the ROOT writers on the Facade::Grouping bus
/// (Grouping::set_nu_census), the same way the TrackFitting is.  It is NOT a
/// point-cloud array, so the saved pctree is unchanged.
///
/// The census is published BEFORE TaggerCheckNeutrino's "no candidate" return,
/// so an event with no T_tagger row still says why: every in-window bundle
/// has a row with a reason code, and in-window clusters with no matched flash
/// (silently dropped by the selection) are counted.
///
/// Units: times in microseconds, lengths in cm, charge in PE.

#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace WireCell::Clus::PR {

    struct NuBundleCensus {
        // The selection's own configuration, so a writer never needs a second
        // copy of it (TaggerCheckNeutrino beam_window_low/high, flash_pair_dt_us).
        double beam_window_low_us{0};
        double beam_window_high_us{0};
        // The same window in internal WCT units: the exact numbers the
        // selection compares cluster_t0 against, so a writer applying the same
        // test cannot flip at a window edge through a unit round trip.
        double beam_window_low{0};
        double beam_window_high{0};
        double flash_pair_dt_us{0};
        int beam_gate{0};     // 1 = beam window configured (low < high)
        int per_bundle{0};    // 1 = nu_per_bundle selection; 0 = legacy single winner (bundles[] then empty)

        // ---- event counters ------------------------------------------- //
        int n_main{0};                 // clusters flagged main_cluster
        int n_in_window_main{0};       // ... with cluster_t0 in the beam window
        int n_in_window_demoted{0};    // demoted mains with cluster_t0 in the beam window
        int n_in_window_nogid{0};      // in-window mains + demoted mains with matched_flash_gid < 0 (dropped by the selection)
        int n_bundles{0};              // in-window bundles (distinct gids) examined
        int n_candidates{0};           // bundles that yielded a neutrino candidate (= T_tagger rows)

        /// Why a bundle did or did not yield a candidate.
        enum Reason {
            kSelectedMain = 0,      // a main activity was selected
            kSelectedDemoted = 1,   // the demoted-main fallback was selected
            kAllCosmic = 2,         // every examined activity was TGM/STM/lm_flag>0 vetoed
            kLengthFloor = 3,       // at least one activity failed only the nu_per_bundle_min_length floor
            kStmOnly = 4,           // rejected by nu_per_bundle_stm_only (and nothing failed the floor)
            kNoEligible = 5,        // no activity examined (e.g. only demoted mains and the fallback is off)
            kDedupFlashGroup = 6,   // sbnd_xin/docs/109 rev 3: a candidate was selected, then dropped
                                    // because another bundle of the SAME flash_group (one physical
                                    // flash seen by both TPCs) kept a longer one (nu_dedup_flash_group)
        };

        struct Bundle {
            int gid{-1};
            int flash_tpc{-1};          // physical drift side of the flash (opflash "apa"); -1 = flash not found
            float flash_time_us{0};
            float flash_pe{0};
            int flash_group{-1};        // shared by different-TPC flashes within flash_pair_dt_us (min gid of the group)
            int n_main{0};              // in-window mains in the bundle
            int n_demoted{0};           // in-window demoted mains in the bundle
            int n_companion{0};         // associated clusters kept as companions (selected bundles only)
            int n_companion_dropped{0}; // companions dropped as cosmic (skip_cosmic_companions)
            int n_rej_cosmic{0};        // activities pick() vetoed as cosmic
            int n_rej_stm_only{0};      // activities pick() rejected under nu_per_bundle_stm_only
            int n_rej_floor{0};         // activities pick() rejected by the length floor
            int reason{kNoEligible};
            int nu_index{-1};           // T_tagger/T_kine row; -1 = no candidate
            int sel_cluster_id{-1};     // activity the selection chose
            int final_cluster_id{-1};   // main cluster the row was written with (after any vertex-driven move)
            float sel_length_cm{0};
        };
        std::vector<Bundle> bundles;    // ascending gid

        struct Flash {
            int gid{-1};
            int tpc{-1};
            float time_us{0};
            float pe{0};
            int in_window{0};           // time in [beam_window_low_us, beam_window_high_us)
            int flash_group{-1};
            int n_matched_clusters{0};  // clusters (any role) whose matched_flash_gid == gid
            int n_matched_main{0};      // ... flagged main_cluster
            int nu_index{-1};           // candidate row using this flash; -1 = none
        };
        std::vector<Flash> flashes;     // ascending gid

        /// Index into flashes/bundles by gid, or -1.
        int flash_index(int gid) const
        {
            for (std::size_t i = 0; i < flashes.size(); ++i)
                if (flashes[i].gid == gid) return static_cast<int>(i);
            return -1;
        }
        int bundle_index(int gid) const
        {
            for (std::size_t i = 0; i < bundles.size(); ++i)
                if (bundles[i].gid == gid) return static_cast<int>(i);
            return -1;
        }
    };

    /// Group flashes: two flashes on DIFFERENT TPCs whose times differ by less
    /// than dt_us are one physical flash seen by both TPCs' light detectors
    /// (doc 108 sec 4.1: the one-flash peak of that |dt| ends near 0.05 us).
    /// Transitive (union-find).  Returns one group id per input flash = the
    /// smallest gid in its group.  Inputs are parallel arrays; dt_us <= 0
    /// gives every flash its own group.  Order-independent.
    std::vector<int> group_flashes(const std::vector<int>& gid, const std::vector<int>& tpc,
                                   const std::vector<double>& time_us, double dt_us);

    /// sbnd_xin/docs/109 rev 3.  Collapse neutrino candidates that come from
    /// one physical flash.  `gid` lists the candidates' bundle gids IN THE
    /// ORDER THE SELECTION RANKED THEM (longest selected activity first), and
    /// `gid_group` maps a gid to its flash group (group_flashes above).
    /// Returns the indices to KEEP: the first candidate of each group, i.e. the
    /// longest, which is also the one the row ordering would have put in slot 0.
    /// A gid absent from the map keeps its own gid as its group, so a flash that
    /// could not be grouped can never be mistaken for a duplicate.
    std::vector<std::size_t> dedup_flash_groups(const std::vector<int>& gid,
                                                const std::map<int, int>& gid_group);

    /// sbnd_xin/docs/109 rev 4 (nu_bundle_flash_group).  The pure parts of
    /// "one candidate per physical flash, but only when the two halves touch".
    /// A flash group (group_flashes) says two bundles saw the SAME light; it
    /// cannot say whether that light came from ONE interaction whose particle
    /// crossed the cathode (one candidate, double counted today) or from TWO
    /// interactions, one per drift volume (two candidates, both real).  The
    /// physical contact of their charge is what separates the two.  Measured on
    /// the colleague's r472 s36 e40 (doc 109 sec 9): the two bundles touch at
    /// 0.4 cm -- at the VERTEX, not at the cathode, because the Q/L matching
    /// had already put the muon's near-side segment into the far flash's
    /// bundle as an associated cluster.  So the contact test is a distance,
    /// and the cathode window is an optional tightening (xcut > 0), off by
    /// default.

    /// The bundle pairs worth testing: same flash group, DIFFERENT TPC, both
    /// gids in `gids`.  Ascending (min gid, max gid) order, so the caller's
    /// union-find is deterministic.  A gid absent from gid_group is its own
    /// group and is never paired; a gid absent from gid_tpc has tpc -1 and is
    /// never paired either (an unknown side cannot be "the other side").
    std::vector<std::pair<int, int>> cathode_pair_candidates(const std::vector<int>& gids,
                                                             const std::map<int, int>& gid_group,
                                                             const std::map<int, int>& gid_tpc);

    /// The contact predicate on one closest-approach pair: distance `d`
    /// between the two clusters' closest points, whose drift coordinates are
    /// `xa` and `xb`.  Contact iff d < gap AND, only when xcut > 0, both
    /// points lie within xcut of cathode_x.  All in the caller's units.  NaN
    /// never contacts.
    bool bundle_contact(double d, double xa, double xb, double cathode_x, double xcut, double gap);

    /// Union-find over `gids` with the pairs the caller found in contact.
    /// Returns gid -> root, the SMALLEST gid of its merged bundle; a gid in no
    /// contacting pair maps to itself (the legacy bundle, byte-for-byte).
    std::map<int, int> merge_bundles(const std::vector<int>& gids,
                                     const std::vector<std::pair<int, int>>& contacts);

}  // namespace WireCell::Clus::PR

#endif  // WIRECELL_CLUS_NUBUNDLECENSUS_H
