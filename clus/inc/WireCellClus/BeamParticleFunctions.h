/** Free functions behind the CheckBeamParticle visitor (doc pdvd/120).

    The three selection rules of the beam-particle PR stage -- which flash
    bundle is the beam bundle, which cluster of that bundle carries the beam
    particle, and which end of that cluster is its entry -- operate on plain
    numbers so they can be unit-tested (clus/test/doctest_beam_particle_
    functions.cxx) without a Grouping, a fitter or a detector.  The visitor
    (clus/src/CheckBeamParticle.cxx) fills the rows from the facade and reads
    the verdict back.

    Iteration is deterministic by construction: bundles are keyed by their
    int gid, candidates are visited in the caller's (ident-sorted) order and
    every tie breaks on a value, never on an address.
 */
#ifndef WIRECELL_CLUS_BEAMPARTICLEFUNCTIONS
#define WIRECELL_CLUS_BEAMPARTICLEFUNCTIONS

#include "WireCellUtil/Point.h"

#include <string>
#include <vector>

namespace WireCell::Clus::PR {

    /// One row per Flags::main_cluster cluster of the grouping.
    struct BeamBundleRow {
        int cluster_id{-1};
        int gid{-1};          ///< scalar "matched_flash_gid"; -1 = unmatched
        double t0{0};         ///< Cluster::get_cluster_t0(), the RAW flash time (internal units)
        double length{0};     ///< Cluster::get_length() (internal units)
    };

    /// The grouping's "opflash" row for one gid (Grouping::flash_by_gid).
    struct BeamFlashInfo {
        int gid{-1};
        double pe{-1};        ///< Flash::value(); -1 when the row is invalid
        bool valid{false};
    };

    struct BeamBundlePick {
        int gid{-1};          ///< the beam bundle; -1 = none
        int n_in_window{0};   ///< mains with lo <= t0 < hi
        int n_gids{0};        ///< distinct gids among them
        std::string why;      ///< "brightest" | "longest" | "none" | "window-off"
    };

    /// The beam bundle: among the mains with lo <= t0 < hi (half-open, the
    /// TaggerCheckNeutrino beam_window convention) grouped by gid, the gid
    /// whose flash is the brightest (the rule the Bee op_beam label uses,
    /// MultiAlgBlobClustering::fill_bee_flashes).  When no in-window gid has
    /// a valid flash row, or the brightest two tie, the gid holding the
    /// longest main wins; a final tie goes to the smallest gid.  gid -1 when
    /// nothing is in the window or lo >= hi (the gate is off => this stage
    /// never runs on "every bundle").
    BeamBundlePick beam_particle_pick_bundle(const std::vector<BeamBundleRow>& mains,
                                             const std::vector<BeamFlashInfo>& flashes,
                                             double lo, double hi);

    /// One bundle member as a main-cluster candidate.
    struct BeamMainCand {
        int cluster_id{-1};
        double dist{1e9};     ///< |closest point - nominal entry| (internal units)
        double length{0};
    };

    /// Index into `cands` of the beam particle's cluster: the smallest dist
    /// among candidates with length >= min_length; ties -> the longer, then
    /// the smaller cluster_id.  -1 when empty or every candidate is filtered.
    int beam_particle_pick_main(const std::vector<BeamMainCand>& cands, double min_length);

    struct BeamEntryChoice {
        WireCell::Point entry;
        WireCell::Point exit;
        double dist{1e9};     ///< |entry - nominal|
        double cos_beam{0};   ///< unit(exit - entry) . unit(beam_dir); 0 when degenerate
        bool tie_by_dir{false};
    };

    /// Which of the two axis-extreme points a and b is the entry: the one
    /// nearer `nominal`; when the two distances differ by less than tie_tol
    /// the end whose (other - this) direction has the larger dot with
    /// beam_dir (the particle travels INTO the detector along beam_dir).
    BeamEntryChoice beam_particle_choose_entry(const WireCell::Point& a, const WireCell::Point& b,
                                               const WireCell::Point& nominal,
                                               const WireCell::Vector& beam_dir, double tie_tol);

}  // namespace WireCell::Clus::PR

#endif
