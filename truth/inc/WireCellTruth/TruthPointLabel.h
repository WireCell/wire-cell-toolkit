/** TruthPointLabel: the per-point truth labels of a live point tree.

    The algorithm of TruthPointLabeler, separated from the pgraph node
    so that it can be unit tested on a hand-built tree and reused by a
    clus visitor (header + one .cxx in WireCellTruth).

    Given the live point tree of PointTreeBuilding (cluster nodes with
    blob nodes carrying the "3d" and "scalar" local point clouds), the
    ICluster it was sampled from, the drifted depos, the association
    result of DepoBlobAssoc (run with keep_depos when refinement is
    wanted), the truth table and the PFMerger result, add to EVERY blob
    node a local point cloud named Schema::truth_pc_name with the keys
    of Schema::truth_pc_keys(), same major size and order as "3d".

    The blob nodes are visited in the order PointTreeBuilding::sample_live
    created them (blob_clusters() over the cluster graph, b-vertices in
    set order) and each node is cross-checked against its blob by the
    "scalar" keys wpid, slice_index_min and the per-plane wire index
    bounds; a blob without a node (sample_live skipped it) is passed
    over and a node without a blob throws.

    Per point (doc 15 sec 2.1):
      tid            track key of the nearest contributing depo when refined,
                     else the blob's largest charge fraction; tid_shared when the
                     blob's largest fraction is below 0.5 (not refined), tid_ghost
                     when the blob has no true charge
      orig_tid       the pre-merge G4 id of that key (sed_ids), else tid
      pid_frac       the blob charge fraction of tid
      true_q_frac    blob true charge / blob solved charge (0 when unsolved)
      instance_valid 0 for ghost, shared and LE points
      pdg, cls, iid, pf_node   from the table and PFMerger (iid = interaction row)
      d_traj, dir_x/y/z, s_t, end, dedx_true   from the DepoPolyline of tid
                     (d_traj capped at d_cap, s_t = arc fraction, end = 1 near
                     the start, 2 near the end, within end_window; dedx in
                     energy per cm)
      v_heat_vtx/end/kink/conv  Gaussian heat (sigma heat_sigma) around the
                     starts of surviving primary nodes, this track's end, this
                     track's kinks, and the starts of "conv" electrons

    wcp-porting-img pr_ml/docs/19 sec 3.6 (T6).
 */

#ifndef WIRECELLTRUTH_TRUTHPOINTLABEL
#define WIRECELLTRUTH_TRUTHPOINTLABEL

#include "WireCellTruth/DepoBlobAssoc.h"
#include "WireCellTruth/DepoPolyline.h"
#include "WireCellTruth/PFMerger.h"
#include "WireCellTruth/TruthTable.h"
#include "WireCellUtil/PointTree.h"
#include "WireCellUtil/Units.h"

#include <functional>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace WireCell::Truth {

    struct PointLabelParams {
        /// Per-point tid from the nearest contributing depo (needs keep_depos).
        bool refine{true};
        double d_cap{25 * units::mm};         // tau_off, ~5 pitches of 4.79 mm
        double end_window{10 * units::mm};    // tau_on, ~2 pitches
        double kink_window{3 * units::cm};
        double kink_angle{15 * units::deg};
        double heat_sigma{1 * units::cm};
        /// Blob-level fraction below which the blob is "shared".
        double shared_frac{0.5};
        /// Do the shared test on the blob's true charge summed per PF
        /// node instead of per G4 track (pr_ml doc 25 sec 2.8: inside a
        /// shower the charge is split among many tracks of one node).
        /// Off = per G4 track, as before.  Recorded in meta() only when on.
        bool shared_by_pf_node{false};
        /// The end flag from the point's own PF node (pr_ml doc 25 sec 3.2):
        /// 1 / 2 when the point lies within end_radius (3-D) of the node's
        /// start / end point (track, delta and Michel nodes: both; other
        /// classes: start only).  Off = the legacy flag (arc position on the
        /// nearest polyline of the key, no distance requirement).  Recorded
        /// in meta() only when on.
        bool end_by_pf_node{false};
        double end_radius{10 * units::mm};
        /// Also store, per point, the 3-D distance to its own PF node's
        /// start / end (truth keys d_node_start, d_node_end; -1 when the
        /// node has no such end in this face or the point has no node), so
        /// a training can choose its own radius (pr_ml doc 25 sec 3.2).
        /// Off = the key set of truth_pc_keys().  Recorded in meta() only
        /// when on.  Independent of end_by_pf_node (both use the node ends).
        bool end_distance{false};
        /// The heat targets from the PF nodes (pr_ml doc 25 sec 4.2):
        /// v_heat_end from the point's own node end; v_heat_kink from the
        /// kinks of its own track-like node's root particle; v_heat_vtx
        /// from true interaction vertices only (neutrino: the primaries'
        /// start; beam / gun: the primary's start inside the face's
        /// sensitive box, else its end; cosmic and radiological: none);
        /// v_heat_conv from the starts of photon-rooted shower nodes.
        /// Adds the truth keys d_vertex / d_conversion (distance to the
        /// nearest vertex / conversion point, -1 when none).  Off = the
        /// legacy key-polyline heat.  Recorded in meta() only when on.
        bool heat_by_pf_node{false};
        /// With heat_by_pf_node: kinks only from track nodes whose root is
        /// a muon or a charged hadron (not delta / Michel electrons, which
        /// scatter, nor neutrons).  Off = every track-like node (track,
        /// delta, Michel).  Recorded in meta() only when on.
        bool kink_tracks_only{false};
        /// Order a particle's depos of equal Geant4 time along the local
        /// direction (order_time_ties) in every polyline and when picking a
        /// PF node's start / end (pr_ml doc 25 sec 4.2: on PDHD / PDVD / SBND
        /// a particle's depos tie in time, and the time-ordered polylines
        /// zigzag 2-4x their length).  Off = time order, ties in file order.
        /// Recorded in meta() only when on.
        bool time_tie_order{false};
        /// With refine: a point may take only a track key whose PF node
        /// (a key without a node counts as its own node) holds at least
        /// this fraction of the blob's true charge; the nearest such depo
        /// gives the point its tid, and when no key qualifies the point
        /// keeps the unrefined label (the blob's top key, or tid_shared).
        /// pr_ml doc 25 sec 7.4: the nearest depo of a keV Compton electron
        /// holding 0.1 % of a blob handed its points that electron's node
        /// and class (2-6 % of ProtoDUNE points, 18 % of an FD-HD nue
        /// event's).  0 = every key (legacy).  Recorded in meta() only when
        /// on.
        double refine_min_frac{0.0};
        /// The tick used by the sampler for slice_index_min (cross-check).
        double tick{0.5 * units::us};

        Configuration meta() const
        {
            Configuration m;
            m["refine_by_nearest_depo"] = refine;
            m["d_cap_mm"] = d_cap / units::mm;
            m["end_window_mm"] = end_window / units::mm;
            m["kink_window_cm"] = kink_window / units::cm;
            m["kink_angle_deg"] = kink_angle / units::deg;
            m["heat_sigma_cm"] = heat_sigma / units::cm;
            m["shared_frac"] = shared_frac;
            if (shared_by_pf_node) m["shared_by_pf_node"] = true;
            if (end_by_pf_node) {
                m["end_by_pf_node"] = true;
                m["end_radius_mm"] = end_radius / units::mm;
            }
            if (end_distance) m["end_distance"] = true;
            if (heat_by_pf_node) m["heat_by_pf_node"] = true;
            if (kink_tracks_only) m["kink_tracks_only"] = true;
            if (time_tie_order) m["time_tie_order"] = true;
            if (refine_min_frac > 0) m["refine_min_frac"] = refine_min_frac;
            return m;
        }
    };

    struct PointLabelSummary {
        size_t nclusters{0}, nnodes{0}, npoints{0};
        size_t nblobs_skipped{0};   // blobs without a node
        size_t nghost{0}, nshared{0}, nrefined_differ{0};
        size_t npolylines{0};
        std::set<int> void_roots;
    };

    /// A PF node's start / end points in the reconstructed frame of one
    /// face (end_by_pf_node); has_* false when the point is not in the
    /// face's drift volume or the node has no such end.
    struct NodeEnds {
        bool has_start{false}, has_end{false};
        Point start, end;
    };
    using node_ends_t = std::map<int, NodeEnds>;   // PF node index -> ends

    /// heat_by_pf_node: whether a PF node's root particle gives kink keypoints.
    inline bool node_has_kinks(const PFNode& n, bool tracks_only)
    {
        if (!tracks_only) return n.cls == Cls::track || n.cls == Cls::delta || n.cls == Cls::michel;
        if (n.cls != Cls::track) return false;
        switch (std::abs(n.pdg)) {
        case 13: case 211: case 321: case 2212: case 3222: case 3112: return true;   // muon, charged hadrons
        default: return false;
        }
    }

    /// refine_min_frac: the track keys of a blob (by_track: key -> (true
    /// charge, energy)) that nearest-depo refinement may give a point --
    /// those whose PF node holds at least min_frac of the blob's true
    /// charge q_total; node_of(key) is the key's PF node, < 0 for none (the
    /// key is then its own node, as in the shared_by_pf_node test).
    /// min_frac <= 0: every key.
    template <class NodeOf>
    std::set<int> refine_eligible_tids(const std::map<int, std::pair<double, double>>& by_track, double q_total,
                                       NodeOf node_of, double min_frac)
    {
        std::set<int> out;
        auto key_of = [&](int tid) {
            const int pfn = node_of(tid);
            return pfn >= 0 ? std::make_pair(0, pfn) : std::make_pair(1, tid);
        };
        std::map<std::pair<int, int>, double> bynode;
        for (const auto& [tid, qe] : by_track) bynode[key_of(tid)] += qe.first;
        for (const auto& [tid, qe] : by_track) {
            if (min_frac <= 0 || (q_total > 0 && bynode[key_of(tid)] >= min_frac * q_total)) out.insert(tid);
        }
        return out;
    }

    /// heat_by_pf_node keypoints in the reconstructed frame of one face.
    struct HeatKeys {
        std::vector<Point> vertices, conversions;
        std::map<int, std::vector<Point>> kinks;   // PF node index -> kinks of its root particle
    };

    /// Reconstructed-frame position of a drifted depo (x from its
    /// arrival time, as PointTreeBuilding places blobs).
    using depo_position_f = std::function<Point(const IDepo&)>;

    /// The polylines of every track key, built from the depos (one
    /// vertex per depo: the reco-frame position, ordered by the prior's
    /// time, energy of the key depo).  Keys follow the association's
    /// id_from_prior / id_map.
    std::map<int, std::shared_ptr<DepoPolyline>>
    build_polylines(const IDepo::vector& depos, const DepoBlobAssocParams& assoc,
                    const depo_position_f& position, bool tie_by_direction = false);

    /// Per-particle polylines of every track key (pr_ml doc 23 sec 7.4):
    /// the depos of a key grouped by their Geant4 particle -- orig_ids[depo
    /// id] when given (Stage A index mode: the sed_ids orig_trackid
    /// column), else the depo id before id_map and folding -- one polyline
    /// per particle, in particle-id order.  A key folded from several
    /// particles is then a set of polylines instead of one time-ordered
    /// zigzag across them.
    using particle_polylines_t = std::map<int, std::vector<std::shared_ptr<DepoPolyline>>>;
    particle_polylines_t
    build_particle_polylines(const IDepo::vector& depos, const DepoBlobAssocParams& assoc,
                             const depo_position_f& position,
                             const std::vector<int>* orig_ids = nullptr,
                             bool tie_by_direction = false);

    /// Label the tree in place.  Throws ValueError when a node has no
    /// matching blob.  With particles, the per-point trajectory fields
    /// (d_traj, dir, s_t, dedx, end) come from the nearest particle
    /// polyline of the key and the end / kink heat from the keypoints of
    /// all of them; without, from the key's single polyline (the legacy
    /// path).  The vertex and conversion keypoints use polylines either way.
    PointLabelSummary
    label_live_tree(PointCloud::Tree::Points::node_t& root,
                    const ICluster::pointer& cluster,
                    const IDepo::vector& depos,
                    const DepoBlobAssocResult& res,
                    const TruthTable& table,
                    const PFResult& pf,
                    const std::map<int, std::shared_ptr<DepoPolyline>>& polylines,
                    const PointLabelParams& params,
                    const depo_position_f& position,
                    Log::logptr_t log = nullptr,
                    const particle_polylines_t* particles = nullptr,
                    const node_ends_t* node_ends = nullptr,
                    const HeatKeys* heat_keys = nullptr);
}

#endif
