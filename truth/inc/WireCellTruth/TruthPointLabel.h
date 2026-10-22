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

    /// Reconstructed-frame position of a drifted depo (x from its
    /// arrival time, as PointTreeBuilding places blobs).
    using depo_position_f = std::function<Point(const IDepo&)>;

    /// The polylines of every track key, built from the depos (one
    /// vertex per depo: the reco-frame position, ordered by the prior's
    /// time, energy of the key depo).  Keys follow the association's
    /// id_from_prior / id_map.
    std::map<int, std::shared_ptr<DepoPolyline>>
    build_polylines(const IDepo::vector& depos, const DepoBlobAssocParams& assoc,
                    const depo_position_f& position);

    /// Label the tree in place.  Throws ValueError when a node has no
    /// matching blob.
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
                    Log::logptr_t log = nullptr);
}

#endif
