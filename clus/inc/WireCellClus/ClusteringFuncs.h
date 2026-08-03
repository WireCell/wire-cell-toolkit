/**
   This header provides various free functions used in clustering.

   Some implementations may be found in clustering_*.cxx.  The rest are in
   ClusteringFuncs.cxx.

 */

#ifndef WIRECELLCLUS_CLUSTERINGFUNCS
#define WIRECELLCLUS_CLUSTERINGFUNCS

#include "WireCellClus/MultiAlgBlobClustering.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/ClusteringFuncsMixins.h"
#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/Graphs.h"

#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"



#include <string>
#include <fstream>
#include <set>

namespace WireCell::Clus::Facade {
    using namespace WireCell::PointCloud::Tree;

    /// Some clustering functions define and react to flags defined on cluster
    /// facades.  We name these flags with strings in this namespace to assure
    /// consistency, add documentation/comments and avoid typos.  Note,
    /// merge_clusters() will call fresh.from(other) to forward flags.
    namespace Flags {

        /// Indicates the cluster is a "live" cluster and connected to a "dead"
        /// cluster.
        inline const std::string live_dead = "live_dead";

        /// Indicates the cluster has a flash coincident with beam timing
        inline const std::string beam_flash = "beam_flash";
        
        /// Indicates the cluster is tagged as through-going muon (TGM)
        inline const std::string tgm = "tgm";
        
        /// Indicates the cluster is tagged as low energy
        inline const std::string low_energy = "low_energy";
        
        /// Indicates the cluster is tagged as light mismatch (LM)
        inline const std::string light_mismatch = "light_mismatch";
        
        /// Indicates the cluster is tagged as fully contained
        inline const std::string fully_contained = "fully_contained";
        
        /// Indicates the cluster is tagged as short track muon (STM)
        inline const std::string short_track_muon = "short_track_muon";
        
        /// Indicates the cluster has full detector dead region
        inline const std::string full_detector_dead = "full_detector_dead";

        // main cluster
        inline const std::string main_cluster = "main_cluster";

        // associated cluster
        inline const std::string associated_cluster = "associated_cluster"; 

        /// This flag is set by ClusteringTaggerCheckSTM algorithm when specific STM conditions are met
        inline const std::string STM = "STM";

        inline const std::string TGM = "TGM";

        /// Set by TaggerCheckFC when cluster_fc_check() finds the main cluster
        /// fully contained inside the fiducial volume.  Tagger-computed sibling
        /// of STM/TGM above -- NOT the same as the lowercase "fully_contained"
        /// flag, which is an imported uBooNE verdict copied out of a
        /// tagger_info point cloud by ClusteringTaggerFlagTransfer.
        inline const std::string FC = "FC";

    }

    struct ClusterLess {
        bool operator()(const Cluster* a, const Cluster* b) const {
            return cluster_less(a, b);
        }
    };

    using cluster_set_t = std::set<const Cluster*>;

    // Each vertex of a cluster connectivity graph represents the index of a
    // cluster in some collection.
    using cluster_connectivity_graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, int>;

    using cluster_vector_t = std::vector<Cluster*>;


    // This function will produce a new cluster in the grouping corresponding to
    // each connected component in the cluster_connectivity_graph_t that as two
    // or more clusters.  Any cluster in a single-cluster component is simply
    // left in place in the grouping.
    //
    // Each new cluster will be given the children (blob nodes) of the clusters
    // in the associated connected component.  These now empty clusters will be
    // removed from the grouping and discarded.
    //
    // If both aname and pcname are given then a representation of the previous
    // clustering of blob nodes will be stored in the new cluster.  This
    // connected component (cc) array is in child-node-order and its integer
    // value counts which original cluster donated the blob to the new cluster.
    //
    // The fresh.from(other) is called to transfer flags, scope and possibly
    // other bits of information.
    //
    // Pointers to the newly created cluster node facades are returned.  These
    // are loaned.  As usual, the cluster node owns the facade and these nodes
    // are in turn owned by the grouping node.
    // When orig_id_aname is non-empty, each merged cluster also gets a per-blob
    // int array (stored under orig_id_aname in PC "pcname") holding the original
    // ident() of the sub-cluster each blob came from, so a downstream consumer
    // (e.g. the Bee writer) can recover the pre-merge cluster identity of every
    // blob.  This is independent of the aname/parent_id ("perblob") array above.
    //
    // SCOPE INVARIANT -- the recorded ids are meaningful only WITHIN one
    // cluster.  Equal values in two DIFFERENT clusters name different objects
    // and must never be joined on.  MultiAlgBlobClustering re-runs
    // Grouping::enumerate_idents() after every visitor, so the ident() values
    // captured here belong to the numbering epoch in force at THIS call, while
    // the save-time fill-in for a cluster that was never merged writes the
    // ident of the LATEST epoch.  Both epochs are dense 1..N and nothing in the
    // array says which one a row belongs to, so the two overlap: measured 31%
    // of distinct "real_cluster_id" values name two clusters on the SBND d52ron
    // 30-event set (all of them recorded-vs-fill-in; recorded-vs-recorded
    // cannot collide, since idents are unique at the instant of this call).
    // See sbnd_xin/docs/53_unmerge-vs-cathode-crossers.md.  Set
    // MultiAlgBlobClustering's "real_cluster_id_global" to re-stamp the array
    // into one globally unique epoch at save time.
    // flags_from_longest (default false = historical behavior): Cluster::from()
    // is called once per member and flags_from() copies the donor's flag VALUES
    // (including zeros), so the merged cluster's flags are those of whichever
    // member happened to be visited last -- arbitrary.  A cluster that was the
    // matched main of its bundle can therefore lose flag_main_cluster to a tiny
    // co-merged fragment.  With this true, the flags are re-applied after the
    // merge from the same representative member that donates the flash (the
    // longest flash-bearing one; longest overall when no member has a flash),
    // mirroring the cluster_t0/flash/matched_flash_gid fixup below so that the
    // merged cluster's flags and flash describe one coherent donor.
    // orig_main_aname: when set (with orig_id_aname and pcname), additionally
    // save a per-blob marker (1/0) of the merged cluster's REPRESENTATIVE
    // member -- the flash/flags donor, i.e. the one main the merged cluster
    // reports -- so its points stay identifiable inside the merged cluster
    // (read back by TaggerCheckTGM main_component_mode="real").  The
    // "main_cluster" flag cannot serve: a flash group can merge several
    // bundles, each of whose mains carries it.
    // Independently of all of the above, a fixed registry of provenance PAIRS
    // (currently just "assoc_cluster_id"/"assoc_cluster_main") is CARRIED across
    // the merge whenever the members have them: the main array is concatenated
    // verbatim and the id array is rebased per member into a fresh dense range,
    // with a member lacking the arrays contributing one fresh id and main=1.
    // See the long comment at the top of merge_clusters() for why each of those
    // three choices is what it is; doc 52 4b has the full argument.  Absent
    // arrays => no-op, so this costs nothing when the writer knob is off.
    std::vector<Cluster*> merge_clusters(cluster_connectivity_graph_t& g, //
                                         Grouping& grouping,
                                         const std::string& aname="",
                                         const std::string& pcname="perblob",
                                         const std::string& orig_id_aname="",
                                         bool flags_from_longest=false,
                                         const std::string& orig_main_aname="");

    // Assign each cluster an integer "flash-time group" id.  Clusters whose
    // matched flash time (cluster_t0) differ by less than `window` share a group
    // id.  A cluster without a valid matched flash (scalar "flash" < 0) gets a
    // unique singleton id so it can never share a group with another cluster.
    // Every cluster in `clusters` is present as a key in the returned map.  This
    // is used by the combined-stage (post QL-matching) clustering functions to
    // restrict merging to clusters coincident in flash time.
    std::map<const Cluster*, int> assign_flash_t0_groups(
        const std::vector<Cluster*>& clusters, double window);

    // Diagnostic-only (env WCT_PROV_CHECK=1, else a fast no-op returning 0):
    // validate the per-blob provenance arrays of every cluster child of the
    // given grouping NODE -- per-key row counts vs blob count, and spatial
    // compactness of each assoc_cluster_main==0 group -- logging any violation
    // under logger "clus.provchk" tagged with `tag`.  Returns the number of
    // problems found.  Takes the raw node (not the facade) so it can run
    // where facades are not set up (PointTreeMerging, serialization
    // round-trips).  See clus/src/prov_check.cxx and doc 52 §11.
    size_t check_perblob_provenance(
        const WireCell::PointCloud::Tree::Points::node_t& groot,
        const std::string& tag);

    // Restore the "perblob" row-i == child-i invariant after a
    // separate-then-merge-back (or total take-back) round trip.  The round
    // trip permutes the blob children -- kept (cc<0) blobs stay in place,
    // split-off groups are re-appended at the END in ascending group-id order
    // (both FacadeParent::merge() and an ascending-key take_children loop do
    // this) -- while the cluster's node-local per-blob Dataset keeps its
    // original row order.  This applies the SAME permutation to the whole
    // Dataset (all arrays at once, via Dataset::subset), exactly as
    // QLMatching::recompose_cluster_groups() does (doc 52 §12.4).  `cc` is
    // the per-blob group array that was given to separate(), in PRE-split
    // children order.  Fail-open no-op (returns false) unless the Dataset
    // exists and both its major size and the cluster's current child count
    // equal cc.size() -- so at any pipeline stage where the array was never
    // written this is byte-identical to not calling it.
    bool realign_perblob_after_regroup(Cluster& cluster,
                                       const std::vector<int>& cc,
                                       const std::string& pcname = "perblob");



    /**
     * Extract geometry information from a grouping
     * @param grouping The input Grouping object
     * @param dv Detector geometry provider
     * @return Tuple of (drift_direction, angle_u, angle_v, angle_w)
     */
    std::tuple<geo_point_t, double, double, double> extract_geometry_params(
        const Grouping& grouping,
        IDetectorVolumes::pointer dv);

    /**
     * Validate that a set of wpids forms ONE drift volume: x-aligned APAs
     * (identical drift-x fiducial metadata) viewed through the SAME face.
     * allow_mixed_faces waives the same-face requirement (NOT the identical
     * FV_x one) for detectors where both faces of an anode share one drift
     * volume (PDVD: faces = y-halves of a CRP).  Raises ValueError naming
     * `who` on violation.
     */
    void validate_drift_group(
        const std::set<WirePlaneId>& wpids,
        IDetectorVolumes::pointer dv,
        bool allow_mixed_faces,
        const std::string& who);

    std::vector<std::pair<geo_point_t, const Blob*>> get_strategic_points(const Cluster& cluster);

    //helper function ..
    double Find_Closest_Points(const Cluster& cluster1,
			       const Cluster& cluster2,
			       double length_1,
			       double length_2,
			       double length_cut,
			       geo_point_t& p1_save, // output
			       geo_point_t& p2_save,  // output
                   bool flag_print = false
			       );
			       
    





    // Fiducial-volume bounds selected for the scope (set of drift volumes) that a
    // clustering pass operates on.  Populated by select_scope_fv().  All values are
    // in WireCell internal units.
    struct ScopeFV {
        double xmin{0}, xmax{0}, ymin{0}, ymax{0}, zmin{0}, zmax{0};
        double xmin_margin{0}, xmax_margin{0}, ymin_margin{0}, ymax_margin{0}, zmin_margin{0}, zmax_margin{0};
        geo_point_t vertical_dir{0, 1, 0}, beam_dir{0, 0, 1};
    };

    // Select the fiducial volume appropriate to the scope `dv` was configured for.
    // The scope is taken from the dv's configured drift volumes (dv->wpident_faces()),
    // NOT from which TPCs happen to have live activity in a given event:
    //   - dv spans >1 APA (or none)  -> the global "overall" (cryostat) FV.  This
    //     reproduces the legacy dv->metadata(WirePlaneId(0)) reads bit-for-bit, so
    //     all-APA stages are unchanged regardless of per-event activity.
    //   - dv spans a single APA      -> the union (outermost envelope) of that APA's
    //     configured per-(APA,face) FV blocks (full APA even if a face is quiet).
    // Any FV field missing from a per-face block falls back to the "overall" value.
    // vertical_dir / beam_dir are detector-global and are always read from "overall".
    //
    // common_face_x (default false, bit-identical): in the multi-APA branch, when
    // ALL configured faces carry identical FV_xmin/FV_xmax metadata (a drift-side
    // group: several APAs imaging one common drift side, e.g. PDHD group02 or a
    // PDVD CRP drift group), use that common x-range (and its margins) instead of
    // the cryostat overall x.  This makes the no-T0 "out-of-time" apparent-x test
    // reflect the group's drift volume rather than the union of both drift sides.
    ScopeFV select_scope_fv(IDetectorVolumes::pointer dv, bool common_face_x = false);

    // These Judge*() functions are used by multiple clustering methods.  They
    // are defined in clustering_separate.cxx.

    // time_slice_length is length span for a slice
    // guard_main_angle (deg): the long/thin/drift-aligned protection guard
    // (toolkit addition e28db401; not in the prototype) additionally requires
    // the cluster MAIN axis within this angle of the drift axis.  Default <0
    // keeps the guard unconditional (bit-identical).  The guard's angle1 tests
    // the 2nd axis against perp-to-drift, which wide isochronous complexes
    // satisfy trivially -- without the main-axis requirement it vetoes exactly
    // the multi-track over-clusters separation exists for.
    bool JudgeSeparateDec_1(const Cluster* cluster, const geo_point_t& drift_dir, const double length,
                            double guard_main_angle = -1);
    /// @attention contains hard-coded distance cuts
    /// @param boundary_points return the boundary points
    /// @param independent_points return the independent points
    /// @param fv scope-appropriate fiducial volume (see select_scope_fv)
    /// @param max_hull_points cap passed to Cluster::get_hull (<0 = use Constants::MaxHullPoints)
    /// @param far_point_x_cut drift-x deviation that promotes a boundary point to a
    ///        "far" point in the two-endpoint test.  Default 140 cm reproduces the
    ///        prototype expression `fabs(dir_3.x()/units::cm) > 14*units::cm` (the
    ///        cm-number compared against 14*cm internal = 140), which is effectively
    ///        dead; detectors may set the evidently intended 14 cm.
    /// @param far_point_mid_dis cap on the distance from the midpoint of the two
    ///        endpoints to the cluster, above which far points are discarded
    ///        (default 25 cm = prototype).  Two forking/diverging tracks can hold
    ///        their endpoint-midpoint in empty space between the prongs; detectors
    ///        may raise it to keep the far-point evidence for such topologies.
    bool JudgeSeparateDec_2(const Cluster* cluster, const geo_point_t& drift_dir,
                               std::vector<geo_point_t>& boundary_points, std::vector<geo_point_t>& independent_points,
                               const double cluster_length, const ScopeFV& fv, int max_hull_points = -1,
                               double far_point_x_cut = 140 * units::cm,
                               double far_point_mid_dis = 25 * units::cm);
    



    // this is used only by ClusteringSeparate but keep it public for symmetry with Separate_2.
    std::vector<Cluster *> Separate_1(const bool use_ctpc, Cluster *cluster,
                                      std::vector<geo_point_t> &boundary_points,
                                      std::vector<geo_point_t> &independent_points,
                                      double length, geo_point_t dir_cosmic, geo_point_t dir_beam, 
                                      IDetectorVolumes::pointer dv, 
                                      IPCTransformSet::pointer pcts,
                                      const Tree::Scope& scope);

    // This is used by multiple clustering methods.
    std::vector<int> Separate_2(Cluster *cluster, 
                                const Tree::Scope& scope,
                                const double dis_cut =  5*units::cm);

    // Function to compute wire plane parameters for clustering algorithms
    void compute_wireplane_params(
        const std::set<WirePlaneId>& wpids,
        const IDetectorVolumes::pointer dv,
        std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>>& wpid_params,
        std::map<WirePlaneId, std::pair<geo_point_t, double>>& wpid_U_dir,
        std::map<WirePlaneId, std::pair<geo_point_t, double>>& wpid_V_dir,
        std::map<WirePlaneId, std::pair<geo_point_t, double>>& wpid_W_dir,
        std::set<int>& apas);

    // Calculate PCA direction for a set of points around a center point
    geo_vector_t calc_pca_dir(const geo_point_t& center, const std::vector<geo_point_t>& points);

    /// Result of a cluster fully-contained (FC) boundary check.
    /// Shared between TaggerCheckNeutrino (which only needs is_fc) and
    /// TaggerCheckSTM (which also needs the exit endpoint data to drive STM analysis).
    struct FCCheckResult {
        /// True if every cluster endpoint lies inside the fiducial volume.
        bool is_fc{false};
        /// Candidate exit points (empty when is_fc == true).
        std::vector<geo_point_t> exit_wcps;
        /// Which steiner boundary endpoints (0=first, 1=second) are exit candidates.
        std::set<int> exit_boundary_set;
        /// The two steiner-graph boundary points (from round-1, flag_cosmic=true).
        /// These are the reference endpoints used by TaggerCheckSTM for path tracking.
        geo_point_t boundary_first{};
        geo_point_t boundary_second{};
    };

    /// Perform the two-round cluster boundary check to determine whether the
    /// cluster is fully contained inside the fiducial volume.
    ///
    /// The logic replicates the FC check originally embedded in
    /// TaggerCheckSTM::check_stm_conditions (round 1 with flag_cosmic=true,
    /// round 2 with flag_cosmic=false).  Returns a default FCCheckResult
    /// (is_fc=false) if the cluster has no steiner_pc or FiducialUtils.
    ///
    /// Used by:
    ///   - TaggerCheckNeutrino to fill tagger_info.match_isFC
    ///   - TaggerCheckSTM to drive STM / TGM classification
    ///   - TaggerCheckFC to set the "FC" cluster flag
    ///
    /// fiducial (optional, default nullptr): when null, the direct
    /// inside/outside tests use the grouping's FiducialUtils
    /// (DetectorVolumes = union of per-face sensitive volumes, no margin) --
    /// the historical behavior, unchanged for every caller that omits it.
    /// When given, those tests instead run against this IFiducial with
    /// fv_tolerance applied (negative = inset), matching how TaggerCheckTGM
    /// evaluates containment so the FC and TGM verdicts can be compared.
    /// Only the DIRECT fiducial tests are redirected: the dead-region and
    /// signal-processing checks keep using FiducialUtils' per-(apa,face)
    /// logic, exactly as TaggerCheckTGM does.
    FCCheckResult cluster_fc_check(Cluster& cluster, IDetectorVolumes::pointer dv,
                                   IFiducial::pointer fiducial = nullptr,
                                   const std::vector<double>& fv_tolerance = {});


}  // namespace WireCell::Clus::Facade

#endif

