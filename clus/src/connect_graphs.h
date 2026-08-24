#ifndef WIRECELLCLUS_PRIVATE_CONNECT_GRAPHS
#define WIRECELLCLUS_PRIVATE_CONNECT_GRAPHS

#include "WireCellClus/Graphs.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/Facade_Cluster.h"

namespace WireCell::Clus::Graphs {


    // See make_graphs.h for bundling of construction and connecting.

    void connect_graph(const Facade::Cluster& cluster, 
                       Weighted::Graph& graph);

    void connect_graph_with_reference(
        const Facade::Cluster& cluster,
        const Facade::Cluster& ref_cluster,
        Weighted::Graph& graph);

    // doc 79 round 2: parameters of the busy-cluster lazy-walk mode of
    // connect_graph_ctpc, selected only by the "ctpc_fast" flavor
    // (ClusteringDeghost graph_name knob).  A cluster whose
    // connected-component count `num` is at or below busy_num_threshold
    // takes the legacy path bit-for-bit; above it, the per-pair 1 cm
    // is_good_point walk (the O(num^2 x path-cm) cost -- 49.5M steps on
    // the mcp1k 280884 monster, doc 79 round 0) is evaluated lazily via
    // a single-pass Kruskal: only pairs that would bridge two union-find
    // components (plus the < 3 cm direct-edge pairs and pairs carrying a
    // directional candidate) are walked.  Same scheme as RelaxedFastCfg
    // above (doc 78); NOT proven bit-identical under exact distance ties,
    // which is why it sits behind the busy gate and A/B validation.
    struct CtpcFastCfg {
        size_t busy_num_threshold{200};
    };

    void connect_graph_ctpc(const Facade::Cluster& cluster,
                            IDetectorVolumes::pointer dv,
                            Clus::IPCTransformSet::pointer pcts,
                            Weighted::Graph& graph,
                            const CtpcFastCfg* fast = nullptr);

    void connect_graph_ctpc_with_reference(
        const Facade::Cluster& cluster,
        const Facade::Cluster& ref_cluster,
        IDetectorVolumes::pointer dv,
        Clus::IPCTransformSet::pointer pcts,
        Weighted::Graph& graph);

    void connect_graph_closely(const Facade::Cluster& cluster,
                               Weighted::Graph& graph, int num_neighbors = 5);

    void connect_graph_closely_pid(const Facade::Cluster& cluster,
                                   Weighted::Graph& graph);                           

    // doc 78 round 2: parameters of the busy-cluster lazy-walk mode of
    // connect_graph_relaxed, selected only by the "relaxed_fast" flavor
    // (ClusteringExamineBundles graph_name knob).  A cluster whose
    // connected-component count `num` is at or below busy_num_threshold
    // takes the legacy path bit-for-bit; above it, the per-pair path walk
    // (the O(num^2 x path-cm) cost, doc 78 sec 2) is evaluated lazily:
    // only pairs selected by the MST (plus the < 3 cm direct-edge pairs)
    // are walked, iterating kill-and-rebuild until the MST is stable.
    // The result is the MST of the walk-filtered graph computed by a
    // different route -- NOT proven bit-identical to the eager order
    // under Prim tie-breaking, which is why it sits behind the busy gate
    // and per-event final-PR validation (doc 78 sec 6).
    struct RelaxedFastCfg {
        size_t busy_num_threshold{200};
        // Safety valve: if the kill-and-rebuild loop exceeds this many
        // iterations, evaluate every remaining selected pair eagerly.
        size_t max_fixpoint_iters{200};
    };

    // ne' overclustering protection
    void connect_graph_relaxed(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts,
        Weighted::Graph& graph,
        const RelaxedFastCfg* fast = nullptr);

    void connect_graph_relaxed_pid(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts,
        Weighted::Graph& graph);

    // doc pr/53 sec 16: stricter fork of connect_graph_relaxed (S1+S2+S3),
    // selected only by ClusteringProtectBundle via graph_name "relaxed_strict".
    // Its kill predicate relaxed_strict_bad() is declared in the public
    // WireCellClus/Graphs.h (doctested there).
    //
    // doc pr/53 round 7 sec 18: image_check adds the S5 3D-image-support
    // OR-kill (relaxed_img_bad(), also in the public header) alongside
    // relaxed_strict_bad(); default false keeps this function byte-for-byte
    // identical to graph_name "relaxed_strict"'s behavior (round 6). The
    // "relaxed_strict_img" flavor (make_graph_relaxed_strict_img) is the
    // only caller that passes true.
    //
    // doc pr/56 round 2: two_d_check adds S6, an additional per-plane 2D
    // wind/tick fired-pixel connectivity OR-kill (two_d_connectivity_bad(),
    // public header) alongside relaxed_strict_bad()/relaxed_img_bad();
    // default false keeps this function byte-for-byte identical to
    // "relaxed_strict_img"'s behavior. The "relaxed_strict_img_2d" flavor
    // (make_graph_relaxed_strict_img_2d) is the only caller that passes
    // true. S6 can only additionally kill an edge S1-S5 already let
    // through -- it never rescues one they killed.
    //
    // doc pr/58: floor_w_override (no-op unless two_d_check is also true)
    // lets an unexcused W-plane-only S6 gap kill an edge even below
    // s6_dis_floor -- W (collection) is far more robust against the
    // induction-plane signal-inefficiency false positives the floor exists
    // to guard against, so a genuine W gap is trusted at any distance.
    // Default false keeps "relaxed_strict_img_2d" byte-for-byte unchanged;
    // the "relaxed_strict_img_2d_wfloor" flavor
    // (make_graph_relaxed_strict_img_2d_wfloor) is the only caller that
    // passes true. Owner-requested investigation, evt 21073 j=6/k=7
    // (0.98cm, W-only gap at the neutrino vertex) is the motivating case.
    //
    // doc pr/57 round 6: two_d_rescue (no-op unless two_d_check is also
    // true) re-examines every S6-killed candidate with
    // Graphs::two_d_rescue_ok() -- a rule fitted against the owner's full
    // hand scan of the S6 separation displays (899 labels; doc pr/57 sec
    // 14) -- and UN-kills the ones it explains as artifacts (dead-W band,
    // prolonged induction signal, direction-consistent long-track breaks,
    // two-view co-location).  It can only remove kills S6 made, never add
    // one, and only for candidates <= 5 cm.  Default false keeps
    // "relaxed_strict_img_2d_wfloor" byte-for-byte unchanged; the
    // "relaxed_strict_img_2d_rescue" flavor
    // (make_graph_relaxed_strict_img_2d_rescue) is the only caller that
    // passes true.
    //
    // doc pr/62: long_check adds S7, an independent additional OR-kill in
    // all three path-check blocks, on exactly the band S6 declines to judge
    // -- candidates at or above s6_dis_cap (30cm), which today reach the MST
    // governed by S1-S3 alone. Graphs::long_corridor_bad() (public header)
    // is the verdict; the evidence is a per-plane bounded flood fill
    // restricted to a narrow corridor around the candidate's own two
    // endpoint lattice cells, which is what makes it O(D) where S6's
    // open-rectangle fill is O(D^2). It can only additionally kill an edge
    // S1-S6 already let through; it never rescues one they killed, and
    // there is deliberately no S7 analogue of two_d_rescue's
    // s6_rescued_dir_pairs emission repair (that repair exists only because
    // a RESCUE could evaporate a dir emission -- S7 never rescues).
    // long_min_planes selects the operating point (default 1, mirroring S6's
    // owner rule; 2 is the conservative fallback -- doc pr/62 has no
    // hand-scan labels in this band yet, unlike S6's, so both are shipped as
    // separate flavors rather than one being guessed as final). Default
    // false/1 keeps this function byte-for-byte identical to
    // "relaxed_strict_img_2d_rescue"'s behavior; "relaxed_strict_img_2d_rescue_long"
    // and "relaxed_strict_img_2d_rescue_long2"
    // (make_graph_relaxed_strict_img_2d_rescue_long{,2}) are the only callers
    // that pass long_check=true. long_check does not require two_d_check:
    // the two tests partition the distance axis at 30cm and cannot both fire
    // on one candidate.
    //
    // doc pr/64 round 4: w_track_excuse (no-op unless two_d_check is also
    // true) re-examines every candidate STILL killed after two_d_rescue with
    // Graphs::two_d_w_track_ok() (public header) -- the W-plane long-track
    // exception, fitted against the same 899-label hand scan at the owner's
    // tightened operating point -- and UN-kills the ones where W is the sole
    // voting plane on a long, thin, globally-collinear track pair (or a
    // dead-W band explains the gap).  Revive-only, like two_d_rescue, but
    // with NO <= 5 cm cap: the exposed population is every killed S6
    // candidate (S6 itself caps at 30 cm).  Default false keeps
    // "relaxed_strict_img_2d_rescue_long" byte-for-byte unchanged; the
    // "relaxed_strict_img_2d_rescue_long_wtrack" flavor
    // (make_graph_relaxed_strict_img_2d_rescue_long_wtrack) is the only
    // caller that passes true.
    void connect_graph_relaxed_strict(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts,
        Weighted::Graph& graph,
        bool image_check = false,
        bool two_d_check = false,
        bool floor_w_override = false,
        bool two_d_rescue = false,
        bool long_check = false,
        int long_min_planes = 1,
        bool w_track_excuse = false);

    bool is_point_good(const Facade::Cluster& cluster, size_t point_index, int ncut = 3);

    std::vector<bool> check_direction(const Facade::Cluster& cluster, Facade::geo_vector_t& v1, int apa, int face, double angle_cut_1 = 12.5, double angle_cut_2 = 10);

    bool check_connectivity(const Facade::Cluster& cluster,  IDetectorVolumes::pointer dv, 
    IPCTransformSet::pointer pcts, std::tuple<int, int, double>& index_index_dis,  std::shared_ptr<Facade::Simple3DPointCloud> pc1, std::vector<size_t> pc1_global_index, std::shared_ptr<Facade::Simple3DPointCloud> pc2, std::vector<size_t> pc2_global_index,
    double step_size = 0.6*units::cm, bool flag_strong_check = false);
}
    
#endif
