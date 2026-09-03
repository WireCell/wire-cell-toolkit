#ifndef WIRECELLCLUS_PRIVATE_MAKE_GRAPHS
#define WIRECELLCLUS_PRIVATE_MAKE_GRAPHS

#include "WireCellClus/Graphs.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellClus/IPCTransform.h"
#include "WireCellClus/Facade_Cluster.h"


namespace WireCell::Clus::Graphs {

    // factory functions wrapping up construction and various connect_graph*
    // functions.

    // just closely connected.
    Weighted::Graph make_graph_closely(
        const Facade::Cluster& cluster);

    // just closely connected.
    Weighted::Graph make_graph_closely_pid(
        const Facade::Cluster& cluster);

    // closely + basic connection
    Weighted::Graph make_graph_basic(
        const Facade::Cluster& cluster);

    // closely_pid + basic connection with reference cluster (empty by default)
    Weighted::Graph make_graph_basic_pid(
        const Facade::Cluster& cluster,
        const Facade::Cluster& ref_cluster = Facade::Cluster{});

    // closely + ctpc connection
    Weighted::Graph make_graph_ctpc(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts);

    // doc 79 round 2: same algorithm with the busy-cluster lazy walk armed
    // (CtpcFastCfg defaults, connect_graphs.h).  Selected only by the
    // "ctpc_fast" flavor (ClusteringDeghost graph_name knob); clusters at
    // or below the busy threshold are bit-identical to make_graph_ctpc.
    Weighted::Graph make_graph_ctpc_fast(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    Weighted::Graph make_graph_ctpc_pid(
        const Facade::Cluster& cluster,
        const Facade::Cluster& ref_cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts);

    // closely + relaxed (overclustering protection)
    Weighted::Graph make_graph_relaxed(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv, 
        IPCTransformSet::pointer pcts);

    // doc 78 round 2: same algorithm with the busy-cluster lazy walk armed
    // (RelaxedFastCfg defaults, connect_graphs.h).  Selected only by the
    // "relaxed_fast" flavor (ClusteringExamineBundles graph_name knob);
    // clusters at or below the busy threshold are bit-identical to
    // make_graph_relaxed.
    Weighted::Graph make_graph_relaxed_fast(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    Weighted::Graph make_graph_relaxed_pid(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    // closely + relaxed_strict (doc pr/53 sec 16: protect_bundle-only
    // stricter overclustering protection)
    Weighted::Graph make_graph_relaxed_strict(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    // closely + relaxed_strict with image_check=true (doc pr/53 round 7 sec
    // 18: protect_bundle-only, adds the S5 3D-image-support OR-kill)
    Weighted::Graph make_graph_relaxed_strict_img(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    // closely + relaxed_strict with image_check=true AND two_d_check=true
    // (doc pr/56 round 2: protect_bundle-only, adds S6 the per-plane 2D
    // wind/tick fired-pixel connectivity OR-kill on top of S1-S5)
    Weighted::Graph make_graph_relaxed_strict_img_2d(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    // closely + relaxed_strict with image_check=true, two_d_check=true AND
    // floor_w_override=true (doc pr/58: protect_bundle-only, owner-requested
    // -- an unexcused W-plane-only S6 gap kills even below s6_dis_floor,
    // since W is far more robust against the induction-plane false
    // positives the floor otherwise guards against)
    Weighted::Graph make_graph_relaxed_strict_img_2d_wfloor(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    // closely + relaxed_strict with image_check, two_d_check,
    // floor_w_override AND two_d_rescue all true (doc pr/57 round 6:
    // protect_bundle-only -- S6-killed candidates that
    // Graphs::two_d_rescue_ok() explains as detector artifacts are
    // un-killed; fitted against the owner's full separation hand scan)
    Weighted::Graph make_graph_relaxed_strict_img_2d_rescue(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    // closely + relaxed_strict with image_check, two_d_check,
    // floor_w_override, two_d_rescue AND long_check all true, long_min_planes=1
    // (doc pr/62: protect_bundle-only -- adds S7, corridor connectivity on
    // candidates at or above the 30cm band S6 skips; owner's S6-style ">= 1
    // non-excused plane" rule)
    Weighted::Graph make_graph_relaxed_strict_img_2d_rescue_long(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    // same as make_graph_relaxed_strict_img_2d_rescue_long but
    // long_min_planes=2 (doc pr/62: the conservative fallback operating
    // point -- S1-S5's "at least 2 views" convention -- shipped alongside
    // the owner-rule variant since this distance band has no hand-scan
    // labels to fit against yet)
    Weighted::Graph make_graph_relaxed_strict_img_2d_rescue_long2(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    // same as make_graph_relaxed_strict_img_2d_rescue_long plus
    // w_track_excuse=true (doc pr/64 round 4: the W-plane long-track
    // exception -- revives S6 kills where W is the sole voting plane on a
    // long, thin, globally-collinear track pair, or a dead-W band explains
    // the gap; Graphs::two_d_w_track_ok is the pure verdict)
    Weighted::Graph make_graph_relaxed_strict_img_2d_rescue_long_wtrack(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

    // doc clus/connect-graph-strict-perf round 2: identical to
    // ..._long_wtrack, plus the doc 78 busy-gated lazy walk (RelaxedFastCfg
    // defaults, busy_num_threshold=200).  Separate flavor NAME rather than a
    // config key, so a job that does not ask for it is byte-identical by
    // construction -- the same wiring choice doc 78 round 2 made for
    // "relaxed_fast".
    Weighted::Graph make_graph_relaxed_strict_img_2d_rescue_long_wtrack_fast(
        const Facade::Cluster& cluster,
        IDetectorVolumes::pointer dv,
        IPCTransformSet::pointer pcts);

}

#endif
