local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import 'pgrapher/common/funcs.jsonnet';
local clus = import "pgrapher/common/clus.jsonnet";


// Top-level function: parameters overridable by the importing entry point.
//   output_dir: directory for mabc-*.zip and trash-*.tar.gz ('' = ./data)
//   time_offset: event-T0 / readout-tick0 compensation applied by the live
//     BlobSampler when converting slice time to drift x
//     (x = xorig + dirx*(t_slice + time_offset)*drift_speed).  Readout tick0
//     sits at -250us relative to the trigger, so -250us would place
//     trigger-time activity at its true x -- but PDVD has no per-event T0
//     determination, so the default is 0: no activity maps behind the anode
//     and live points stay consistent with the (offset-free) dead sampler.
//     Restore a real value once a T0 measurement exists.
//   relax_containment_filter: make the T0Correction scope filter (used by
//     switch_scope in the all-APA stage) accept containment in ANY sensitive
//     volume instead of each point's own (apa,face) volume.  Without a T0,
//     out-of-time activity apparently sits past the cathode in the opposite
//     drift volume; the strict filter would exclude such clusters from ALL
//     all-APA merge passes (see pdvd/docs/clustering-boundary-merge.md).
//   trigger_offset: per-event light-vs-charge time-base offset (the calibrated
//     PDVD analogue of the PDHD opflash offset_us) applied DOWNSTREAM rather
//     than baked into x_raw -- with time_offset=0 the raw imaging x stays
//     offset-free, and T0Correction (x_t0cor scope) adds (cluster_t0 +
//     trigger_offset) so flash-matched charge lands on the light time base.
//     QLMatching folds the same value into its matching geometry.  Default 0
//     => x_t0cor unchanged (bit-identical).
//   trigger_offset_top: optional separate offset for the TOP drift volume
//     (anodes 4-7).  The PDVD TDE/BDE charge-readout crates open their windows
//     up to ~32 us apart, so the light->charge offset is per crate:
//     trigger_offset = bottom (BDE), trigger_offset_top = top (TDE).  null
//     (default) => trigger_offset everywhere (bit-identical).
function (output_dir='', runNo=1, subRunNo=1, eventNo=1, stepped_center_fallback=false,
          time_offset=0 * wc.us, relax_containment_filter=true,
          trigger_offset=0 * wc.us, trigger_offset_top=null,
          drift_speed_b=null, drift_speed_t=null)

// Calibrated from PDVD data (anode->cathode crossing tracks: reconstructed drift
// x-span vs the collection-plane->cathode-surface distance 338.55 cm; cathode
// surface corrected 2.54->3.0 cm, so 1.57->1.568, v proportional to D); was 1.6.
// See pdvd/docs/clus-workflow.md (drift-velocity calibration).
local drift_speed = 1.568 * wc.mm / wc.us;
// Optional per-side calibrated speeds (sec 8.12 of
// pdvd-anode-time-consistency.md): bottom volume = anodes 0-3, top = 4-7.
// null => the common legacy value above; compiled config byte-identical when
// both are null.  These feed BOTH the per-face BlobSampler (x_raw) and the
// dvm blocks below (T0Correction / grouping x_t0cor) from the same resolved
// value, so the two coordinates cannot drift apart.
local drift_speed_bot = if drift_speed_b == null then drift_speed else drift_speed_b;
local drift_speed_top = if drift_speed_t == null then drift_speed else drift_speed_t;

local initial_index = "0";
local index = std.parseInt(initial_index);


local common_coords = ["x", "y", "z"];
local common_corr_coords = ["x_t0cor", "y", "z"];

// Drift-side groups: anodes 0-3 (bottom drift, x0=-341.5cm) and anodes 4-7
// (top drift, x0=+341.5cm).  Each anode's two faces fold into its group
// automatically (routing is per-anode via wpid.apa()).  These define both the
// stage-3 per-drift-group clustering scope (clus_per_group) and the Bee
// display groups (clustering points + dead-area output in the all-TPC dump),
// so each drift side shows as a single Bee instance.
local apa_drift_groups = [
    { name: "group0123", apas: [0, 1, 2, 3] },
    { name: "group4567", apas: [4, 5, 6, 7] },
];


// ProtoDUNE-VD geometry parameters
// 8 anodes total: anodes 0-3 are bottom drift (centerline x=-3415.5mm, drift in +x direction)
//                 anodes 4-7 are top    drift (centerline x=+3415.5mm, drift in -x direction)
// apa_plane = 0.4mm (first-induction U plane, PDHD convention; 2026-07-12,
// was 57.15mm = the fictional DocDB-203 grid plane -- see params.jsonnet and
// pdvd/docs/qlmatch/pdvd-anode-time-consistency.md), cpa_plane = 3385.5mm
// Bottom drift anode face x ~ -3415.1mm, cathode x ~ -30.0mm
// Top    drift anode face x ~  3415.1mm, cathode x ~   30.0mm
local dvm = {
    overall: {
        FV_xmin: -3415.5 * wc.mm,
        FV_xmax:  3415.5 * wc.mm,
        FV_ymin: -3214.0 * wc.mm,  // active -3364.0 mm + 15 cm inset (see clus/docs/clustering-separate-fv.md)
        FV_ymax:  3214.0 * wc.mm,  // active 3364.0 mm - 15 cm inset
        FV_zmin: 150.5 * wc.mm,    // active 0.5 mm + 15 cm inset
        FV_zmax: 2842.5 * wc.mm,   // active 2992.5 mm - 15 cm inset
        FV_xmin_margin: 2 * wc.cm,
        FV_xmax_margin: 2 * wc.cm,
        FV_ymin_margin: 2.5 * wc.cm,
        FV_ymax_margin: 2.5 * wc.cm,
        FV_zmin_margin: 3 * wc.cm,
        FV_zmax_margin: 3 * wc.cm,
        vertical_dir: [0,1,0],
        beam_dir: [0,0,1]
    },
    // Bottom drift (anodes 0-3): anode face at ~-3415.1mm, cathode at ~-30.0mm
    // Both faces share same x-bounds (2-sided CRP, same drift direction)
    a0f0pA: {
        drift_speed: drift_speed_bot,
        tick: 0.5 * wc.us,  // 0.5 mm per tick
        tick_drift: self.drift_speed * self.tick,
        time_offset: time_offset,
        trigger_offset: trigger_offset,
        nticks_live_slice: 4,
        FV_xmin: -3415.1 * wc.mm,
        FV_xmax: -30.0 * wc.mm,
        FV_xmin_margin: 2 * wc.cm,
        FV_xmax_margin: 2 * wc.cm,
    },
    a0f1pA: $.a0f0pA,
    a1f0pA: $.a0f0pA,
    a1f1pA: $.a0f0pA,
    a2f0pA: $.a0f0pA,
    a2f1pA: $.a0f0pA,
    a3f0pA: $.a0f0pA,
    a3f1pA: $.a0f0pA,
    // Top drift (anodes 4-7): cathode at ~30.0mm, anode face at ~3415.1mm
    a4f0pA: $.a0f0pA + {
        drift_speed: drift_speed_top,
        FV_xmin: 30.0 * wc.mm,
        FV_xmax: 3415.1 * wc.mm,
        // Top volume charge is read by the TDE crate (its own window start).
        trigger_offset: if trigger_offset_top == null then trigger_offset
                        else trigger_offset_top,
    },
    a4f1pA: $.a4f0pA,
    a5f0pA: $.a4f0pA,
    a5f1pA: $.a4f0pA,
    a6f0pA: $.a4f0pA,
    a6f1pA: $.a4f0pA,
    a7f0pA: $.a4f0pA,
    a7f1pA: $.a4f0pA,
};

local anodes_name(anodes, face="") =
    std.join("-", [std.toString(a.data.ident) for a in anodes]) + if face == "" then "" else "-" + std.toString(face);


local detector_volumes(anodes, face="") = {
    "type": "DetectorVolumes",
    "name": "dv-apa" + anodes_name(anodes, face),
    "data": {
        "anodes": [wc.tn(anode) for anode in anodes],
        metadata:
            {overall: dvm["overall"]} +
            {
                [ "a" + std.toString(a.data.ident) + "f0pA" ]:
                    dvm[ "a" + std.toString(a.data.ident) + "f0pA" ]
                for a in anodes
            } +
            {
                [ "a" + std.toString(a.data.ident) + "f1pA" ]:
                    dvm[ "a" + std.toString(a.data.ident) + "f1pA" ]
                for a in anodes
            }
    },
    uses: anodes
};


local pctransforms(dv) = {
    type: "PCTransformSet",
    name: dv.name,
    data: {
        detector_volumes: wc.tn(dv),
        relax_containment_filter: relax_containment_filter,
    },
    uses: [dv]
};



local bs_live_face(apa, face, center_fallback=false, speed=drift_speed) = {
    type: "BlobSampler",
    name: "live-%s-%d"%[apa, face],
    data: {
        drift_speed: speed,
        time_offset: time_offset,
        // center_fallback: emit one point at the blob center when the stepped
        // grid yields none (tiny 1-wire blobs); default off -> bit-identical.
        strategy: [{name: "stepped", center_fallback: center_fallback}],
        extra: [".*wire_index", ".*charge_val", ".*charge_unc", "wpid"]
    }
};
local bs_dead_face(apa, face) = {
    type: "BlobSampler",
    name: "dead-%s-%d"%[apa, face],
    data: {
        strategy: ["center"],
        extra: [".*"] // want all the extra
    }
};
// The factory used to give blob samplers to ClusteringRetile ("rt").
local bs_rt_face = bs_live_face;


local clus_per_face (
    anode,
    face,
    dump = true,
    bee_dir = "data",
    runNo = 1,
    subRunNo = 1,
    eventNo = 1,
    stepped_center_fallback = false,
    ) =
{

    local dv = detector_volumes([anode], face),
    local pcts = pctransforms(dv),


    local cluster_scope_filter_live = g.pnode({
        type: "ClusterScopeFilter",
        name: "csf-live-%s-%d"%[anode.name, face],
        data: {
            face_index: face,
        }
    }, nin=1, nout=1, uses=[]),

    local cluster_scope_filter_dead = g.pnode({
        type: "ClusterScopeFilter",
        name: "csf-dead-%s-%d"%[anode.name, face],
        data: {
            face_index: face,
        }
    }, nin=1, nout=1, uses=[]),

    local bsl = bs_live_face(anode.name, face, center_fallback=stepped_center_fallback,
                             speed=if anode.data.ident < 4 then drift_speed_bot
                                   else drift_speed_top),
    local bsd = bs_dead_face(anode.name, face),

    local ptb = g.pnode({
        type: "PointTreeBuilding",
        name: "%s-%d"%[anode.name, face],
        data:  {
            samplers: {
                "3d": wc.tn(bsl),
                "dead": wc.tn(bsd),
            },
            multiplicity: 2,
            tags: ["live", "dead"],
            anode: wc.tn(anode),
            face: face,
            detector_volumes: wc.tn(dv),
        }
    }, nin=2, nout=1, uses=[bsl, bsd, dv]),

    local cluster2pct = g.intern(
        innodes = [cluster_scope_filter_live, cluster_scope_filter_dead],
        centernodes = [],
        outnodes = [ptb],
        edges = [
            g.edge(cluster_scope_filter_live, ptb, 0, 0),
            g.edge(cluster_scope_filter_dead, ptb, 0, 1)
        ]
    ),
    // local cluster2pct = ptb,

    local face_name = "%s-%d"%[anode.name, face],

    local cm = clus.clustering_methods(prefix=face_name,
                                       detector_volumes=dv,
                                       pc_transforms=pcts,
                                       coords=common_coords),
    local cm_pipeline = [
        cm.pointed(),
        // cm.ctpointcloud(),
        cm.live_dead(dead_live_overlap_offset=2),
        cm.extend(flag=4, length_cut=60*wc.cm, num_try=0, length_2_cut=15*wc.cm, num_dead_try=1),
        cm.regular(name="-one", length_cut=60*wc.cm, flag_enable_extend=false),
        cm.regular(name="_two", length_cut=30*wc.cm, flag_enable_extend=true),
        cm.parallel_prolong(length_cut=35*wc.cm),
        cm.close(length_cut=1.2*wc.cm),
        cm.extend_loop(num_try=3),
        // separate moved to the per-drift-group stage (stage 3)
        cm.connect1(),
        // cm.isolated(),
        // cm.retile(cut_time_low=3*wc.us, cut_time_high=5*wc.us, anodes=[anode], samplers=[clus.sampler(bsl, apa=anode.data.ident, face=face)]),
    ],

    local mabc = g.pnode({
        local name = "%s-%d"%[anode.name, face],
        type: "MultiAlgBlobClustering",
        name: name,
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            // grouping2file_prefix: "grouping%s-%d"%[anode.name, face],
            perf: true,
            bee_dir: bee_dir, // "data/0/0", // not used
            bee_zip: "%s/mabc-%s-face%d.zip"%[bee_dir, anode.name, face],
            bee_detector: "protodunevd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            dead_area_version: 2,  // v2 wrapper (tpc=apa) so the dead slab lands on the correct PDVD anode face
            anodes: [wc.tn(anode)],
            face: face,
            detector_volumes: wc.tn(dv),
            bee_points_sets: [  // New configuration for multiple bee points sets
                {
                    name: "clustering",         // Name of the bee points set
                    detector: "protodunevd",         // Detector name
                    algorithm: "clustering",    // Algorithm identifier
                    pcname: "3d",           // Which scope to use
                    coords: ["x", "y", "z"],    // Coordinates to use
                    individual: true            // Output individual APA/Face
                }
            ],
            pipeline: wc.tns(cm_pipeline),
        }
    }, nin=1, nout=1, uses=[dv, anode, pcts]+cm_pipeline),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "clus_per_face-%s-%d"%[anode.name, face],
        data: {
            outname: "%s/trash-%s-face%d.tar.gz"%[bee_dir, anode.name, face],
            prefix: "clustering_", // json, numpy, dummy
            dump_mode: true,
        }
    }, nin=1, nout=0),

    local end = if dump
    then g.pipeline([mabc, sink])
    else g.pipeline([mabc]),

    ret :: g.pipeline([cluster2pct, end], "clus_per_face-%s-%d"%[anode.name, face])
}.ret;

local clus_per_apa (
    anode,
    dump = true,
    bee_dir = "data",
    runNo = 1,
    subRunNo = 1,
    eventNo = 1,
    stepped_center_fallback = false,
    ) =
{
    local cfout_live = g.pnode({
        type:'ClusterFanout',
        name: 'clus_per_apa-cfout_live-%s'%anode.name,
        data: {
            multiplicity: 2
        }}, nin=1, nout=2),

    local cfout_dead = g.pnode({
        type:'ClusterFanout',
        name: 'clus_per_apa-cfout_dead-%s'%anode.name,
        data: {
            multiplicity: 2
        }}, nin=1, nout=2),

    local per_face_pipes = [
        clus_per_face(anode, face=0, dump=false, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo, stepped_center_fallback=stepped_center_fallback),
        clus_per_face(anode, face=1, dump=false, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo, stepped_center_fallback=stepped_center_fallback),
    ],

    local pcmerging = g.pnode({
        type: "PointTreeMerging",
        name: "%s"%[anode.name],
        data:  {
            multiplicity: 2,
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
        }
    }, nin=2, nout=1),

    local dv = detector_volumes([anode]),
    local pcts = pctransforms(dv),

    local cm = clus.clustering_methods(prefix=anode.name,
                                       detector_volumes=dv,
                                       pc_transforms=pcts,
                                       coords=common_coords),
    local cm_pipeline = [
        cm.deghost(),
        cm.protect_overclustering(),
    ],

    local mabc = g.pnode({
        local name = anode.name,
        type: "MultiAlgBlobClustering",
        name: "clus_per_apa-%s"%[name],
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            // grouping2file_prefix: "grouping%s-%d"%[anode.name, face],
            perf: true,
            bee_dir: bee_dir, // "data/0/0", // not used
            bee_zip: "%s/mabc-%s.zip"%[bee_dir, anode.name],
            bee_detector: "protodunevd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            dead_area_version: 2,  // v2 wrapper (tpc=apa) so the dead slab lands on the correct PDVD anode face
            anodes: [wc.tn(anode)],
            detector_volumes: wc.tn(dv),
            bee_points_sets: [
                {
                    name: "clustering",
                    detector: "protodunevd",
                    algorithm: "clustering",
                    pcname: "3d",
                    coords: ["x", "y", "z"],
                    individual: false,
                }
            ],
            pipeline: wc.tns(cm_pipeline),
        }
    }, nin=1, nout=1, uses=[anode, dv, pcts]+cm_pipeline),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "clus_per_apa-%s"%[anode.name],
        data: {
            outname: "%s/trash-%s.tar.gz"%[bee_dir, anode.name],
            prefix: "clustering_", // json, numpy, dummy
            dump_mode: true,
        }
    }, nin=1, nout=0),

    local end = if dump
    then g.pipeline([mabc, sink])
    else g.pipeline([mabc]),

    ret :: g.intern(
        innodes = [cfout_live, cfout_dead],
        centernodes = per_face_pipes + [pcmerging],
        outnodes = [end],
        edges = [
            g.edge(cfout_live, per_face_pipes[0], 0, 0),
            g.edge(cfout_dead, per_face_pipes[0], 0, 1),
            g.edge(cfout_live, per_face_pipes[1], 1, 0),
            g.edge(cfout_dead, per_face_pipes[1], 1, 1),
            g.edge(per_face_pipes[0], pcmerging, 0, 0),
            g.edge(per_face_pipes[1], pcmerging, 0, 1),
            g.edge(pcmerging, end, 0, 0),
        ]
    ),
}.ret;

// Stage 3: per drift-side group (anodes 0-3 bottom drift, anodes 4-7 top
// drift), raw x,y,z coordinates.  Runs the long-range merge family across the
// four CRPs of one drift side, then the topology passes (separate moved here
// from the per-face stage; examine_x_boundary requires all wpids in the group
// to share the same FV_x metadata, true within one drift side).
local clus_per_group (
    anodes,
    group_name,
    dump = true,
    bee_dir = "data",
    runNo = 1,
    subRunNo = 1,
    eventNo = 1,
    ) = {
    local nanodes = std.length(anodes),
    local pcmerging = g.pnode({
        type: "PointTreeMerging",
        name: "clus_per_group-%s"%group_name,
        data:  {
            multiplicity: nanodes,
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            tolerate_missing: true,
        }
    }, nin=nanodes, nout=1),

    local dv = detector_volumes(anodes),
    local pcts = pctransforms(dv),

    local cm = clus.clustering_methods(prefix=group_name,
                                       detector_volumes=dv,
                                       pc_transforms=pcts,
                                       coords=common_coords),
    local cm_pipeline = [
        cm.extend(flag=4, length_cut=60*wc.cm, num_try=0, length_2_cut=15*wc.cm, num_dead_try=1),
        cm.regular(name="1", length_cut=60*wc.cm, flag_enable_extend=false),
        cm.regular(name="2", length_cut=30*wc.cm, flag_enable_extend=true),
        cm.parallel_prolong(length_cut=35*wc.cm),
        cm.close(length_cut=1.2*wc.cm),
        cm.extend_loop(num_try=3),
        // max_hull_points raised from the 10k default: get_hull() bailing silently
        // disabled the separation decision on exactly the clusters that need it
        // (a 102k-pt giant in PDVD 39324 evt 339890 slipped past the earlier 100k).
        cm.separate(use_ctpc=true, max_hull_points=1000000, collinear_recover=true, collinear_interior=true,
                    collinear_member_merge=true,
                    track_repartition=true, band_merge_back=true, band_recarve=true, drift_side_fv_x=true,
                    far_point_x_cut=14*wc.cm, far_point_mid_dis=60*wc.cm, track_recarve=true, dec1_guard_main_angle=45,
                    iso_slab_split=true, tag_family=true, collinear_global_merge=true),
        // MicroBooNE order after separate: connect1 (reconnect dashed-line
        // fragments, e.g. drift-direction tracks split across the group) then
        // deghost (remove ghosts that only the group scope can adjudicate).
        // A PDVD drift group mixes faces by construction (an anode's two
        // faces are the y-halves of one CRP, sharing one drift volume and
        // identical FV_x metadata), so waive the same-face check here and in
        // examine_x_boundary below.
        // empty_view_unique: required at group scope -- see common clus.jsonnet.
        cm.connect1(allow_mixed_faces=true, respect_separate_family=true),
        cm.deghost(allow_mixed_faces=true, empty_view_unique=true),
        cm.examine_x_boundary(allow_mixed_faces=true),
        cm.neutrino(protect_iso_band=true),
        cm.isolated(),
        // examine_bundles() DISABLED at the per-drift-group stage (2026-06-14),
        // matching the PDHD chain.  At this stage (use_flash_t0=false) it only
        // rewrites the "isolated"/"perblob" array (main -1 / associated >=0); it
        // never changes cluster membership.  PDVD runs no Q/L matching, so that
        // perblob array is never consumed downstream -- removing this pass leaves
        // the cluster grouping byte-identical (verified A/B on a PDVD event).  On
        // PDHD the same change fixes a boundary-flag bug where QLMatching flagged
        // only the main sub-cluster; see cfg/.../pdhd/clus.jsonnet and
        // pdhd/docs/clustering-algorithm.md.  Re-enable in all-TPC if needed.
        // cm.examine_bundles(),
    ],

    local mabc = g.pnode({
        type: "MultiAlgBlobClustering",
        name: "clus_per_group-%s"%group_name,
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            perf: true,
            bee_dir: bee_dir, // "data/0/0", // not used
            bee_zip: "%s/mabc-%s.zip"%[bee_dir, group_name],
            bee_detector: "protodunevd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            dead_area_version: 2,  // v2 wrapper (tpc=apa) so the dead slab lands on the correct PDVD anode face
            anodes: [wc.tn(a) for a in anodes],
            detector_volumes: wc.tn(dv),
            bee_points_sets: [
                {
                    name: "clustering",
                    detector: "protodunevd",
                    algorithm: "clustering",
                    pcname: "3d",
                    coords: ["x", "y", "z"],
                    individual: false,
                }
            ],
            pipeline: wc.tns(cm_pipeline),
        }
    }, nin=1, nout=1, uses=anodes+[dv, pcts]+cm_pipeline),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "clus_per_group-%s"%group_name,
        data: {
            outname: "%s/trash-%s.tar.gz"%[bee_dir, group_name],
            prefix: "clustering_", // json, numpy, dummy
            dump_mode: true,
        }
    }, nin=1, nout=0),

    local end = if dump
    then g.pipeline([mabc, sink])
    else g.pipeline([mabc]),

    ret :: g.intern(
        innodes = [pcmerging],
        centernodes = [],
        outnodes = [end],
        edges = [
            g.edge(pcmerging, end, 0, 0),
        ]
    ),
}.ret;

// Stage 4: all-TPC, merging the two drift-side groups.  switch_scope applies
// the per-cluster T0 correction (x_t0cor; with no flash matching on PDVD this
// is the apparent x) and its containment scope filter.
local clus_all_tpc (
    anodes,
    ngroups = 2,
    dump = true,
    bee_dir = "data",
    runNo = 1,
    subRunNo = 1,
    eventNo = 1,
    save_opflash = false,
    premerged = false,   // skip the input PointTreeMerging (joint QLMatching already merged)
    ) = {
    local pcmerging = g.pnode({
        type: "PointTreeMerging",
        name: "clus_all_tpc",
        data:  {
            multiplicity: ngroups,
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            tolerate_missing: true,
        }
    }, nin=ngroups, nout=1),

    local dv = detector_volumes(anodes),
    local pcts = pctransforms(dv),


    local cm_old = clus.clustering_methods(prefix="all",
                                           detector_volumes=dv,
                                           pc_transforms=pcts,
                                           coords=common_coords),


    local cm = clus.clustering_methods(prefix="all",
                                       detector_volumes=dv,
                                       pc_transforms=pcts,
                                       coords=common_corr_coords),

    local cm_pipeline = [
        cm_old.switch_scope(),

        // Cathode-crossing connector (SBND-tuned parameters as placeholder;
        // PDVD central cathode is at x=0, the C++ default cathode_x — the
        // sensitive volumes end at -+30.0mm so cathode_x_cut=5cm spans the
        // |x|<3.0cm gap).  Without Q/L matching (premerged=false)
        // use_flash_t0=false (the flash-coincidence gate would veto every
        // pair).  With Q/L matching upstream (premerged=true) every cluster
        // carries its matched cluster_t0: gate pairs on flash-time coincidence.
        // NB PDVD has ONE all-PD flash, so a genuine crosser's two halves
        // normally share the SAME matched flash — the 1 us window is a
        // formality that also admits a rare split match.
        cm.cathode_connect(cathode_x_cut=5*wc.cm, drift_cut=8*wc.cm,
                           min_length_short=2*wc.cm, short_dir_len=25*wc.cm,
                           conn_short_cut=30.0, use_flash_t0=premerged,
                           flash_t0_window=1*wc.us),
    ],

    local mabc = g.pnode({
        type: "MultiAlgBlobClustering",
        name: "clus_all_tpc",
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            // grouping2file_prefix: "grouping%s-%d"%[anode.name, face],
            perf: true,
            bee_dir: bee_dir, // "data/0/0", // not used
            bee_zip: "%s/mabc-all-apa.zip"%[bee_dir],
            bee_detector: "protodunevd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            // Dump the optical "op" bee instance (measured flash PE + Q/L
            // predicted PE per matched cluster) at the pre-pipeline point.
            // Needs do_qlmatch (the QLMatching "opflash" root PC); no-op
            // otherwise.  Default false => bit-identical.
            save_opflash: save_opflash,
            // One "op" row per flash carrying ALL its matched cluster ids
            // (summed predicted PE) — with PDVD's single all-PD flash a flash
            // is typically matched to many clusters, so per-flash rows are the
            // only readable display (SBND/PDHD parity).
            bee_flash_per_flash: true,
            // No flash-flash grouping: PDVD has ONE all-PD flash per time (no
            // per-side pair to group).  0 = off (no column).
            flash_group_window: 0,
            dead_area_version: 2,  // v2 wrapper (tpc=apa) so the dead slab lands on the correct PDVD anode face
            dead_apa_groups: apa_drift_groups,  // group dead area by drift side -> 2 dead instances
            anodes: [wc.tn(a) for a in anodes],
            detector_volumes: wc.tn(dv),
            bee_points_sets: [  // New configuration for multiple bee points sets
            {
                    name: "clustering",         // Name of the bee points set
                    detector: "protodunevd",         // Detector name
                    algorithm: "clustering",    // Algorithm identifier
                    pcname: "3d",           // Which scope to use
                    // With Q/L matching upstream (premerged) draw the
                    // T0-corrected drift coords (x_t0cor = x_raw - dirx*
                    // (cluster_t0 + trigger_offset)*drift_speed, materialised
                    // by switch_scope from the QLMatching cluster_t0; an
                    // unmatched cluster fails the containment filter and is
                    // dropped here — it still appears in the raw "img" set).
                    // Without matching keep the raw coords (bit-identical).
                    coords: if premerged then common_corr_coords else ["x", "y", "z"],
                    individual: false            // Output individual APA/Face
                },
            {
                    // The hard-coded "img" name is MABC's pre-pipeline hook: this
                    // set is dumped from the live grouping BEFORE the all-TPC
                    // clustering pipeline runs -- i.e. the raw imaged charge (the
                    // per-anode matched clusters).  With no apa_groups it is
                    // dumped whole -> img-global (PDHD/SBND parity).  Its cluster
                    // ids are the SAME enumeration the pre-pipeline "op" dump
                    // writes into each flash's op_cluster_ids (fill_bee_flashes
                    // runs at the same point), so the Bee viewer's flash<->cluster
                    // pairing must be checked against THIS instance -- neither the
                    // post-pipeline clustering-global (re-enumerated) nor a
                    // separate bee-blobs imaging pass (own numbering) matches.
                    // The hook can host only one set, so emitting img-global here
                    // replaces the former per-drift-side stage-3
                    // clustering-group0123/4567 dump (cf. PDHD clus.jsonnet).
                    name: "img",
                    detector: "protodunevd",
                    algorithm: "img",           // -> img-global (raw imaged charge)
                    pcname: "3d",
                    coords: ["x", "y", "z"],    // raw drift coords, matching the op dump point
                    individual: false,
                }
            ],
            pipeline: wc.tns(cm_pipeline),
        },
    }, nin=1, nout=1, uses=anodes+[dv, pcts]+cm_pipeline),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "clus_all_tpc",
        data: {
            outname: "%s/trash-all-apa.tar.gz"%[bee_dir],
            prefix: "clustering_", // json, numpy, dummy
            dump_mode: true,
        }
    }, nin=1, nout=0),
    local end = if dump
    then g.pipeline([mabc, sink])
    else g.pipeline([mabc]),
    // premerged: the upstream joint QLMatching node already emits ONE merged tree, so
    // feed mabc directly (its single input is the all-TPC stage's input).  Else merge the
    // ngroups per-side outputs through PointTreeMerging first (historical path).
    ret :: if premerged then end else g.intern(
        innodes = [pcmerging],
        centernodes = [],
        outnodes = [end],
        edges = [
            g.edge(pcmerging, end, 0, 0),
        ]
    ),
}.ret;


{
    local bee_dir = if output_dir == '' then 'data' else output_dir,
    per_face(anode, face=0, dump=true) :: clus_per_face(anode, face=face, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo, stepped_center_fallback=stepped_center_fallback),
    per_apa(anode, dump=true) :: clus_per_apa(anode, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo, stepped_center_fallback=stepped_center_fallback),
    per_group(anodes, group_name, dump=true) :: clus_per_group(anodes, group_name, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo),
    all_tpc(anodes, ngroups=2, dump=true, save_opflash=false, premerged=false) :: clus_all_tpc(anodes, ngroups=ngroups, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo, save_opflash=save_opflash, premerged=premerged),
    // Expose the DetectorVolumes node builder so the Q/L matching graph can
    // reference the SAME all-anode DV the clustering uses (deterministic by name).
    detector_volumes(anodes, face="") :: detector_volumes(anodes, face),
}
