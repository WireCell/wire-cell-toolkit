local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import 'pgrapher/common/funcs.jsonnet';
local clus = import "pgrapher/common/clus.jsonnet";


// Top-level function: parameters overridable by the importing entry point.
//   output_dir: directory for mabc-*.zip and trash-*.tar.gz ('' = ./data)
//   time_offset: event-T0 / readout-tick0 compensation applied by the live
//     BlobSampler when converting slice time to drift x
//     (x = xorig + dirx*(t_slice + time_offset)*drift_speed) and recorded in
//     the DetectorVolumes metadata used by T0Correction (x_t0cor scope).
//     Readout tick0 sits at -250us relative to the trigger, so -250us would
//     place trigger-time activity at its true x -- but PDHD has no per-event
//     T0 determination, so the default is 0 (250us = 40cm at 1.6 mm/us): no
//     activity maps behind the anode and the points stay consistent with the
//     (offset-free) imaging frame.  Restore a real value once a T0
//     measurement exists.
function (output_dir='', runNo=1, subRunNo=1, eventNo=1,
          time_offset=0 * wc.us)

local drift_speed = 1.6 * wc.mm / wc.us;

local initial_index = "0";
local index = std.parseInt(initial_index);


local common_coords = ["x", "y", "z"];
local common_corr_coords = ["x_t0cor", "y", "z"];

// Drift-side groups: APA0+APA2 (face0, drift -x) and APA1+APA3 (face1, drift
// +x).  These define both the stage-3 per-drift-group clustering scope
// (clus_per_group) and the Bee display groups (clustering points + dead-area
// output in the all-TPC dump), so each pair shows as a single Bee instance.
local apa_drift_groups = [
    { name: "group02", apas: [0, 2] },
    { name: "group13", apas: [1, 3] },
];


local dvm = {
    overall: {
        FV_xmin: -3579.85 * wc.mm,
        FV_xmax: 3579.85 * wc.mm,
        FV_ymin: 76.1 * wc.mm,
        FV_ymax: 6060.0 * wc.mm,
        FV_zmin: 2.34345 * wc.mm,
        FV_zmax: 4622.97 * wc.mm,
        FV_xmin_margin: 2 * wc.cm,
        FV_xmax_margin: 2 * wc.cm,
        FV_ymin_margin: 2.5 * wc.cm,
        FV_ymax_margin: 2.5 * wc.cm,
        FV_zmin_margin: 3 * wc.cm,
        FV_zmax_margin: 3 * wc.cm,
        vertical_dir: [0,1,0],
        beam_dir: [0,0,1]
    },
    a0f0pA: {
        drift_speed: drift_speed,
        tick: 0.5 * wc.us,  // 0.5 mm per tick
        tick_drift: self.drift_speed * self.tick,
        time_offset: time_offset,
        nticks_live_slice: 4,
        FV_xmin: -3579.85 * wc.mm,
        FV_xmax: -25.4 * wc.mm,
        FV_xmin_margin: 2 * wc.cm,
        FV_xmax_margin: 2 * wc.cm,
    },
    a0f1pA: $.a0f0pA + {
        FV_xmin: -3579.85 * wc.mm,
        FV_xmax: -3579.85 * wc.mm,
    },
    a1f0pA: $.a0f0pA + {
        FV_xmin: 3579.85 * wc.mm,
        FV_xmax: 3579.85 * wc.mm,
    },
    a1f1pA: $.a0f0pA + {
        FV_xmin: 25.4 * wc.mm,
        FV_xmax: 3579.85 * wc.mm,
    },
    a2f0pA: $.a0f0pA,
    a2f1pA: $.a0f1pA,
    a3f0pA: $.a1f0pA,
    a3f1pA: $.a1f1pA,
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
    data: { detector_volumes: wc.tn(dv) },
    uses: [dv]
};



local bs_live_face(apa, face) = {
    type: "BlobSampler",
    name: "live-%s-%d"%[apa, face],
    data: {
        drift_speed: drift_speed,
        time_offset: time_offset,
        strategy: ["stepped"],
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

    local bsl = bs_live_face(anode.name, face),
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
            bee_detector: "protodunehd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            dead_area_version: 2,  // v2 wrapper (tpc=apa) so the dead slab lands on the correct PD anode face
            anodes: [wc.tn(anode)],
            face: face,
            detector_volumes: wc.tn(dv),
            bee_points_sets: [  // New configuration for multiple bee points sets
                {
                    name: "clustering",         // Name of the bee points set
                    detector: "protodunehd",         // Detector name
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
        clus_per_face(anode, face=0, dump=false, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo),
        clus_per_face(anode, face=1, dump=false, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo),
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
            bee_detector: "protodunehd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            dead_area_version: 2,  // v2 wrapper (tpc=apa) so the dead slab lands on the correct PD anode face
            anodes: [wc.tn(anode)],
            detector_volumes: wc.tn(dv),
            bee_points_sets: [
                {
                    name: "clustering",
                    detector: "protodunehd",
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

// Stage 3: per drift-side group ({APA0,APA2} drift -x, {APA1,APA3} drift +x),
// raw x,y,z coordinates.  Runs the long-range merge family across the two
// stacked APAs of one drift side, then the topology passes (separate moved
// here from the per-face stage; examine_x_boundary requires all wpids in the
// group to share the same FV_x metadata, true within one drift side).
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
        cm.separate(use_ctpc=true),
        cm.examine_x_boundary(),
        cm.neutrino(),
        cm.isolated(),
        cm.examine_bundles(),
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
            bee_detector: "protodunehd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            dead_area_version: 2,  // v2 wrapper (tpc=apa) so the dead slab lands on the correct PD anode face
            anodes: [wc.tn(a) for a in anodes],
            detector_volumes: wc.tn(dv),
            bee_points_sets: [
                {
                    name: "clustering",
                    detector: "protodunehd",
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
// the per-cluster T0 correction (x_t0cor; with no flash matching on PDHD this
// is the apparent x) and its containment scope filter.
local clus_all_tpc (
    anodes,
    ngroups = 2,
    dump = true,
    bee_dir = "data",
    runNo = 1,
    subRunNo = 1,
    eventNo = 1,
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
        // PDHD central cathode is at x=0, the C++ default cathode_x —
        // dimensions to be confirmed).  use_flash_t0=false because PDHD has
        // no flash matching (the flash-coincidence gate would veto every
        // pair).  Commented out for now; uncomment to enable.
        // cm.cathode_connect(cathode_x_cut=5*wc.cm, drift_cut=8*wc.cm,
        //                    min_length_short=2*wc.cm, short_dir_len=25*wc.cm,
        //                    conn_short_cut=30.0, use_flash_t0=false),
        #cm.retile(cut_time_low=3*wc.us,
        #          cut_time_high=5*wc.us,
        #          anodes=anodes,
        #          samplers=[
        #              clus.sampler(bs_rt_face(0,0), apa=0, face=0),
        #              clus.sampler(bs_rt_face(0,1), apa=0, face=1),
        #              clus.sampler(bs_rt_face(1,0), apa=1, face=0),
        #              clus.sampler(bs_rt_face(1,1), apa=1, face=1),
        #              clus.sampler(bs_rt_face(2,0), apa=2, face=0),
        #              clus.sampler(bs_rt_face(2,1), apa=2, face=1),
        #              clus.sampler(bs_rt_face(3,0), apa=3, face=0),
        #              clus.sampler(bs_rt_face(3,1), apa=3, face=1),
        #          ]),
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
            bee_detector: "protodunehd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: runNo,
            subRunNo: subRunNo,
            eventNo: eventNo,
            save_deadarea: true,
            dead_area_version: 2,  // v2 wrapper (tpc=apa) so the dead slab lands on the correct PD anode face
            dead_apa_groups: apa_drift_groups,  // group dead area by drift side -> 2 dead instances
            anodes: [wc.tn(a) for a in anodes],
            detector_volumes: wc.tn(dv),
            bee_points_sets: [  // New configuration for multiple bee points sets
            //    {
            //        name: "img",                // Name of the bee points set
            //        detector: "protodunehd",         // Detector name
            //        algorithm: "img",           // Algorithm identifier
            //        pcname: "3d",           // Which scope to use
            //        coords: ["x", "y", "z"],    // Coordinates to use
            //        individual: false           // Whether to output as a whole or individual APA/Face
            //    },
            {
                    name: "clustering",         // Name of the bee points set
                    detector: "protodunehd",         // Detector name
                    algorithm: "clustering",    // Algorithm identifier
                    pcname: "3d",           // Which scope to use
                    coords: ["x", "y", "z"],    // Coordinates to use (uncorrected; x_t0cor needs flash-associated t0)
                    individual: false            // Output individual APA/Face
                },
            {
                    // name "img" dumps the live grouping BEFORE the all-TPC
                    // pipeline -> the per-drift-group clustering result (stage
                    // 3), grouped by drift side: clustering-group02 /
                    // clustering-group13.  The "clustering" set above (end
                    // dump) gives clustering-global (full-detector
                    // clustering).
                    name: "img",
                    detector: "protodunehd",
                    algorithm: "clustering",    // -> clustering-group02 / clustering-group13
                    pcname: "3d",
                    coords: ["x", "y", "z"],
                    individual: false,
                    apa_groups: apa_drift_groups,
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
    ret :: g.intern(
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
    per_face(anode, face=0, dump=true) :: clus_per_face(anode, face=face, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo),
    per_apa(anode, dump=true) :: clus_per_apa(anode, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo),
    per_group(anodes, group_name, dump=true) :: clus_per_group(anodes, group_name, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo),
    all_tpc(anodes, ngroups=2, dump=true) :: clus_all_tpc(anodes, ngroups=ngroups, dump=dump, bee_dir=bee_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo),
}
