local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import 'pgrapher/common/funcs.jsonnet';
// local params = import "pgrapher/experiment/dune-vd/params.jsonnet";
// local tools_maker = import 'pgrapher/common/tools.jsonnet';
// local tools = tools_maker(params({}));
// local anodes = tools.anodes;


local time_offset = 0 * wc.us;
local drift_speed = 1.6 * wc.mm / wc.us;
local bee_dir = "data";
local bee_zip = "mabc.zip";

local initial_index = "0";
local initial_runNo = "1";
local initial_subRunNo = "1";
local initial_eventNo = "1";
local index = std.parseInt(initial_index);
local LrunNo = std.parseInt(initial_runNo);
local LsubRunNo = std.parseInt(initial_subRunNo);
local LeventNo  = std.parseInt(initial_eventNo);


local common_coords = ["x", "y", "z"];
    
local geom_helper = {
    type: "SimpleClusGeomHelper",
    name: "uboone",
    data: {
        a0f0: {
            face: 0,
            pitch_u: 3 * wc.mm,
            pitch_v: 3 * wc.mm,
            pitch_w: 3 * wc.mm,
            angle_u: 1.0472,    // 60 degrees
            angle_v: -1.0472,   // -60 degrees
            angle_w: 0,         // 0 degrees
            drift_speed: drift_speed*1,
            tick: 0.5 * wc.us,  // 0.5 mm per tick
            tick_drift: self.drift_speed * self.tick,
            time_offset: time_offset,
            nticks_live_slice: 4,
            FV_xmin: 1 * wc.cm,
            FV_xmax: 255 * wc.cm,
            FV_ymin: -99.5 * wc.cm,
            FV_ymax: 101.5 * wc.cm,
            FV_zmin: 15 * wc.cm,
            FV_zmax: 1022 * wc.cm,
            FV_xmin_margin: 2 * wc.cm,
            FV_xmax_margin: 2 * wc.cm,
            FV_ymin_margin: 2.5 * wc.cm,
            FV_ymax_margin: 2.5 * wc.cm,
            FV_zmin_margin: 3 * wc.cm,
            FV_zmax_margin: 3 * wc.cm
        },
        a1f0: self.a0f0 {},
        a2f0: self.a0f0 {},
        a3f0: self.a0f0 {},
        a4f0: self.a0f0 {},
        a5f0: self.a0f0 {},
        a6f0: self.a0f0 {},
        a7f0: self.a0f0 {},
        a8f0: self.a0f0 {},
        a9f0: self.a0f0 {},
        a10f0: self.a0f0 {},
        a11f0: self.a0f0 {},
        a12f0: self.a0f0 {},
        a13f0: self.a0f0 {},
        a14f0: self.a0f0 {},
        a15f0: self.a0f0 {},
        a16f0: self.a0f0 {},
        a17f0: self.a0f0 {},
        a18f0: self.a0f0 {},
        a19f0: self.a0f0 {},
        a20f0: self.a0f0 {},
        a21f0: self.a0f0 {},
        a22f0: self.a0f0 {},
        a23f0: self.a0f0 {},
        a0f1: self.a0f0 {face: 1},
        a1f1: self.a0f0 {face: 1},
        a2f1: self.a0f0 {face: 1},
        a3f1: self.a0f0 {face: 1},
        a4f1: self.a0f0 {face: 1},
        a5f1: self.a0f0 {face: 1},
        a6f1: self.a0f0 {face: 1},
        a7f1: self.a0f0 {face: 1},
        a8f1: self.a0f0 {face: 1},
        a9f1: self.a0f0 {face: 1},
        a10f1: self.a0f0 {face: 1},
        a11f1: self.a0f0 {face: 1},
        a12f1: self.a0f0 {face: 1},
        a13f1: self.a0f0 {face: 1},
        a14f1: self.a0f0 {face: 1},
        a15f1: self.a0f0 {face: 1},
        a16f1: self.a0f0 {face: 1},
        a17f1: self.a0f0 {face: 1},
        a18f1: self.a0f0 {face: 1},
        a19f1: self.a0f0 {face: 1},
        a20f1: self.a0f0 {face: 1},
        a21f1: self.a0f0 {face: 1},
        a22f1: self.a0f0 {face: 1},
        a23f1: self.a0f0 {face: 1},
    }
};

local clus_per_face (
    anode,
    face,
    dump = true,
    ) =
{

    // Note, the "sampler" must be unique to the "sampling".
    local bs_live = {
        type: "BlobSampler",
        name: "%s-%d"%[anode.name, face],
        data: {
            drift_speed: drift_speed,
            time_offset: time_offset,
            strategy: [
                // "center",
                // "corner",
                // "edge",
                // "bounds",
                "stepped",
                // {name:"grid", step:1, planes:[0,1]},
                // {name:"grid", step:1, planes:[1,2]},
                // {name:"grid", step:1, planes:[2,0]},
                // {name:"grid", step:2, planes:[0,1]},
                // {name:"grid", step:2, planes:[1,2]},
                // {name:"grid", step:2, planes:[2,0]},
            ],
            // extra: [".*"] // want all the extra
            extra: [".*wire_index", "wpid"] //
            // extra: [] //
        }},
    local bs_dead = {
        type: "BlobSampler",
        name: "%s-%d"%[anode.name, face],
        data: {
            strategy: [
                "center",
            ],
            extra: [".*"] // want all the extra
        }},


    local detector_volumes = 
    {
        "type": "DetectorVolumes",
        "name": "dv-%s-%d"%[anode.name, face],
        "data": {
            "anodes": [wc.tn(anode)],
        }
    },

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

    local ptb = g.pnode({
        type: "PointTreeBuilding",
        name: "%s-%d"%[anode.name, face],
        data:  {
            samplers: {
                "3d": wc.tn(bs_live),
                "dead": wc.tn(bs_dead),
            },
            multiplicity: 2,
            tags: ["live", "dead"],
            anode: wc.tn(anode),
            face: face,
            geom_helper: wc.tn(geom_helper),
            detector_volumes: wc.tn(detector_volumes),
        }
    }, nin=2, nout=1, uses=[bs_live, bs_dead, detector_volumes]),

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

    local mabc = g.pnode({
        type: "MultiAlgBlobClustering",
        name: "%s-%d"%[anode.name, face],
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            // grouping2file_prefix: "grouping%s-%d"%[anode.name, face],
            perf: true,
            bee_dir: bee_dir, // "data/0/0", // not used
            bee_zip: "mabc-%s-face%d.zip"%[anode.name, face],
            bee_detector: "sbnd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: LrunNo,
            subRunNo: LsubRunNo,
            eventNo: LeventNo,
            save_deadarea: true, 
            anodes: [wc.tn(anode)],
            face: face,
            geom_helper: wc.tn(geom_helper),
            detector_volumes: wc.tn(detector_volumes),
            func_cfgs: [
                // {name: "clustering_test", detector_volumes: wc.tn(detector_volumes)},
                // {name: "clustering_ctpointcloud", detector_volumes: wc.tn(detector_volumes)},
                // {name: "clustering_switch_scope", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: ["x", "y", "z"], correction_name: "T0Correction"},
                {name: "clustering_live_dead", dead_live_overlap_offset: 2, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_extend", flag: 4, length_cut: 60 * wc.cm, num_try: 0, length_2_cut: 15 * wc.cm, num_dead_try: 1, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_regular", length_cut: 60*wc.cm, flag_enable_extend: false, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_regular", length_cut: 30*wc.cm, flag_enable_extend: true, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_parallel_prolong", length_cut: 35*wc.cm, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_close", length_cut: 1.2*wc.cm, pc_name: "3d", coords: common_coords},
                // {name: "clustering_extend_loop", num_try: 3, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_separate", use_ctpc: true, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_connect1", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_deghost", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_examine_x_boundary", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_protect_overclustering", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_neutrino", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_isolated", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
            ],
        }
    }, nin=1, nout=1, uses=[geom_helper]),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "clus_per_face-%s-%d"%[anode.name, face],
        data: {
            outname: "trash-%s-face%d.tar.gz"%[anode.name, face],
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
        clus_per_face(anode, face=0, dump=false),
        clus_per_face(anode, face=1, dump=false),
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

    local detector_volumes = 
    {
        "type": "DetectorVolumes",
        "name": "dv-%s"%[anode.name],
        "data": {
            "anodes": [wc.tn(anode)],
        }
    },

    local mabc = g.pnode({
        type: "MultiAlgBlobClustering",
        name: "clus_per_apa-%s"%[anode.name],
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            // grouping2file_prefix: "grouping%s-%d"%[anode.name, face],
            perf: true,
            bee_dir: bee_dir, // "data/0/0", // not used
            bee_zip: "mabc-%s.zip"%[anode.name],
            bee_detector: "sbnd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: LrunNo,
            subRunNo: LsubRunNo,
            eventNo: LeventNo,
            save_deadarea: true,
            anodes: [wc.tn(anode)],
            // face: face,
            geom_helper: wc.tn(geom_helper),
            detector_volumes: wc.tn(detector_volumes),
            func_cfgs: [
                // {name: "clustering_test", detector_volumes: wc.tn(detector_volumes)},
                // {name: "clustering_ctpointcloud", detector_volumes: wc.tn(detector_volumes)},
                // {name: "clustering_switch_scope", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: ["x", "y", "z"], correction_name: "T0Correction"},
                // {name: "clustering_live_dead", dead_live_overlap_offset: 2, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_extend", flag: 4, length_cut: 60 * wc.cm, num_try: 0, length_2_cut: 15 * wc.cm, num_dead_try: 1, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_regular", length_cut: 60*wc.cm, flag_enable_extend: false, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_regular", length_cut: 30*wc.cm, flag_enable_extend: true, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_parallel_prolong", length_cut: 35*wc.cm, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_close", length_cut: 1.2*wc.cm, pc_name: "3d", coords: common_coords},
                // {name: "clustering_extend_loop", num_try: 3, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_separate", use_ctpc: true, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_connect1", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_deghost", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_examine_x_boundary", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_protect_overclustering", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_neutrino", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_isolated", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
            ],
        }
    }, nin=1, nout=1, uses=[geom_helper, detector_volumes]),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "clus_per_apa-%s"%[anode.name],
        data: {
            outname: "trash-%s.tar.gz"%[anode.name],
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

local clus_all_apa (
    anodes,
    dump = true,
    ) = {
    local nanodes = std.length(anodes),
    local pcmerging = g.pnode({
        type: "PointTreeMerging",
        name: "clus_all_apa",
        data:  {
            multiplicity: nanodes,
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
        }
    }, nin=nanodes, nout=1),
    local detector_volumes = 
    {
        "type": "DetectorVolumes",
        "name": "clus_all_apa",
        "data": {
            "anodes": [wc.tn(anode) for anode in anodes],
        }
    },
    local mabc = g.pnode({
        type: "MultiAlgBlobClustering",
        name: "clus_all_apa",
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            // grouping2file_prefix: "grouping%s-%d"%[anode.name, face],
            perf: true,
            bee_dir: bee_dir, // "data/0/0", // not used
            bee_zip: "mabc-all-apa.zip",
            bee_detector: "sbnd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: LrunNo,
            subRunNo: LsubRunNo,
            eventNo: LeventNo,
            save_deadarea: true, 
            anodes: [wc.tn(a) for a in anodes],
            // face: face,
            geom_helper: wc.tn(geom_helper),
            detector_volumes: wc.tn(detector_volumes),
            bee_points_sets: [  // New configuration for multiple bee points sets
                {
                    name: "img",                // Name of the bee points set
                    detector: "uboone",         // Detector name
                    algorithm: "img",           // Algorithm identifier
                    pcname: "3d",           // Which scope to use
                    coords: ["x", "y", "z"],    // Coordinates to use
                    individual: false           // Whether to output as a whole or individual APA/Face
                },
                {
                    name: "clustering",         // Name of the bee points set
                    detector: "uboone",         // Detector name
                    algorithm: "clustering",    // Algorithm identifier
                    pcname: "3d",           // Which scope to use
                    coords: ["x_t0cor", "y", "z"],    // Coordinates to use
                    individual: true            // Output individual APA/Face
                }
            ],
            func_cfgs: [
                // {name: "clustering_test", detector_volumes: wc.tn(detector_volumes)},
                // {name: "clustering_ctpointcloud", detector_volumes: wc.tn(detector_volumes)},
                // {name: "clustering_switch_scope", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: ["x", "y", "z"], correction_name: "T0Correction"},
                // {name: "clustering_live_dead", dead_live_overlap_offset: 2, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_extend", flag: 4, length_cut: 60 * wc.cm, num_try: 0, length_2_cut: 15 * wc.cm, num_dead_try: 1, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_regular", length_cut: 60*wc.cm, flag_enable_extend: false, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_regular", length_cut: 30*wc.cm, flag_enable_extend: true, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_parallel_prolong", length_cut: 35*wc.cm, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_close", length_cut: 1.2*wc.cm, pc_name: "3d", coords: common_coords},
                // {name: "clustering_extend_loop", num_try: 3, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_separate", use_ctpc: true, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_connect1", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_deghost", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_examine_x_boundary", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_protect_overclustering", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_neutrino", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
                // {name: "clustering_isolated", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
            ],
        },
    }, nin=1, nout=1, uses=[geom_helper, detector_volumes]),

    local sink = g.pnode({
        type: "TensorFileSink",
        name: "clus_all_apa",
        data: {
            outname: "trash-all-apa.tar.gz",
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


function () {
    per_face(anode, face=0, dump=true) :: clus_per_face(anode, face=face, dump=dump),
    per_apa(anode, dump=true) :: clus_per_apa(anode, dump=dump),
    all_apa(anodes, dump=true) :: clus_all_apa(anodes, dump=dump),
}