local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import 'pgrapher/common/funcs.jsonnet';
local clus = import "pgrapher/common/clus.jsonnet";


local time_offset = -250 * wc.us;
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
local common_corr_coords = ["x_t0cor", "y", "z"];


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

// WARNING: wcp-porting-img had the two blob samplers (live and dead) with the
// same name!  These two functions make them distinct.  But in order to
// reproduce the results, we use just bs_face().  Once things are working, the
// code below should delete bs_face() and use one of bs_{live,dead}_face().


// Note, the "sampler" must be unique to the "sampling".
local bs_live_face(apa, face) = {
    type: "BlobSampler",
    name: "live-%s-%d"%[apa, face],
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
        extra: [".*wire_index", "wpid"]
    }
};
local bs_dead_face(apa, face) = {
    type: "BlobSampler",
    name: "dead-%s-%d"%[apa, face],
    data: {
        strategy: [
            "center",
        ],
        extra: [".*"] // want all the extra
    }
};

local bs_face(apa, face) = {
    type: "BlobSampler",
    name: "apa%s-%d"%[apa, face],
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
        extra: [".*wire_index", "wpid"]
    }
};




local clus_per_face (
    anode,
    face,
    dump = true,
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

    // local bsl = bs_live_face(anode.name, face),
    // local bsd = bs_dead_face(anode.name, face),

    // WARNING, likely a bug: these two are the same!
    local bsl = bs_face(anode.data.ident, face),
    local bsd = bs_face(anode.data.ident, face),

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
        // cm.ctpointcloud(),
        cm.live_dead(dead_live_overlap_offset=2),
        cm.extend(flag=4, length_cut=60*wc.cm, num_try=0, length_2_cut=15*wc.cm, num_dead_try=1),
        cm.regular(name="-one", length_cut=60*wc.cm, flag_enable_extend=false),
        cm.regular(name="_two", length_cut=30*wc.cm, flag_enable_extend=true),
        cm.parallel_prolong(length_cut=35*wc.cm),
        cm.close(length_cut=1.2*wc.cm),
        cm.extend_loop(num_try=3),
        cm.separate(use_ctpc=true),
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
            clustering_methods: wc.tns(cm_pipeline),
        }
    }, nin=1, nout=1, uses=[dv, anode, pcts]+cm_pipeline),

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
            bee_zip: "mabc-%s.zip"%[anode.name],
            bee_detector: "sbnd",
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: LrunNo,
            subRunNo: LsubRunNo,
            eventNo: LeventNo,
            save_deadarea: true,
            anodes: [wc.tn(anode)],
            detector_volumes: wc.tn(dv),
            clustering_methods: wc.tns(cm_pipeline),
        }
    }, nin=1, nout=1, uses=[anode, dv, pcts]+cm_pipeline),

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
        // cm_old.examine_x_boundary(),
        cm_old.switch_scope(),

        cm.extend(flag=4, length_cut=60*wc.cm, num_try=0, length_2_cut=15*wc.cm, num_dead_try= 1),
        cm.regular(name="1", length_cut=60*wc.cm, flag_enable_extend=false),
        cm.regular(name="2", length_cut=30*wc.cm, flag_enable_extend=true),
        cm.parallel_prolong(length_cut=35*wc.cm),
        cm.close(length_cut=1.2*wc.cm),
        cm.extend_loop(num_try=3),
        cm.separate(use_ctpc=true),
        cm.neutrino(),
        cm.isolated(),
        cm.examine_bundles(),
        cm.retile(cut_time_low=3*wc.us,
                  cut_time_high=5*wc.us,
                  anodes=anodes, 
                  samplers=[
                      clus.sampler(bs_face(0,0), apa=0, face=0),
                      clus.sampler(bs_face(0,1), apa=0, face=1),
                      clus.sampler(bs_face(1,0), apa=1, face=0),
                      clus.sampler(bs_face(1,1), apa=1, face=1),
                      clus.sampler(bs_face(2,0), apa=2, face=0),
                      clus.sampler(bs_face(2,1), apa=2, face=1),
                      clus.sampler(bs_face(3,0), apa=3, face=0),
                      clus.sampler(bs_face(3,1), apa=3, face=1),
                  ]),
    ],

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
                    coords: ["x_t0cor", "y", "z"],    // Coordinates to use
                    individual: false            // Output individual APA/Face
                }
            ],
            clustering_methods: wc.tns(cm_pipeline),
        },
    }, nin=1, nout=1, uses=anodes+[dv, pcts]+cm_pipeline),

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
