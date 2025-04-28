local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import 'pgrapher/common/funcs.jsonnet';
local params = import "pgrapher/experiment/uboone/simparams.jsonnet";
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);
local anodes = tools.anodes;

local cluster_source(fname) = g.pnode({
    type: "ClusterFileSource",
    name: fname,
    data: {
        inname: fname,
        anodes: [wc.tn(a) for a in anodes],
    }
}, nin=0, nout=1, uses=anodes);

function (
    active_clusters = "active-clusters-anode0.npz",
    masked_clusters = "masked-clusters-anode0.npz",
    output = "tensor-apa-uboone.tar.gz",
    bee_dir = "data",
    bee_zip = "mabc.zip",
    initial_index = "0",
    initial_runNo = "1",
    initial_subRunNo = "1",
    initial_eventNo = "1")

    local index = std.parseInt(initial_index);

    local LrunNo = std.parseInt(initial_runNo);
    local LsubRunNo = std.parseInt(initial_subRunNo);
    local LeventNo  = std.parseInt(initial_eventNo);

    local active = cluster_source(active_clusters);
    local masked = cluster_source(masked_clusters);

    // Note, the "sampler" must be unique to the "sampling".
    local bs_live = {
        type: "BlobSampler",
        name: "bs_live",
        data: {
            drift_speed: 1.101 * wc.mm / wc.us,
            time_offset: -1600 * wc.us + 6 * wc.mm/self.drift_speed,
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
        }};
    local bs_dead = {
        type: "BlobSampler",
        name: "bs_dead",
        data: {
            strategy: [
                "center",
            ],
            extra: [".*"] // want all the extra
        }};

    local detector_volumes = 
    {
        type: "DetectorVolumes",
        name: "",
        data: {
            anodes: [wc.tn(a) for a in anodes],
            metadata:
                {overall: {
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
                    FV_zmax_margin: 3 * wc.cm,
                    vertical_dir: [0,1,0],
                    beam_dir: [0,0,1]
                }} +
                {
                    [ "a" + std.toString(a.data.ident) + "f0pA" ]: {
                        drift_speed: 1.101 * wc.mm / wc.us,
                        tick: 0.5 * wc.us,  // 0.5 mm per tick
                        tick_drift: self.drift_speed * self.tick,
                        time_offset: -1600 * wc.us + 6 * wc.mm/self.drift_speed,
                        nticks_live_slice: 4,
                        FV_xmin: 1 * wc.cm,
                        FV_xmax: 255 * wc.cm,
                        FV_xmin_margin: 2 * wc.cm,
                        FV_xmax_margin: 2 * wc.cm,
                    } for a in anodes
                }
        },
    };

    local pctransforms = {
        type: "PCTransformSet",
        name: "",
        data: { detector_volumes: wc.tn(detector_volumes) },
        uses: [detector_volumes]
    };


    local ptb = g.pnode({
        type: "PointTreeBuilding",
        name: "",
        data:  {
            samplers: {
                "3d": wc.tn(bs_live),
                "dead": wc.tn(bs_dead),
            },
            multiplicity: 2,
            tags: ["live", "dead"],
            anode: wc.tn(anodes[0]),
            face: 0,
            detector_volumes: wc.tn(detector_volumes),
        }
    }, nin=2, nout=1, uses=[bs_live, bs_dead, detector_volumes]);

    local front_end = g.intern(
        innodes = [active, masked],
        outnodes = [ptb],
        edges = [
            g.edge(active, ptb, 0, 0),
            g.edge(masked, ptb, 0, 1),
        ],
        name = "front-end");

    local common_coords = ["x_t0cor", "y", "z"];
    // local common_coords = ["x", "y", "z"];
    local mabc = g.pnode({
        type: "MultiAlgBlobClustering",
        name: "",
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            perf: true,
            bee_dir: bee_dir, // "data/0/0",
            bee_zip: bee_zip,
            initial_index: index,   // New RSE configuration
            use_config_rse: true,  // Enable use of configured RSE
            runNo: LrunNo,
            subRunNo: LsubRunNo,
            eventNo: LeventNo,
            save_deadarea: true, 
            anodes: [wc.tn(a) for a in anodes],
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
               // {name: "clustering_test", detector_volumes: wc.tn(detector_volumes), pc_transforms: wc.tn(pctransforms)},
               // {name: "clustering_ctpointcloud", detector_volumes: wc.tn(detector_volumes), pc_transforms: wc.tn(pctransforms)},
               {name: "clustering_switch_scope", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: ["x", "y", "z"], correction_name: "T0Correction", pc_transforms: wc.tn(pctransforms)},
               {name: "clustering_live_dead", dead_live_overlap_offset: 2, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
               {name: "clustering_extend", flag: 4, length_cut: 60 * wc.cm, num_try: 0, length_2_cut: 15 * wc.cm, num_dead_try: 1, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
               {name: "clustering_regular", length_cut: 60*wc.cm, flag_enable_extend: false, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
               {name: "clustering_regular", length_cut: 30*wc.cm, flag_enable_extend: true, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
               {name: "clustering_parallel_prolong", length_cut: 35*wc.cm, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
               {name: "clustering_close", length_cut: 1.2*wc.cm, pc_name: "3d", coords: common_coords},
               {name: "clustering_extend_loop", num_try: 3, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
               {name: "clustering_separate", use_ctpc: true, detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords, pc_transforms: wc.tn(pctransforms)},
               {name: "clustering_connect1", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
               {name: "clustering_deghost", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords, pc_transforms: wc.tn(pctransforms)},
               {name: "clustering_examine_x_boundary", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
               {name: "clustering_protect_overclustering", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords, pc_transforms: wc.tn(pctransforms)},
               {name: "clustering_neutrino", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
               {name: "clustering_isolated", detector_volumes: wc.tn(detector_volumes), pc_name: "3d", coords: common_coords},
            ],
        }
    }, nin=1, nout=1, uses=[detector_volumes, pctransforms]);

    local sink = g.pnode({
        type: "TensorFileSink",
        name: output,
        data: {
            outname: output,
            prefix: "clustering_", // json, numpy, dummy
            dump_mode: true,
        }
    }, nin=1, nout=0);

    local graph = g.pipeline([front_end, mabc, sink]);

    local app = {
        type: 'Pgrapher', //Pgrapher, TbbFlow
        data: {
            edges: g.edges(graph),
        },
    };
    local cmdline = {
        type: "wire-cell",
        data: {
            plugins: ["WireCellGen", "WireCellPgraph", /*"WireCellTbb",*/
                      "WireCellSio", "WireCellSigProc", "WireCellRoot", "WireCellImg", "WireCellClus"],
        apps: [wc.tn(app)]
        },
    };

    [cmdline] + g.uses(graph) + [app]
