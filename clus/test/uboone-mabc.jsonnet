// This job inputs Uboone ROOT files (TC, etc), runs MABC, dumbs to bee.
//
// Use like:
//
// wire-cell -l stderr -L debug \
//   -A infiles=nuselEval_5384_137_6852.root \
//      clus/test/uboone-mabc.jsonnet

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

//// Experimental new style 
// local detector = "uboone";
// local params = high.params(detector);
// local mid = high.api(detector, params);
// local anode = mid.anodes()[0];
//// Old style:
local params = import "pgrapher/experiment/uboone/simparams.jsonnet";
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);
local anode = tools.anodes[0];



// The TDM datapath to find point trees.  This needs coordination between a few
// nodes.  This provides the default but the construction functions allow
// override.
local pointtree_datapath = "pointtrees/%d/uboone";

// This object holds a bunch of functions that construct parts of the graph.
// We use this object below to build the full graph.
local ub = {
    anode: anode,

    bs_live : {
        type: "BlobSampler",
        name: "live",
        data: {
            time_offset: -1600 * wc.us,
            drift_speed: 1.101 * wc.mm / wc.us,
            strategy: [
                "stepped",
            ],
            extra: [".*wire_index"] //
        }
    },

    bs_dead : {
        type: "BlobSampler",
        name: "dead",
        data: {
            strategy: [
                "center",
            ],
            extra: [".*"] // want all the extra
        }
    },
    
    UbooneBlobSource(fname, kind /*live or dead*/, views /* uvw, uv, vw, wu */) :: pg.pnode({
        type: 'UbooneBlobSource',
        name: kind+'-'+views,
        data: {
            input: fname,
            anode: wc.tn(anode),
            kind: kind,
            views: views,
        }
    }, nin=0, nout=1, uses=[anode]),

    UbooneClusterSource(fname, sampler=$.bs_live, datapath=pointtree_datapath, optical=true) :: pg.pnode({
        type: 'UbooneClusterSource',
        name: sampler.name,
        data: {
            input: fname,       // file name or list
            datapath: datapath + '/live', // see issue #375
            sampler: wc.tn(sampler),
        } + if optical then {
            light: "light", flash: "flash", flashlight: "flashlight"
        } else {}
    }, nin=1, nout=1, uses=[sampler]),
        
    ClusterFlashDump(datapath=pointtree_datapath) :: pg.pnode({
        type: 'ClusterFlashDump',
        name: "",
        data: {
            datapath: datapath + '/live', // see issue #375
        },
    }, nin=1, nout=0),

    BlobSetMerge(kind, multiplicity) :: pg.pnode({
        type: "BlobSetMerge",
        name: kind,
        data: { multiplicity: multiplicity, },
    }, nin=multiplicity, nout=1),

    // Make one UbooneBlobSource for view and multiplex their output
    multiplex_blob_views(iname, kind, views) ::
        local nviews = std.length(views);
        local srcs = [ $.UbooneBlobSource(iname, kind, view), for view in views ];
        local bsm = $.BlobSetMerge(kind, nviews);
        pg.intern(innodes = srcs, outnodes=[bsm],
                  edges = [ pg.edge(srcs[ind], bsm, 0, ind),
                            for ind in std.range(0, nviews-1) ]),

    BlobClustering(name) :: pg.pnode({
        type: 'BlobClustering',
        name: name,
        data: {
            policy: "uboone",
        },
    }, nin=1, nout=1),

    ClusterFileSource(fname) :: pg.pnode({
        type: 'ClusterFileSource',
        name: fname,
        data: {
            inname: fname,
            anodes: [wc.tn(anode)],
        },
    }, nin=0, nout=1, uses=[anode]),

    ClusterFileSink(fname) :: pg.pnode({
        type: 'ClusterFileSink',
        name: fname,
        data: {
            format: "numpy",
            outname: fname,
        },
    }, nin=1, nout=0),

// generators of the live pipeline elements
    ProjectionDeghosting(name) :: pg.pnode({
        type: 'ProjectionDeghosting',
        name: name,
        data: {},
    }, nin=1, nout=1),

    InSliceDeghosting(name, round /*1,2,3*/) :: pg.pnode({
        type: "InSliceDeghosting",
        name: name,
        data:  {
            config_round: round,
        }
    }, nin=1, nout=1),

    BlobGrouping(name) :: pg.pnode({
        type: "BlobGrouping",
        name: name,
        data:  { }
    }, nin=1, nout=1),

    ChargeSolving(name, weighting /* uniform, uboone */) :: pg.pnode({
        type: "ChargeSolving",
        name: name,
        data:  {
            weighting_strategies: [weighting],
        }
    }, nin=1, nout=1),

    LocalGeomClustering(name) :: pg.pnode({
        type: "LocalGeomClustering",
        name: name,
        data:  { },
    }, nin=1, nout=1),
        

    GlobalGeomClustering(name, policy="uboone") :: pg.pnode({
        type: "GlobalGeomClustering",
        name: name,
        data:  {
            clustering_policy: policy,
        },
    }, nin=1, nout=1),


    SimpleClusGeomHelper() :: {
        type: "SimpleClusGeomHelper",
        name: "uboone",
        data: {
            a0f0: {
                pitch_u: 3 * wc.mm,
                pitch_v: 3 * wc.mm,
                pitch_w: 3 * wc.mm,
                angle_u: 1.0472,    // 60 degrees
                angle_v: -1.0472,   // -60 degrees
                angle_w: 0,         // 0 degrees
                drift_speed: 1.101 * wc.mm / wc.us,
                tick: 0.5 * wc.us,  // 0.5 mm per tick
                tick_drift: self.drift_speed * self.tick,
                time_offset: -1600 * wc.us + 6 * wc.mm/self.drift_speed,
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
        }
    },

    PointTreeBuilding(geom_helper = $.SimpleClusGeomHelper()) :: pg.pnode({
        type: "PointTreeBuilding",
        name: "",
        data:  {
            samplers: {
                "3d": wc.tn($.bs_live),
                "dead": wc.tn($.bs_dead),
            },
            multiplicity: 2,
            tags: ["live", "dead"],
            anode: wc.tn(anode),
            face: 0,
            geom_helper: wc.tn(geom_helper),
        }
    }, nin=2, nout=1, uses=[$.bs_live, $.bs_dead, geom_helper]),

    point_tree_source(livefn, deadfn) ::
        local livesrc = $.ClusterFileSource(livefn);
        local deadsrc = $.ClusterFileSource(deadfn);
        local ptb = $.PointTreeBuilding();
        pg.intern(innodes=[livesrc, deadsrc], outnodes=[ptb],
                  edges=[ pg.edge(livesrc, ptb, 0, 0),
                          pg.edge(deadsrc, ptb, 0, 1) ]
                 ),

    BeeBlobSink(fname, sampler) :: pg.pnode({
        type: "BeeBlobSink",
        name: fname,
        data: {
            geom: "uboone",
            type: "wcp",
            outname: fname,
            samplers: wc.tn(sampler)
        },
    }, nin=1, nout=0, uses=[sampler]),

    BeeBlobTap(fname) ::
        local sink = $.BeeBlobSink(fname);
        local fan = pg.pnode({
            type:'BlobSetFanout',
            name:fname,
            data: { multiplicity: 2 },
        }, nin=1, nout=2);
        pg.intern(innodes=[fan], centernodes=[sink],
                  edges=[ pg.edge(fan, sink, 1, 0) ]),

    MultiAlgBlobClustering(beezip, datapath=pointtree_datapath, geom_helper = $.SimpleClusGeomHelper()) :: pg.pnode({
        type: "MultiAlgBlobClustering",
        name: "",
        data:  {
            inpath: pointtree_datapath,
            outpath: pointtree_datapath,
            perf: true,
            bee_zip: beezip,
            initial_index: 0,
            use_config_rse: true,  // Enable use of configured RSE
            runNo: 1,
            subRunNo: 1,
            eventNo: 1,
            save_deadarea: true, 
            anode: wc.tn(anode),
            face: 0,
            geom_helper: wc.tn(geom_helper),
            func_cfgs: [
                {name: "clustering_ctpointcloud"},
                // {name: "clustering_live_dead", dead_live_overlap_offset: 2},
                // {name: "clustering_extend", flag: 4, length_cut: 60 * wc.cm, num_try: 0, length_2_cut: 15 * wc.cm, num_dead_try: 1},
                // {name: "clustering_regular", length_cut: 60*wc.cm, flag_enable_extend: false},
                // {name: "clustering_regular", length_cut: 30*wc.cm, flag_enable_extend: true},
                // {name: "clustering_parallel_prolong", length_cut: 35*wc.cm},
                // {name: "clustering_close", length_cut: 1.2*wc.cm},
                // {name: "clustering_extend_loop", num_try: 3},
                // {name: "clustering_separate", use_ctpc: true},
                // {name: "clustering_connect1"},
                // {name: "clustering_deghost"},
                // {name: "clustering_examine_x_boundary"},
                // {name: "clustering_protect_overclustering"},
                // {name: "clustering_neutrino"},
                // {name: "clustering_isolated"},
            ],
        }
    }, nin=1, nout=1, uses=[geom_helper]),

    TensorFileSink(fname) :: pg.pnode({
        type: "TensorFileSink",
        name: fname,
        data: {
            outname: fname,
            prefix: "clustering_",
            dump_mode: true,
        }
    }, nin=1, nout=0),

    main(graph, app='Pgrapher', extra_plugins = []) ::
        local uses = pg.uses(graph);
        local plugins = [
            "WireCellSio", "WireCellAux",
            "WireCellGen", "WireCellSigProc", "WireCellImg", "WireCellClus",
            "WireCellRoot",
            "WireCellApps"] + {
                'TbbFlow': ["WireCellTbb"],
                'Pgrapher': ["WireCellPgraph"],
            }[app];

        local appcfg = {
            type: app,
            data: {
                edges: pg.edges(graph)
            },
        };
        local cmdline = {
            type: "wire-cell",
            data: {
                plugins: plugins,
                apps: [appcfg.type]
            }
        };
        [cmdline] + pg.uses(graph) + [appcfg],
};

    


local graph(infiles, beezip, datapath=pointtree_datapath) = pg.pipeline([
    ub.multiplex_blob_views(infiles, 'live', ["uvw","uv","vw","wu"]),
    ub.UbooneClusterSource(infiles, datapath=datapath),
    ub.MultiAlgBlobClustering(beezip, datapath=datapath),
    ub.ClusterFlashDump(datapath=datapath)
]);

local extra_plugins = ["WireCellAux", "WireCellRoot", "WireCellClus"];

function(infiles="uboone.root", beezip="bee.zip")
    local g = graph(wc.listify(infiles), beezip);
    ub.main(g, "Pgrapher", extra_plugins)
    
