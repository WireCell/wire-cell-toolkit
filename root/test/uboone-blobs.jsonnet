// Clustering testing with Uboone and blobs originating from WCP.
//
// This provides config for various "kinds" of subgraphs to be run separately
// and be connected together via files.:
// 
// - live :: load WCP "live" blobs from root file, produce npz
// - live :: load WCP "dead" blobs from root file, produce npz
// - clus :: read back live and dead npz files produce bee and tensor files
//
// Best to run full graph via uboone-blobs.smake

local high = import "layers/high.jsonnet";
local wc = high.wc;
local pg = high.pg;
local detector = "uboone";
local params = high.params(detector);
local mid = high.api(detector, params);
local anodes = mid.anodes();
local anode = anodes[0];

// live/dead symmetries
local UbooneBlobSource(fname, kind /*live or dead*/, views /* uvw, uv, vw, wu */) = pg.pnode({
    type: 'UbooneBlobSource',
    name: kind+'-'+views,
    data: {
        input: fname,
        anode: wc.tn(anode),
        kind: kind,
        views: views,
    }
}, nin=0, nout=1, uses=[anode]);
local multi_source = function(iname, kind, views)
    local nviews = std.length(views);
    local srcs = [ UbooneBlobSource(iname, kind, view), for view in views ];
    local bsm = pg.pnode({
        type: "BlobSetMerge",
        name: kind,
        data: { multiplicity: nviews, },
    }, nin=4, nout=1);
    pg.intern(innodes = srcs, outnodes=[bsm],
              edges = [ pg.edge(srcs[ind], bsm, 0, ind),
                        for ind in std.range(0, nviews-1) ]);

local BlobClustering(name) = pg.pnode({
    type: 'BlobClustering',
    name: name,
    data: {
        policy: "uboone",
    },
}, nin=1, nout=1);
local ClusterFileSource(fname) = pg.pnode({
    type: 'ClusterFileSource',
    name: fname,
    data: {
        inname: fname,
        anodes: [wc.tn(a) for a in anodes],
    },
}, nin=0, nout=1, uses=anodes);
local ClusterFileSink(fname) = pg.pnode({
    type: 'ClusterFileSink',
    name: fname,
    data: {
        format: "numpy",
        outname: fname,
    },
}, nin=1, nout=0);


// generators of the live pipeline elements
local ProjectionDeghosting(name) = pg.pnode({
    type: 'ProjectionDeghosting',
    name: name,
    data: {},
}, nin=1, nout=1);
local InSliceDeghosting(name, round /*1,2,3*/) = pg.pnode({
    type: "InSliceDeghosting",
    name: name,
    data:  {
        config_round: round,
    }
}, nin=1, nout=1);
local BlobGrouping(name) = pg.pnode({
    type: "BlobGrouping",
    name: name,
    data:  { }
}, nin=1, nout=1);
local ChargeSolving(name, weighting /* uniform, uboone */) = pg.pnode({
    type: "ChargeSolving",
    name: name,
    data:  {
        weighting_strategies: [weighting],
    }
}, nin=1, nout=1);
local LocalGeomClustering(name) = pg.pnode({
    type: "LocalGeomClustering",
    name: name,
    data:  { },
}, nin=1, nout=1);
local GlobalGeomClustering(name, policy="uboone") = pg.pnode({
    type: "GlobalGeomClustering",
    name: name,
    data:  {
        clustering_policy: policy,
    },
}, nin=1, nout=1);

local bs_live = {
    type: "BlobSampler",
    name: "live",
    data: {
        time_offset: -1600 * wc.us,
        drift_speed: 1.101 * wc.mm / wc.us,
        strategy: [
            "stepped",
        ],
        extra: [".*wire_index"] //
    }};
local bs_dead = {
    type: "BlobSampler",
    name: "dead",
    data: {
        strategy: [
            "center",
        ],
        extra: [".*"] // want all the extra
    }};
local PointTreeBuilding() = pg.pnode({
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
    }
}, nin=2, nout=1, uses=[bs_live, bs_dead]);
local point_tree_source = function(livefn, deadfn)
    local livesrc = ClusterFileSource(livefn);
    local deadsrc = ClusterFileSource(deadfn);
    local ptb = PointTreeBuilding();
    pg.intern(innodes=[livesrc, deadsrc], outnodes=[ptb],
              edges=[ pg.edge(livesrc, ptb, 0, 0),
                      pg.edge(deadsrc, ptb, 0, 1) ]
             );

    
local BeeBlobSink(fname, sampler) = pg.pnode({
    type: "BeeBlobSink",
    name: fname,
    data: {
        geom: "uboone",
        type: "wcp",
        outname: fname,
        samplers: wc.tn(sampler)
    },
}, nin=1, nout=0, uses=[sampler]);
local BeeBlobTap = function(fname)
    local sink = BeeBlobSink(fname);
    local fan = pg.pnode({
        type:'BlobSetFanout',
        name:fname,
        data: { multiplicity: 2 },
    }, nin=1, nout=2);
    pg.intern(innodes=[fan], centernodes=[sink],
              edges=[ pg.edge(fan, sink, 1, 0) ]);


local MultiAlgBlobClustering(beezip) = pg.pnode({
    type: "MultiAlgBlobClustering",
        name: beezip,
        data:  {
            inpath: "pointtrees/%d",
            outpath: "pointtrees/%d",
            perf: true,
            bee_zip: beezip,
            save_deadarea: true, 
            dead_live_overlap_offset: 2,
            anode: wc.tn(anode),
        }
    }, nin=1, nout=1, uses=[anode]);

local TensorFileSink(fname) = pg.pnode({
    type: "TensorFileSink",
        name: fname,
        data: {
            outname: fname,
            prefix: "clustering_",
            dump_mode: true,
        }
    }, nin=1, nout=0);

local default_files = {
    root: "result_5384_130_6501.root",
    live: "live-clus.npz",
    dead: "dead-clus.npz",
    tens: "live-dead.npz",
    beez: "live-dead.bee"
};

local live(iname, oname) = pg.pipeline([
    multi_source(iname, "live", ["uvw","uv","vw","wu"]),
    BlobClustering("live"),
    // BlobGrouping("0"),

    // "standard":
    // ProjectionDeghosting("1"),
    // BlobGrouping("1"), ChargeSolving("1a","uniform"), LocalGeomClustering("1"), ChargeSolving("1b","uboone"),
    // InSliceDeghosting("1",1),
    // ProjectionDeghosting("2"),
    // BlobGrouping("2"), ChargeSolving("2a","uniform"), LocalGeomClustering("2"), ChargeSolving("2b","uboone"),
    // InSliceDeghosting("2",2),
    // BlobGrouping("3"), ChargeSolving("3a","uniform"), LocalGeomClustering("3"), ChargeSolving("3b","uboone"),
    // InSliceDeghosting("3",3),
    GlobalGeomClustering(""),
    ClusterFileSink(oname),
]);

local dead(iname, oname) = pg.pipeline([
    multi_source(iname, "dead", ["uv","vw","wu"]),
    BlobClustering("dead"),
    GlobalGeomClustering("", "dead_clus"),
    ClusterFileSink(oname),
]);

local clus = function(iname, oname)
    local live_dead = std.split(iname, ",");
    local beezip_tensor = std.split(oname, ",");
    pg.pipeline([
        point_tree_source(live_dead[0], live_dead[1]),
        MultiAlgBlobClustering(beezip_tensor[0]),
        TensorFileSink(beezip_tensor[1])]);
    
local extra_plugins = ["WireCellRoot","WireCellClus"];

function(iname, oname, kind /*clus or live or dead*/)
    if kind == "live" then
        high.main(live(iname, oname), "Pgrapher", extra_plugins)
    else if kind=="dead" then
        high.main(dead(iname, oname), "Pgrapher", extra_plugins)
    else if kind=="clus" then
        high.main(clus(iname, oname), "Pgrapher", extra_plugins)
    else
        error "uknown kind: "+kind

