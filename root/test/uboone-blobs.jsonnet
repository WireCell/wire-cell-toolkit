// This loads the Uboone ROOT file with Trun, TC and TDC TTrees to produce
// "live" and "dead" blob sets with two UbooneBlobSource nodes.  It then runs a
// version of wct-uboone-img.jsonnet which can the be followed with
// wct-uboone-clustering from https://github.com/HaiwangYu/wcp-porting-img/.

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
local BlobClustering(name) = pg.pnode({
    type: 'BlobClustering',
    name: name,
    data: {
        policy: "uboone",
    },
}, nin=1, nout=1);
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
local GlobalGeomClustering(name) = pg.pnode({
    type: "GlobalGeomClustering",
    name: name,
    data:  { },
}, nin=1, nout=1);

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
    

local live_sampler = {
    type: "BlobSampler",
    name: "live",
    data: {
        time_offset: -1600 * wc.us,
        drift_speed: 1.101 * wc.mm / wc.us,
        strategy: [
            "center",
            "stepped",
        ],
    }};
local dead_sampler = {
    type: "BlobSampler",
    name: "dead",
    data: {
        strategy: [
            "center",
        ],
    }};
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


local live(iname, oname) = pg.pipeline([
    multi_source(iname, "live", ["uvw","uv","vw","wu"]),
    BeeBlobSink(oname, live_sampler),

    // BeeBlobTap("live.zip"),

    // BlobClustering("live"),
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
    // GlobalGeomClustering(""),
    // ClusterFileSink(oname),
]);


local dead(iname, oname) = pg.pipeline([
    multi_source(iname, "dead", ["uv","vw","wu"]),
    BlobClustering("dead"), ClusterFileSink(oname),
]);

local extra_plugins = ["WireCellRoot","WireCellClus"];

function(iname, oname, kind /*live or dead*/)
    if kind == "live"
    then high.main(live(iname, oname), "Pgrapher", extra_plugins)
    else high.main(dead(iname, oname), "Pgrapher", extra_plugins)

