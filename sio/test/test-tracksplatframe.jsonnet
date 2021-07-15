// Test a TrackDepos -> DepoSplat -> FrameSaver graph using
// schema-based configuration.

local g = import "pgraph.jsonnet";
local wc = import "wirecell.jsonnet";

local td = import "schema/gen/TrackDepos.jsonnet";
local h = wc.mergeObjects([
    import "schema/util/Base.jsonnet",
    import "schema/gen/TrackDepos.jsonnet",
    import "schema/gen/Ductor.jsonnet",
    import "schema/sio/NumpyFrameSaver.jsonnet"
]);

local stubby = h.util.Base.Ray(tail=wc.point(1000.0, 3.0, 100.0, wc.mm),
                               head=wc.point(1100.0, 3.0, 200.0, wc.mm));
local tracklist = [
    h.gen.TrackDepos.Track(time=1*wc.ms, charge=-5000, ray=stubby),
];

local depos = g.pnode({
    type: 'TrackDepos',
    data: h.gen.TrackDepos.Config(tracks=tracklist),
}, nin=0, nout=1);

local ductor = g.pnode({
    type: 'DepoSplat',
    name: "deposplat",
    data: h.gen.Ductor.Config(),
}, nin=1, nout=1);

local saver = g.pnode({
    type: 'NumpyFrameSaver',
    name: "saver",
    data: h.sio.NumpyFrameSaver.Config(),
}, nin=1, nout=1);

local sink = g.pnode({
    type: "DumpFrames",
    name: "sink",
    data: {},
}, nin=1, nout=0);


local graph = g.pipeline([depos, ductor, saver, sink]);


local app = {
    type: "Pgrapher",
    data: {
        edges: g.edges(graph),
    },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc"],
        apps: ["Pgrapher"]
    }
};

[cmdline] + g.uses(graph) + [app]

//depos
//ductor


