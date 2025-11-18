// This is a main entry point for configuring a wire-cell CLI job to
// simulate protoDUNE-SP.  It is simplest signal-only simulation with
// one set of nominal field response function.  It excludes noise.
// The kinematics are a mixture of Ar39 "blips" and some ideal,
// straight-line MIP tracks.
//
// Output is a Python numpy .npz file.

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";

local io = import "fileio.jsonnet";
local params = import "pgrapher/experiment/pdsp/simparams.jsonnet";
local tools_maker = import "pgrapher/common/tools.jsonnet";

local tools = tools_maker(params);

local sim_maker = import "pgrapher/experiment/pdsp/sim.jsonnet";
local sim = sim_maker(params, tools);

local stubby = {
    tail: wc.point(1000.0, 3.0, 100.0, wc.mm),
    head: wc.point(1100.0, 3.0, 200.0, wc.mm),
};

// Something close to APA 0 (smallest Y,Z)
local close0 = {
    tail: wc.point(-3.000, 3.0, 1.000, wc.m),
    head: wc.point(-3.100, 3.0, 1.100, wc.m),
};


local tracklist = [
    // {
    //     time: 1*wc.ms,
    //     charge: -5000,         
    //     ray: stubby,
    // },
    {
        time: 0*wc.ms,
        charge: -5000,
        ray: close0,
    },
   // {
   //     time: 0,
   //     charge: -5000,         
   //     ray: params.det.bounds,
   // },
];
local output = "wct-pdsp-sim-ideal-sn.npz";

    
//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
local depos = sim.tracks(tracklist);


local deposio = io.depo_file_source(output);
local drifter = sim.drifter;
local bagger = sim.make_bagger();
local sn_pipes = sim.splusn_pipelines;
local sn_graph = f.fanpipe('DepoSetFanout', sn_pipes, 'FrameFanin', "sn");


local sink = io.frame_file_sink(output);

local graph = g.pipeline([depos, drifter, bagger, sn_graph, sink]);

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


// Finally, the configuration sequence which is emitted.
[cmdline] + g.uses(graph) + [app]

