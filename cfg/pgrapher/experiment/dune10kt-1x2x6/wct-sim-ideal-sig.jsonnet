// This is a main entry point for configuring a wire-cell CLI job to
// simulate protoDUNE-SP.  It is simplest signal-only simulation with
// one set of nominal field response function.  It excludes noise.
// The kinematics are a mixture of Ar39 "blips" and some ideal,
// straight-line MIP tracks.
//
// Output is a Python numpy .npz file.

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";

local params_maker = import 'pgrapher/experiment/dune10kt-1x2x6/simparams.jsonnet';
local params = params_maker({});
local tools_maker = import "pgrapher/common/tools.jsonnet";

local tools = tools_maker(params);

local sim_maker = import "pgrapher/experiment/dune10kt-1x2x6/sim.jsonnet";
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
local output = "wct-sim-ideal-sig.tar.bz2";

    
//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
local depos = sim.tracks(tracklist);


local drifter = sim.drifter;
local bagger = sim.make_bagger();
local signal = sim.signal;

// Save the simulated raw frame.  Replaces the removed
// pgrapher/common/fileio.jsonnet npz helper -- follows the inline
// FrameFileSink pattern of the pdhd_sim configs (output: tar.bz2).
local frame_sink = g.pnode({
    type: "FrameFileSink",
    name: "framesink",
    data: {
        outname: output,
        tags: [],          // empty => save all traces
        digitize: false,
        masks: false,
    },
}, nin=1, nout=0);

local graph = g.pipeline([depos, drifter, bagger, signal, frame_sink]);

local app = {
    type: "Pgrapher",
    data: {
        edges: g.edges(graph),
    },
};

// The command-line node.  Replaces the removed pgrapher/ui/cli/nodes.jsonnet
// helper -- inlined, following the pdhd_sim configs.
local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellTbb"],
        apps: ["Pgrapher"],
    },
};

// Finally, the configuration sequence which is emitted.

[cmdline] + g.uses(graph) + [app]

