// This is a main entry point to configure WCT via wire-cell CLI to
// run simulation, (an essentially empty) noise filtering and signal
// processing.
// 
// Simulation is signal and noise with Ar39 and some ideal line MIP
// tracks as initial kinematics.

// Output is a Python numpy .npz file.

local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";

local cli = import "pgrapher/ui/cli/nodes.jsonnet";

local io = import "pgrapher/common/fileio.jsonnet";
local params_maker = import 'pgrapher/experiment/dune10kt-1x2x2/simparams.jsonnet';
local params = params_maker({});
local tools_maker = import "pgrapher/common/tools.jsonnet";
local sim_maker = import "pgrapher/experiment/dune10kt-1x2x2/sim.jsonnet";
// Fixme: currently, no noise filter.  Need to at least add a "null" NF to produce thresholds.
// Or, maybe better, move that into OSP.  W/out it, behavior is undefined.
// local nf = ...
local sp_maker = import "pgrapher/experiment/dune10kt-1x2x2/sp.jsonnet";

local tools = tools_maker(params);

local sim = sim_maker(params, tools);

local sp = sp_maker(params, tools);


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

local output = "wct-dune10kt-1x2x2-sim-ideal-sn-nf-sp.npz";
    
//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
local depos = sim.tracks(tracklist);


local deposio = io.numpy.depos(output);
local drifter = sim.drifter;
local bagger = sim.make_bagger();

// signal plus noise pipelines
local sn_pipes = sim.splusn_pipelines;
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];
local sn_sp = [g.pipeline([sn_pipes[n], sp_pipes[n]], "sn_sp_pipe_%d" % n)
               for n in std.range(0, std.length(tools.anodes)-1)];
local snsp_graph = f.fanpipe('DepoSetFanout', sn_sp, 'FrameFanin', "snsp");

local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;


local graph = g.pipeline([depos, deposio, drifter, bagger, snsp_graph, frameio, sink]);

local app = {
    type: "Pgrapher",
    data: {
        edges: g.edges(graph),
    },
};

// Finally, the configuration sequence which is emitted.


[cli.cmdline] + g.uses(graph) + [app]

