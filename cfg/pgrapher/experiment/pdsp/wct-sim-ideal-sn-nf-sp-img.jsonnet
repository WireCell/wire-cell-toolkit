// This is a main entry point to configure WCT via wire-cell CLI to
// run simulation, noise filtering, signal
// processing, and imaging.


local wc = import "wirecell.jsonnet";
local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";

local params = import "pgrapher/experiment/pdsp/simparams.jsonnet";
local tools_maker = import "pgrapher/common/tools.jsonnet";
local sim_maker = import "pgrapher/experiment/pdsp/sim.jsonnet";
// img.jsonnet exports a zero-arg function; call it to get the maker
// object, whose public method is per_anode(anode, outname_prefix).
local img = (import "pgrapher/experiment/pdsp/img.jsonnet")();

local sp_maker = import "pgrapher/experiment/pdsp/sp.jsonnet";

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
    {
        time: 1*wc.ms,
        charge: -5000,         
        ray: stubby,
    },
    {
        time: 0*wc.ms,
        charge: -5000,
        ray: close0,
    },
   {
       time: 0,
       charge: -5000,         
       ray: params.det.bounds,
   },
];

    
//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
local depos = sim.tracks(tracklist);


local drifter = sim.drifter;
local bagger = sim.make_bagger();

// signal plus noise pipelines
// fixme, make this called more like sp and img....
local sn_pipes = sim.splusn_pipelines;

// This makes a per-APA pipeline.  To avoid unneeded merge/split it
// conseptually spans both simulation parts and sigproc+imaging.  To
// "glue" this into data we need something a bit different.
local make_a_pipe = function(ind) {
    local anode = tools.anodes[ind],
    local aname = anode.name,
    ret : g.pipeline([
        sn_pipes[ind],
        sp.make_sigproc(anode, 'sigproc-'+aname),
        // per_anode = pre_proc -> slicing -> tiling -> solving -> dump.
        // Replaces the stale per-step img.slicing/tiling/solving/dump calls,
        // which used an API img.jsonnet no longer exposes.
        img.per_anode(anode),
    ], "sn-sp-img-" + aname)
}.ret;

local perapa_pipes = [ make_a_pipe(ind) for ind in std.range(0, std.length(tools.anodes) - 1)];

local snspimg = f.fansink('DepoSetFanout', perapa_pipes, "snspimg");

local graph = g.pipeline([depos, drifter, bagger, snspimg]);


local app = {
    type: "Pgrapher",
    data: {
        edges: g.edges(graph),
    },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellImg"],
        apps: ["Pgrapher"]
    },
};

// Finally, the configuration sequence which is emitted.

[cmdline] + g.uses(graph) + [app]

