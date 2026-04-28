local g = import "pgraph.jsonnet";
local f = import "pgrapher/experiment/dune-vd/funcs.jsonnet";
local wc = import "wirecell.jsonnet";

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params = import 'pgrapher/experiment/dune-vd/params.jsonnet';

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/dune-vd/sim.jsonnet';
local sim = sim_maker(params, tools);

// Deposit and drifter ///////////////////////////////////////////////////////////////////////////////

// APA 35
local stubby = {
  tail: wc.point(290, 340.0,   755.0, wc.cm),
  head: wc.point(50,  508.101,  905.0, wc.cm),
};

local tracklist = [
  {
    time: 0,
    charge: -5000,
    ray: stubby, // params.det.bounds,
  },
];


local depos = g.pnode({
        type: 'TrackDepos',
        data: {
            step_size: 1.0*wc.mm,
            tracks: tracklist
        },
}, nin=0, nout=1);


local depos = sim.tracks(tracklist);
local drifter = sim.drifter;
local bagger = sim.make_bagger();


// Parallel part //////////////////////////////////////////////////////////////////////////////


// local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;


local parallel_pipes = [ g.pipeline([ sn_pipes[n] ], 'parallel_pipe_%d' % n)  for n in std.range(0, std.length(tools.anodes) - 1)];

local outtags = ['orig%d' % n for n in std.range(0, std.length(tools.anodes) - 1)];
local parallel_graph = f.multifanpipe('DepoSetFanout', parallel_pipes, 'FrameFanin', 6, 'sn_mag', outtags);


// Frame save  ////////////////////////////////////////////////////////////////////////////

local npzsink = g.pnode({
    type: 'NumpyFrameSaver',
    name: 'npzsink',
    data: {
        filename:  'dune-vd-sim-check.npz',
        digitize:  true, 
        frame_tags: ['orig%d' % n for n in std.range(0, std.length(tools.anodes) - 1)]
    },
  }, nin=1, nout=1);

local dumpframes = g.pnode({
        type: "DumpFrames",
        name: 'dumpframe',
    }, nin=1, nout=0);


// Final pipeline //////////////////////////////////////////////////////////////////////////////

local graph = g.pipeline([depos, drifter, bagger, parallel_graph, npzsink, dumpframes], "main");

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot"],
        apps: ["Pgrapher"]
    }
};

[cmdline] + g.uses(graph) + [app]