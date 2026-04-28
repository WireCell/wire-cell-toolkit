local g = import "pgraph.jsonnet";
local f = import "pgrapher/experiment/dune-vd/funcs.jsonnet";
local wc = import "wirecell.jsonnet";

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local param_maker = import 'pgrapher/experiment/dune-vd/params.jsonnet';
local params = param_maker(10*wc.cm) {
};

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/dune-vd/sim.jsonnet';
local sim = sim_maker(params, tools);

// Deposit and drifter ///////////////////////////////////////////////////////////////////////////////

// APA 35
// local stubby = {
//   tail: wc.point(290, 340.0,   755.0, wc.cm),
//   head: wc.point(50,  508.101,  905.0, wc.cm),
// };

// max diagonal
// local stubby = {
//   tail: wc.point(290, -508,   0.3, wc.cm),
//   head: wc.point(50,  508.101,  905.0, wc.cm),
// };

// horizontal line
local stubby = {
  tail: wc.point(140.507, -350,   0.3, wc.cm),
  head: wc.point(140.507,  -350,  700, wc.cm),
};

local tracklist = [
  {
    time: 0*wc.ms,
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

local drifter = sim.drifter;
local bagger = sim.make_bagger();

// Magnify definition /////////////////////////////////////////////////////////////////////////

// Origin Magnify

local origmagnify = [ 
  g.pnode({
    type: 'MagnifySink',
    name: 'origmag%d' % n,
    data: {
        output_filename: 'dune-vd-sim-check.root',
        root_file_mode: 'UPDATE',
        frames: ['orig%d' % n ],
        trace_has_tag: false,
        anode: wc.tn(tools.anodes[n]), 
    },
  }, nin=1, nout=1) for n in std.range(0, std.length(tools.anodes) - 1)];


local origmagnify_pipe = [g.pipeline([origmagnify[n]], name='origmagnifypipes%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];


// Parallel part //////////////////////////////////////////////////////////////////////////////


// local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

local sp_maker = import 'pgrapher/experiment/dune-vd/sp.jsonnet';
local sp = sp_maker(params, tools, { use_roi_debug_mode: false,} );
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local spmagnify = [ 
g.pnode({
    type: 'MagnifySink',
    name: 'spmag%d' % n,
    data: {
        output_filename: 'dune-vd-sim-check.root',
        root_file_mode: 'UPDATE',
        frames: ['gauss%d' % n ],
        trace_has_tag: false,
        anode: wc.tn(tools.anodes[n]), 
    },
  }, nin=1, nout=1) for n in std.range(0, std.length(tools.anodes) - 1)];
local spmagnify_pipe = [g.pipeline([spmagnify[n]], name='spmagnifypipes%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local magoutput = 'sim-check.root';
local magnify = import 'pgrapher/experiment/pdsp/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);

local parallel_pipes = [
  g.pipeline([ 
                sn_pipes[n],
                // origmagnify_pipe[n],
                sinks.orig_pipe[n],
                sp_pipes[n],
                // spmagnify_pipe[n],
                sinks.decon_pipe[n],
                sinks.debug_pipe[n], // use_roi_debug_mode=true in sp.jsonnet
          ], 
          'parallel_pipe_%d' % n) 
  for n in std.range(0, std.length(tools.anodes) - 1)];

local outtags = ['orig%d' % n for n in std.range(0, std.length(tools.anodes) - 1)];
local parallel_graph = f.multifanpipe('DepoSetFanout', parallel_pipes, 'FrameFanin', 6, 'sn_mag', outtags);


// Only one sink ////////////////////////////////////////////////////////////////////////////


local sink = sim.frame_sink;


// Final pipeline //////////////////////////////////////////////////////////////////////////////

local graph = g.pipeline([depos, drifter, bagger, parallel_graph, sink], "main");

local app = {
  type: 'Pgrapher', //Pgrapher, TbbFlow
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot", "WireCellTbb"],
        apps: ["Pgrapher"]
    }
};

[cmdline] + g.uses(graph) + [app]
