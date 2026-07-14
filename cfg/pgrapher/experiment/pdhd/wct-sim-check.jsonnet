# usage: wire-cell -l stdout wct-sim-check.jsonnet

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local base = import 'pgrapher/experiment/pdhd/simparams.jsonnet';
local params = base {
  lar: super.lar {
        // Longitudinal diffusion constant
        DL :  6.2 * wc.cm2/wc.s,
        // Transverse diffusion constant
        DT : 16.3 * wc.cm2/wc.s,
        lifetime : 50*wc.ms,
        drift_speed : 1.565*wc.mm/wc.us,
  },
};

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/pdhd/sim.jsonnet';
local sim = sim_maker(params, tools);

local beam_dir = [-0.178177, -0.196387, 0.959408];
local beam_center = [-27.173, 421.445, 0];

local track0 = {
  head: wc.point(beam_center[0], beam_center[1], beam_center[2], wc.cm),
  tail: wc.point(beam_center[0] + 100*beam_dir[0], beam_center[1] + 100*beam_dir[1], beam_center[2] + 100*beam_dir[2], wc.cm),

  // head: wc.point(-260,300, 50,wc.cm), // apa1
  // tail: wc.point(-260,300, 200,wc.cm), // apa1 w
  // tail: wc.point(-260,300 - 0.58364 * 300, 50 + 0.812013 * 300 ,wc.cm), // apa1 u
  // tail: wc.point(-260,300 + 0.58364 * 300, 50 + 0.812013 * 300 ,wc.cm), // apa1 v

  // head: wc.point(260,300, 50,wc.cm), // apa2
  // tail: wc.point(260,300, 200,wc.cm), // apa2 w
  // tail: wc.point(260,300 + 0.58364 * 300, 50 + 0.812013 * 300,wc.cm), // apa2 u
  // tail: wc.point(260,300 - 0.58364 * 300, 50 + 0.812013 * 300,wc.cm), // apa2 v

};


local tracklist = [

  {
    time: 0 * wc.us,
    charge: -500, // negative means # electrons per step (see below configuration) 
    ray: track0, // params.det.bounds,
  },

];

local depos = sim.tracks(tracklist, step=0.1 * wc.mm); // MIP <=> 5000e/mm

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes-1);
local anode_idents = [anode.data.ident for anode in tools.anodes];

// local output = 'wct-sim-ideal-sig.npz';
// local deposio = io.numpy.depos(output);
local drifter = sim.drifter;
local bagger = sim.make_bagger();
// signal plus noise pipelines
local sn_pipes = sim.splusn_pipelines;
// local analog_pipes = sim.analog_pipelines;

local perfect = import 'pgrapher/experiment/pdhd/chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];

local nf_maker = import 'pgrapher/experiment/pdhd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp_override = {
    sparse: true,
    use_roi_debug_mode: false,
    use_multi_plane_protection: true,
    process_planes: [0, 1, 2]
};

local sp_maker = import 'pgrapher/experiment/pdhd/sp.jsonnet';
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local magoutput = 'protodunehd-sim-check.root';
local magnify = import 'pgrapher/experiment/pdhd/magnify-sinks.jsonnet';
local magnifyio = magnify(tools, magoutput);

local parallel_pipes = [
  g.pipeline([
               sn_pipes[n],
               // magnifyio.orig_pipe[n],
               nf_pipes[n],
               magnifyio.raw_pipe[n],
               sp_pipes[n],
               // magnifyio.debug_pipe[n],
               magnifyio.decon_pipe[n],
             ],
             'parallel_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];
local outtags = ['raw%d' % n for n in std.range(0, std.length(tools.anodes) - 1)];
local parallel_graph = f.fanpipe('DepoSetFanout', parallel_pipes, 'FrameFanin', 'sn_mag_nf', outtags);

//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;
local graph = g.pipeline([depos, drifter, bagger, parallel_graph, sink]);


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


// Finally, the configuration sequence which is emitted.

[cmdline] + g.uses(graph) + [app]
