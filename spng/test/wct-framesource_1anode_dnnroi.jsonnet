function(
  input_file='tensor_frames.npz',
) {# usage: wire-cell -l stdout wct-sim-check.jsonnet

local g = import 'pgraph.jsonnet',
local f = import 'pgrapher/common/funcs.jsonnet',
local wc = import 'wirecell.jsonnet',

local io = import 'pgrapher/common/fileio.jsonnet',
local fileio = import 'layers/high/fileio.jsonnet',
local tools_maker = import 'pgrapher/common/tools.jsonnet',

local base = import 'perfect_pdhd/simparams.jsonnet',
local sim_maker = import 'perfect_pdhd/sim.jsonnet',
local perfect = import 'perfect_pdhd/chndb-base.jsonnet',
local nf_maker = import 'perfect_pdhd/nf.jsonnet',
local sp_maker = import 'perfect_pdhd/sp.jsonnet',
local magoutput = 'protodunehd-sim-check.root',

local params = base {
  lar: super.lar {
        // Longitudinal diffusion constant
        DL :  6.2 * wc.cm2/wc.s,
        // Transverse diffusion constant
        DT : 16.3 * wc.cm2/wc.s,
        lifetime : 50*wc.ms,
        drift_speed : 1.565*wc.mm/wc.us,
  },
},

local tools = tools_maker(params),
local nanodes = std.length(tools.anodes),

// local outtags = ['raw%d' % n for n in std.range(0, std.length(tools.anodes) - 1)],


local frame_input = fileio.frame_tensor_file_source(input_file),


local selectors = [
  g.pnode({
      type: 'ChannelSelector',
      name: 'chanselect-'+tools.anodes[n].name,
      data: {
          anode: wc.tn(tools.anodes[n]),
          tags: ['raw'],           // ?? what do?
          channels: std.range(2560 * n, 2560 * (n + 1) - 1)
      },
  }, nin=1, nout=1) for n in std.range(0, nanodes-1)
],

// local simple_pipes = [
//   g.pipeline([
//                nf_pipes[n],
//                selectors[n],
//              ],
//              'parallel_pipe_%d' % n)
//   for n in std.range(0, std.length(tools.anodes) - 1)
// ],

local fanout_apa_rules =
[
    {
        frame: {
            //'.*': 'number%d' % n,
            //'.*': 'gauss%d' % n,
            //'.*': 'framefanout%d ' % n,
            // '.*': 'orig%d' % n,
        },
        trace: {
            // fake doing Nmult SP pipelines
            //orig: ['wiener', 'gauss'],
            // gauss: 'gauss%d' % n, //uncommented Ewerton 2023-09-27
            // wiener: 'wiener%d' % n, //created Ewerton 2023-09-27
            // ['orig%d'%n]: 'orig',
            ['raw%d'%n]: 'raw',
            // ['orig%d'%n]: 'orig',
        },
    }
    for n in std.range(0, std.length(tools.anodes) - 1)
],

// local parallel_graph = f.fanpipe('FrameFanout', parallel_pipes, 'FrameFanin', 'sn_mag_nf', outtags, fanout_apa_rules),
// // local fanout_graph = g.fan.fanout('FrameFanout', simple_pipes, 'sn_mag_nf', fanout_apa_rules),
// local fanout_graph = g.fan.fanout('FrameFanout', simple_pipes, 'sn_mag_nf', fanout_apa_rules),



local torch_maker = import 'torch_1anode_dnnroi.jsonnet',
local torch_nodes = torch_maker(
  tools,
),
local spng_stacked = torch_nodes.stacked_spng,


// local fanout = g.pnode({
//     type: "FrameFanout",
//     name: "sn_mag_nf",
//     data: {
//         multiplicity: 4,
//         tag_rules: fanout_apa_rules,
//     },
// }, nin=1, nout=4),

// local load_to_fanout = g.intern(
//   innodes=[frame_input],
//   outnodes=[fanout],
//   edges = [
//     g.edge(frame_input, fanout),
//   ]
// ),

local graph = 
    g.intern(
      innodes=[frame_input],
      centernodes=selectors,
      outnodes=[spng_stacked],
      edges = [
        g.edge(frame_input, selectors[0], 0),
        g.edge(selectors[0], spng_stacked, 0, 0),
      ]
    ),

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
},

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: [
          "WireCellGen", "WireCellPgraph", "WireCellSio", 
          "WireCellSigProc",
          "WireCellRoot",
          "WireCellSpng"],
        apps: ["Pgrapher"]
    }
},


// Finally, the configuration sequence which is emitted.
ret: [cmdline] + g.uses(graph) + [app]
}.ret
