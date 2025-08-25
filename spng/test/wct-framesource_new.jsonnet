function(
  input_file='tensor_frames.npz',
  skip_noise='False',
  spng_flag = 1,
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

local sim = sim_maker(params, tools),

local nanodes = std.length(tools.anodes),
local anode_iota = std.range(0, nanodes-1),
local anode_idents = [anode.data.ident for anode in tools.anodes],

local output = 'wct-sim-ideal-sigproc.npz',
local drifter = sim.drifter,
local bagger = sim.make_bagger(),
local sn_pipes = sim.splusn_pipelines,

local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota],

local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)],

local sp_override = {
    sparse: true,
    use_roi_debug_mode: true,
    save_negtive_charge: true,
    use_multi_plane_protection: true,
    process_planes: [0, 1, 2],
    debug_no_frer: false,
    debug_no_wire_filter: false,
},

local sp = sp_maker(params, tools, sp_override),
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes],

// local magnify = import 'perfect_pdhd/magnify-sinks.jsonnet',
local magnify = import 'magnify-sinks.jsonnet',
local magnifyio = magnify(tools, magoutput),

local outtags = ['raw%d' % n for n in std.range(0, std.length(tools.anodes) - 1)],


local frame_input = fileio.frame_tensor_file_source(input_file),
// local frame_input = fileio.frame_file_source('frames.tar', tags=outtags),

local parallel_pipes = [
  g.pipeline([

               nf_pipes[n],
               sp_pipes[n],
               magnifyio.decon_pipe[n],
             ],
             'parallel_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
],


local simple_pipes = [
  g.pipeline([
               nf_pipes[n],
             ],
             'parallel_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
],

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
            ['orig%d'%n]: 'orig',
            // ['orig%d'%n]: 'orig',
        },
    }
    for n in std.range(0, std.length(tools.anodes) - 1)
],

local parallel_graph = f.fanpipe('FrameFanout', parallel_pipes, 'FrameFanin', 'sn_mag_nf', outtags, fanout_apa_rules),
// local fanout_graph = g.fan.fanout('FrameFanout', simple_pipes, 'sn_mag_nf', fanout_apa_rules),
local fanout_graph = g.fan.fanout('FrameFanout', simple_pipes, 'sn_mag_nf', fanout_apa_rules),



local torch_maker = import 'torch2.jsonnet',
local torch_nodes = torch_maker(
  tools,
  // debug_force_cpu=false,
  // apply_gaus=(std.extVar("ApplyGaus") == 1),
  // do_roi_filters=(std.extVar("ROI") == 1),
  // do_collate_apa=(std.extVar("CollateAPAs") == 1),
  // do_run_roi=(std.extVar("RunROI") == 1),
  // do_tiling=(std.extVar("DoTiling") == 1),
),
local spng_decons = torch_nodes.spng_decons,
local spng_stacked = torch_nodes.stacked_spng,

local load_to_fanout = g.intern(
  innodes=[frame_input],
  outnodes=[fanout_graph],
  edges = [
    g.edge(frame_input, fanout_graph),
  ]
),

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

local sink = sim.frame_sink,

//local spng_flag = std.extVar("SPNG"),

local graph = if (spng_flag == 0) then
    g.pipeline([frame_input, parallel_graph, sink])

  else if (spng_flag == 1) then
    g.intern(
      innodes=[load_to_fanout],
      outnodes=spng_decons,
      centernodes=[],
      edges = [
        g.edge(load_to_fanout, spng_decons[0], 0),
        g.edge(load_to_fanout, spng_decons[1], 1),
        g.edge(load_to_fanout, spng_decons[2], 2),
        g.edge(load_to_fanout, spng_decons[3], 3),
      ]
    )

  else if (spng_flag == 2) then
    g.intern(
      innodes=[load_to_fanout],
      outnodes=[spng_stacked],
      edges = [
        g.edge(load_to_fanout, spng_stacked, 0, 0),
        g.edge(load_to_fanout, spng_stacked, 1, 1),
        g.edge(load_to_fanout, spng_stacked, 2, 2),
        g.edge(load_to_fanout, spng_stacked, 3, 3),
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
          "WireCellSigProc", "WireCellRoot", "WireCellSpng"],
        apps: ["Pgrapher"]
    }
},


// Finally, the configuration sequence which is emitted.
ret: [cmdline] + g.uses(graph) + [app]
}.ret
