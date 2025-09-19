function(
  input_file='tensor_frames.npz',
  do_spng=true,
  ts_model_file="/nfs/data/1/abashyal/spng/spng_dev_050525/Pytorch-UNet/ts-model-2.3/unet-l23-cosmic500-e50.ts",
  device='gpu',
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
          tags: ['raw0'],           // ?? what do?
          channels: std.range(2560 * n, 2560 * (n + 1) - 1)
      },
  }, nin=1, nout=1) for n in std.range(0, nanodes-1)
],

local sp_override = {
    sparse: false,
    use_roi_debug_mode: true,
    use_roi_refinement: true,
    save_negtive_charge: false,
    use_multi_plane_protection: true,
    process_planes: [0, 1, 2],
    debug_no_frer: false,
    debug_no_wire_filter: false,
    break_roi_loop1_tag: "",
    break_roi_loop2_tag: "",
    shrink_roi_tag: "",
    extend_roi_tag: "",
    // decon_charge_tag: "",
    do_not_mp_protect_traditional: true, // do_not_mp_protect_traditional to 
                                         // make a clear ref, defualt is false
    mp_tick_resolution: 4,
},

local sp = sp_maker(params, tools, sp_override),
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes],
local dnnroi = import 'perfect_pdhd/dnnroi.jsonnet',
local ts = {
    type: "TorchService",
    name: "dnnroi",
    data: {
        model: ts_model_file,//"ts-model/unet-cosmic390-newwc-depofluxsplat-pdhd.ts",
        device: (if device == 'gpu' then 'gpucpu' else 'cpu'), // "gpucpu",
        concurrency: 1,
    },
},

local toolkit_pipe = g.pipeline([
  sp_pipes[0],
  dnnroi(tools.anodes[0], ts, output_scale=1.0, nticks=params.daq.nticks, nchunks=1)
]),


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
  ts_model_file=ts_model_file,
  debug_force_cpu=(device=='cpu'),
),
local spng_stacked = torch_nodes.stacked_spng,

local frame_output = fileio.frame_tensor_file_sink('toolkit_dnnroi_output.tar', mode='tagged'),
local sim = sim_maker(params, tools),
local sink = sim.frame_sink,

local sp_graph = g.intern(
    innodes=[frame_input],
    centernodes=selectors + [toolkit_pipe],
    outnodes=[frame_output],
    edges = [
      g.edge(frame_input, selectors[0]),
      g.edge(selectors[0], toolkit_pipe),
      g.edge(toolkit_pipe, frame_output),
    ]
),

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

local target_graph = (if do_spng then graph else sp_graph),

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(target_graph),
  },
},

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: [
          "WireCellGen",
          "WireCellPgraph",
          "WireCellSio", 
          "WireCellSigProc",
          "WireCellRoot",
          "WireCellSpng"] + (if do_spng then [] else ["WireCellPytorch"]),
        apps: ["Pgrapher"]
    }
},


// Finally, the configuration sequence which is emitted.
ret: [cmdline] + g.uses(target_graph) + [app]
}.ret
