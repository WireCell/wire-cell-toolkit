# usage: wire-cell -l stdout wct-sim-check.jsonnet
function(outname="tensor_frames.npz", skip_noise="False", do_depo="False") {
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

  local beam_dir = [-0.178177, -0.196387, 0.959408],
  local beam_center = [-27.173, 421.445, 0],
  local beam_center1 = [27.173, 421.445, 0],

  local track0 = {
    head: wc.point(beam_center[0], beam_center[1], beam_center[2], wc.cm),
    tail: wc.point(beam_center[0] + 100*beam_dir[0], beam_center[1] + 100*beam_dir[1], beam_center[2] + 100*beam_dir[2], wc.cm),
  },

  local track1 = {
    head: wc.point(beam_center1[0], beam_center1[1], beam_center1[2], wc.cm),
    tail: wc.point(beam_center1[0] + 100*beam_dir[0], beam_center1[1] + 100*beam_dir[1], beam_center1[2] + 100*beam_dir[2], wc.cm),
  },

  local depo_track = {
    head: wc.point(beam_center[0], beam_center[1], 50.0, wc.cm),
    tail: wc.point(beam_center[0], beam_center[1], 50.1, wc.cm),
  },
  local tracklist = [

    {
      time: 0 * wc.us,
      charge: -500, // negative means # electrons per step (see below configuration) 
      ray: (if do_depo=="False" then track0 else depo_track), // params.det.bounds,
    },

    // {
    //   time: 0 * wc.us,
    //   charge: -500, // negative means # electrons per step (see below configuration) 
    //   ray: track1, // params.det.bounds,
    // },
  ],

  local depos = sim.tracks(tracklist, step=0.1 * wc.mm), // MIP <=> 5000e/mm


  local nanodes = std.length(tools.anodes),
  local anode_iota = std.range(0, nanodes-1),
  local anode_idents = [anode.data.ident for anode in tools.anodes],

  local drifter = sim.drifter,
  local bagger = sim.make_bagger(),
  local sn_pipes = sim.splusn_pipelines,
  local signal_pipes = sim.signal_pipelines,

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
      process_planes: [0, 1, 2]
  },

  local sp = sp_maker(params, tools, sp_override),
  local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes],


  


  local parallel_pipes = [
    g.pipeline([
                (if skip_noise == "False" then sn_pipes[n] else signal_pipes[n]),
                // sn_pipes[n]
                nf_pipes[n],
              ],
              'parallel_pipe_%d' % n)
    for n in std.range(0, std.length(tools.anodes) - 1)
  ],
  // local outtags = ['orig%d' % n for n in std.range(0, std.length(tools.anodes) - 1)],
  local outtags = ['raw%d' % n for n in std.range(0, std.length(tools.anodes) - 1)],
  local parallel_graph = f.fanpipe('DepoSetFanout', parallel_pipes, 'FrameFanin', 'sn_mag_nf', outtags),


  // local outname =
  //   (try std.extVar("OutName")
    //  catch "tensor_frames.npz"),

  local frame_output = fileio.frame_tensor_file_sink(outname, mode='tagged'),

  local graph = g.pipeline([depos, drifter, bagger, parallel_graph, frame_output]),


  local app = {
    type: 'Pgrapher',
    data: {
      edges: g.edges(graph),
    },
  },

  local cmdline = {
      type: "wire-cell",
      data: {
          plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellRoot"],
          apps: ["Pgrapher"]
      }
  },


  // Finally, the configuration sequence which is emitted.

  ret: [cmdline] + g.uses(graph) + [app]
}.ret