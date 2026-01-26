# usage: wire-cell -l stdout wct-sim-check.jsonnet
// , nanodes=4
function(input, outname="frame", style='osp', device='cpu', verbosity=0) {
  local g = import 'pgraph.jsonnet',
  local f = import 'pgrapher/common/funcs.jsonnet',
  local wc = import 'wirecell.jsonnet',

  // local io = import 'pgrapher/common/fileio.jsonnet',
  local fileio = import 'layers/high/fileio.jsonnet',

  local tools_maker = import 'pgrapher/common/tools.jsonnet',


  local base = import 'perfect_pdhd/simparams.jsonnet',
  local sp_maker = import 'perfect_pdhd/sp.jsonnet',
  local util = import 'perfect_pdhd/funcs.jsonnet',

  local tpcids = [0],
  local control_mod = import "spng/control.jsonnet",

  local controls = control_mod(device=device, verbosity=wc.intify(verbosity)),
  local detconf = import "spng/detconf.jsonnet",
  local det = detconf.get('pdhd', tpcids),


  local torch_stuff = import 'torch_components.jsonnet',
  local torch_nodes = torch_stuff(det.tpcs[0], controls, outname='sigproc-spng-%s.pt'%outname),
  local spng_graph = torch_nodes.subgraph,

  // local newsim = import 'spng/sim.jsonnet';
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
  local source = fileio.frame_file_source(input),

  local nanodes = std.length(tools.anodes),
  local anode_iota = std.range(0, nanodes-1),
  local anode_idents = [anode.data.ident for anode in tools.anodes],

  local sp_override = {
      sparse: true,
      use_roi_debug_mode: true,
      save_negtive_charge: true,
      use_multi_plane_protection: true,
      process_planes: [0, 1, 2]
  },
  local sp = sp_maker(params, tools, sp_override),
  local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes],


local hio_sp = g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_sp',
      data: {
        anode: wc.tn(tools.anodes[0]),
        trace_tags: [
          'loose_lf0',
          'mp3_roi0',
          'mp2_roi0',
          'gauss0'],
        filename:"sigproc-%s-%s-0.h5"% [style, outname],
        chunk: [0, 0], // ncol, nrow
        gzip: 2,
        high_throughput: true,
      },
    }, nin=1, nout=1),

local sig_reframer = g.pnode({
        type: 'Reframer',
        name: 'reframer-sig-'+tools.anodes[0].name,
        uses: [tools.anodes[0]],
        data: {
            frame_tag: 'sigproc0',
            tags: [
              'loose_lf0',
              'mp3_roi0',
              'mp2_roi0',
              'gauss0'
            ],
            anode: wc.tn(tools.anodes[0]),
            nticks: 6000,
        },
    }, nin=1, nout=1),

local reco_fork = if style == 'osp' then g.pipeline([
              sp_pipes[0],
              sig_reframer,
              hio_sp,
              g.pnode({ type: 'DumpFrames', name: 'reco_fork' }, nin=1, nout=0),
             ],
             'reco_fork') else g.pipeline([
              spng_graph,
              g.pnode({ type: 'DumpFrames', name: 'reco_fork' }, nin=1, nout=0),
             ]),
local graph = g.pipeline([source, reco_fork], 'graph'),

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
},
local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellSpng", "WireCellHio"],
        apps: ["Pgrapher"]
    }
},


  // Finally, the configuration sequence which is emitted.
  ret: [cmdline] + g.uses(graph) + [app]
}.ret