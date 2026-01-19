# usage: wire-cell -l stdout wct-sim-check.jsonnet
// , nanodes=4
function(input, outname="frame.npz") {
  local g = import 'pgraph.jsonnet',
  local f = import 'pgrapher/common/funcs.jsonnet',
  local wc = import 'wirecell.jsonnet',

  // local io = import 'pgrapher/common/fileio.jsonnet',
  local fileio = import 'layers/high/fileio.jsonnet',

  local tools_maker = import 'pgrapher/common/tools.jsonnet',


  local base = import 'perfect_pdhd/simparams.jsonnet',
  local sim_maker = import 'perfect_pdhd/sim.jsonnet',
  local perfect = import 'perfect_pdhd/chndb-base.jsonnet',
  local nf_maker = import 'perfect_pdhd/nf.jsonnet',
  local sp_maker = import 'perfect_pdhd/sp.jsonnet',
  local util = import 'perfect_pdhd/funcs.jsonnet',

  // local newsim = import 'spng/sim.jsonnet';
  local frame_sink = sim.frame_sink,
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

  // local source = fileio.depo_file_source(input),
  // local source = g.pnode({
  //     type: 'NumpyDepoLoader',
  //     data: { filename: input }
  // }, nin=0, nout=1),
  local source = g.pnode({
      type: 'NumpyDepoSetLoader',
      data: { filename: input }
  }, nin=0, nout=1),

  local nanodes = std.length(tools.anodes),
  local anode_iota = std.range(0, nanodes-1),
  local anode_idents = [anode.data.ident for anode in tools.anodes],
  local deposplats = util.splatonly(params, tools, tools.anodes[0]),

  local drifter = sim.drifter,
  local setdrifter = g.pnode({
            type: 'DepoSetDrifter',
            data: {
                drifter: "Drifter"
            }
        }, nin=1, nout=1,
        uses=[drifter]),

  local sn_pipes = sim.splusn_pipelines,

  local chndb = [{
    type: 'OmniChannelNoiseDB',
    name: 'ocndbperfect%d' % n,
    data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
    uses: [tools.anodes[n], tools.field, tools.dft],
  } for n in anode_iota],
  local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)],

  //TBD Remove
  local sp_override = {
      sparse: true,
      use_roi_debug_mode: true,
      save_negtive_charge: true,
      use_multi_plane_protection: true,
      process_planes: [0, 1, 2]
  },
  local sp = sp_maker(params, tools, sp_override),
  local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes],



  local depo_fanout = g.pnode({
    type:'DepoSetFanout',
    name:'deposet_fanout',
    data:{
        multiplicity:2,
        tags: [], 
    }}, nin=1, nout=2),


  local parallel_pipes = [
    g.pipeline([
                sn_pipes[n],
                nf_pipes[n],
                // sp_pipes[n],
              ],
              'parallel_pipe_%d' % n)
    for n in std.range(0, std.length(tools.anodes) - 1)
  ],
  local outtags = ['raw%d' % n for n in std.range(0, std.length(tools.anodes) - 1)],
  local parallel_graph = f.fanpipe('DepoSetFanout', parallel_pipes, 'FrameFanin', 'sn_mag_nf', outtags),


  // local drifted_frame_output = fileio.frame_file_sink(
  //   'drifted-%s' % outname, tags=[]
  // ),
  local splat_frame_output = fileio.frame_file_sink(
    outname, tags=[]
  ),


  local splat_reframers = [
      g.pnode({
          type: 'Reframer',
          name: 'reframer-splat-'+a.name,
          data: {
              anode: wc.tn(a),
              nticks: 6000,
          },
      }, nin=1, nout=1) for a in tools.anodes],

  local sig_reframers = [
      g.pnode({
          type: 'Reframer',
          name: 'reframer-sig-'+a.name,
          data: {
              anode: wc.tn(a),
              nticks: 6000,
          },
      }, nin=1, nout=1) for a in tools.anodes],
local dump = g.pnode({type: 'DumpDepos', name: 'dump_depos', data: {}}, nin=1, nout=0),
local sink = g.pnode({
    name: "deposink",
    type: "DepoFileSink",
    data: {
        outname: 'test.tar.bz2',
    }
}, nin=1, nout=0),
local graph = g.pipeline([source, setdrifter, deposplats, splat_reframers[0], splat_frame_output], 'graph'),
  // local graph = g.intern(
  //   innodes=[source],
  //   outnodes=[
  //     splat_frame_output,
  //     drifted_frame_output,
  //     // frame_sink
  //   ],
  //   centernodes=[
  //     setdrifter,
  //     parallel_graph,
  //     depo_fanout,
  //     deposplats,
  //   ]+ splat_reframers + sig_reframers,
  //   edges=[
  //     g.edge(source, setdrifter),
  //     // g.edge(setdrifter, parallel_graph),
  //     g.edge(setdrifter, depo_fanout),
  //     g.edge(depo_fanout, parallel_graph, 0, 0),
  //     g.edge(depo_fanout, deposplats, 1, 0),
  //     g.edge(deposplats, splat_reframers[0]),
  //     g.edge(parallel_graph, drifted_frame_output),
  //     g.edge(splat_reframers[0], splat_frame_output),
  //   ],
  // ),


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