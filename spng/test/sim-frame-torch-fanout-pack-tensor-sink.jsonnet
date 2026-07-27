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
              //  magnifyio.raw_pipe[n],
              //  sp_pipes[n],
              //  // magnifyio.debug_pipe[n],
              //  magnifyio.decon_pipe[n],
             ],
             'parallel_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];
local outtags = ['raw%d' % n for n in std.range(0, std.length(tools.anodes) - 1)];

local parallel_graph = f.fanpipe('DepoSetFanout', parallel_pipes, 'FrameFanin', 'sn_mag_nf', outtags);
local depo_sn_nf_fanout = g.fan.fanout('DepoSetFanout', parallel_pipes, 'sn_mag_nf');

local tf_maker = import 'torch.jsonnet';
local tf = tf_maker(tools);
local tf_fans = [tf.make_fanout(a) for a in tools.anodes];

// local frame_to_torch_fanout = g.pnode({
//   type: "FrameToTorchFanout",
//   data: {
//     planes: [800, 800, 960],
//     expected_nticks: 6000,
//     multiplicity: 3,

//     channel_ranges: [
//       [2, [[0,480], [2080, 2560]]],
//       [0, [[480, 1280]]],
//       [1, [[1280, 2080]]],
//     ]

//   }
// }, nin=1, nout=3);

local tensor_sink = g.pnode({
    type: 'TensorFileSink',
    name: 'tfsink',
    data: {
        outname: 'testout2.tar',
        prefix: ''
    },
}, nin=1, nout=0);


local tensor_sinks = [
   g.pnode({
    type: 'TensorFileSink',
    name: 'tfsink%d' % n,
    data: {
        outname: 'testout_fan%d.tar' % n,
        prefix: ''
    },
  //  }, nin=1, nout=0) for n in std.range(0, std.length(tools.anodes) - 1)
   }, nin=1, nout=0) for n in std.range(0, 15)
];

//local frameio = io.numpy.frames(output);
// local sink = frame_sink;
// local graph = g.pipeline([depos, drifter, bagger, parallel_graph, sink]);

local torch_packer = g.pnode({
  type: 'TorchPacker',
  name: 'torchpacker',
  data: {
    'multiplicity': 3
  },
}, nin=3, nout=1);

local torch_packers = [g.pnode({
  type: 'TorchPacker',
  name: 'torchpacker%d' % n,
  data: {
    'multiplicity': 3
  },
}, nin=3, nout=1) for n in std.range(0, std.length(tools.anodes)-1)];

local torch_to_tensor = g.pnode({
  type: 'TorchToTensor',
  name: 'torchtotensor',
  data: {},
}, nin=1, nout=1);

local torch_to_tensors = [g.pnode({
  type: 'TorchToTensor',
  name: 'torchtotensor%d' % n,
  data: {},
// }, nin=1, nout=1) for n in std.range(0, std.length(tools.anodes) - 1)];
}, nin=1, nout=1) for n in std.range(0, 15)];

// local spng_sigproc = g.pnode({
//   type: 'SPNGSigProc',
//   name: 'spngsigproc',
//   data: {}
// }, nin=1, nout=1);

// local e0 = g.edge(frame_to_torch_fanout, torch_packer, 0, 0);
// local e1 = g.edge(frame_to_torch_fanout, torch_packer, 1, 1);
// local e2 = g.edge(frame_to_torch_fanout, torch_packer, 2, 2);


// local connector = g.intern(
//   innodes = [frame_to_torch_fanout],
//   outnodes = [torch_packer],
//   edges = [e0, e1, e2],
// );

#4 x 3-->1
// local tfs_to_packers = g.intern(
//   innodes = tf_fans,
//   outnodes = torch_packers,
//   edges =  [
//     g.edge(tf_fans[0], torch_packers[0], 0, 0),
//     g.edge(tf_fans[0], torch_packers[0], 1, 1),
//     g.edge(tf_fans[0], torch_packers[0], 2, 2),
//     g.edge(tf_fans[1], torch_packers[1], 0, 0),
//     g.edge(tf_fans[1], torch_packers[1], 1, 1),
//     g.edge(tf_fans[1], torch_packers[1], 2, 2),
//     g.edge(tf_fans[2], torch_packers[2], 0, 0),
//     g.edge(tf_fans[2], torch_packers[2], 1, 1),
//     g.edge(tf_fans[2], torch_packers[2], 2, 2),
//     g.edge(tf_fans[3], torch_packers[3], 0, 0),
//     g.edge(tf_fans[3], torch_packers[3], 1, 1),
//     g.edge(tf_fans[3], torch_packers[3], 2, 2)
//   ]
// );

// local tfs_to_packers_to_sinks = g.intern(
//     innodes=[tfs_to_packers],
//     centernodes=torch_to_tensors,
//     outnodes=tensor_sinks,
//     edges=[
//       g.edge(tfs_to_packers, torch_to_tensors[0], 0),
//       g.edge(tfs_to_packers, torch_to_tensors[1], 1),
//       g.edge(tfs_to_packers, torch_to_tensors[2], 2),
//       g.edge(tfs_to_packers, torch_to_tensors[3], 3),
//       g.edge(torch_to_tensors[0], tensor_sinks[0]),
//       g.edge(torch_to_tensors[1], tensor_sinks[1]),
//       g.edge(torch_to_tensors[2], tensor_sinks[2]),
//       g.edge(torch_to_tensors[3], tensor_sinks[3]),
//     ],
// );

local tfs_to_tensors = g.intern(
    innodes=tf_fans,
    // centernodes=[],
    outnodes=torch_to_tensors,
    edges=[
      g.edge(tf_fans[0], torch_to_tensors[0], 0),
      g.edge(tf_fans[0], torch_to_tensors[1], 1),
      g.edge(tf_fans[0], torch_to_tensors[2], 2),
      g.edge(tf_fans[0], torch_to_tensors[3], 3),

      g.edge(tf_fans[1], torch_to_tensors[4], 0),
      g.edge(tf_fans[1], torch_to_tensors[5], 1),
      g.edge(tf_fans[1], torch_to_tensors[6], 2),
      g.edge(tf_fans[1], torch_to_tensors[7], 3),

      g.edge(tf_fans[2], torch_to_tensors[8], 0),
      g.edge(tf_fans[2], torch_to_tensors[9], 1),
      g.edge(tf_fans[2], torch_to_tensors[10], 2),
      g.edge(tf_fans[2], torch_to_tensors[11], 3),

      g.edge(tf_fans[3], torch_to_tensors[12], 0),
      g.edge(tf_fans[3], torch_to_tensors[13], 1),
      g.edge(tf_fans[3], torch_to_tensors[14], 2),
      g.edge(tf_fans[3], torch_to_tensors[15], 3),

      // g.edge(torch_to_tensors[0], tensor_sinks[0]),
      // g.edge(torch_to_tensors[1], tensor_sinks[1]),
      // g.edge(torch_to_tensors[2], tensor_sinks[2]),
      // g.edge(torch_to_tensors[3], tensor_sinks[3]),
    ],
);

local tfs_to_tensors_to_sinks = g.intern(
    innodes=[tfs_to_tensors],
    // centernodes=[],
    outnodes=tensor_sinks,
    edges=[
      g.edge(tfs_to_tensors, tensor_sinks[0], 0),
      g.edge(tfs_to_tensors, tensor_sinks[1], 1),
      g.edge(tfs_to_tensors, tensor_sinks[2], 2),
      g.edge(tfs_to_tensors, tensor_sinks[3], 3),
      
      g.edge(tfs_to_tensors, tensor_sinks[4], 4),
      g.edge(tfs_to_tensors, tensor_sinks[5], 5),
      g.edge(tfs_to_tensors, tensor_sinks[6], 6),
      g.edge(tfs_to_tensors, tensor_sinks[7], 7),
      
      g.edge(tfs_to_tensors, tensor_sinks[8], 8),
      g.edge(tfs_to_tensors, tensor_sinks[9], 9),
      g.edge(tfs_to_tensors, tensor_sinks[10], 10),
      g.edge(tfs_to_tensors, tensor_sinks[11], 11),

      g.edge(tfs_to_tensors, tensor_sinks[12], 12),
      g.edge(tfs_to_tensors, tensor_sinks[13], 13),
      g.edge(tfs_to_tensors, tensor_sinks[14], 14),
      g.edge(tfs_to_tensors, tensor_sinks[15], 15),

    ],
);


// local graph = g.pipeline([
//   depos, drifter, bagger, parallel_graph,
//   // frame_to_torch_fanout,
//   // // // spng_sigproc,
//   // torch_packer,
//   connector,
//   torch_to_tensor,
//   tensor_sink,
// ]);



local depos_to_fanout = g.intern(
  innodes=[depos],
  centernodes=[drifter, bagger], 
  outnodes=[depo_sn_nf_fanout],
  edges = [
    g.edge(depos, drifter),
    g.edge(drifter, bagger),
    g.edge(bagger, depo_sn_nf_fanout),
  ]
);

local new_graph = g.intern(
  innodes=[depos_to_fanout],
  // outnodes=[to_sinks],
  // outnodes=[tfs_to_packers_to_sinks],
  outnodes=[tfs_to_tensors_to_sinks],
  centernodes=[],
  edges = [
    g.edge(depos_to_fanout, tfs_to_tensors_to_sinks, 0, 0),
    g.edge(depos_to_fanout, tfs_to_tensors_to_sinks, 1, 1),
    g.edge(depos_to_fanout, tfs_to_tensors_to_sinks, 2, 2),
    g.edge(depos_to_fanout, tfs_to_tensors_to_sinks, 3, 3),
    // g.edge(depos_to_fanout, to_sinks[0], 0, 0),
    // g.edge(depos_to_fanout, to_sinks[1], 0, 0),
    // g.edge(depos_to_fanout, to_sinks[2], 0, 0),
    // g.edge(depos_to_fanout, to_sinks[3], 0, 0),
  ]
);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(new_graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph",
                  "WireCellSio", "WireCellSigProc",
                  "WireCellSpng"],
        apps: ["Pgrapher"]
    }
};

// Finally, the configuration sequence which is emitted.
[cmdline] + g.uses(new_graph) + [app]