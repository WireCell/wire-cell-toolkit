# usage: wire-cell -l stdout wct-sim-check.jsonnet

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local fileio = import 'layers/high/fileio.jsonnet';
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

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes-1);
local anode_idents = [anode.data.ident for anode in tools.anodes];

local output = 'wct-sim-ideal-sigproc.npz';
local drifter = sim.drifter;
local bagger = sim.make_bagger();
local sn_pipes = sim.splusn_pipelines;

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
    use_roi_debug_mode: true,
    save_negtive_charge: true,
    use_multi_plane_protection: true,
    process_planes: [0, 1, 2]
};

// local sp_maker = import 'pgrapher/experiment/pdhd/sp.jsonnet';
local sp_maker = 
  if std.extVar('OLDPDHD') == 0 
  then import 'newsp.jsonnet'
  else import 'oldsp.jsonnet';
// local sp_maker = import 'oldsp.jsonnet';
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local magoutput = 'protodunehd-sim-check.root';
// local magnify = import 'pgrapher/experiment/pdhd/magnify-sinks.jsonnet';
local magnify = import 'magnify-sinks.jsonnet';
local magnifyio = magnify(tools, magoutput);

local outtags = ['raw%d' % n for n in std.range(0, std.length(tools.anodes) - 1)];


local frame_input = fileio.frame_tensor_file_source('tensor_frames.npz');
// local frame_input = fileio.frame_file_source('frames.tar', tags=outtags);

local parallel_pipes = [
  g.pipeline([
              //  sn_pipes[n],
               // magnifyio.orig_pipe[n],
               nf_pipes[n],
              //  magnifyio.raw_pipe[n],
               sp_pipes[n],
               // magnifyio.debug_pipe[n],
               magnifyio.decon_pipe[n],
             ],
             'parallel_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];


local simple_pipes = [
  g.pipeline([
              //  sn_pipes[n],
               // magnifyio.orig_pipe[n],
               nf_pipes[n],
              //  magnifyio.raw_pipe[n],
              //  sp_pipes[n],
               // magnifyio.debug_pipe[n],
              //  magnifyio.decon_pipe[n],
             ],
             'parallel_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

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
            ['orig%d'%n]: if std.extVar('OLDPDHD') == 1 then 'orig%d'%n else 'orig',
        },
    }
    for n in std.range(0, std.length(tools.anodes) - 1)
];

local parallel_graph = f.fanpipe('FrameFanout', parallel_pipes, 'FrameFanin', 'sn_mag_nf', outtags, fanout_apa_rules);
local fanout_graph = g.fan.fanout('FrameFanout', simple_pipes, 'sn_mag_nf', fanout_apa_rules);


local tf_maker = import 'torch.jsonnet';
local tf = tf_maker(tools);
local tf_fans = [tf.make_fanout(a) for a in tools.anodes];


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

local torch_to_tensors = [g.pnode({
  type: 'TorchToTensor',
  name: 'torchtotensor%d' % n,
  data: {},
// }, nin=1, nout=1) for n in std.range(0, std.length(tools.anodes) - 1)];
}, nin=1, nout=1) for n in std.range(0, 15)];


local nchans = [800, 800, 480];

// One FRER per field file.
local torch_frers = [[{
    type: "TorchFRERSpectrum",
    name: "torch_frer%d_plane%d" % [anode.data.ident, iplane],
    uses: [
      if std.extVar('OLDPDHD') == 0 then tools.fields[anode.data.ident] else tools.fields[0],
      tools.elec_resp
    ],
    data: {
      field_response: wc.tn(if std.extVar('OLDPDHD') == 0 then tools.fields[anode.data.ident] else tools.fields[0]),#"FieldResponse:field%d"% anode.data.ident,
      fr_plane_id: iplane,
      ADC_mV: 11702142857.142859,
      inter_gain: 1.0,
      default_nchans : nchans[iplane],
      default_nticks: 6000,
      default_period: 500.0, #512.0,
      extra_scale: 1.0,
      anode_num: anode.data.ident,
    }
  } for iplane in std.range(0,2)] for anode in tools.anodes];


local wire_filters = [
  {
    type: 'HfFilter',
    name: 'Wire_ind',
    data: {
      max_freq: 1,  // warning: units
      power: 2,
      flag: false,
      sigma: 1.0 / wc.sqrtpi * 0.75,  // caller should provide
    }
  },
  {
    type: 'HfFilter',
    name: 'Wire_col',
    data: {
      max_freq: 1,  // warning: units
      power: 2,
      flag: false,
      sigma: 1.0 / wc.sqrtpi * 10.0,  // caller should provide
    }
  },
];

local torch_wire_filters = [
  {
    type: "Torch1DSpectrum",
    name: "torch_1dspec_ind",
    uses: [wire_filters[0]],
    data: {
      spectra: [
        wc.tn(wire_filters[0]),
      ]
    },
  },
  {
    type: "Torch1DSpectrum",
    name: "torch_1dspec_col",
    uses: [wire_filters[1]],
    data: {
      spectra: [
        wc.tn(wire_filters[1]),
      ]
    },
  },
];

local torch_decons = [
  ([g.pnode({
    type: 'SPNGDecon',
    name: 'spng_decon_apa%d_plane%d' % [anode.data.ident, iplane],
    data: {
      frer_spectrum: wc.tn(torch_frers[anode.data.ident][iplane]),
      wire_filter: wc.tn(if iplane == 2 then torch_wire_filters[1] else torch_wire_filters[0]), #put in if statement
      coarse_time_offset: 1000,
    },
  },
  nin=1, nout=1,
  uses=[torch_frers[anode.data.ident][iplane], if iplane == 2 then torch_wire_filters[1] else torch_wire_filters[0]])
  for iplane in std.range(0, 2)] + [ //Duplicate the collection plane
  g.pnode({
  type: 'SPNGDecon',
  name: 'spng_decon_apa%d_plane%d_opp' % [anode.data.ident, 2],
  data: {
    frer_spectrum: wc.tn(torch_frers[anode.data.ident][2]),
    wire_filter: wc.tn(torch_wire_filters[1]),
    coarse_time_offset: 1000,
  },
  }, nin=1, nout=1, uses=[torch_frers[anode.data.ident][2], torch_wire_filters[1]])
])
 for anode in tools.anodes];

local tfs_to_tensors = g.intern(
    innodes=tf_fans,
    centernodes=std.flattenArrays(torch_decons),
    outnodes=torch_to_tensors,
    edges=[
      g.edge(tf_fans[0], torch_decons[0][0], 0),
      g.edge(tf_fans[0], torch_decons[0][1], 1),
      g.edge(tf_fans[0], torch_decons[0][2], 2),
      g.edge(tf_fans[0], torch_decons[0][3], 3),
      g.edge(tf_fans[1], torch_decons[1][0], 0),
      g.edge(tf_fans[1], torch_decons[1][1], 1),
      g.edge(tf_fans[1], torch_decons[1][2], 2),
      g.edge(tf_fans[1], torch_decons[1][3], 3),
      g.edge(tf_fans[2], torch_decons[2][0], 0),
      g.edge(tf_fans[2], torch_decons[2][1], 1),
      g.edge(tf_fans[2], torch_decons[2][2], 2),
      g.edge(tf_fans[2], torch_decons[2][3], 3),
      g.edge(tf_fans[3], torch_decons[3][0], 0),
      g.edge(tf_fans[3], torch_decons[3][1], 1),
      g.edge(tf_fans[3], torch_decons[3][2], 2),
      g.edge(tf_fans[3], torch_decons[3][3], 3),

      g.edge(torch_decons[0][0], torch_to_tensors[0]),
      g.edge(torch_decons[0][1], torch_to_tensors[1]),
      g.edge(torch_decons[0][2], torch_to_tensors[2]),
      g.edge(torch_decons[0][3], torch_to_tensors[3]),
      g.edge(torch_decons[1][0], torch_to_tensors[4]),
      g.edge(torch_decons[1][1], torch_to_tensors[5]),
      g.edge(torch_decons[1][2], torch_to_tensors[6]),
      g.edge(torch_decons[1][3], torch_to_tensors[7]),
      g.edge(torch_decons[2][0], torch_to_tensors[8]),
      g.edge(torch_decons[2][1], torch_to_tensors[9]),
      g.edge(torch_decons[2][2], torch_to_tensors[10]),
      g.edge(torch_decons[2][3], torch_to_tensors[11]),
      g.edge(torch_decons[3][0], torch_to_tensors[12]),
      g.edge(torch_decons[3][1], torch_to_tensors[13]),
      g.edge(torch_decons[3][2], torch_to_tensors[14]),
      g.edge(torch_decons[3][3], torch_to_tensors[15]),

    ],
);

local tfs_to_tensors_to_sinks = g.intern(
    innodes=[tfs_to_tensors],
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

local load_to_fanout = g.intern(
  innodes=[frame_input],
  outnodes=[fanout_graph],
  edges = [
    g.edge(frame_input, fanout_graph),
  ]
);

local sink = sim.frame_sink;


local graph = if std.extVar("SPNG") == 0 then
  g.pipeline([frame_input, parallel_graph, sink])
  else g.intern(
    innodes=[load_to_fanout],
    outnodes=[tfs_to_tensors_to_sinks],
    centernodes=[],
    edges = [
      g.edge(load_to_fanout, tfs_to_tensors_to_sinks, 0, 0),
      g.edge(load_to_fanout, tfs_to_tensors_to_sinks, 1, 1),
      g.edge(load_to_fanout, tfs_to_tensors_to_sinks, 2, 2),
      g.edge(load_to_fanout, tfs_to_tensors_to_sinks, 3, 3),
    ]
  );



local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: [
          "WireCellGen", "WireCellPgraph", "WireCellSio", 
          "WireCellSigProc", "WireCellRoot", "WireCellSpng"],
        apps: ["Pgrapher"]
    }
};


// Finally, the configuration sequence which is emitted.
[cmdline] + g.uses(graph) + [app]
