local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local util = import 'pgrapher/experiment/protodunevd/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local base = import 'pgrapher/experiment/protodunevd/simparams.jsonnet';
local params = base {
  daq: super.daq {
    nticks: std.extVar('nticks'),
  },
  lar: super.lar {
    // Longitudinal diffusion constant
    DL: std.extVar('DL') * wc.cm2 / wc.ns,
    // Transverse diffusion constant
    DT: std.extVar('DT') * wc.cm2 / wc.ns,
    // Electron lifetime
    lifetime: std.extVar('lifetime') * wc.us,
    // Electron drift speed
    // drift_speed: std.extVar('driftSpeed') * wc.mm / wc.us,
    drift_speed: util.drift_velocity(std.extVar('efield'), std.extVar('temperature')) * wc.mm / wc.us,
  },
};


local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/protodunevd/sim.jsonnet';
local sim = sim_maker(params, tools);

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes - 1);


local output = 'wct-sim-ideal-sig.npz';
local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);
local wcls_input = {
  depos: wcls.input.depos(name='electron', art_tag='IonAndScint'),
  // depos: wcls.input.depos(name='electron'),  // default art_tag="blopper"
};

// Collect all the wc/ls output converters for use below.  Note the
// "name" MUST match what is used in theh "outputers" parameter in the
// FHiCL that loads this file.
local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};

local resolution = params.adc.resolution;
local fullscale = params.adc.fullscale[1] - params.adc.fullscale[0];
local ADC_mV_ratio = ((1 << resolution) - 1 ) / fullscale;

local wcls_output = {
  // signal waveform from simulation
  sim_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'simsignals',
    data: {
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['sig'],
      frame_scale: [ADC_mV_ratio], // convert mV to ADC
      nticks: -1,
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  // ADC output from simulation
  // sim_digits: wcls.output.digits(name="simdigits", tags=["orig"]),
  sim_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'simdigits',
    data: {
      // anode: wc.tn(tools.anode),
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['daq'],
      // nticks: params.daq.nticks,
      nticks: -1,
      // chanmaskmaps: ['bad'],
      pedestal_mean: 'native',
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  // The noise filtered "ADC" values.  These are truncated for
  // art::Event but left as floats for the WCT SP.  Note, the tag
  // "raw" is somewhat historical as the output is not equivalent to
  // "raw data".
  nf_digits: wcls.output.digits(name='nfdigits', tags=['raw']),

  // The output of signal processing.  Note, there are two signal
  // sets each created with its own filter.  The "gauss" one is best
  // for charge reconstruction, the "wiener" is best for S/N
  // separation.  Both are used in downstream WC code.
  sp_signals: wcls.output.signals(name='spsignals', tags=['gauss', 'wiener']),

  // save "threshold" from normal decon for each channel noise
  // used in imaging
  sp_thresholds: wcls.output.thresholds(name='spthresholds', tags=['threshold']),
};

//local deposio = io.numpy.depos(output);
local drifter = sim.drifter;
local bagger = sim.make_bagger();

// signal plus noise pipelines
local signal_pipes = sim.analog_pipelines;
// local signal_pipes = sim.signal_pipelines;
// local sn_pipes = sim.splusn_pipelines;

local perfect = import 'pgrapher/experiment/protodunevd/chndb-perfect.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n){dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];

local nf_maker = import 'pgrapher/experiment/protodunevd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local sp_maker = import 'pgrapher/experiment/protodunevd/sp.jsonnet';
local sp = sp_maker(params, tools);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local rng = tools.random;
local wcls_simchannel_sink = g.pnode({
  type: 'wclsSimChannelSink',
  name: 'postdrift',
  data: {
    artlabel: 'simpleSC',  // where to save in art::Event
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
    rng: wc.tn(rng),
    tick: params.daq.tick,
    response_plane: params.det.response_plane,
    start_time: -0.25 * wc.ms,
    readout_time: params.daq.readout_time,
    nsigma: 3.0,
    drift_speed: params.lar.drift_speed,
    u_to_rp: params.det.response_plane,
    v_to_rp: params.det.response_plane,
    y_to_rp: params.det.response_plane,
    u_time_offset: 0.0 * wc.us,
    v_time_offset: 0.0 * wc.us,
    y_time_offset: 0.0 * wc.us,
    g4_ref_time: -250 * wc.us,
    use_energy: true,
  },
}, nin=1, nout=1, uses=tools.anodes);

local magnify = import 'pgrapher/experiment/protodunevd/magnify-sinks.jsonnet';
local magnifyio = magnify(tools, "protodunevd-sim-check.root");

local multipass = [
  g.pipeline([
               signal_pipes[n],
               // magnifyio.orig_pipe[n],
             ],
             'multipass%d' % n)
  for n in anode_iota
];
local outtags = ['orig%d' % n for n in anode_iota];
local bi_manifold = f.fanpipe('DepoSetFanout', multipass, 'FrameFanin', 'sn_mag_nf', outtags);
// local bi_manifold = f.fanpipe('DepoFanout', multipass, 'FrameFanin', 'sn_mag_nf', outtags);

local retagger = g.pnode({
  type: 'Retagger',
  data: {
    // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
    tag_rules: [{
      // Retagger also handles "frame" and "trace" like fanin/fanout
      // merge separately all traces like gaussN to gauss.
      frame: {
        '.*': 'orig',
      },
      merge: {
        // 'orig\\d+': 'daq',
        'orig\\d': 'sig',
      },
    }],
  },
}, nin=1, nout=1);

local make_noise_model = function(anode, csdb=null) {
    type: "EmpiricalNoiseModel",
    name: "empericalnoise%s"% anode.name,
    data: {
        anode: wc.tn(anode),
        dft: wc.tn(tools.dft),
        chanstat: if std.type(csdb) == "null" then "" else wc.tn(csdb),
        spectra_file: params.files.noises[ if anode.data.ident<4 then 0 else 1 ],
        nsamples: params.daq.nticks,
        period: params.daq.tick,
        wire_length_scale: 1.0*wc.cm, // optimization binning
    },
    uses: [anode, tools.dft] + if std.type(csdb) == "null" then [] else [csdb],
};
local noise_models = [make_noise_model(anode) for anode in tools.anodes];

local add_noise = function(model) g.pnode({
    type: "AddNoise",
    name: "addnoise%s"%[model.name],
    data: {
        rng: wc.tn(tools.random),
        dft: wc.tn(tools.dft),
        model: wc.tn(model),
        nsamples: params.daq.nticks,
        replacement_percentage: 1.0, // random optimization
    }}, nin=1, nout=1, uses=[tools.random, tools.dft, model]);
local noises = [add_noise(model) for model in noise_models];

local digitizer(anode, name="", tag="") = g.pnode({
  type: "Digitizer",
  name: name,
  data : params.adc {
      anode: wc.tn(anode),
      frame_tag: tag,
      baselines: if anode.data.ident < 4
                 then [1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt]
                 else [1.0*wc.volt, 1.0*wc.volt, 1.0*wc.volt],
      fullscale: if anode.data.ident < 4
                 then [0.2*wc.volt, 1.6*wc.volt]
                 else [0.0*wc.volt, 2.0*wc.volt],
  }
}, nin=1, nout=1, uses=[anode]);
local digitizers_noise = [
    digitizer(tools.anodes[n], name="digitizer-noise%d"%n, tag="splusn%d"%n)
    for n in std.range(0,nanodes-1)];

local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: util.anode_channels(n),
      //tags: ['orig%d' % n], // traces tag
      tags: ['sig'],  // traces tag
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local multipass_noise = [
  g.pipeline([
               chsel_pipes[n],
               noises[n],
               digitizers_noise[n],
               // magnifyio.orig_pipe[n],
             ],
             'multipass-noise%d' % n)
  for n in anode_iota
];
local outtags_noise = ['splusn%d' % n for n in anode_iota];
local bi_manifold_noise = f.fanpipe('FrameFanout', multipass_noise, 'FrameFanin', 'noisedigit', outtags_noise);
local retagger_noise = g.pnode({
  type: 'Retagger',
  name: 'retagger-noise',
  data: {
    // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
    tag_rules: [{
      // Retagger also handles "frame" and "trace" like fanin/fanout
      // merge separately all traces like gaussN to gauss.
      frame: {
        '.*': 'splusn',
      },
      merge: {
        'splusn\\d': 'daq',
      },
    }],
  },
}, nin=1, nout=1);


local sink = sim.frame_sink;

local graph = g.pipeline([wcls_input.depos, drifter, wcls_simchannel_sink, bagger, bi_manifold, retagger, wcls_output.sim_signals, bi_manifold_noise, retagger_noise, wcls_output.sim_digits, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};


// Finally, the configuration sequence which is emitted.

g.uses(graph) + [app]
