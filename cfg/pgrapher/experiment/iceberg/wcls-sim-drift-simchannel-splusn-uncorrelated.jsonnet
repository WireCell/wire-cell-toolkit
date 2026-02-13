// ============================================================================
// wcls-sim-drift-simchannel-splusn-uncorrelated.jsonnet
// ============================================================================
//
// Assumptions about the stored noise model (WCT base units):
//   - freq_ghz : GHz  (== 1/ns in WCT base time units)
//   - avg_mag  : MV   (WCT base voltage)
//
// UncorrelatedAddNoise injects noise in MV into frames in MV.
//
// ============================================================================

local g  = import 'pgraph.jsonnet';
local f  = import 'pgrapher/common/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io           = import 'pgrapher/common/fileio.jsonnet';
local tools_maker  = import 'pgrapher/common/tools.jsonnet';
local params_maker = import 'pgrapher/experiment/iceberg/simparams.jsonnet';

// ----------------------------------------------------------------------------
// (1) FHiCL-provided externals -> params/tools/sim
// ----------------------------------------------------------------------------
local fcl_params = {
  G4RefTime: std.extVar('G4RefTime') * wc.us,
};

local params = params_maker(fcl_params) {
  lar: super.lar {
    DL: std.extVar('DL') * wc.cm2 / wc.s,
    DT: std.extVar('DT') * wc.cm2 / wc.s,
    lifetime: std.extVar('lifetime') * wc.ms,
    drift_speed: std.extVar('driftSpeed') * wc.mm / wc.us,
  },
};

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/iceberg/sim.jsonnet';
local sim = sim_maker(params, tools);

// Anode bookkeeping
local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes - 1);

// ----------------------------------------------------------------------------
// (2) Noise model path (hard-coded)
// ----------------------------------------------------------------------------
local noise_model = "/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/Iceberg/noise_simulation/ICEBERG_Uncorrelated_noise_model.json.bz2";

// ----------------------------------------------------------------------------
// (3) WCLS input + MegaAnode plane definition
// ----------------------------------------------------------------------------
local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);

local wcls_input = {
  depos: wcls.input.depos(
    name='',
    art_tag='largeant:LArG4DetectorServicevolTPCActive'
  ),
};

local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};

// ----------------------------------------------------------------------------
// (4) Digitization scale (ADC per WCT-voltage-unit)
// ----------------------------------------------------------------------------
// params.adc.fullscale is in WCT voltage units (MV base).
local resolution = params.adc.resolution;
local fullscale  = params.adc.fullscale[1] - params.adc.fullscale[0];
local ADC_per_Vwct = ((1 << resolution) - 1) / fullscale;

// ----------------------------------------------------------------------------
// (5) Outputs (frame savers + optional SP outputs)
// ----------------------------------------------------------------------------
local wcls_output = {

  // Signal waveform from simulation (converted to ADC using frame_scale)
  sim_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'simsignals',
    data: {
      anode: wc.tn(mega_anode),
      digitize: true,
      frame_tags: ['sig'],
      frame_scale: [ADC_per_Vwct],
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  // Digitized output ("daq")
  sim_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'simdigits',
    data: {
      anode: wc.tn(mega_anode),
      digitize: true,
      frame_tags: ['daq'],
      pedestal_mean: 'native',
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  // Optional outputs (kept as-is)
  nf_digits: wcls.output.digits(name='nfdigits', tags=['raw']),
  sp_signals: wcls.output.signals(name='spsignals', tags=['gauss', 'wiener']),
  sp_thresholds: wcls.output.thresholds(name='spthresholds', tags=['threshold']),
};

// ----------------------------------------------------------------------------
// (6) Core sim chain (drift -> simchannel sink -> bagger -> signal pipelines)
// ----------------------------------------------------------------------------
local drifter = sim.drifter;
local bagger  = sim.make_bagger();
local signal_pipes = sim.signal_pipelines;

local rng = tools.random;

// SimChannel sink node
local wcls_simchannel_sink = g.pnode({
  type: 'wclsSimChannelSink',
  name: 'postdrift',
  data: {
    artlabel: 'simpleSC',
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
    rng: wc.tn(rng),
    tick: 0.5 * wc.us,
    start_time: -0.25 * wc.ms,
    readout_time: self.tick * 6000,
    nsigma: 3.0,
    drift_speed: params.lar.drift_speed,
    u_to_rp: 100 * wc.mm,
    v_to_rp: 100 * wc.mm,
    y_to_rp: 100 * wc.mm,
    u_time_offset: 0.0 * wc.us,
    v_time_offset: 0.0 * wc.us,
    y_time_offset: 0.0 * wc.us,
    g4_ref_time: fcl_params.G4RefTime,
    use_energy: true,
  },
}, nin=1, nout=1, uses=tools.anodes);

// Fan out across anodes and simulate signals
local multipass_signal = [
  g.pipeline([ signal_pipes[n] ], 'multipass%d' % n)
  for n in anode_iota
];

local outtags_signal = ['orig%d' % n for n in anode_iota];
local bi_manifold_signal =
  f.fanpipe('DepoSetFanout', multipass_signal, 'FrameFanin', 'sn_mag_nf', outtags_signal);

// Merge anode tags -> single "sig"
local retagger_signal = g.pnode({
  type: 'Retagger',
  name: 'retagger-signal',
  data: {
    tag_rules: [{
      frame: { '.*': 'orig' },
      merge: { 'orig\\d+': 'sig' },
    }],
  },
}, nin=1, nout=1);

// ----------------------------------------------------------------------------
// (7) Noise injection branch (UNCORRELATED)
//     select -> add noise -> digitize -> merge
// ----------------------------------------------------------------------------

// Select per-anode channel ranges from the 'sig' frame
local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: std.range(1280 * n, 1280 * (n + 1) - 1),
      tags: ['sig'],
    },
  }, nin=1, nout=1)
  for n in anode_iota
];

// One UncorrelatedAddNoise per lane
local noise_adders = [
  g.pnode({
    type: 'UncorrelatedAddNoise',
    name: 'uncorrnoise%d' % n,
    data: {
      rng: wc.tn(tools.random),
      model_file: noise_model,

      // These are in WCT base units already:
      nsamples: params.daq.nticks,
      dt: params.daq.tick,       // ns (WCT base time units)

      // Keep at 1.0 unless you need a quick normalization tweak
      ifft_scale: 1.0,
    },
  }, nin=1, nout=1, uses=[tools.random])
  for n in anode_iota
];

// Digitize noisy frames into splusnN
local digitizers_noise = [
  sim.digitizer(tools.anodes[n], name='digitizer-noise%d' % n, tag='splusn%d' % n)
  for n in anode_iota
];

// Per-lane pipeline: select -> add noise -> digitize
local multipass_noise = [
  g.pipeline([
    chsel_pipes[n],
    noise_adders[n],
    digitizers_noise[n],
  ], 'multipass-noise%d' % n)
  for n in anode_iota
];

// Fanout/fanin to collect the noisy digitized lanes
local outtags_noise = ['splusn%d' % n for n in anode_iota];
local bi_manifold_noise =
  f.fanpipe('FrameFanout', multipass_noise, 'FrameFanin', 'noisedigit', outtags_noise);

// Merge lane tags -> final "daq"
local retagger_noise = g.pnode({
  type: 'Retagger',
  name: 'retagger-noise',
  data: {
    tag_rules: [{
      frame: { '.*': 'splusn' },
      merge: { 'splusn\\d+': 'daq' },
    }],
  },
}, nin=1, nout=1);

// ----------------------------------------------------------------------------
// (8) Sink + app wrapper
// ----------------------------------------------------------------------------
local sink = sim.frame_sink;

local graph = g.pipeline([
  // inputs + drift
  wcls_input.depos,
  drifter,
  wcls_simchannel_sink,
  bagger,

  // signal path
  bi_manifold_signal,
  retagger_signal,
  wcls_output.sim_signals,

  // noise+digitize path
  bi_manifold_noise,
  retagger_noise,
  wcls_output.sim_digits,

  // final sink
  sink
]);

local app = {
  type: 'Pgrapher',
  data: { edges: g.edges(graph) },
};

g.uses(graph) + [app]