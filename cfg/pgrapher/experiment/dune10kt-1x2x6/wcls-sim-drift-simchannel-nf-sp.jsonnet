// This is a main entry point for configuring a wire-cell CLI job to
// simulate protoDUNE-SP.  It is simplest signal-only simulation with
// one set of nominal field response function.  It excludes noise.
// The kinematics are a mixture of Ar39 "blips" and some ideal,
// straight-line MIP tracks.
//
// Output is a Python numpy .npz file.

local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/experiment/dune-vd/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';
local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params_maker = import 'pgrapher/experiment/dune10kt-1x2x6/simparams.jsonnet';


local fcl_params = {
    G4RefTime: std.extVar('G4RefTime') * wc.us,
    apa_sparsity: std.extVar('apa_sparsity'),
    use_dnnroi: std.extVar('use_dnnroi'),
    dnn_model: std.extVar('dnn_model'), // "ts-model/CP49_mobilenetv3.ts",
    inference_service: "TorchService", // "TorchService" or "TritonService"
    nchunks: 1,
    dnn_output_sparcify: true,
};
local params = params_maker(fcl_params) {
  lar: super.lar {
    // Longitudinal diffusion constant
    DL: std.extVar('DL') * wc.cm2 / wc.ns,
    // Transverse diffusion constant
    DT: std.extVar('DT') * wc.cm2 / wc.ns,
    // Electron lifetime
    lifetime: std.extVar('lifetime') * wc.us,
    // Electron drift speed, assumes a certain applied E-field
    drift_speed: std.extVar('driftSpeed') * wc.mm / wc.us,
  },
};


local tools_all = tools_maker(params);
local tools = 
if fcl_params.apa_sparsity == "apa0" then
tools_all {
  anodes: [tools_all.anodes[0]]
}
else
  tools_all;

local sim_maker = import 'pgrapher/experiment/dune10kt-1x2x6/sim.jsonnet';
local sim = sim_maker(params, tools);

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes - 1);


local inputTag = 'IonAndScint';
local output = 'wct-sim-ideal-sig.npz';


//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
// local depos = sim.tracks(tracklist, step=1.0 * wc.mm);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);
local wcls_input = {
  depos: wcls.input.depos(name='', art_tag=inputTag),
  deposet: g.pnode({
    type: 'wclsSimDepoSetSource',
    name: "",
    data: {
      model: "",
      scale: -1, //scale is -1 to correct a sign error in the SimDepoSource converter.
      art_tag: inputTag, //name of upstream art producer of depos "label:instance:processName"
      id_is_track: false,
      assn_art_tag: "",
    },
  }, nin=0, nout=1),
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
local wcls_output = {
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
      // chanmaskmaps: ['bad'],
      pedestal_mean: 'native',
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  // The noise filtered "ADC" values.  These are truncated for
  // art::Event but left as floats for the WCT SP.  Note, the tag
  // "raw" is somewhat historical as the output is not equivalent to
  // "raw data".
  nf_digits: wcls.output.digits(name='nfdigits', tags=['raw']),

  // this wcls.output.signals one only use one APA?
  // sp_signals: wcls.output.signals(name='spsignals', tags=['gauss', 'wiener']),
  sp_signals: g.pnode({
  type: 'wclsFrameSaver',
  name: 'spsignals',
  data: {
    anode: wc.tn(mega_anode),
    digitize: false,  // true means save as RawDigit, else recob::Wire
    frame_tags: ['gauss', 'wiener','dnnsp'],
    frame_scale: [0.005, 0.005, 0.005],
    chanmaskmaps: [],
    nticks: params.daq.nticks,
  },
  }, nin=1, nout=1, uses=[mega_anode]),

  // save "threshold" from normal decon for each channel noise
  // used in imaging
  sp_thresholds: wcls.output.thresholds(name='spthresholds', tags=['threshold']),
};

//local deposio = io.numpy.depos(output);
local drifter = sim.drifter;
local setdrifter = g.pnode({
            type: 'DepoSetDrifter',
            data: {
                drifter: "Drifter"
            }
        }, nin=1, nout=1,
        uses=[drifter]);

local bagger = sim.make_bagger();
// local bagger = g.pnode({
//   type: 'DepoBagger',
//   name: 'bagger',
//   data: {
//     gate: [-250 * wc.us, 2750 * wc.us],  // fixed
//   },
// }, nin=1, nout=1);

// signal plus noise pipelines
//local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

local perfect = import 'pgrapher/experiment/dune10kt-1x2x6/chndb-perfect.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n),
  uses: [tools.anodes[n], tools.field],  // pnode extension
} for n in anode_iota];

//local chndb_maker = import 'pgrapher/experiment/pdsp/chndb.jsonnet';
//local noise_epoch = "perfect";
//local noise_epoch = "after";
//local chndb_pipes = [chndb_maker(params, tools.anodes[n], tools.fields[n]).wct(noise_epoch)
//                for n in std.range(0, std.length(tools.anodes)-1)];
local nf_maker = import 'pgrapher/experiment/dune10kt-1x2x6/nf.jsonnet';
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb_pipes[n]) for n in std.range(0, std.length(tools.anodes)-1)];
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

// local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"
local sp_maker = import 'pgrapher/experiment/dune10kt-1x2x6/sp.jsonnet';
local sp_override = if fcl_params.use_dnnroi then
{
    sparse: false,
    use_roi_debug_mode: true,
    m_save_negative_charge: false,
    use_multi_plane_protection: true,
    mp_tick_resolution: 10,
    tight_lf_tag: "",
    cleanup_roi_tag: "",
    break_roi_loop1_tag: "",
    break_roi_loop2_tag: "",
    shrink_roi_tag: "",
    extend_roi_tag: "",
}
else
{
    sparse: false,
    use_roi_debug_mode: false,
    m_save_negative_charge: false,
    use_multi_plane_protection: false,
    mp_tick_resolution: 10,
};
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local dnnroi = import 'pgrapher/experiment/dune-vd/dnnroi.jsonnet';
local ts =
if fcl_params.inference_service == "TorchService" then
{
    type: "TorchService",
    name: "dnnroi",
    data: {
        model: fcl_params.dnn_model,
        device: "cpu",
        concurrency: 1,
    },
}
else if fcl_params.inference_service == "TritonService" then
{
    type: "TritonService",
    name: "dnnroi",
    data: {
        url: "ailab01.fnal.gov:8001",
        model: fcl_params.dnn_model, // "dnn"
    },
}
else error("unsupported inference_service: " + fcl_params.inference_service);

local rng = tools.random;
local wcls_simchannel_sink = g.pnode({
  type: 'wclsSimChannelSink',
  name: 'postdrift',
  data: {
    artlabel: 'simpleSC',  // where to save in art::Event
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
    rng: wc.tn(rng),
    tick: 0.5 * wc.us,
    start_time: -0.25 * wc.ms,
    readout_time: self.tick * 6000,
    nsigma: 3.0,
    drift_speed: params.lar.drift_speed,
    u_to_rp: 100 * wc.mm,  // 90.58 * wc.mm,
    v_to_rp: 100 * wc.mm,  // 95.29 * wc.mm,
    y_to_rp: 100 * wc.mm,
    u_time_offset: 0.0 * wc.us,
    v_time_offset: 0.0 * wc.us,
    y_time_offset: 0.0 * wc.us,
    g4_ref_time: fcl_params.G4RefTime, // -250 * wc.us,
    use_energy: true,
  },
}, nin=1, nout=1, uses=tools.anodes);

local wcls_depoflux_writer = g.pnode({
  type: 'wclsDepoFluxWriter',
  name: 'postdrift',
  data: {
    anodes: [wc.tn(anode) for anode in tools.anodes],
    field_response: wc.tn(tools.field),
    tick: 0.5 * wc.us,
    window_start: 0.0 * wc.ms,
    window_duration: self.tick * 6000,
    nsigma: 3.0,
    reference_time: fcl_params.G4RefTime,
    energy: 1, # equivalent to use_energy = true
    simchan_label: 'simpleSC',
    sed_label: 'IonAndScint',
    sparse: false,
  },
}, nin=1, nout=1, uses=tools.anodes + [tools.field]);

local magoutput = 'protodune-data-check.root';
local magnify = import 'pgrapher/experiment/dune-vd/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);

local multipass = [
  g.pipeline([
               sn_pipes[n],
               // sinks.orig_pipe[n],
               // nf_pipes[n],
               sp_pipes[n],
               // sinks.decon_pipe[n],
             ] + if fcl_params.use_dnnroi then [
               dnnroi(
                tools.anodes[n],
                ts,
                output_scale=1.0,
                nchunks=fcl_params.nchunks,
                sparcify=fcl_params.dnn_output_sparcify,),
                // sinks.dnnroi_pipe[n],
             ] else [],
             'multipass%d' % n)
  for n in anode_iota
];
local outtags = [];
local tag_rules = {
    frame: {
        '.*': 'framefanin',
    },
    trace: {['gauss%d' % anode.data.ident]: ['gauss%d' % anode.data.ident] for anode in tools.anodes}
        + {['wiener%d' % anode.data.ident]: ['wiener%d' % anode.data.ident] for anode in tools.anodes}
        + {['threshold%d' % anode.data.ident]: ['threshold%d' % anode.data.ident] for anode in tools.anodes}
        + {['dnnsp%d' % anode.data.ident]: ['dnnsp%d' % anode.data.ident] for anode in tools.anodes},
};


local make_switch_pipe = function(d2f, anode ) {
    local ds_filter = g.pnode({
        type: "DepoSetFilter",
        name: "ds-filter-switch-%d" % anode.data.ident,
        data: {anode: wc.tn(anode)},
        }, nin=1, nout=1, uses=[anode]),
    local dorb = g.pnode({
        type: "DeposOrBust",
        name: "dorb-switch-%d" % anode.data.ident,
        }, nin=1, nout=2),
    local frame_sync = g.pnode({
        type: "FrameSync",
        name: "frame-sync-switch-%d" % anode.data.ident,
        }, nin=2, nout=1),
    ret1: g.intern(
        innodes=[ds_filter],
        outnodes=[frame_sync],
        centernodes=[dorb, d2f],
        edges=
            [g.edge(ds_filter, dorb, 0, 0),
            g.edge(dorb, d2f, 0, 0),
            g.edge(d2f, frame_sync, 0, 0),
            g.edge(dorb, frame_sync, 1, 1)]),
    ret2: g.pipeline([ds_filter, d2f]),
}.ret1;

local switch_pipes = [
    make_switch_pipe(multipass[n], tools.anodes[n]),
    for n in std.range(0, std.length(tools.anodes) - 1)
];

// local bi_manifold = 
local bi_manifold =
if fcl_params.apa_sparsity == "apa0" then
    f.multifanpipe('DepoSetFanout', multipass, 'FrameFanin', [1,1], [1,1], [1,1], [1,1], 'sn_mag_nf', outtags, tag_rules)
else if fcl_params.apa_sparsity == "hydra" then
    f.multifanpipe('DepoSetFanout', switch_pipes, 'FrameFanin', [1,2], [2,6], [1,2], [2,6], 'sn_mag_nf', outtags, tag_rules)
else if fcl_params.apa_sparsity == "all" then
    f.multifanpipe('DepoSetFanout', multipass, 'FrameFanin', [1,2], [2,6], [1,2], [2,6], 'sn_mag_nf', outtags, tag_rules)
else
    error("unsupported apa_sparsity: " + fcl_params.apa_sparsity);

local retagger = g.pnode({
  type: 'Retagger',
  data: {
    // Note: retagger keeps tag_rules an array to be like frame fanin/fanout.
    tag_rules: [{
      // Retagger also handles "frame" and "trace" like fanin/fanout
      // merge separately all traces like gaussN to gauss.
      frame: {
        '.*': 'retagger',
      },
      merge: {
        'gauss\\d+': 'gauss',
        'wiener\\d+': 'wiener',
        'dnnsp\\d+': 'dnnsp',
      },
    }],
  },
}, nin=1, nout=1);

//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;

local graph = g.pipeline([wcls_input.deposet, setdrifter, wcls_depoflux_writer, bi_manifold, retagger, wcls_output.sp_signals, sink]);

local app = {
  type: std.extVar('engine'), //Pgrapher, TbbFlow
  data: {
    edges: g.edges(graph),
  },
};


// Finally, the configuration sequence which is emitted.

g.uses(graph) + [app]
