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


//local depos = g.join_sources(g.pnode({type:"DepoMerger", name:"BlipTrackJoiner"}, nin=2, nout=1),
//                             [sim.ar39(), sim.tracks(tracklist)]);
// local depos = sim.tracks(tracklist, step=1.0 * wc.mm);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);
local wcls_input = {
  // depos: wcls.input.depos(name="", art_tag="IonAndScint"),
  // depos: wcls.input.depos(name='electron', art_tag='IonAndScint'),  // default art_tag="blopper"
  deposet: g.pnode({
    type: 'wclsSimDepoSetSource',
    name: "",
    data: {
      model: "",
      scale: -1, //scale is -1 to correct a sign error in the SimDepoSource converter.
      art_tag: "IonAndScint", //name of upstream art producer of depos "label:instance:processName"
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
      // Three options for nticks:
      // - If nonzero, force number of ticks in output waveforms.
      // - If zero, use whatever input data has. (default)
      // - If -1, use value as per LS's detector properties service.
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
local setdrifter = g.pnode({
            type: 'DepoSetDrifter',
            data: {
                drifter: "Drifter"
            }
        }, nin=1, nout=1, uses=[drifter]);
local bagger = sim.make_bagger();

// signal plus noise pipelines
//local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

local perfect = import 'pgrapher/experiment/protodunevd/chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n) {dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];


// local nf_maker = import 'pgrapher/experiment/protodunevd/nf.jsonnet';
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb_pipes[n]) for n in std.range(0, std.length(tools.anodes)-1)];
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local sp_maker = import 'pgrapher/experiment/protodunevd/sp.jsonnet';
local sp = sp_maker(params, tools);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local rng = tools.random;
local wcls_simchannel_sink = g.pnode({
  type: 'wclsDepoSetSimChannelSink',
  name: 'postdrift',
  data: {
    artlabel: 'simpleSC',  // where to save in art::Event
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
    rng: wc.tn(rng),
    tick: params.daq.tick,
    start_time: -0.25 * wc.ms, 
    readout_time: params.daq.readout_time,
    nsigma: 3.0,
    drift_speed: params.lar.drift_speed,
    u_to_rp: 189.5 * wc.mm,
    v_to_rp: 189.5 * wc.mm,
    y_to_rp: 189.5 * wc.mm,
    u_time_offset: 0.0 * wc.us,
    v_time_offset: 0.0 * wc.us,
    y_time_offset: 0.0 * wc.us,
    g4_ref_time: -250 * wc.us,
    use_energy: true,
  },
}, nin=1, nout=1, uses=tools.anodes);

// local magoutput = 'protodune-data-check.root';
// local magnify = import 'pgrapher/experiment/protodunevd/magnify-sinks.jsonnet';
// local magio = magnify(tools, magoutput);

local multipass = [
  g.pipeline([
               sn_pipes[n],
               // magio.orig_pipe[n],
               // nf_pipes[n],
               // sp_pipes[n],
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
      // merge separately all traces like origN to orig.
      frame: {
        '.*': 'orig',
      },
      merge: {
        'orig\\d': 'daq',
      },
    }],
  },
}, nin=1, nout=1);

//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;

local plainbagger = g.pnode({
        type:'DepoBagger',
        name:'plainbagger',
        data: {
            gate: [0, 0],
        },
    }, nin=1, nout=1);

local deposet_rotate = g.pnode({
        type:'DepoSetRotate',
        name:'deposet_rotate',
        data: {
            rotate: true,
            transpose: [1,0,2],
            scale: [-1,1,1],
        },
    }, nin=1, nout=1);
local deposet_rotate_rev = g.pnode({
        type:'DepoSetRotate',
        name:'deposet_rotate_rev',
        data: {
            rotate: true,
            transpose: [1,0,2],
            scale: [1,-1,1],
        },
    }, nin=1, nout=1);

// FIXME: need a "bagger" to gate depos
local graph = g.pipeline([wcls_input.deposet, deposet_rotate_rev, setdrifter, wcls_simchannel_sink, bi_manifold, retagger, wcls_output.sim_digits, sink]);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};


// Finally, the configuration sequence which is emitted.

g.uses(graph) + [app]
