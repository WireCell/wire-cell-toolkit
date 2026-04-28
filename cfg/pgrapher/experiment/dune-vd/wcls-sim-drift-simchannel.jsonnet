// This is a main entry point for configuring a wire-cell CLI job to
// simulate protoDUNE-SP.  It is simplest signal-only simulation with
// one set of nominal field response function.  It excludes noise.
// The kinematics are a mixture of Ar39 "blips" and some ideal,
// straight-line MIP tracks.
//
// Output is a Python numpy .npz file.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local dnnroi = import 'pgrapher/experiment/dune-vd/dnnroi.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';
local params_maker = import 'pgrapher/experiment/dune-vd/params.jsonnet';
local response_plane = std.extVar('response_plane')*wc.cm;
local fcl_params = {
    G4RefTime: std.extVar('G4RefTime') * wc.us,
    response_plane: std.extVar('response_plane')*wc.cm,
    nticks: std.extVar('nticks'),
    ncrm: std.extVar('ncrm'),
    use_dnnroi: std.extVar('use_dnnroi'),
    use_hydra: std.extVar('use_hydra'),
    save_rawdigits: std.extVar('save_rawdigits'),
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
  files: super.files {
      wires: std.extVar('files_wires'),
      fields: [ std.extVar('files_fields'), ],
      noise: std.extVar('files_noise'),
  },
};

local tools = tools_maker(params);

local sim_maker = import 'pgrapher/experiment/dune-vd/sim.jsonnet';
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
  depos: wcls.input.depos(name='', art_tag='IonAndScint'),
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
// local bagger = g.pnode({
//   type: 'DepoBagger',
//   name: 'bagger',
//   data: {
//     gate: [-250 * wc.us, 2750 * wc.us],  // fixed
//   },
// }, nin=1, nout=1);

// signal plus noise pipelines
// local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

local perfect = import 'pgrapher/experiment/dune10kt-1x2x6/chndb-perfect.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: perfect(params, tools.anodes[n], tools.field, n) {dft:wc.tn(tools.dft)},
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in anode_iota];

//local chndb_maker = import 'pgrapher/experiment/pdsp/chndb.jsonnet';
//local noise_epoch = "perfect";
//local noise_epoch = "after";
//local chndb_pipes = [chndb_maker(params, tools.anodes[n], tools.fields[n]).wct(noise_epoch)
//                for n in std.range(0, std.length(tools.anodes)-1)];
local nf_maker = import 'pgrapher/experiment/dune-vd/nf.jsonnet';
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb_pipes[n]) for n in std.range(0, std.length(tools.anodes)-1)];
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local sp_maker = import 'pgrapher/experiment/dune-vd/sp.jsonnet';
local sp_override = if fcl_params.use_dnnroi then
{
    sparse: true,
    use_roi_debug_mode: true,
    use_multi_plane_protection: true,
    process_planes: [0, 1, 2]
} else {
    sparse: true,
};
local sp = sp_maker(params, tools, sp_override);
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

local ts = {
    type: "TorchService",
    name: "dnnroi",
    data: {
        model: "ts-model/unet-l23-cosmic500-e50.ts",
        device: "gpucpu",
        concurrency: 1,
    },
};

local rng = tools.random;
local wcls_simchannel_sink = g.pnode({
  type: 'wclsSimChannelSink',
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
    u_to_rp: response_plane,  // 90.58 * wc.mm,
    v_to_rp: response_plane,  // 95.29 * wc.mm,
    y_to_rp: response_plane,
    u_time_offset: 0.0 * wc.us,
    v_time_offset: 0.0 * wc.us,
    y_time_offset: 0.0 * wc.us,
    g4_ref_time: fcl_params.G4RefTime,
    use_energy: true,
    response_plane: response_plane,
  },
}, nin=1, nout=1, uses=tools.anodes);

local magoutput = 'dune-vd-sim-check.root';
local magnify = import 'pgrapher/experiment/pdsp/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);

local full_sim_pipes = [
  g.pipeline([
                sn_pipes[n],
                // sinks.orig_pipe[n],
             ],
             'apa_dense_pipes%d' % n)
  for n in anode_iota
];

local full_sp_pipes = [
  g.pipeline([
                sp_pipes[n],
                // sinks.decon_pipe[n],
                // sinks.debug_pipe[n], // use_roi_debug_mode=true in sp.jsonnet
             ] + if fcl_params.use_dnnroi then [
                 dnnroi(tools.anodes[n], ts, output_scale=1.2),
                //  sinks.dnnroi_pipe[n],
             ] else [],
             'full_sp_pipes%d' % n)
  for n in anode_iota
];

local multipass = [
    g.pipeline([
        full_sim_pipes[n],
        full_sp_pipes[n]],
        'multipass%d' % n)
    for n in anode_iota
];


local make_switch_pipe = function(sim, sp, anode ) {
    local d2f = g.pipeline([sim, sp]),
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
    // direct pipe
    local ret1 = g.pipeline([ds_filter, d2f]),
    // hydra shortcut
    local ret2 =  g.intern(
        innodes=[ds_filter],
        outnodes=[frame_sync],
        centernodes=[dorb, d2f],
        edges=
            [g.edge(ds_filter, dorb, 0, 0),
            g.edge(dorb, d2f, 0, 0),
            g.edge(d2f, frame_sync, 0, 0),
            g.edge(dorb, frame_sync, 1, 1)]),
    // special case to tapout and sync rawdigits
    local fout_bust = g.pnode({
        type: "FrameFanout",
        name: "fout-switch-bust-%d" % anode.data.ident,
        data: {multiplicity: 2},
        }, nin=1, nout=2),
    local fout_rawdigits = g.pnode({
        type: "FrameFanout",
        name: "fout-switch-rawdigits-%d" % anode.data.ident,
        data: {multiplicity: 2},
        }, nin=1, nout=2),
    local frame_sync_rawdigits = g.pnode({
        type: "FrameSync",
        name: "frame-sync-switch-rawdigits-%d" % anode.data.ident,
        }, nin=2, nout=1),
    
    local dump_rawdigits = g.pnode(
        { type: 'DumpFrames', name:"switch-dump-rawdigits-%d" % anode.data.ident, data:{} },
        nin=1, nout=0),
    // hydra shortcut with rawdigits
    local ret3 = g.intern(
        innodes=[ds_filter],
        outnodes=[frame_sync],
        centernodes=[dorb, sim, sp, fout_bust, fout_rawdigits, frame_sync_rawdigits, dump_rawdigits],
        edges=
            [g.edge(ds_filter, dorb, 0, 0),
            g.edge(dorb, sim, 0, 0),
            g.edge(dorb, fout_bust, 1, 0),
            g.edge(sim, fout_rawdigits, 0, 0),
            g.edge(fout_rawdigits, sp, 0, 0),
            g.edge(sp, frame_sync, 0, 0),
            g.edge(fout_bust, frame_sync, 0, 1),
            g.edge(fout_rawdigits, frame_sync_rawdigits, 1, 0),
            g.edge(fout_bust, frame_sync_rawdigits, 1, 1),
            g.edge(frame_sync_rawdigits, dump_rawdigits, 0, 0),
            ]
    ),
    ret: if fcl_params.save_rawdigits then ret3 else ret2,
}.ret;

local switch_pipes = [
    make_switch_pipe(full_sim_pipes[n], full_sp_pipes[n], tools.anodes[n]),
    for n in std.range(0, std.length(tools.anodes) - 1)
];

local process_pipes = if fcl_params.use_hydra then switch_pipes else multipass;

// assert (fcl_params.ncrm == 36 || fcl_params.ncrm == 112) : "only ncrm == 36 or 112 are configured";
local f = import 'pgrapher/experiment/dune-vd/funcs.jsonnet';
local outtags = ['orig%d' % n for n in anode_iota];
// local bi_manifold = f.multifanpipe('DepoSetFanout', multipass, 'FrameFanin', 6, 'sn_mag', outtags);
local bi_manifold =
    if fcl_params.ncrm == 36
    then f.multifanpipe('DepoSetFanout', process_pipes, 'FrameFanin', [1,6], [6,6], [1,6], [6,6], 'sn_mag', outtags)
    else if fcl_params.ncrm == 48
    then f.multifanpipe('DepoSetFanout', process_pipes, 'FrameFanin', [1,8], [8,6], [1,8], [8,6], 'sn_mag', outtags)
    else if fcl_params.ncrm == 112
    then f.multifanpipe('DepoSetFanout', process_pipes, 'FrameFanin', [1,8,16], [8,2,7], [1,8,16], [8,2,7], 'sn_mag', outtags);

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
        'orig\\d+': 'daq',
      },
    }],
  },
}, nin=1, nout=1);

//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;

local graph = g.pipeline([wcls_input.depos, drifter, wcls_simchannel_sink, bagger, bi_manifold, retagger, wcls_output.sim_digits, sink]);

local app = {
    type: 'TbbFlow', //Pgrapher, TbbFlow
    data: {
        edges: g.edges(graph),
    },
};


// Finally, the configuration sequence which is emitted.

g.uses(graph) + [app]
