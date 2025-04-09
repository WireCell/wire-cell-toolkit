
local g = import 'pgraph.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local util = import 'pgrapher/experiment/pdhd/funcs.jsonnet';
local wc = import 'wirecell.jsonnet';

local io = import 'pgrapher/common/fileio.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';

local params = import 'pgrapher/experiment/pdhd/simparams.jsonnet';

// local tools = tools_maker(params);
local btools = tools_maker(params);
local tools = btools {
    anodes : [btools.anodes[0], ],
};

local sim_maker = import 'pgrapher/experiment/pdhd/sim.jsonnet';
local sim = sim_maker(params, tools);

local nanodes = std.length(tools.anodes);
local anode_iota = std.range(0, nanodes-1);

local wcls_maker = import "pgrapher/ui/wcls/nodes.jsonnet";
local wcls = wcls_maker(params, tools);
local wcls_input = {
    depos: wcls.input.depos(name="", art_tag="IonAndScint"),
};

local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};
local wcls_output = {
    // ADC output from simulation
    sim_digits: g.pnode({
      type: 'wclsFrameSaver',
      name: 'simdigits',
      data: {
        anode: wc.tn(mega_anode),
        digitize: true,  // true means save as RawDigit, else recob::Wire
        frame_tags: ['daq'],
        // nticks: params.daq.nticks,
        // chanmaskmaps: ['bad'],
      },
    }, nin=1, nout=1, uses=[mega_anode]),

    // The noise filtered "ADC" values.  These are truncated for
    // art::Event but left as floats for the WCT SP.  Note, the tag
    // "raw" is somewhat historical as the output is not equivalent to
    // "raw data".
    nf_digits: wcls.output.digits(name="nfdigits", tags=["raw"]),

    // The output of signal processing.  Note, there are two signal
    // sets each created with its own filter.  The "gauss" one is best
    // for charge reconstruction, the "wiener" is best for S/N
    // separation.  Both are used in downstream WC code.
    sp_signals: wcls.output.signals(name="spsignals", tags=["gauss", "wiener"]),

    // save "threshold" from normal decon for each channel noise
    // used in imaging
    sp_thresholds: wcls.output.thresholds(name="spthresholds", tags=["threshold"]),
};

//local deposio = io.numpy.depos(output);
local drifter = sim.drifter;
// local bagger = sim.make_bagger();
local bagger = [sim.make_bagger("bagger%d"%n) for n in anode_iota];

// signal plus noise pipelines
//local sn_pipes = sim.signal_pipelines;
local sn_pipes = sim.splusn_pipelines;

// local perfect = import 'pgrapher/experiment/pdhd/chndb-perfect.jsonnet';
local base = import 'pgrapher/experiment/pdhd/chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  // data: perfect(params, tools.anodes[n], tools.field, n),
  data: base(params, tools.anodes[n], tools.field, n),
  uses: [tools.anodes[n], tools.field],  // pnode extension
} for n in anode_iota];

//local chndb_maker = import 'pgrapher/experiment/pdhd/chndb.jsonnet';
//local noise_epoch = "perfect";
//local noise_epoch = "after";
//local chndb_pipes = [chndb_maker(params, tools.anodes[n], tools.fields[n]).wct(noise_epoch)
//                for n in std.range(0, std.length(tools.anodes)-1)];
local nf_maker = import 'pgrapher/experiment/pdhd/nf.jsonnet';
// local nf_pipes = [nf_maker(params, tools.anodes[n], chndb_pipes[n]) for n in std.range(0, std.length(tools.anodes)-1)];
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in anode_iota];

local sp_maker = import 'pgrapher/experiment/pdhd/sp.jsonnet';
local sp = sp_maker(params, tools, { sparse: true, use_roi_debug_mode: true, use_multi_plane_protection: true, mp_tick_resolution: 4, });
local sp_pipes = [sp.make_sigproc(a) for a in tools.anodes];

// local deposplats = [sim.make_ductor('splat%d'%n, tools.anodes[n], tools.pirs[0], 'DepoSplat', 'deposplat%d'%n) for n in anode_iota] ;
local deposplats = [util.splat(params, tools, tools.anodes[n]) for n in anode_iota] ;

local rng = tools.random;
local wcls_simchannel_sink = g.pnode({
    type: 'wclsSimChannelSink',
    name: 'postdrift',
    data: {
      artlabel: "simpleSC",    // where to save in art::Event
      anodes_tn: [wc.tn(anode) for anode in tools.anodes],  
      rng: wc.tn(rng),
      tick: 0.5*wc.us,
      start_time: -0.25*wc.ms,
      readout_time: self.tick*6000,
      nsigma: 3.0,
      drift_speed: params.lar.drift_speed,
      u_to_rp: 90.58*wc.mm,
      v_to_rp: 95.29*wc.mm,
      y_to_rp: 100*wc.mm,
      u_time_offset: 0.0*wc.us,
      v_time_offset: 0.0*wc.us,
      y_time_offset: 0.0*wc.us,
      use_energy: true,
    },
}, nin=1, nout=1, uses=tools.anodes);

local magoutput = 'g4-rec-0.root';
local magnify = import 'pgrapher/experiment/pdhd/magnify-sinks.jsonnet';
local sinks = magnify(tools, magoutput);

local hio_truth = [g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_truth%d' % n,
      data: {
        anode: wc.tn(tools.anodes[n]),
        trace_tags: ['deposplat%d'%n],
        filename: "g4-tru-%d.h5" % n,
        chunk: [0, 0], // ncol, nrow
        gzip: 2,
        high_throughput: true,
      },  
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
    ];

local hio_orig = [g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_orig%d' % n,
      data: {
        anode: wc.tn(tools.anodes[n]),
        trace_tags: ['orig%d'%n],
        filename: "g4-rec-%d.h5" % n,
        chunk: [0, 0], // ncol, nrow
        gzip: 2,
        high_throughput: true,
      },  
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
    ];

local hio_sp = [g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_sp%d' % n,
      data: {
        anode: wc.tn(tools.anodes[n]),
        trace_tags: ['loose_lf%d' % n
        , 'tight_lf%d' % n
        , 'cleanup_roi%d' % n
        , 'break_roi_1st%d' % n
        , 'break_roi_2nd%d' % n
        , 'shrink_roi%d' % n
        , 'extend_roi%d' % n
        , 'mp3_roi%d' % n
        , 'mp2_roi%d' % n
        , 'decon_charge%d' % n
        , 'gauss%d' % n],
        filename: "g4-rec-%d.h5" % n,
        chunk: [0, 0], // ncol, nrow
        gzip: 2,
        high_throughput: true,
      },  
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
    ];

local hio_dnn = [g.pnode({
      type: 'HDF5FrameTap',
      name: 'hio_dnn%d' % n,
      data: {
        anode: wc.tn(tools.anodes[n]),
        // trace_tags: ['dnn_sp%d' % n],
        trace_tags: ['dnnsp%d' % n],
        filename: "g4-rec-%d.h5" % n,
        chunk: [0, 0], // ncol, nrow
        gzip: 2,
        high_throughput: true,
      },  
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
    ];

local rio_orig = [g.pnode({
      type: 'ExampleROOTAna',
      name: 'rio_orig_apa%d' % n,
      data: {
        output_filename: "g4-rec-%d.root" % n,
        anode: wc.tn(tools.anodes[n]),
      },  
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
    ];

local rio_nf = [g.pnode({
      type: 'ExampleROOTAna',
      name: 'rio_nf_apa%d' % n,
      data: {
        output_filename: "g4-rec-%d.root" % n,
        anode: wc.tn(tools.anodes[n]),
      },  
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
    ];

local rio_sp = [g.pnode({
      type: 'ExampleROOTAna',
      name: 'rio_sp_apa%d' % n,
      data: {
        output_filename: "g4-rec-%d.root" % n,
        anode: wc.tn(tools.anodes[n]),
      },  
    }, nin=1, nout=1),
    for n in std.range(0, std.length(tools.anodes) - 1)
    ];

// Note: better switch to layers
local dnnroi = import 'pgrapher/experiment/pdhd/dnnroi.jsonnet';
local ts = {
    type: "TorchService",
    name: "dnnroi",
    data: {
        // model: "ts-model/unet-l23-cosmic500-e50.ts",
        model: "ts-model/CP49.ts",
        device: "cpu", // "gpucpu",
        concurrency: 1,
    },
};


local reco_fork = [
  g.pipeline([
              // wcls_simchannel_sink[n],
              bagger[n],
              sn_pipes[n],
              // hio_orig[n],
              // nf_pipes[n],
              // rio_nf[n],
              sp_pipes[n],
              hio_sp[n],

              // dnn_roi_finding[n],

              dnnroi(tools.anodes[n], ts, output_scale=1.2),

              hio_dnn[n],
              // rio_sp[n],
              g.pnode({ type: 'DumpFrames', name: 'reco_fork%d'%n }, nin=1, nout=0),
              // perapa_img_pipelines[n],
             ],
             'reco_fork%d' % n)
  for n in anode_iota
];

local truth_fork = [
  g.pipeline([
               deposplats[n],
               hio_truth[n],
               g.pnode({ type: 'DumpFrames', name: 'truth_fork%d'%n  }, nin=1, nout=0)
             ],
             'truth_fork%d' % n)
  for n in anode_iota
];

local depo_fanout = [g.pnode({
    type:'DepoFanout',
    name:'depo_fanout-%d'%n,
    data:{
        multiplicity:2,
        tags: [],
    }}, nin=1, nout=2) for n in anode_iota];
local frame_fanin = [g.pnode({
    type: 'FrameFanin',
    name: 'frame_fanin-%d'%n,
    data: {
        multiplicity: 2, 
        tags: [],
    }}, nin=2, nout=1) for n in anode_iota];

local frame_sink = g.pnode({ type: 'DumpFrames' }, nin=1, nout=0);

local multipass = [g.intern(innodes=[depo_fanout[n]], centernodes=[truth_fork[n], reco_fork[n]], outnodes=[],
                     edges = [
                       g.edge(depo_fanout[n], truth_fork[n],  0, 0),
                       g.edge(depo_fanout[n], reco_fork[n],   1, 0)]) for n in anode_iota];

// local multipass = [reco_fork[n] for n in anode_iota];

local outtags = ['orig%d' % n for n in anode_iota];
// local bi_manifold = f.fanpipe('DepoFanout', multipass, 'FrameFanin', 'sn_mag_nf', outtags);


local depo_fanout_1st = g.pnode({
    type:'DepoFanout',
    name:'depo_fanout_1st',
    data:{
        multiplicity:nanodes,
        tags: [],
    }}, nin=1, nout=nanodes);
local bi_manifold = g.intern(innodes=[depo_fanout_1st], centernodes=multipass, outnodes=[],
                      edges = [
                        g.edge(depo_fanout_1st, multipass[n],  n, 0) for n in anode_iota
                      ],
);

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
        'orig\\d': 'daq',
      },
    }],
  },
}, nin=1, nout=1);

//local frameio = io.numpy.frames(output);
local sink = sim.frame_sink;

// g4 sim as input
local graph = g.intern(innodes=[wcls_input.depos], centernodes=[drifter, depo_fanout_1st]+multipass, outnodes=[],
                      edges = 
                      [
                        g.edge(wcls_input.depos, drifter, 0, 0),
                        g.edge(drifter, depo_fanout_1st, 0, 0),
                      ] +
                      [g.edge(depo_fanout_1st, multipass[n],  n, 0) for n in anode_iota],
                      );
local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};
g.uses(graph) + [app]
