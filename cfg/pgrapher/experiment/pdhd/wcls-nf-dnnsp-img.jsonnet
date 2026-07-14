
local reality = std.extVar('reality');
local sigoutform = std.extVar('signal_output_form');  // eg "sparse" or "dense"
local save_tradsp = true;

// Wire L1SP-after-DNN-ROI (hybrid mode) inside the per-anode pipeline.
// Baked 'true' matches the deployed PDHD chain documented in
// experiments/stage_a_pu_round4/deploy_round4.md.  Flip to false to
// revert to the pre-L1SP DNN-only behaviour.
local use_l1sp_dnn = true;

local wc = import 'wirecell.jsonnet';
local f = import "pgrapher/common/funcs.jsonnet";
local g = import 'pgraph.jsonnet';

local raw_input_label = std.extVar('raw_input_label');  // eg "daq"


local data_params = import 'params.jsonnet';
local simu_params = import 'simparams.jsonnet';
local params = if reality == 'data' then data_params else simu_params;


local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);

//local chndb_maker = import "pgrapher/experiment/pdsp/chndb.jsonnet";

local sp_maker = import 'pgrapher/experiment/pdhd/sp.jsonnet';

local use_resampler = (reality == 'data');

// Collect the WC/LS input converters for use below.  Make sure the
// "name" argument matches what is used in the FHiCL that loads this
// file.  In particular if there is no ":" in the inputer then name
// must be the emtpy string.
local wcls_input = {
  adc_digits: g.pnode({
    type: 'wclsRawFrameSource',
    name: '',
    data: {
      art_tag: raw_input_label,
      frame_tags: ['orig'],  // this is a WCT designator
      // nticks: params.daq.nticks,
      tick: 512*wc.ns,
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
  // The noise filtered "ADC" values.  These are truncated for
  // art::Event but left as floats for the WCT SP.  Note, the tag
  // "raw" is somewhat historical as the output is not equivalent to
  // "raw data".
  nf_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'nfsaver',
    data: {
      // anode: wc.tn(tools.anode),
      anode: wc.tn(mega_anode),
      digitize: true,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['raw'],
      // nticks: params.daq.nticks,
      chanmaskmaps: ['bad'],
    },
  }, nin=1, nout=1, uses=[mega_anode]),


  // The output of signal processing.  Note, there are two signal
  // sets each created with its own filter.  The "gauss" one is best
  // for charge reconstruction, the "wiener" is best for S/N
  // separation.  Both are used in downstream WC code.
  sp_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'spsaver',
    data: {
      // anode: wc.tn(tools.anode),
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['gauss', 'wiener'],
      frame_scale: [0.001, 0.001],
      // nticks: params.daq.nticks,
      chanmaskmaps: [],
      nticks: -1,
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  dnn_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'dnnsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: false,  // true means save as RawDigit, else recob::Wire
      frame_tags: ['dnnsp'],
      frame_scale: [0.001],
      nticks: -1,

    },
  }, nin=1, nout=1, uses=[mega_anode]),
};

// local perfect = import 'chndb-perfect.jsonnet';
local base = import 'chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  // data: perfect(params, tools.anodes[n], tools.field, n) { dft:wc.tn(tools.dft) },
  data: base(params, tools.anodes[n], tools.field, n) { dft:wc.tn(tools.dft) },
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in std.range(0, std.length(tools.anodes) - 1)];

local nf_maker = import 'pgrapher/experiment/pdhd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], n, name='nf%d' % n) for n in std.range(0, std.length(tools.anodes) - 1)];

local sp_override = { // assume all tages sets in base sp.jsonnet
    sparse: sigoutform == 'sparse',
    // wiener_tag: "",
    // gauss_tag: "",
    use_roi_refinement: true,
    use_roi_debug_mode: true,
    save_negtive_charge: false, // no negative charge in gauss
    tight_lf_tag: "",
    // loose_lf_tag: "",
    cleanup_roi_tag: "",
    break_roi_loop1_tag: "",
    break_roi_loop2_tag: "",
    shrink_roi_tag: "",
    extend_roi_tag: "",
    // decon_charge_tag: "",
    use_multi_plane_protection: true,
    do_not_mp_protect_traditional: true, // do_not_mp_protect_traditional to 
                                         // make a clear ref, defualt is false
    mp_tick_resolution: 10,
};
//local sp = sp_maker(params, tools, { sparse: sigoutform == 'sparse' });
local sp = sp_maker(params, tools, sp_override);
// L1SP inside sp.make_sigproc must stay OFF in the DNN-ROI chain — the
// embedded FrameMerger drops the SP debug tags (loose_lf*, mp2_roi*,
// mp3_roi*, decon_charge*) that DNN-ROI's 6-channel input model needs.
// L1SP-after-DNN is wired separately via the l1sp_after_dnnroi envelope
// when use_l1sp_dnn=true (see below).
local sp_pipes = [sp.make_sigproc(a, l1sp_pd_mode='') for a in tools.anodes];

local img = import 'pgrapher/experiment/pdhd/img.jsonnet';
local img_maker = img({use_dnn_img: true});
local img_pipes = [img_maker.per_anode(a) for a in tools.anodes];

//local util = import 'pgrapher/experiment/pdhd/funcs.jsonnet';
local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % n,
    data: {
      channels: std.range(2560 * n, 2560 * (n + 1) - 1),
      // tags: ['orig%d' % n], // traces tag
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local dnnroi = import 'pgrapher/experiment/pdhd/dnnroi_pp.jsonnet';
local ts = {
    type: "TorchService",
    name: "dnnroi",
    data: {
        // Deployed 6-ch FP32 KD model (Dice 0.9107).  See
        // experiments/stage_a_pu_round4/deploy_round4.md.
        model: "dnnroi/pdhd/pipe_distill_transformer_6ch.ts",
        device: "cpu", // "gpucpu",
        concurrency: 1,
    },
};

// L1SP DNN tagger TorchService — round-4 hybrid deploy.  Loaded only
// when use_l1sp_dnn=true; null otherwise so the heuristic-only path
// doesn't pull libtorch a second time.
local l1sp_ts = if use_l1sp_dnn then {
    type: "TorchService",
    name: "l1sp_dnn_pdhd",
    data: {
        model: "l1sp/pdhd/l1sp_dnn_pdhd_v1.ts",
        device: "cpu",
        concurrency: 1,
    },
} else null;

local l1sp_dnn_maker = import 'pgrapher/experiment/pdhd/l1sp_after_dnnroi.jsonnet';

// Per-anode DNN-ROI inner subgraph, bound so the L1SP envelope can take
// it as a positional argument (alongside the SP pipe).
local dnnroi_inner = [
    dnnroi(tools.anodes[n], ts, output_scale=1.0,
           nticks=params.daq.nticks, nchunks=1)
    for n in std.range(0, std.length(tools.anodes) - 1)
];

// PDHD round-4 hybrid-mode deploy: loose-heur (300, 10, 0.20) +
// DNN threshold 0.5 + L1SP DNN model l1sp/pdhd/l1sp_dnn_pdhd_v1.ts.
// Numbers track experiments/stage_a_pu_round4/deploy_round4.md.
// Threshold raised 0.10 → 0.35 → 0.5 (2026-05-25); see
// l1sp_dl_tagger/docs/13-l1sp-deploy-thresh0p5.md for the analysis.
// Also enables DNN veto on adjacency-promoted ROIs.
local l1sp_envelope = if use_l1sp_dnn then [
    l1sp_dnn_maker(tools.anodes[n], sp_pipes[n], dnnroi_inner[n],
                   tools, params,
                   l1sp_pd_dump_mode='hybrid',
                   l1sp_pd_torch_service=l1sp_ts,
                   l1sp_pd_dnn_threshold=0.5,
                   l1sp_pd_adj_dnn_veto=true,
                   l1sp_pd_gmax_min=300.0,
                   l1sp_pd_min_length=10,
                   l1sp_pd_energy_frac_thr=0.20)
    for n in std.range(0, std.length(tools.anodes) - 1)
] else [];

local resamplers_config = import 'pgrapher/common/resamplers.jsonnet';
local load_resamplers = resamplers_config(g, wc, tools);
local resamplers = load_resamplers.resamplers;

local magoutput = 'protodunehd-data-check.root';
local magnify = import 'pgrapher/experiment/pdhd/magnify-sinks.jsonnet';
local magio = magnify(tools, magoutput);

local use_magnify = std.extVar("use_magnify");

// Per-anode pipeline body.  When use_l1sp_dnn=true the envelope wraps
// SP + DNN-ROI + L1SP into one subgraph that emits L1SP-corrected
// gauss%d / wiener%d / raw%d (the dnnsp%d tag is renamed inside the
// envelope, so magio.dnnsp_pipe is dropped; magio.decon_pipe after the
// envelope captures the L1SP-corrected gauss/wiener instead).
local body(n) =
  if use_l1sp_dnn then (
    if use_magnify == 'true' then [
      magio.orig_pipe[n],
      nf_pipes[n],
      magio.raw_pipe[n],
      l1sp_envelope[n],
      magio.decon_pipe[n],
      img_pipes[n],
    ]
    else [
      nf_pipes[n],
      l1sp_envelope[n],
      img_pipes[n],
    ]
  ) else (
    if use_magnify == 'true' then [
      magio.orig_pipe[n],
      nf_pipes[n],
      magio.raw_pipe[n],
      sp_pipes[n],
      magio.decon_pipe[n],
      dnnroi_inner[n],
      magio.dnnsp_pipe[n],
      img_pipes[n],
    ]
    else [
      nf_pipes[n],
      sp_pipes[n],
      dnnroi_inner[n],
      img_pipes[n],
    ]
  );

local nfsp_pipes = [
  g.pipeline(
    [chsel_pipes[n]] +
    (if use_resampler then [resamplers[n]] else []) +
    body(n),
    'nfsp_pipe_%d' % n)
    for n in std.range(0, std.length(tools.anodes) - 1)
];

// local fanpipe = util.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'sn_mag_nf');
local fanout_tag_rules = [ 
          {
            frame: {
              '.*': 'orig%d' % tools.anodes[n].data.ident,
            },
            trace: {
              // fake doing Nmult SP pipelines
              //orig: ['wiener', 'gauss'],
              //'.*': 'orig',
            },
          }
          for n in std.range(0, std.length(tools.anodes) - 1)
        ];

local anode_ident = [tools.anodes[n].data.ident for n in std.range(0, std.length(tools.anodes) - 1)];
local fanin_tag_rules = [
          {
            frame: {
              //['number%d' % n]: ['output%d' % n, 'output'],
              '.*': 'framefanin',
            },
            trace: {
              ['gauss%d'%ind]:'gauss%d'%ind,
              ['wiener%d'%ind]:'wiener%d'%ind,
              ['threshold%d'%ind]:'threshold%d'%ind,
              // ['tight_lf%d'%ind]:'tight_lf%d'%ind,
              ['loose_lf%d'%ind]:'loose_lf%d'%ind,
            },

          }
          for ind in anode_ident
        ];
// local fanpipe = util.fanpipe('FrameFanout', nfsp_pipes, 'FrameFanin', 'nfsp', [], fanout_tag_rules, fanin_tag_rules);

local nanodes = std.length(tools.anodes);
local fanpipe = f.multifanout('FrameFanout', nfsp_pipes, [1,nanodes], [nanodes,1], 'sn_mag', fanin_tag_rules);

local retagger = g.pnode({
  type: 'Retagger',
  name: 'dnnout',
  data: {
    tag_rules: [{
      frame: {'.*': 'dnnretagger',},
      merge: {'dnnsp\\d': 'dnnsp',},
    }],
  },
}, nin=1, nout=1);

local sink = g.pnode({ type: 'DumpFrames' }, nin=1, nout=0);
// local graph = g.pipeline([wcls_input.adc_digits, fanpipe, retagger, wcls_output.dnn_signals, sink]);
local graph = g.pipeline([wcls_input.adc_digits, fanpipe], retagger);

// Build an incomplete subgraph ending to be spliced for saving out frames 
local ofanin = g.pnode({ 
      type: 'FrameFanin',
      name:"outfanin",
      data:{
          multiplicity: std.length(tools.anodes),
          tag_rules: [
            {
              frame: {'.*': 'outfanin',},
              trace: {
                ['gauss%d' % n]: ['gauss%d' % n],
                ['wiener%d' % n]: ['wiener%d' % n],
                // ['threshold%d' % n]: ['threshold%d' % n],
              },
            }
            for n in std.range(0, std.length(tools.anodes) - 1)
          ],
      } 
      }, nin=std.length(tools.anodes), nout=1);
local osink = g.pnode({ type: 'DumpFrames', name:"outsink", data:{} }, nin=1, nout=0);
// local outsgr = g.intern(innodes=[ofanin,], centernodes = [osink],
//                       edges=[ g.edge(ofanin, osink) ], name="outsgr");
local outretagger = g.pnode({
  type: 'Retagger',
  name: 'spout',
  data: {
    tag_rules: [{
      frame: {'.*': 'spretagger',},
      merge: {
        'gauss\\d': 'gauss',
        'wiener\\d': 'wiener',
        // 'threshold\\d': 'threshold',
      },
    }],
  },
}, nin=1, nout=1);
local outgr = g.pipeline([ofanin, outretagger, wcls_output.sp_signals, osink]);

local edge_selector(e) = std.startsWith(e.tail.node, "OmnibusSigProc:");
local fanout_factory(n,e) = { type:'FrameFanout', name:"splice%d"%n, data:{multiplicity: 2} }; // "2-wire" splice

local spliced_graph =
  if save_tradsp then
    g.splice(graph, outgr, edge_selector, fanout_factory)
  else
    graph;

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(spliced_graph),
  },
};

// Finally, the configuration sequence
g.uses(spliced_graph) + [app]
