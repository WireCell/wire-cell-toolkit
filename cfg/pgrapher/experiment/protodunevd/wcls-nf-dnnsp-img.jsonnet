// LArSoft (WC/LS) entry point for the deployed PDVD DNN-ROI + L1SP-DNN
// hybrid processing chain.  Mirrors pgrapher/experiment/pdhd/wcls-nf-dnnsp-img.jsonnet
// (PDHD round-4 deploy) with PDVD-specific defaults from
// experiments/stage_a_pu_round2_pdvd/deploy_round2.md:
//   - DNN-ROI model dnnroi/pdvd/pipe_distill_nestedunet_6ch.ts
//   - L1SP DNN model l1sp/pdvd/l1sp_dnn_pdvd_v1.ts
//   - Hybrid mode with loose-heur (gmax>=300, min_length>=5, energy_frac>=0.20)
//   - Per-CRP DNN thresholds (bottom apa<4 = 0.16, top apa>=4 = 0.46)
//   - Track-veto auto-disabled in hybrid mode (the DNN is the FP suppressor;
//     see deploy_round2.md "Deviations from PDHD")
// All 8 anodes (4 bottom-CRP + 4 top-CRP) are processed.

local reality = std.extVar('reality');
local sigoutform = std.extVar('signal_output_form');  // 'sparse' or 'dense'
local save_tradsp = true;

// Wire L1SP-after-DNN-ROI (hybrid mode) inside the per-anode pipeline.
// Baked 'true' matches the deployed PDVD chain.  Flip to false to
// revert to pre-L1SP DNN-only behaviour.
local use_l1sp_dnn = true;

local wc = import 'wirecell.jsonnet';
local f = import 'pgrapher/common/funcs.jsonnet';
local g = import 'pgraph.jsonnet';

local raw_input_label = std.extVar('raw_input_label');  // eg "daq"

local data_params = import 'params.jsonnet';
local simu_params = import 'simparams.jsonnet';
local params = if reality == 'data' then data_params else simu_params;

local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools = tools_maker(params);

local wcls_maker = import 'pgrapher/ui/wcls/nodes.jsonnet';
local wcls = wcls_maker(params, tools);

local sp_maker = import 'pgrapher/experiment/protodunevd/sp.jsonnet';

local wcls_input = {
  adc_digits: g.pnode({
    type: 'wclsRawFrameSource',
    name: '',
    data: {
      art_tag: raw_input_label,
      frame_tags: ['orig'],
    },
  }, nin=0, nout=1),
};

local mega_anode = {
  type: 'MegaAnodePlane',
  name: 'meganodes',
  data: {
    anodes_tn: [wc.tn(anode) for anode in tools.anodes],
  },
};

local wcls_output = {
  nf_digits: g.pnode({
    type: 'wclsFrameSaver',
    name: 'nfsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: true,
      frame_tags: ['raw'],
      chanmaskmaps: ['bad'],
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  // L1SP-corrected gauss/wiener under the standard 'gauss'/'wiener'
  // frame_tags.  PDVD historical frame_scale 0.005 matches the existing
  // wcls-nf-sp-img.jsonnet.
  sp_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'spsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: false,
      frame_tags: ['gauss', 'wiener'],
      frame_scale: [0.005, 0.005],
      chanmaskmaps: [],
      nticks: -1,
    },
  }, nin=1, nout=1, uses=[mega_anode]),

  dnn_signals: g.pnode({
    type: 'wclsFrameSaver',
    name: 'dnnsaver',
    data: {
      anode: wc.tn(mega_anode),
      digitize: false,
      frame_tags: ['dnnsp'],
      frame_scale: [0.005],
      nticks: -1,
    },
  }, nin=1, nout=1, uses=[mega_anode]),
};

local base = import 'chndb-base.jsonnet';
local chndb = [{
  type: 'OmniChannelNoiseDB',
  name: 'ocndbperfect%d' % n,
  data: base(params, tools.anodes[n], tools.field, tools.anodes[n].data.ident) { dft: wc.tn(tools.dft) },
  uses: [tools.anodes[n], tools.field, tools.dft],
} for n in std.range(0, std.length(tools.anodes) - 1)];

local nf_maker = import 'pgrapher/experiment/protodunevd/nf.jsonnet';
local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n],
                           tools.anodes[n].data.ident,
                           name='nf%d' % tools.anodes[n].data.ident)
                  for n in std.range(0, std.length(tools.anodes) - 1)];

local sp_override = {
  sparse: sigoutform == 'sparse',
  use_roi_refinement: true,
  use_roi_debug_mode: true,
  save_negtive_charge: false,
  tight_lf_tag: '',
  cleanup_roi_tag: '',
  break_roi_loop1_tag: '',
  break_roi_loop2_tag: '',
  shrink_roi_tag: '',
  extend_roi_tag: '',
  use_multi_plane_protection: true,
  do_not_mp_protect_traditional: true,
  mp_tick_resolution: 10,
};
local sp = sp_maker(params, tools, sp_override);
// L1SP inside sp.make_sigproc must stay OFF in the DNN-ROI chain — the
// embedded FrameMerger drops the SP debug tags (loose_lf*, mp2_roi*,
// mp3_roi*, decon_charge*) that DNN-ROI's 6-channel input model needs.
// L1SP-after-DNN is wired separately via the l1sp_after_dnnroi envelope
// when use_l1sp_dnn=true.
local sp_pipes = [sp.make_sigproc(a, l1sp_pd_mode='') for a in tools.anodes];

local img = import 'pgrapher/experiment/protodunevd/img.jsonnet';
local img_maker = img();
local img_pipes = [img_maker.per_anode(a) for a in tools.anodes];

local util = import 'pgrapher/experiment/protodunevd/funcs.jsonnet';
local chsel_pipes = [
  g.pnode({
    type: 'ChannelSelector',
    name: 'chsel%d' % tools.anodes[n].data.ident,
    data: {
      channels: util.anode_channels(tools.anodes[n].data.ident),
    },
  }, nin=1, nout=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local dnnroi = import 'pgrapher/experiment/protodunevd/dnnroi_pp.jsonnet';
local ts = {
  type: 'TorchService',
  name: 'dnnroi',
  data: {
    // Deployed 6-ch FP32 KD model (Dice 0.7816).  PDVD-bottom CRP-
    // trained; top CRP runs the same model.
    model: 'dnnroi/pdvd/pipe_distill_nestedunet_6ch.ts',
    device: 'cpu',
    concurrency: 1,
  },
};

// L1SP DNN tagger TorchService — round-2 hybrid deploy.  Loaded only
// when use_l1sp_dnn=true; null otherwise so the heuristic-only path
// doesn't pull libtorch a second time.
local l1sp_ts = if use_l1sp_dnn then {
  type: 'TorchService',
  name: 'l1sp_dnn_pdvd',
  data: {
    model: 'l1sp/pdvd/l1sp_dnn_pdvd_v1.ts',
    device: 'cpu',
    concurrency: 1,
  },
} else null;

local l1sp_dnn_maker = import 'pgrapher/experiment/protodunevd/l1sp_after_dnnroi.jsonnet';

// Per-anode DNN-ROI inner subgraph, bound so the L1SP envelope can take
// it as a positional argument (alongside the SP pipe).
local dnnroi_inner = [
  dnnroi(tools.anodes[n], ts, output_scale=1.0,
         nticks=params.daq.nticks, nchunks=1)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

// PDVD round-2 hybrid-mode deploy.  Per-CRP DNN thresholds unified at
// 0.5 (bottom and top) and PDVD-tuned loose-heur (300, 5, 0.20) from
// experiments/stage_a_pu_round2_pdvd/deploy_round2.md.  Bottom
// threshold raised 0.16 → 0.35 → 0.5; top threshold raised 0.46 → 0.5
// on 2026-05-25 to align with the PDHD round-4 cleanup and to push
// the operating point further along precision (see doc 13).  The
// envelope auto-disables PDVD's track-veto when
// l1sp_pd_dump_mode=='hybrid' (see deploy_round2.md "Deviations
// from PDHD").  Also enables DNN veto on adjacency-promoted ROIs.
local l1sp_envelope = if use_l1sp_dnn then [
  l1sp_dnn_maker(tools.anodes[n], sp_pipes[n], dnnroi_inner[n],
                 tools, params,
                 l1sp_pd_dump_mode='hybrid',
                 l1sp_pd_torch_service=l1sp_ts,
                 l1sp_pd_dnn_threshold_bottom=0.5,
                 l1sp_pd_dnn_threshold_top=0.5,
                 l1sp_pd_adj_dnn_veto=true,
                 l1sp_pd_gmax_min=300.0,
                 l1sp_pd_min_length=5,
                 l1sp_pd_energy_frac_thr=0.20)
  for n in std.range(0, std.length(tools.anodes) - 1)
] else [];

local magoutput = 'protodune-data-check.root';
local magnify = import 'pgrapher/experiment/protodunevd/magnify-sinks.jsonnet';
local magio = magnify(tools, magoutput);

local use_magnify = std.extVar('use_magnify');

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
    [chsel_pipes[n]] + body(n),
    'nfsp_pipe_%d' % n)
  for n in std.range(0, std.length(tools.anodes) - 1)
];

local anode_ident = [tools.anodes[n].data.ident for n in std.range(0, std.length(tools.anodes) - 1)];
local fanin_tag_rules = [
  {
    frame: { '.*': 'framefanin' },
    trace: {
      ['gauss%d' % ind]: 'gauss%d' % ind,
      ['wiener%d' % ind]: 'wiener%d' % ind,
      ['threshold%d' % ind]: 'threshold%d' % ind,
      ['loose_lf%d' % ind]: 'loose_lf%d' % ind,
    },
  }
  for ind in anode_ident
];

local nanodes = std.length(tools.anodes);
local fanpipe = f.multifanout('FrameFanout', nfsp_pipes, [1, nanodes], [nanodes, 1], 'sn_mag', fanin_tag_rules);

local retagger = g.pnode({
  type: 'Retagger',
  name: 'dnnout',
  data: {
    tag_rules: [{
      frame: { '.*': 'dnnretagger' },
      merge: { 'dnnsp\\d': 'dnnsp' },
    }],
  },
}, nin=1, nout=1);

local sink = g.pnode({ type: 'DumpFrames' }, nin=1, nout=0);
local graph = g.pipeline([wcls_input.adc_digits, fanpipe], retagger);

// Subgraph for saving the trad SP frame.  Spliced in below at the
// OmnibusSigProc output edge so the SP gauss/wiener (pre-DNN, pre-L1SP)
// is preserved alongside the L1SP-corrected output.
local ofanin = g.pnode({
  type: 'FrameFanin',
  name: 'outfanin',
  data: {
    multiplicity: std.length(tools.anodes),
    tag_rules: [
      {
        frame: { '.*': 'outfanin' },
        trace: {
          ['gauss%d' % n]: ['gauss%d' % n],
          ['wiener%d' % n]: ['wiener%d' % n],
        },
      }
      for n in std.range(0, std.length(tools.anodes) - 1)
    ],
  },
}, nin=std.length(tools.anodes), nout=1);

local osink = g.pnode({ type: 'DumpFrames', name: 'outsink', data: {} }, nin=1, nout=0);
local outretagger = g.pnode({
  type: 'Retagger',
  name: 'spout',
  data: {
    tag_rules: [{
      frame: { '.*': 'spretagger' },
      merge: {
        'gauss\\d': 'gauss',
        'wiener\\d': 'wiener',
      },
    }],
  },
}, nin=1, nout=1);
local outgr = g.pipeline([ofanin, outretagger, wcls_output.sp_signals, osink]);

local edge_selector(e) = std.startsWith(e.tail.node, 'OmnibusSigProc:');
local fanout_factory(n, e) = { type: 'FrameFanout', name: 'splice%d' % n, data: { multiplicity: 2 } };

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

g.uses(spliced_graph) + [app]
