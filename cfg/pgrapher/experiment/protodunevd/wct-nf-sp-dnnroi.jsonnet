// DNN-ROI variant of wct-nf-sp.jsonnet for ProtoDUNE-VD.
//
// Inserts a per-anode DNN-ROI subgraph (DNN_ROI_SP model exported as
// wire-cell-data/dnnroi/pdvd/*.ts) after the standard SP pipeline.  The
// PDVD models are 4-channel and stacked-multi-plane: the two induction
// planes U+V are fed to DNNROIFindingMultiPlane in one call, the W
// collection plane is passed through from SP gauss.
//
// All other behaviour mirrors wct-nf-sp.jsonnet.  L1SP-after-DNN is not
// wired (PDVD has no L1SP-after-DNN config); L1SP inside sp.jsonnet stays
// OFF in the DNN chain so the DNN-ROI debug input tags survive.
//
// Run example (single anode, single event):
//   wire-cell -l stdout -L debug \
//     --tla-str orig_prefix="protodune-orig-frames" \
//     --tla-str sp_prefix="protodune-sp-dnnroi-frames" \
//     --tla-code anode_indices='[0]' \
//     --tla-str dnnroi_device="cpu" \
//     -c pgrapher/experiment/protodunevd/wct-nf-sp-dnnroi.jsonnet
//
// To skip the DNN step and behave like wct-nf-sp.jsonnet pass
//   --tla-code use_dnnroi='false'.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

local params = import 'pgrapher/experiment/protodunevd/params.jsonnet';

local tools_maker = import 'pgrapher/common/tools.jsonnet';
local tools_all = tools_maker(params);

function(
  orig_prefix   = 'protodune-orig-frames',  // input prefix; reads {prefix}-anode{N}.tar.bz2
  sp_prefix     = 'protodune-sp-dnnroi-frames',  // output prefix for post-DNN frames
  reality       = 'data',                   // 'data' enables the 512->500 ns Resampler on bottom anodes (n<4)
  sigoutform    = 'dense',                  // 'sparse' or 'dense'
  anode_indices = std.range(0, std.length(tools_all.anodes) - 1),
  use_freqmask  = true,
  debug_dump_path = '',
  debug_dump_groups = [],
  shield_dump_path = '',

  // DNN-ROI specific
  use_dnnroi    = true,
  // Default = FP32 best KD (6-ch).  Resolved via WIRECELL_PATH.
  //   FP32 best KD:    pipe_distill_nestedunet_6ch.ts        (Dice 0.7816)
  //   INT8 best QAT:   pipe_qat_nestedunet_6ch_ep0_int8.ts   (Dice 0.7797, CPU only)
  // The shell wrapper run_nf_sp_dnnroi_evt.sh -P fp32|int8 sets this via tla-str.
  dnnroi_model  = 'dnnroi/pdvd/pipe_distill_nestedunet_6ch.ts',
  // dnnroi_model = 'dnnroi/pdvd/pipe_qat_nestedunet_6ch_ep0_int8.ts',  // INT8: requires dnnroi_device='cpu'
  dnnroi_device = 'cpu',                    // 'cpu' or 'gpu'.  INT8 graph is CPU-only.
  dnnroi_nchan  = 6,                        // PDVD 6-ch deployment (DAGMan 287)
  dnnroi_concurrency = 1,
  dnnroi_nticks = 6000,
  dnnroi_tick_per_slice = 4,                // training rebin=4
  dnnroi_output_scale = 1.0,
  dnnroi_mask_thresh = 0.2,
  dnnroi_nchunks = 1,
  dnnroi_debugfile = '',                    // if non-empty, C++ node dumps per-call .pt

  // L1SP-after-DNN-ROI envelope (wraps SP + DNN-ROI + L1SP into one per-anode
  // subgraph).  When use_l1sp_dnn=false (default), the file behaves exactly as
  // before — L1SP stays off in the DNN-ROI chain.  When true, L1SP runs on
  // the DNN-relabeled gauss/wiener tags using the PDVD-tuned config in
  // pgrapher/experiment/protodunevd/l1sp_after_dnnroi.jsonnet.
  use_l1sp_dnn = false,
  l1sp_pd_mode = '',                        // 'process' | 'dump' | 'dnn' | 'hybrid' | '' (off)
  l1sp_pd_dump_path = '',                   // dump-mode NPZ dir
  l1sp_pd_wf_dump_path = '',                // process-mode per-ROI waveform dump dir
  l1sp_pd_dump_all_rois = false,            // dump every ROI, not just triggered
  l1sp_pd_adj_enable = true,
  l1sp_pd_adj_max_hops = 3,
  l1sp_pd_planes = null,                    // null -> envelope default [0,1] (U+V)
  // ── L1SP DNN tagger (mode == 'dnn' || 'hybrid').  Model file is
  //    resolved via WIRECELL_PATH and ships with wire-cell-data/l1sp/pdvd/.
  //    See experiments/stage_a_pu_round2_pdvd/deploy_round2.md for the
  //    round-2 derivation of the per-CRP thresholds.  Defaults updated
  //    2026-05-25: bottom 0.16 → 0.35 → 0.5, top 0.46 → 0.5 — both
  //    CRPs now share a single deployed threshold.  See doc 13 for
  //    motivation.  Single-threshold fallback also 0.5.
  l1sp_pd_dnn_model        = 'l1sp/pdvd/l1sp_dnn_pdvd_v1.ts',
  l1sp_pd_dnn_device       = 'cpu',
  l1sp_pd_dnn_concurrency  = 1,
  l1sp_pd_dnn_threshold    = 0.5,             // single-threshold fallback
  l1sp_pd_dnn_threshold_bottom = 0.5,         // per-CRP override (apa < 4)
  l1sp_pd_dnn_threshold_top    = 0.5,         // per-CRP override (apa >= 4)
  l1sp_pd_dnn_window_ticks = 256,
  l1sp_pd_dnn_debug_path   = '',
  // Run DNN veto on adjacency-promoted ROIs too (default true since
  // 2026-05-25).  Without it the heuristic-driven adjacency chain
  // bypasses the DNN entirely.
  l1sp_pd_adj_dnn_veto     = true,
  // ── Loose-heur pre-filter overrides (DNN-chain only).  Defaults match
  //    the PDVD-deployed C++ values.  Use the runner's --loose-heur to
  //    relax to (300, 10, 0.20) for hybrid mode; see deploy_round2.md
  //    §"Phase A".
  l1sp_pd_gmax_min         = 1500.0,
  l1sp_pd_min_length       = 30,
  l1sp_pd_energy_frac_thr  = 0.66,
  // One charge-scale constant for the TOP (TDE) electronics, SP only: passed
  // to both the top OmnibusSigProc postgain and the top L1SP gain_scale
  // (toolkit protodunevd/sp.jsonnet make_sigproc + l1sp_after_dnnroi.jsonnet).
  // PRODUCTION 0.889 since doc pdvd/100 sec 7.5 (owner 2026-09-13; measured in
  // doc pdvd/99, s = 0.889).  1.0 = the pre-flip SP (compiled config byte-identical
  // to the job before the knob); the runner's --top-gain-scale overrides it.  The
  // toolkit sp.jsonnet default stays 1.0, so no other importer moves.
  top_gain_scale           = 0.889,
  // Wire-geometry file for this job only (e.g. 'protodunevd-wires-larsoft-v5.json.bz2', the file the July
  // SP frames were made with).  '' = params.files.wires (production) => compiled config byte-identical.
  // The runner's --wires sets it.  doc pdvd/99 sec 4.3 (SP depends on the per-plane channel order).
  wires_file               = '',
)

  local tools = if wires_file == '' then tools_all
                else tools_maker(params { files+: { wires: wires_file } });
  local use_resampler = (reality == 'data');

  local base = import 'pgrapher/experiment/protodunevd/chndb-base.jsonnet';
  local chndb = [{
    type: 'OmniChannelNoiseDB',
    name: 'ocndbperfect%d' % n,
    data: base(params, tools.anodes[n], tools.field, tools.anodes[n].data.ident, use_freqmask=use_freqmask) { dft: wc.tn(tools.dft) },
    uses: [tools.anodes[n], tools.field, tools.dft],
  } for n in std.range(0, std.length(tools.anodes) - 1)];

  local nf_maker = import 'pgrapher/experiment/protodunevd/nf.jsonnet';
  local nf_pipes = [nf_maker(params, tools.anodes[n], chndb[n], tools.anodes[n].data.ident, name='nf%d' % tools.anodes[n].data.ident,
                             debug_dump_path=debug_dump_path, debug_dump_groups=debug_dump_groups,
                             shield_dump_path=shield_dump_path)
                    for n in std.range(0, std.length(tools.anodes) - 1)];

  // DNN-ROI input tags (loose_lf*, mp2_roi*, mp3_roi*, decon_charge*) require
  // OmnibusSigProc to be in debug + multi-plane-protection mode.
  local sp_maker = import 'pgrapher/experiment/protodunevd/sp.jsonnet';
  local sp_override = { sparse: sigoutform == 'sparse' }
                      + (if use_dnnroi
                         then { use_roi_debug_mode: true, use_multi_plane_protection: true }
                         else {});
  local sp = sp_maker(params, tools, sp_override);
  // L1SP OFF (l1sp_pd_mode='') so the DNN-ROI debug input tags survive.
  // ctoffset = +4 us (= the toolkit sp.jsonnet default; per-side knob left at
  // its legacy value).  REVERTED 2026-07-13 from the -5.5 us of §8.11 back to
  // +4 us: the raw-vs-deconvolved W-plane rising-edge alignment (§8.13 of
  // docs/qlmatch/pdvd-anode-time-consistency.md) is the primary, per-hit clock
  // reference, and it favours +4 us -- the deconvolved (gauss) peak sits on the
  // raw collection rising edge (charge arrival) at +4 us, whereas -5.5 us places
  // it ~9.5 us (19 ticks) too early.  The §8.11 -5.5 us was an endpoint-geometry
  // shift that folded the ~+1.5 cm near-anode DETECTION LOSS into the clock (its
  // own "competing reading"); that reading is now adopted -- the +1.5 cm is real
  // physics, not a clock error, so ctoffset is restored to its FR-consistent +4 us.
  // (Velocity / §8.12 convention re-evaluation deferred.)
  local sp_pipes = [sp.make_sigproc(a, l1sp_pd_mode='',
                                    ctoffset_b=4*wc.microsecond,
                                    ctoffset_t=4*wc.microsecond,
                                    top_gain_scale=top_gain_scale)
                    for a in tools.anodes];

  // TorchService instance shared by all per-anode DNN-ROI nodes.
  local ts = {
    type: 'TorchService',
    name: 'dnnroi_pdvd',
    data: {
      model: dnnroi_model,
      device: dnnroi_device,
      concurrency: dnnroi_concurrency,
    },
  };

  local dnnroi_maker = import 'pgrapher/experiment/protodunevd/dnnroi_pp.jsonnet';
  // Per-anode debug-file basename; when empty, the C++ node skips the dump.
  local _per_anode_dbg(n) =
    if dnnroi_debugfile == '' then ''
    else '%s_anode%d' % [dnnroi_debugfile, n];
  local dnnroi_inner_pipes = [dnnroi_maker(tools.anodes[n], ts,
                                           nticks=dnnroi_nticks,
                                           tick_per_slice=dnnroi_tick_per_slice,
                                           output_scale=dnnroi_output_scale,
                                           mask_thresh=dnnroi_mask_thresh,
                                           nchunks=dnnroi_nchunks,
                                           nchan=dnnroi_nchan,
                                           debugfile=_per_anode_dbg(tools.anodes[n].data.ident))
                              for n in std.range(0, std.length(tools.anodes) - 1)];

  local resamplers_config = import 'pgrapher/common/resamplers.jsonnet';
  local load_resamplers = resamplers_config(g, wc, tools);
  local resamplers = load_resamplers.resamplers;

  // (Tick-count alignment for the DNN-ROI model is handled inside the C++
  // node via dnnroi_pp.jsonnet's `tick_pad_multiple=128` — see
  // DNNROIFinding.cxx (pad_mult branch).  No Reframer needed in this driver.)

  // L1SP-after-DNN envelope.  When use_l1sp_dnn=true && use_dnnroi=true the
  // envelope wraps SP + DNN-ROI + L1SP into a single per-anode subgraph (the
  // FrameSplitter sits before SP since OmnibusSigProc drops raw%d).
  local l1sp_dnn_maker = import 'pgrapher/experiment/protodunevd/l1sp_after_dnnroi.jsonnet';

  // TorchService for the L1SP DNN tagger.  Built when mode is 'dnn' or
  // 'hybrid'; null otherwise so the heuristic-only path doesn't pull
  // libtorch a second time.
  local l1sp_torch_service =
    if (l1sp_pd_mode == 'dnn' || l1sp_pd_mode == 'hybrid') then {
      type: 'TorchService',
      name: 'l1sp_dnn_pdvd',
      data: {
        model:       l1sp_pd_dnn_model,
        device:      l1sp_pd_dnn_device,
        concurrency: l1sp_pd_dnn_concurrency,
      },
    } else null;

  // Final frame sink: write the post-DNN frame with standard SP trace tags
  // gauss%d (= DNN-ROI output, U+V from the model + W passthrough) and
  // wiener%d (= SP Wiener, carrying the per-channel threshold summary).
  // When use_l1sp_dnn=true, also include raw%d so the L1SP-modified frame's
  // raw channel survives to the archive.
  // When use_l1sp_dnn=false, the post-DNN frame only carries the DNN-ROI
  // output (no wiener traces — they live upstream in SP and are not
  // propagated by the inner subgraph) so we save just gauss%d.  The
  // dnnsp_to_gauss_retag below renames the inner subgraph's dnnsp%d trace
  // tag to gauss%d so this sink config matches what downstream Magnify
  // expects.
  local final_frame_sink = function(n)
    g.pnode({
      type: 'FrameFileSink',
      name: 'dnnroiframesink%d' % n,
      data: {
        outname: '%s-anode%d.tar.bz2' % [sp_prefix, n],
        tags: if use_l1sp_dnn
              then ['gauss%d' % n, 'wiener%d' % n, 'raw%d' % n]
              else ['gauss%d' % n],
        digitize: false,
        masks: true,
      },
    }, nin=1, nout=0);

  // Rename the inner DNN-ROI subgraph's unified trace tag dnnsp%d → gauss%d
  // for the L1SP-off path so Magnify's expectation (gauss%d = DNN-ROI
  // output) is met without changing the canonical dnnroi_pp.jsonnet
  // (which other chains, including the L1SP-after-DNN envelope, depend on).
  local dnn_to_gauss_retag = function(n)
    g.pnode({
      type: 'Retagger',
      name: 'dnnsp_to_gauss%d' % n,
      data: {
        tag_rules: [{
          frame: { ['dnnsp%d' % n]: 'gauss%d' % n },
          merge: { ['dnnsp%d' % n]: 'gauss%d' % n },
        }],
      },
    }, nin=1, nout=1);

  local per_anode_graph(n) =
    local anode = tools.anodes[n];
    local aid = anode.data.ident;

    local src = g.pnode({
      type: 'FrameFileSource',
      name: 'origframesrc%d' % aid,
      data: {
        inname: '%s-anode%d.tar.bz2' % [orig_prefix, aid],
        tags: ['orig'],
      } + (
        // Top-CRP (n>=4) orig frames are mislabeled 512 ns by the upstream
        // extraction; the top is physically 500 ns.  Re-stamp the tick at
        // read time so NF/SP run on the correct grid.  See wct-nf-sp.jsonnet.
        if use_resampler && n >= 4 then { tick: 500 * wc.ns } else {}
      ),
    }, nin=0, nout=1);

    local sp_dnn_l1sp_segment =
      if use_dnnroi && use_l1sp_dnn
      then [l1sp_dnn_maker(anode, sp_pipes[n], dnnroi_inner_pipes[n],
                           tools, params,
                           l1sp_pd_adj_enable=l1sp_pd_adj_enable,
                           l1sp_pd_adj_max_hops=l1sp_pd_adj_max_hops,
                           l1sp_pd_planes=l1sp_pd_planes,
                           l1sp_pd_dump_mode=l1sp_pd_mode,
                           l1sp_pd_dump_path=l1sp_pd_dump_path,
                           l1sp_pd_wf_dump_path=l1sp_pd_wf_dump_path,
                           l1sp_pd_dump_all_rois=l1sp_pd_dump_all_rois,
                           l1sp_pd_torch_service=l1sp_torch_service,
                           l1sp_pd_dnn_threshold=l1sp_pd_dnn_threshold,
                           l1sp_pd_dnn_threshold_bottom=l1sp_pd_dnn_threshold_bottom,
                           l1sp_pd_dnn_threshold_top=l1sp_pd_dnn_threshold_top,
                           l1sp_pd_adj_dnn_veto=l1sp_pd_adj_dnn_veto,
                           l1sp_pd_dnn_window_ticks=l1sp_pd_dnn_window_ticks,
                           l1sp_pd_dnn_debug_path=l1sp_pd_dnn_debug_path,
                           l1sp_pd_gmax_min=l1sp_pd_gmax_min,
                           l1sp_pd_min_length=l1sp_pd_min_length,
                           l1sp_pd_energy_frac_thr=l1sp_pd_energy_frac_thr,
                           top_gain_scale=top_gain_scale)]
      else [sp_pipes[n]]
           + (if use_dnnroi
              then [dnnroi_inner_pipes[n], dnn_to_gauss_retag(tools.anodes[n].data.ident)]
              else []);

    g.pipeline(
      [src]
      + (if use_resampler && n < 4 then [resamplers[n]] else [])
      + [nf_pipes[n]]
      + sp_dnn_l1sp_segment
      + [final_frame_sink(aid)],
      'nfspdnn_pipe_%d' % n);

  local graphs = [per_anode_graph(n) for n in anode_indices];

  local all_edges = std.foldl(function(acc, gr) acc + g.edges(gr), graphs, []);
  local all_uses  = std.foldl(function(acc, gr) acc + g.uses(gr),  graphs, []);

  local app = {
    type: 'Pgrapher',
    data: { edges: all_edges },
  };

  local cmdline = {
    type: 'wire-cell',
    data: {
      plugins: [
        'WireCellGen',
        'WireCellPgraph',
        'WireCellSio',
        'WireCellSigProc',
        'WireCellAux',
        'WireCellPytorch',
      ],
      apps: ['Pgrapher'],
    },
  };

  [cmdline] + all_uses + [app]
