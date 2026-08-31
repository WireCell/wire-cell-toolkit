// Per-anode envelope that wires L1SPFilterPD AFTER a DNN-ROI subgraph,
// feeding the DNN output to L1SP as the "signal" channel.  PDVD port of
// pdhd/l1sp_after_dnnroi.jsonnet — same envelope shape, PDVD-tuned
// L1SPFilterPD config block copied from protodunevd/sp.jsonnet.
//
// Architecture (mirrors the PDHD version):
//
//   post-NF frame (has raw%d) --FrameSplitter--> [port 0] sp -> dnn -> Retagger(dnnsp%d -> gauss%d & wiener%d)
//                                                [port 1] -------------------------- (raw%d preserved) ---+
//                                                  [port 0 = DNN gauss/wiener] -----+                     |
//                                                                                   v                     v
//                                                                              FrameMerger ('replace')
//                                                                              mergemap: gauss, wiener, raw
//                                                                                                |
//                                                                                                v
//                                                                                         L1SPFilterPD
//                                                                                                |
//                                                                                                v
//                                                                                         FrameMerger ('replace')
//                                                                                                |
//                                                                                                v
//                                                                              { raw, gauss (L1SP), wiener (L1SP) }
//
// L1SPFilterPD config below is a verbatim copy of the PDVD block in
// protodunevd/sp.jsonnet:165-243 (per-region kernels, gain_scale,
// gauss_filter suffix, PDVD-tuned trigger-gate overrides and opt-in
// track-veto, peak/mean thresholds).  Keep the two in sync if either is
// re-tuned.

local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

function(anode, sp_pipe, dnnroi_pipe, tools, params,
         l1sp_pd_adj_enable=true,
         l1sp_pd_adj_max_hops=3,
         l1sp_pd_planes=null,
         l1sp_pd_dump_mode='process',
         l1sp_pd_dump_path='',
         l1sp_pd_wf_dump_path='',
         l1sp_pd_dump_all_rois=false,
         // ── DNN-mode opts (mode=='dnn' or mode=='hybrid').  The PDVD L1SP
         //    DNN model lives at wire-cell-data/l1sp/pdvd/l1sp_dnn_pdvd_v1.ts;
         //    deploy_round2.md documents the round-2-derived thresholds.
         l1sp_pd_torch_service=null,
         l1sp_pd_dnn_threshold=0.94,         // single-threshold fallback
         l1sp_pd_dnn_threshold_bottom=null,  // per-CRP override (apa<4)
         l1sp_pd_dnn_threshold_top=null,     // per-CRP override (apa>=4)
         // When true, the DNN also vetoes adjacency-promoted ROIs
         // (i.e. ROIs whose heuristic flag_l1 == 0 but flag_l1_adj == 1).
         // Default true since 2026-05-25 — set to false to recover the
         // original behaviour where adj-promoted ROIs bypassed the DNN.
         l1sp_pd_adj_dnn_veto=true,
         l1sp_pd_dnn_window_ticks=256,
         l1sp_pd_dnn_debug_path='',
         // ── Loose-heur overrides (DNN-chain only) ─────────────────────
         // Defaults match the C++/PDVD-deployed values; loosen via the
         // runner's --loose-heur preset when running heuristic on DNN ROIs.
         l1sp_pd_gmax_min=1500.0,
         l1sp_pd_min_length=30,
         l1sp_pd_energy_frac_thr=0.66)

  local n = anode.data.ident;
  local sfx = if n < 4 then '_b' else '_t';
  local l1sp_planes = if l1sp_pd_planes != null then l1sp_pd_planes else [0, 1];
  // Per-CRP threshold resolution: prefer the explicit bottom/top TLA
  // when set, otherwise fall back to the single-threshold scalar.
  local per_crp_thr = if n < 4
                      then (if l1sp_pd_dnn_threshold_bottom != null
                            then l1sp_pd_dnn_threshold_bottom
                            else l1sp_pd_dnn_threshold)
                      else (if l1sp_pd_dnn_threshold_top != null
                            then l1sp_pd_dnn_threshold_top
                            else l1sp_pd_dnn_threshold);
  // In hybrid mode the DNN is the FP suppressor; double-suppressing
  // with the PDVD-specific track veto starves the DNN of candidates,
  // so disable track veto only when mode=='hybrid'.  In all other
  // modes (process/dump/dnn) keep the deployed PDVD veto on.
  local track_veto_enabled = !(l1sp_pd_dump_mode == 'hybrid');
  // Per-region gain reference: bottom (ident<4) follows params.elec.gain
  // at the 7.8 mV/fC PDVD-bottom reference; top (ident>=4) uses
  // JsonElecResponse and is gain-invariant.  Same convention as
  // protodunevd/sp.jsonnet:165-167.
  local gain_scale = if n < 4
                     then params.elec.gain / (7.8 * wc.mV / wc.fC)
                     else 1.0;
  local kernels_file = if n < 4
                       then 'pdvd_bottom_l1sp_kernels.json.bz2'
                       else 'pdvd_top_l1sp_kernels.json.bz2';

  local pre_split = pg.pnode({
    type: 'FrameSplitter', name: 'predannsplit%d' % n,
  }, nin=1, nout=2);

  // protodunevd/dnnroi_pp.jsonnet emits a single dnnsp%d trace tag (its
  // final Retagger merges dnnsp%du/v/w → dnnsp%d).  L1SP downstream
  // expects gauss%d (signal) and wiener%d (carrying the per-channel
  // threshold trace_summary).  A trace/merge Retagger relabels:
  //   trace: dnnsp%d → gauss%d   (rename, carries summary)
  //   merge: dnnsp%d → wiener%d  (alias the same data + summary)
  // Mirrors the PDHD envelope's dnn_relabel.
  local dnn_relabel = pg.pnode({
    type: 'Retagger',
    name: 'dnnsp_to_gauss%d' % n,
    data: {
      tag_rules: [{
        trace: {
          ['dnnsp%d' % n]: 'gauss%d' % n,
        },
        merge: {
          ['dnnsp%d' % n]: 'wiener%d' % n,
        },
      }],
    },
  }, nin=1, nout=1);

  local sp_dnn_chain = pg.pipeline(
    [sp_pipe, dnnroi_pipe, dnn_relabel],
    'sp_dnn_chain_%d' % n);

  local rawsigmerge = pg.pnode({
    type: 'FrameMerger', name: 'dnnrawsigmerge%d' % n,
    data: {
      rule: 'replace',
      mergemap: [
        ['gauss%d'  % n, 'gauss%d'  % n, 'gauss%d'  % n],
        ['wiener%d' % n, 'wiener%d' % n, 'wiener%d' % n],
        ['raw%d'    % n, 'raw%d'    % n, 'raw%d'    % n],
      ],
    },
  }, nin=2, nout=1);

  // L1SPFilterPD — config block copied from protodunevd/sp.jsonnet:178-243.
  local l1sp_node = pg.pnode({
    type: 'L1SPFilterPD',
    name: 'l1sppd_dnn%d' % n,
    data: {
      dft: wc.tn(tools.dft),
      anode: wc.tn(anode),
      kernels_file: kernels_file,
      adctag: 'raw%d' % n,
      sigtag: 'gauss%d' % n,
      outtag: 'gauss%d' % n,
      process_planes: l1sp_planes,
      kernels_scale:       gain_scale,
      l1_raw_asym_eps:     20.0 * gain_scale,
      raw_ROI_th_adclimit: 10.0 * gain_scale,
      adc_sum_threshold:  160.0 * gain_scale,
      gauss_filter: 'HfFilter:Gaus_wide' + sfx,
      l1_adj_enable:  l1sp_pd_adj_enable,
      l1_adj_max_hops: l1sp_pd_adj_max_hops,
      // Loose-heur overrides — exposed as TLAs so the runner's
      // --loose-heur preset can relax the three pre-filters when running
      // heuristic on DNN ROIs (which are typically shorter / lower-peak
      // than trad ROIs for the same physical signal).
      l1_gmax_min:        l1sp_pd_gmax_min,
      l1_min_length:      l1sp_pd_min_length,
      l1_energy_frac_thr: l1sp_pd_energy_frac_thr,
      mode: l1sp_pd_dump_mode,   // 'process' | 'dump' | 'dnn' | 'hybrid'
      dump_mode: l1sp_pd_dump_mode == 'dump',
      dump_path: l1sp_pd_dump_path,
      dump_tag: 'apa%d' % n,
      waveform_dump_path: l1sp_pd_wf_dump_path,
      dump_all_rois: l1sp_pd_dump_all_rois,
      // DNN-mode plumbing — required for mode='dnn' and mode='hybrid'.
      forward: if (l1sp_pd_dump_mode == 'dnn' || l1sp_pd_dump_mode == 'hybrid')
                  && l1sp_pd_torch_service != null
               then wc.tn(l1sp_pd_torch_service)
               else '',
      dnn_threshold:    per_crp_thr,
      l1_adj_dnn_veto:  l1sp_pd_adj_dnn_veto,
      dnn_window_ticks: l1sp_pd_dnn_window_ticks,
      dnn_debug_path:   l1sp_pd_dnn_debug_path,
    } + {
      // PDVD-tuned trigger-gate overrides (verbatim from sp.jsonnet:226-242,
      // except track-veto, which is disabled when mode=='hybrid' so the DNN
      // sees the full loose-heur candidate set rather than the track-veto-
      // residual; see deploy_round2.md §"Deviations from PDHD" for rationale.)
      l1_len_long_mod:         180,
      l1_len_fill_shape:       90,
      l1_fill_shape_fill_thr:  0.30,
      l1_fill_shape_fwhm_thr:  0.25,
      l1_pdvd_track_veto_enable: track_veto_enabled,
      l1_pdvd_track_high_asym:   0.85,
      l1_pdvd_track_long_cl:     170,
      l1_pdvd_track_med_cl:      100,
      l1_pdvd_track_med_fill:    0.40,
      l1_pdvd_track_med_fwhm:    0.40,
      peak_threshold: 1000,
      mean_threshold: 500,
    },
  }, nin=1, nout=1,
     uses=[tools.dft, anode] +
          (if (l1sp_pd_dump_mode == 'dnn' || l1sp_pd_dump_mode == 'hybrid')
              && l1sp_pd_torch_service != null
           then [l1sp_pd_torch_service] else []));

  local final_merger = pg.pnode({
    type: 'FrameMerger', name: 'dnnl1spfinal%d' % n,
    data: {
      rule: 'replace',
      mergemap: [
        ['gauss%d' % n, 'gauss%d'  % n, 'gauss%d'  % n],
        ['gauss%d' % n, 'wiener%d' % n, 'wiener%d' % n],
        ['raw%d'   % n, 'raw%d'    % n, 'raw%d'    % n],
      ],
    },
  }, nin=2, nout=1);

  local post_merge_split = pg.pnode({
    type: 'FrameSplitter', name: 'dnnpostmergesplit%d' % n,
  }, nin=1, nout=2);

  pg.intern(
    innodes=[pre_split],
    centernodes=[sp_dnn_chain, rawsigmerge, post_merge_split, l1sp_node, final_merger],
    edges=[
      pg.edge(pre_split,        sp_dnn_chain,     0, 0),
      pg.edge(sp_dnn_chain,     rawsigmerge,      0, 0),
      pg.edge(pre_split,        rawsigmerge,      1, 1),
      pg.edge(rawsigmerge,      post_merge_split, 0, 0),
      pg.edge(post_merge_split, l1sp_node,        0, 0),
      pg.edge(post_merge_split, final_merger,     1, 1),
      pg.edge(l1sp_node,        final_merger,     0, 0),
    ],
    oports=[final_merger.oports[0]],
    name='dnnroi_l1sp_%d' % n,
  )
