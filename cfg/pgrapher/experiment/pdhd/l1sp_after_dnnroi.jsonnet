// Per-anode envelope that wires L1SPFilterPD AFTER a DNN-ROI subgraph,
// feeding the DNN output to L1SP as the "signal" channel.
//
// The envelope ENCOMPASSES the SP node + DNN-ROI subgraph because
// OmnibusSigProc drops the raw%d trace from its output, so the
// FrameSplitter that preserves raw must sit BEFORE SP (matching the
// pattern in sp.jsonnet's standard L1SP block).
//
// Architecture (mirrors sp.jsonnet:215-228 but with DNN-ROI inserted
// between SP and the merger):
//
//   post-NF frame (has raw%d) --FrameSplitter--> [port 0] sp_pipe -> dnnroi_pipe
//                                                                     -> Retagger(dnnsp%d -> gauss%d & wiener%d)
//                                                [port 1] -------------------------- (raw%d preserved) ---+
//                                                                                                         |
//                                                  [port 0 = DNN-relabeled gauss/wiener] -----+           |
//                                                                                             v           v
//                                                                                       FrameMerger ('replace')
//                                                                                       mergemap: gauss, wiener, raw
//                                                                                             |
//                                                                                             v
//                                                                                       L1SPFilterPD
//                                                                                       (adctag=raw, sigtag=gauss, outtag=gauss)
//                                                                                             |
//                                                                                             v
//                                                                                       FrameMerger ('replace')
//                                                                                       mergemap: L1SP-gauss -> gauss & wiener; raw passthrough
//                                                                                             |
//                                                                                             v
//                                                                                          output frame
//                                                                                  { raw, gauss (L1SP-DNN), wiener (L1SP-DNN) }
//
// L1SPFilterPD config block is copy-paste from sp.jsonnet:130-174
// (the canonical PDHD tuning).  Keep the two in sync if either is
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
         // ── DNN-mode opts (l1sp_pd_dump_mode == 'dnn') ────────────────
         // ``l1sp_pd_torch_service`` is the pre-constructed TorchService
         // jsonnet object that L1SPFilterPD will call via ITensorForward.
         // Must be non-null when l1sp_pd_dump_mode == 'dnn'.
         l1sp_pd_torch_service=null,
         l1sp_pd_dnn_threshold=0.9945,
         // When true, the DNN also vetoes adjacency-promoted ROIs
         // (i.e. ROIs whose heuristic flag_l1 == 0 but flag_l1_adj == 1).
         // Default true since 2026-05-25 — set to false to recover the
         // original behaviour where adj-promoted ROIs bypassed the DNN.
         l1sp_pd_adj_dnn_veto=true,
         l1sp_pd_dnn_window_ticks=256,
         l1sp_pd_dnn_debug_path='',
         // ── Loose-heur overrides (DNN-chain only) ─────────────────────
         // Defaults match the C++ defaults (trad-chain values); loosen
         // when running heuristic on DNN ROIs.
         l1sp_pd_gmax_min=1500.0,
         l1sp_pd_min_length=30,
         l1sp_pd_energy_frac_thr=0.66)

  local n = anode.data.ident;
  local l1sp_planes = if l1sp_pd_planes != null then l1sp_pd_planes
                      else if n == 0 then [0] else [0, 1];
  local gain_scale = params.elec.gain / (14.0 * wc.mV / wc.fC);

  // Splitter at the head: input is the post-NF frame (still carries
  // raw%d).  Port 0 → SP → DNN → relabel.  Port 1 → raw side
  // (unchanged copy carrying raw%d).
  local pre_split = pg.pnode({
    type: 'FrameSplitter', name: 'predannsplit%d' % n,
  }, nin=1, nout=2);

  // Inline the SP + optional SP frame tap + DNN-ROI + Retagger chain
  // on port 0 of pre_split.  Using pg.pipeline composes them into a
  // single pnode whose input port is sp_pipe's input and whose output
  // is the last node's output.
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

  // Merge DNN-relabeled side (port 0: gauss%d, wiener%d) with the
  // raw-preserving side (port 1: raw%d) so L1SP sees raw + gauss in one
  // frame.  'replace' rule: port-0 wins per channel.
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

  // L1SPFilterPD — identical config block to sp.jsonnet:130-174.
  local l1sp_node = pg.pnode({
    type: 'L1SPFilterPD',
    name: 'l1sppd_dnn%d' % n,
    data: {
      dft: wc.tn(tools.dft),
      anode: wc.tn(anode),
      kernels_file: 'pdhd_l1sp_kernels.json.bz2',
      adctag: 'raw%d' % n,
      sigtag: 'gauss%d' % n,
      outtag: 'gauss%d' % n,
      process_planes: l1sp_planes,
      kernels_scale:       gain_scale,
      l1_raw_asym_eps:     20.0 * gain_scale,
      raw_ROI_th_adclimit: 10.0 * gain_scale,
      adc_sum_threshold:  160.0 * gain_scale,
      gauss_filter: 'HfFilter:Gaus_wide',
      l1_len_very_long:  140,
      l1_asym_very_long: 0.35,
      l1_adj_enable: l1sp_pd_adj_enable,
      l1_adj_max_hops: l1sp_pd_adj_max_hops,
      // Phase-C loose-heur overrides: DNN ROIs are typically shorter than
      // trad ROIs for the same signal, so default pre-filters (gmax>=1500,
      // min_length>=30, energy_frac>=0.66) reject many candidates that the
      // trad chain catches. Since the DNN L1SP refines downstream, the
      // heuristic here can be intentionally loose -- false positives at
      // this stage are vetoed by DNN.
      l1_gmax_min:           l1sp_pd_gmax_min,
      l1_min_length:         l1sp_pd_min_length,
      l1_energy_frac_thr:    l1sp_pd_energy_frac_thr,
      mode: l1sp_pd_dump_mode,   // 'process' | 'dump' | 'dnn' | 'hybrid'
      dump_mode: l1sp_pd_dump_mode == 'dump',
      dump_path: l1sp_pd_dump_path,
      dump_tag: 'apa%d' % n,
      waveform_dump_path: l1sp_pd_wf_dump_path,
      dump_all_rois: l1sp_pd_dump_all_rois,
      // DNN-mode plumbing — required for mode='dnn' and mode='hybrid'.
      // Ignored by L1SPFilterPD when mode is process/dump.
      forward: if (l1sp_pd_dump_mode == 'dnn' || l1sp_pd_dump_mode == 'hybrid')
                  && l1sp_pd_torch_service != null
               then wc.tn(l1sp_pd_torch_service)
               else '',
      dnn_threshold:    l1sp_pd_dnn_threshold,
      l1_adj_dnn_veto:  l1sp_pd_adj_dnn_veto,
      dnn_window_ticks: l1sp_pd_dnn_window_ticks,
      dnn_debug_path:   l1sp_pd_dnn_debug_path,
    },
  }, nin=1, nout=1,
     uses=[tools.dft, anode] +
          (if (l1sp_pd_dump_mode == 'dnn' || l1sp_pd_dump_mode == 'hybrid')
              && l1sp_pd_torch_service != null
           then [l1sp_pd_torch_service] else []));

  // Final merger: L1SP-modified gauss replaces gauss AND wiener;
  // raw passes through from the rawsigmerge.
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

  // Splitter after the rawsigmerge so that the same combined frame can
  // feed BOTH the L1SP node (which will modify gauss) AND the final
  // merger port 1 (which carries the original gauss/wiener/raw for
  // channels L1SP did not touch).  Mirrors sigsplit in sp.jsonnet:188.
  local post_merge_split = pg.pnode({
    type: 'FrameSplitter', name: 'dnnpostmergesplit%d' % n,
  }, nin=1, nout=2);

  pg.intern(
    innodes=[pre_split],
    centernodes=[sp_dnn_chain, rawsigmerge, post_merge_split, l1sp_node, final_merger],
    edges=[
      // Port 0 of split → SP+DNN chain → merger port 0
      pg.edge(pre_split,        sp_dnn_chain,     0, 0),
      pg.edge(sp_dnn_chain,     rawsigmerge,      0, 0),
      // Port 1 of split → merger port 1 (carries raw)
      pg.edge(pre_split,        rawsigmerge,      1, 1),
      // Merger → post-merge split → L1SP / final-merger.1
      pg.edge(rawsigmerge,      post_merge_split, 0, 0),
      pg.edge(post_merge_split, l1sp_node,        0, 0),
      pg.edge(post_merge_split, final_merger,     1, 1),
      pg.edge(l1sp_node,        final_merger,     0, 0),
    ],
    oports=[final_merger.oports[0]],
    name='dnnroi_l1sp_%d' % n,
  )
