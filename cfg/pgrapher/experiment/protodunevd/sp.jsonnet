// This provides signal processing related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

local spfilt = import 'pgrapher/experiment/protodunevd/sp-filters.jsonnet';

function(params, tools, override = {}) {

  local pc = tools.perchanresp_nameuses,

  // pDSP needs a per-anode sigproc
  //
  // l1sp_pd_mode: 'dump' (default; ROI tagger ON, LASSO writeback OFF — for
  //               tagger validation prior to kernel generation)
  //               / 'process' (full L1SP fit and replacement; requires
  //               kernels_file populated below)
  //               / '' (OFF, bypass L1SP entirely — bare OmnibusSigProc output)
  // l1sp_pd_dump_path: directory to write per-event NPZ files when mode='dump'
  // l1sp_pd_wf_dump_path: directory to write per-triggered-ROI waveform NPZ
  //               files when mode='process' (raw/decon/lasso/smeared)
  // l1sp_pd_planes: plane indices processed by L1SPFilterPD (default [0,1] = U+V)
  // l1sp_pd_adj_enable / l1sp_pd_adj_max_hops: cross-channel adjacency expansion
  //               knobs, mirror the PDHD defaults (see sigproc/docs/l1sp/L1SPFilterPD.md).
  make_sigproc(anode, name=null,
               l1sp_pd_mode='dump',
               l1sp_pd_dump_path='',
               l1sp_pd_wf_dump_path='',
               l1sp_pd_planes=[0, 1],
               l1sp_pd_adj_enable=true,
               l1sp_pd_adj_max_hops=3)::
    // Top (_t) vs bottom (_b) anode filter suffix.  Bottom = ident 0..3,
    // top = ident 4..7.  See sp-filters.jsonnet for the registered names.
    local sfx = if anode.data.ident < 4 then '_b' else '_t';
    local sp_node = g.pnode({
      type: 'OmnibusSigProc',
      name:
        if std.type(name) == 'null'
        then anode.name + 'sigproc%d' % anode.data.ident
        else name,

      data: {
      /**  
       *  Default SP parameters (till May 2019)
       */
      // anode: wc.tn(anode),
      // field_response: wc.tn(tools.field),
      // per_chan_resp: pc.name,
      // fft_flag: 0,  // 1 is faster but higher memory, 0 is slightly slower but lower memory
      // postgain: 1,  // default 1.2
      // ADC_mV: 4096 / (1400.0 * wc.mV),  // default 4096/2000
      // r_fake_signal_low_th: 400,  // default 500
      // r_fake_signal_high_th: 800,  // default 1000
      // r_fake_signal_low_th_ind_factor: 1.5,  // default 1
      // r_fake_signal_high_th_ind_factor: 1.5,  // default 1
      // troi_col_th_factor: 5.0,  // default 5
      // troi_ind_th_factor: 3.5,  // default 3
      // r_th_factor: 3.5,  // default 3

      /**  
       *  Optimized SP parameters (May 2019)
       *  Associated tuning in sp-filters.jsonnet
       */

      local resolution = params.adc.resolution,
      local fullscale = if anode.data.ident < 4
                        then params.adc.fullscale[1] - params.adc.fullscale[0]
                        else 2.0*wc.volt,
      local ADC_mV_ratio = ((1 << resolution) - 1 ) / fullscale,

      anode: wc.tn(anode),
      dft: wc.tn(tools.dft),
      field_response: wc.tn(tools.field),
      // elecresponse: wc.tn(tools.elec_resp),
      elecresponse: if anode.data.ident < 4
                    then wc.tn(tools.elec_resps[0])
                    else wc.tn(tools.elec_resps[1]),

      // Per-anode-side filter type-name overrides.  All values are equal
      // between top and bottom for now; the split is structural.
      ROI_tight_lf_filter:   'ROI_tight_lf'   + sfx,
      ROI_tighter_lf_filter: 'ROI_tighter_lf' + sfx,
      ROI_loose_lf_filter:   'ROI_loose_lf'   + sfx,
      Gaus_wide_filter:      'Gaus_wide'      + sfx,
      Wiener_tight_filters: ['Wiener_tight_U' + sfx,
                             'Wiener_tight_V' + sfx,
                             'Wiener_tight_W' + sfx],
      Wiener_wide_filters:  ['Wiener_wide_U'  + sfx,
                             'Wiener_wide_V'  + sfx,
                             'Wiener_wide_W'  + sfx],
      // Default Wire_filters layout is [ind, ind, col]; preserve it.
      Wire_filters:         ['Wire_ind'       + sfx,
                             'Wire_ind'       + sfx,
                             'Wire_col'       + sfx],
      ftoffset: 0.0, // default 0.0
      // ctoffset: 1.0*wc.microsecond, // default -8.0
      ctoffset: 4*wc.microsecond, //consistent with FR: protodunevd_FR_imbalance3p_260501.json.bz2
      per_chan_resp: pc.name,
      fft_flag: 0,  // 1 is faster but higher memory, 0 is slightly slower but lower memory
      postgain: 1.0,  // default 1.2
      ADC_mV: ADC_mV_ratio, // 4096 / (1400.0 * wc.mV), 
      troi_col_th_factor: 5.0,  // default 5
      troi_ind_th_factor: 3.0,  // default 3
      lroi_rebin: 6, // default 6
      lroi_th_factor: 3.5, // default 3.5
      lroi_th_factor1: 0.7, // default 0.7
      lroi_jump_one_bin: 1, // default 0

      r_th_factor: 3.0,  // default 3
      r_fake_signal_low_th: 375,  // default 500
      r_fake_signal_high_th: 750,  // default 1000
      r_fake_signal_low_th_ind_factor: 1.0,  // default 1
      r_fake_signal_high_th_ind_factor: 1.0,  // default 1      
      r_th_peak: 3.0, // default 3.0
      r_sep_peak: 6.0, // default 6.0
      r_low_peak_sep_threshold_pre: 1200, // default 1200


      // frame tags
      wiener_tag: 'wiener%d' % anode.data.ident,
      wiener_threshold_tag: 'threshold%d' % anode.data.ident,
      decon_charge_tag: 'decon_charge%d' % anode.data.ident,
      gauss_tag: 'gauss%d' % anode.data.ident,

      use_roi_debug_mode: false,
      tight_lf_tag: 'tight_lf%d' % anode.data.ident,
      loose_lf_tag: 'loose_lf%d' % anode.data.ident,
      cleanup_roi_tag: 'cleanup_roi%d' % anode.data.ident,
      break_roi_loop1_tag: 'break_roi_1st%d' % anode.data.ident,
      break_roi_loop2_tag: 'break_roi_2nd%d' % anode.data.ident,
      shrink_roi_tag: 'shrink_roi%d' % anode.data.ident,
      extend_roi_tag: 'extend_roi%d' % anode.data.ident,

      use_multi_plane_protection: true,
      mp3_roi_tag: 'mp3_roi%d' % anode.data.ident,
      mp2_roi_tag: 'mp2_roi%d' % anode.data.ident,
      
      isWrapped: false,
      // process_planes: [0, 2],

      } + override,
    }, nin=1, nout=1, uses=[anode, tools.dft, tools.field, tools.elec_resps[0], tools.elec_resps[1] ] + pc.uses + spfilt);

    // PDVD top electronics not yet validated for the L1SP fit; cap 'process'
    // to 'dump' on top anodes (ident >= 4) so callers can request process
    // mode globally without accidentally enabling LASSO writeback on top.
    local _eff_mode = if anode.data.ident >= 4 && l1sp_pd_mode == 'process'
                      then 'dump'
                      else l1sp_pd_mode;
    if _eff_mode == '' then sp_node
    else
      local n = anode.data.ident;
      // Per-region (top vs bottom) L1SP parameters.  Bottom = anodes 0..3
      // (params.elec[0], 7.8 mV/fC reference); top = anodes 4..7
      // (JsonElecResponse, fixed reference, no runtime gain knob).
      //
      // Raw-ADC thresholds and kernel amplitudes are tuned at the per-region
      // reference electronics and scale with runtime FE gain via
      // gain_scale.  Deconvolved-domain thresholds (l1_gmax_min, asym ratios,
      // lengths, energy fractions) operate on gain-normalised signals and
      // are gain-invariant — same convention as chndb-base.
      local gain_scale = if anode.data.ident < 4
                         then params.elec.gain / (7.8 * wc.mV / wc.fC)
                         else 1.0;
      // Per-region kernel JSON, generated offline via
      //   wirecell-sigproc gen-l1sp-kernels -d pdvd-bottom  pdvd_bottom_l1sp_kernels.json.bz2
      //   wirecell-sigproc gen-l1sp-kernels -d pdvd-top     pdvd_top_l1sp_kernels.json.bz2
      // Both files live in wire-cell-data/ (resolved via WIRECELL_PATH).
      // The path is set unconditionally; in dump mode init_resp() is
      // guarded by !m_dump_mode (L1SPFilterPD.cxx) so the file is not
      // actually loaded until the user opts into process mode (-w on
      // run_nf_sp_evt.sh).  Default behaviour therefore stays unchanged:
      // tagger runs, no LASSO writeback, SP output bit-identical to a
      // no-L1SP run.
      local kernels_file = if anode.data.ident < 4
                           then 'pdvd_bottom_l1sp_kernels.json.bz2'
                           else 'pdvd_top_l1sp_kernels.json.bz2';
      local l1sp_node = g.pnode({
        type: 'L1SPFilterPD',
        name: 'l1sppd%d' % n,
        data: {
          dft: wc.tn(tools.dft),
          anode: wc.tn(anode),
          kernels_file: kernels_file,
          adctag: 'raw%d' % n,
          sigtag: 'gauss%d' % n,
          outtag: 'gauss%d' % n,
          process_planes: l1sp_pd_planes,
          // ADC-domain thresholds and kernel amplitudes, scaled to runtime FE gain.
          // The PDHD-tuned numerical defaults (20/10/160 at the 14 mV/fC reference)
          // are reused as a starting point; PDVD calibration may retune them
          // once dump-mode tagger validation is complete.
          kernels_scale:       gain_scale,
          l1_raw_asym_eps:     20.0 * gain_scale,
          raw_ROI_th_adclimit: 10.0 * gain_scale,
          adc_sum_threshold:  160.0 * gain_scale,
          // Derive time-domain smearing kernel by IFFT of the SP Gaus_wide
          // filter so both are driven by the same sigma.  Use the per-side
          // (_b/_t) instance to match this anode's OmnibusSigProc.
          gauss_filter: 'HfFilter:Gaus_wide' + sfx,
          // Cross-channel adjacency expansion (default ON, hops=3 — PDHD
          // defaults).  PDVD-side calibration may justify different values.
          l1_adj_enable: l1sp_pd_adj_enable,
          l1_adj_max_hops: l1sp_pd_adj_max_hops,
          // PDHD's "very-long" arm (l1_len_very_long=140, l1_asym_very_long=0.35)
          // is left at the C++ default (OFF) here; revisit after PDVD tagger
          // validation if long-but-moderate-asym artifacts are observed.
          dump_mode: _eff_mode == 'dump',
          dump_path: l1sp_pd_dump_path,
          dump_tag: 'apa%d' % n,
          waveform_dump_path: l1sp_pd_wf_dump_path,
        // ── PDVD-tuned trigger-gate overrides (BOTTOM anodes only) ──────
        // Tuned against pdvd/sp_plot/handscan_039324_anode0.csv (run
        // 39324 events 0-5, bottom anode 0).  See pdvd/sp_plot/
        // eval_l1sp_trigger_pdvd.py for the validation harness.
        // Evaluated per-ROI on the C++ NPZ dumps: these four overrides
        // drop the per-ROI FP count from 70 → 62 with no loss of GT-
        // positive hits (27/28).  l1_asym_mod is left at C++ default
        // 0.40 even though the Python cluster-level sweep prefers 0.50:
        // at the per-ROI gate, raising it kills ch 386 of the evt-0
        // U 386-388 cluster (per-ROI asym=0.43; cluster max-asym=0.56,
        // but C++ checks each ROI individually).  Top anodes (n >= 4)
        // are left at the PDHD-inherited C++ defaults pending their own
        // hand-scan validation — toggle on by extending the conditional.
        } + (if n < 4 then {
          l1_len_long_mod:         180,     // C++ default 100
          l1_len_fill_shape:       90,      // C++ default 50
          l1_fill_shape_fill_thr:  0.30,    // C++ default 0.38
          l1_fill_shape_fwhm_thr:  0.25,    // C++ default 0.30
          // PDVD-only opt-in track veto (added in this commit).
          // Rejects sub-windows that look like real prolonged tracks rather
          // than L1SP unipolar lobes — the residual FP class against
          // pdvd/sp_plot/handscan_039324_anode0.csv (evt 0 U 82-87, 89-91,
          // 96; evt 3 U 319-331).  PDHD never sets this so its trigger
          // remains bit-identical to the pre-veto behaviour.
          l1_pdvd_track_veto_enable: true,
          l1_pdvd_track_high_asym:   0.85,
          l1_pdvd_track_long_cl:     170,
          l1_pdvd_track_med_cl:      100,
          l1_pdvd_track_med_fill:    0.40,
          l1_pdvd_track_med_fwhm:    0.40,
          peak_threshold: 1000,
          mean_threshold: 500,
        } else {}),
      }, nin=1, nout=1, uses=[tools.dft, anode]);
      // L1SPFilterPD needs both raw{n} and gauss{n} in the same frame.
      // OmnibusSigProc drops raw traces from its output, so we split the
      // input frame, run SP on one copy, then merge raw+gauss for L1SP.
      // After L1SP, a final merger replaces BOTH the gauss and wiener
      // traces of the original sp output with the L1SP-modified gauss
      // (the L1SP fit is the canonical deconvolved signal post-L1, so
      // both deconvolved tags should reflect it).  In dump mode the L1SP
      // node passes the frame through unchanged, so the merger output is
      // bit-identical to a no-L1SP run.
      local rawsplit     = g.pnode({type: 'FrameSplitter', name: 'rawsplit%d' % n}, nin=1, nout=2);
      local sigsplit     = g.pnode({type: 'FrameSplitter', name: 'sigsplit%d' % n}, nin=1, nout=2);
      local rawsigmerge  = g.pnode({
        type: 'FrameMerger', name: 'rawsigmerge%d' % n,
        data: {
          rule: 'replace',
          mergemap: [
            ['raw%d' % n, 'raw%d' % n, 'raw%d' % n],
            ['gauss%d' % n, 'gauss%d' % n, 'gauss%d' % n],
          ],
        },
      }, nin=2, nout=1);
      local final_merger = g.pnode({
        type: 'FrameMerger', name: 'l1spfinal%d' % n,
        data: {
          rule: 'replace',
          mergemap: [
            ['gauss%d' % n, 'gauss%d'  % n, 'gauss%d'  % n],   // L1SP gauss -> output gauss
            ['gauss%d' % n, 'wiener%d' % n, 'wiener%d' % n],   // L1SP gauss -> output wiener
          ],
        },
      }, nin=2, nout=1);
      g.intern(
        innodes=[rawsplit],
        centernodes=[sp_node, sigsplit, rawsigmerge, l1sp_node, final_merger],
        edges=[
          g.edge(rawsplit,     sp_node,       0, 0),
          g.edge(sp_node,      sigsplit,      0, 0),
          g.edge(sigsplit,     rawsigmerge,   1, 0),
          g.edge(rawsplit,     rawsigmerge,   1, 1),
          g.edge(rawsigmerge,  l1sp_node,     0, 0),
          g.edge(l1sp_node,    final_merger,  0, 0),
          g.edge(sigsplit,     final_merger,  0, 1),
        ],
        oports=[final_merger.oports[0]],
        name='sigproc_l1sppd_%d' % n
      ),

}
