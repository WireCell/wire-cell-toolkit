// This provides signal processing related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// BIG FAT FIXME: we are taking from uboone.  If PDSP needs tuning do
// four things: 0) read this comment, 1) cp this file into pdsp/, 2)
// fix the import and 3) delete this comment.
local spfilt = import 'pgrapher/experiment/pdhd/sp-filters.jsonnet';

function(params, tools, override = {}) {

  local pc = tools.perchanresp_nameuses,
  local fltr = tools.fltrespuses,

  local resolution = params.adc.resolution,
  local fullscale = params.adc.fullscale[1] - params.adc.fullscale[0],
  local ADC_mV_ratio = ((1 << resolution) - 1 ) / fullscale,

  // pDSP needs a per-anode sigproc
  //
  // l1sp_pd_mode: 'process' (default, ON; replaces gauss/wiener with L1SP fit)
  //               / 'dump' (calibration dump of per-ROI asymmetry quantities to NPZ)
  //               / '' (OFF, bypass L1SP entirely — bare OmnibusSigProc output)
  // l1sp_pd_dump_path: directory to write per-event NPZ files when mode='dump'
  // l1sp_pd_planes: plane indices processed by L1SPFilterPD.
  //   null (default): APA0 → [0] (U only; V anomalous), APA1-3 → [0,1] (U+V).
  //   explicit list: overrides the per-APA default.
  make_sigproc(anode, name=null,
               l1sp_pd_mode='process',
               l1sp_pd_dump_path='',
               l1sp_pd_wf_dump_path='',
               l1sp_pd_planes=null,
               // Cross-channel adjacency expansion (default ON; see
               // sigproc/docs/l1sp/L1SPFilterPD.md).  Pass false to recover
               // the pre-2026-05-02 behaviour (no neighbour-driven promotion).
               l1sp_pd_adj_enable=true,
               // Iterative-expansion hop cap (default 3 ⇔ ±3 channels from any
               // original donor).  Set 1 to recover the pre-2026-05-03
               // single-hop behaviour (donors must be originally-triggered).
               l1sp_pd_adj_max_hops=3)::
    local l1sp_planes = if l1sp_pd_planes != null then l1sp_pd_planes
                        else if anode.data.ident == 0 then [0] else [0, 1];
    local sp_node = g.pnode({
      type: 'OmnibusSigProc',
      name:
        if std.type(name) == 'null'
        then anode.name + 'sigproc%d' % anode.data.ident
        else name,

      data: {
        anode: wc.tn(anode),
      dft: wc.tn(tools.dft),
      field_response: wc.tn(tools.fields[anode.data.ident]),
      filter_responses_tn: if anode.data.ident == 0
                           then ["FilterResponse:plane0",
                                 "FilterResponse:plane2", "FilterResponse:plane1"]
                           else [ ],
      elecresponse: wc.tn(tools.elec_resp),
      ftoffset: 0.0, // default 0.0
      ctoffset: 1.0*wc.microsecond, // default -8.0
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

      r_th_factor: if anode.data.ident==0 then 2.5 else 3.0,  // default 3
      r_fake_signal_low_th: 375,  // default 500
      r_fake_signal_high_th: 750,  // default 1000
      r_fake_signal_low_th_ind_factor: 1.0,  // default 1
      r_fake_signal_high_th_ind_factor: 1.0,  // default 1
      r_th_peak: 3.0, // default 3.0
      r_sep_peak: 6.0, // default 6.0
      r_low_peak_sep_threshold_pre: 1200, // default 1200


      // frame tags
      wiener_tag: 'wiener%d' % anode.data.ident,
      // wiener_threshold_tag: 'threshold%d' % anode.data.ident, // deprecated
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

      use_multi_plane_protection: false,
      mp3_roi_tag: 'mp3_roi%d' % anode.data.ident,
      mp2_roi_tag: 'mp2_roi%d' % anode.data.ident,
      // mp_th1: if anode.data.ident==0 then 200 else 1000,
      // mp_th2: if anode.data.ident==0 then 100 else 500,
      
      isWrapped: false,
      // process_planes: if anode.data.ident==0 then [0, 1] else [0, 1, 2],
      plane2layer: if anode.data.ident==0 then [0,2,1] else [0,1,2],

      Wiener_tight_filters: if anode.data.ident==0
                            then ["Wiener_tight_U_APA1", "Wiener_tight_W_APA1", "Wiener_tight_V_APA1"] // ind, ind, col
                            else ["Wiener_tight_U", "Wiener_tight_V", "Wiener_tight_W"],

      } + override,
    }, nin=1, nout=1, uses=[anode, tools.dft, tools.field, tools.fields[1], tools.fields[2], tools.fields[3], tools.elec_resp] + pc.uses + fltr.uses + spfilt);

    if l1sp_pd_mode == '' then sp_node
    else
      local n = anode.data.ident;
      // Raw-ADC thresholds are tuned at the 14 mV/fC reference and scale
      // with FE gain.  Deconvolved-domain thresholds (l1_gmax_min, asym
      // ratios, lengths, energy fractions) operate on gain-normalised
      // signals and are gain-invariant.  Same convention as chndb-base.
      local gain_scale = params.elec.gain / (14.0 * wc.mV / wc.fC);
      local l1sp_node = g.pnode({
        type: 'L1SPFilterPD',
        name: 'l1sppd%d' % n,
        data: {
          dft: wc.tn(tools.dft),
          anode: wc.tn(anode),
          kernels_file: "pdhd_l1sp_kernels.json.bz2",
          adctag: 'raw%d' % n,
          sigtag: 'gauss%d' % n,
          outtag: 'gauss%d' % n,
          process_planes: l1sp_planes,
          // ADC-domain thresholds and kernel amplitudes, scaled to runtime FE gain.
          // kernels_file is generated at 14 mV/fC; gain_scale corrects to actual gain.
          kernels_scale:       gain_scale,
          l1_raw_asym_eps:     20.0 * gain_scale,
          raw_ROI_th_adclimit: 10.0 * gain_scale,
          adc_sum_threshold:  160.0 * gain_scale,
          // Derive time-domain smearing kernel by IFFT of the SP Gaus_wide
          // filter so both are driven by the same sigma.
          gauss_filter: 'HfFilter:Gaus_wide',
          // Cross-channel adjacency expansion (default ON).  Threshold knobs
          // (overlap_pad, gap_max, len_ratio, loose_*) keep the values baked
          // into L1SPFilterPD.h unless overridden in the data block.
          l1_adj_enable: l1sp_pd_adj_enable,
          l1_adj_max_hops: l1sp_pd_adj_max_hops,
          dump_mode: l1sp_pd_mode == 'dump',
          dump_path: l1sp_pd_dump_path,
          dump_tag: 'apa%d' % n,
          waveform_dump_path: l1sp_pd_wf_dump_path,
        },
      }, nin=1, nout=1, uses=[tools.dft, anode]);
      // L1SPFilterPD needs both raw{n} and gauss{n} in the same frame.
      // OmnibusSigProc drops raw traces from its output, so we split the
      // input frame, run SP on one copy, then merge raw+gauss for L1SP.
      // After L1SP, a final merger replaces BOTH the gauss and wiener
      // traces of the original sp output with the L1SP-modified gauss
      // (the L1SP fit is the canonical deconvolved signal post-L1, so
      // both deconvolved tags should reflect it). With FrameMerger's
      // 'replace' rule (cxx:93-113) the *first* input wins per channel,
      // so L1SP feeds input 0 and its modified gauss traces are emitted
      // under both the gauss and wiener output tags. The untouched sp
      // copy on sigsplit port 0 (input 1) only contributes channels not
      // covered by L1SP, which in practice is none.
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
            ['gauss%d' % n, 'gauss%d'  % n, 'gauss%d'  % n],   // L1SP gauss → output gauss
            ['gauss%d' % n, 'wiener%d' % n, 'wiener%d' % n],   // L1SP gauss → output wiener
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
