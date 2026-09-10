// This provides signal processing related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// BIG FAT FIXME: we are taking from uboone.  If SBND needs tuning do
// four things: 0) read this comment, 1) cp this file into sbnd/, 2)
// fix the import and 3) delete this comment.
local spfilt = import 'pgrapher/experiment/sbnd/sp-filters.jsonnet';

function(params, tools, override = {}) {

  local pc = tools.perchanresp_nameuses,

  local resolution = params.adc.resolution,
  local fullscale = params.adc.fullscale[1] - params.adc.fullscale[0],
  local ADC_mV_ratio = ((1 << resolution) - 1 ) / fullscale,

  // SBND needs a per-anode sigproc
  make_sigproc(anode, name=null):: g.pnode({
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
      anode: wc.tn(anode),
      dft: wc.tn(tools.dft),
      field_response: wc.tn(tools.field),
      elecresponse: wc.tn(tools.elec_resp),
      // Rebaselining is done in NF with a more sophisticated method.  This key
      // MUST be present: OmnibusSigProc only reads it `if (config.isMember(...))`
      // and its C++ default is {0, 1, 2}, so omitting it rebaselines ALL planes
      // rather than none.
      rebase_planes: [],
      ftoffset: 0.0, // default 0.0
      ctoffset: 1.0*wc.microsecond, // default -8.0
      per_chan_resp: pc.name,
      fft_flag: 0,  // 1 is faster but higher memory, 0 is slightly slower but lower memory
      postgain: 1.0,  // default 1.2
      ADC_mV: ADC_mV_ratio, // 4096 / (1400.0 * wc.mV),

      // Tight-ROI thresholds.  These used to be selected here from the
      // 'enableLowROIThresholds' extVar, which made this file unusable outside
      // art/wcls (standalone toolkit jobs define no extVars) and, worse, failed
      // silently: jsonnet ignores an extVar nobody reads, so a copy of this file
      // with the switch dropped ran at the high thresholds while the fcl asked
      // for the low ones.  The switch now lives in the wcls jsonnets, which pass
      // the high values through `override`; see wcls-nf-sp.jsonnet.
      //
      // The defaults below are the LOW values, i.e. what every production fcl
      // asks for, so a caller that forgets the override lands on production
      // behaviour instead of silently coarser ROIs.  C++ defaults are 5 and 3.
      troi_col_th_factor: 3.0,
      troi_ind_th_factor: 1.8,

      // Prolonged-W-signal fix, part 1: MAD-based cal_RMS in ROI finding
      // (C++ default false).  A long track-along-drift W signal occupying
      // >~16% of the readout corrupts ROI_formation::cal_RMS's legacy
      // (16,50,84)-percentile noise estimate, pushing the tight-ROI
      // threshold (5*rms+1) above the signal's own median so only the
      // tallest dE/dx peaks form ROIs.  Measured on SBND MC run 270/6/46
      // ch 10038: decon RMS 2036 vs 62-95 on a normal W channel, i.e.
      // signal/RMS 3.6 vs ~82, and the ROI collapsed to 9 ticks.  MAD stays
      // robust to 50% occupancy.  Generic estimator, all planes/anodes.
      // See ai-helper issue #10.
      roi_mad_rms: true,

      // Prolonged-W-signal fix, part 2: disable BreakROI on the collection
      // plane.  With part 1 in place the long multi-peak W ROI survives to
      // refinement, where BreakROI would subtract a valley-to-valley linear
      // "baseline" that is actually real track charge (collection decon has
      // no LF filter, so its baseline needs no such fix), re-fragmenting the
      // signal into per-peak islands.  SBND uses the standard [U, V, W] slot
      // order on both anodes (no filter_responses_tn remap), so W is slot 2;
      // U and V keep the production 2 break loops.  Override with
      // [2, 2, 2] to recover the pre-tune behaviour.
      r_break_roi_loop_planes: [2, 2, 0],

      lroi_rebin: 6, // default 6
      lroi_th_factor: 3.5, // default 3.5
      lroi_th_factor1: 0.7, // default 0.7
      lroi_jump_one_bin: 1, // default 0

      r_th_factor: 3.0,  // default 3
      r_fake_signal_low_th: 400,  // default 500
      r_fake_signal_high_th: 800,  // default 1000
      r_fake_signal_low_th_ind_factor: 1.0,  // default 1
      r_fake_signal_high_th_ind_factor: 1.0,  // default 1      
      r_th_peak: 3.0, // default 3.0
      r_sep_peak: 6.0, // default 6.0
      r_low_peak_sep_threshold_pre: 1200, // default 1200


      // frame tags
      wiener_tag: 'wiener%d' % anode.data.ident,
      // 'wiener_threshold_tag' dropped: OmnibusSigProc logs it as obsolete and
      // carries thresholds in the summary of the 'wiener' tagged traces, which
      // is what magnify-sinks.jsonnet's threshold sink already reads.
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
      
      isWrapped: false,
      // process_planes: [0, 2],

    } + override,
  }, nin=1, nout=1, uses=[anode, tools.dft, tools.field, tools.elec_resp] + pc.uses + spfilt),

}
