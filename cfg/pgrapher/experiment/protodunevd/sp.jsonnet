// This provides signal processing related pnodes,

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

local spfilt = import 'pgrapher/experiment/protodunevd/sp-filters.jsonnet';

function(params, tools, override = {}) {

  local pc = tools.perchanresp_nameuses,

  // pDSP needs a per-anode sigproc
  //
  // l1sp_pd_mode: '' (default, OFF) / 'process' (process triggered ROIs, still stubbed)
  //               / 'dump' (calibration dump of per-ROI asymmetry quantities to NPZ)
  // l1sp_pd_dump_path: directory to write per-event NPZ files when mode='dump'
  // l1sp_pd_planes: plane indices processed by L1SPFilterPD (default [0,1] = U+V)
  make_sigproc(anode, name=null,
               l1sp_pd_mode='',
               l1sp_pd_dump_path='',
               l1sp_pd_planes=[0, 1])::
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

    if l1sp_pd_mode == '' then sp_node
    else
      local n = anode.data.ident;
      local l1sp_node = g.pnode({
        type: 'L1SPFilterPD',
        name: 'l1sppd%d' % n,
        data: {
          dft: wc.tn(tools.dft),
          anode: wc.tn(anode),
          fields: wc.tn(tools.field),
          adctag: 'raw%d' % n,
          sigtag: 'gauss%d' % n,
          outtag: 'gauss%d' % n,
          process_planes: l1sp_pd_planes,
          // Derive time-domain smearing kernel by IFFT of the SP Gaus_wide
          // filter so both are driven by the same sigma.  Use the per-side
          // (_b/_t) instance to match this anode's OmnibusSigProc.
          gauss_filter: 'HfFilter:Gaus_wide' + sfx,
          dump_mode: l1sp_pd_mode == 'dump',
          dump_path: l1sp_pd_dump_path,
          dump_tag: 'apa%d' % n,
        },
      }, nin=1, nout=1, uses=[tools.dft, anode, tools.field]);
      // L1SPFilterPD needs both raw{n} and gauss{n} in the same frame.
      // OmnibusSigProc drops raw traces from its output, so we split the
      // input frame, run SP on one copy, then merge raw+gauss for L1SP.
      // The final output (sigsplit port 0 = gauss+wiener) is bit-identical
      // to a run without L1SP; L1SP's output is discarded via l1sp_sink.
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
      local l1sp_sink    = g.pnode({type: 'DumpFrames', name: 'l1spsnk%d' % n}, nin=1, nout=0);
      g.intern(
        innodes=[rawsplit],
        centernodes=[sp_node, sigsplit, rawsigmerge, l1sp_node, l1sp_sink],
        edges=[
          g.edge(rawsplit,    sp_node,      0, 0),
          g.edge(sp_node,     sigsplit,     0, 0),
          g.edge(sigsplit,    rawsigmerge,  1, 0),
          g.edge(rawsplit,    rawsigmerge,  1, 1),
          g.edge(rawsigmerge, l1sp_node,   0, 0),
          g.edge(l1sp_node,   l1sp_sink,   0, 0),
        ],
        oports=[sigsplit.oports[0]],
        name='sigproc_l1sppd_%d' % n
      ),

}
