// PDVD OSP "mirror" config: the SPNG-side OmnibusSigProc used by adc-to-osp and
// by the check-dunereco-config faithfulness checker.  Modeled on the SPNG PDHD
// mirror (spng/detconfs/pdhd/sp.jsonnet); values follow the official oracle
// reference/dunereco/dunereco/DUNEWireCell/protodunevd/sp.jsonnet.
//
// PDVD split: anode ident 0..3 = bottom drift ("_b" filters, ColdElec 7.8mV/fC),
// ident 4..7 = top drift ("_t" filters, JSON elec response).  The suffix and the
// per-half electronics response (tpc.er) are provided by the detconf per TPC, so
// this file selects filter names by (anode.data.ident < 4).

local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// OSP hard-wires filter instance names; keep this include so they are built.
local spfilt = import 'sp-filters.jsonnet';

function(tpc)

    local anode = tpc.anode;
    local n = anode.data.ident;
    local sfx = if n < 4 then '_b' else '_t';

    // Per-job dump prefix injected onto tpc by detector.subset(); fall back to
    // the legacy fixed name when absent.
    local dump_prefix = if std.objectHasAll(tpc, 'sp_dump_prefix')
                        then tpc.sp_dump_prefix else 'osp_dump';

    pg.pnode({
        type: 'OmnibusSigProc',
        name: tpc.name,
        data: {
            anode: wc.tn(anode),
            dft: 'FftwDFT',
            dump_2d_spectra: false,
            dump_2d_prefix: dump_prefix,
            do_not_mp_protect_traditional: true,
            field_response: wc.tn(tpc.fr),
            filter_responses_tn: [],
            elecresponse: wc.tn(tpc.er),  // per-half (bottom ColdElec / top JSON)
            ftoffset: 0.0,
            ctoffset: 4.0 * wc.microsecond,  // PDVD (cf. PDHD 1us); see oracle sp.jsonnet
            per_chan_resp: '',

            ADC_mV: 1.0 / tpc.adc.lsb_voltage,
            postgain: 1.0,
            fft_flag: 0,
            troi_col_th_factor: 5.0,
            troi_ind_th_factor: 3.0,
            lroi_rebin: 6,
            lroi_th_factor: 3.5,
            lroi_th_factor1: 0.7,
            lroi_jump_one_bin: 1,

            // PDVD's sp-filters define ONLY the "_b"/"_t" suffixed instances (no
            // bare names), so the OSP C++ MUST be told every filter it uses by
            // name -- the compiled-in defaults ("Gaus_wide", "ROI_loose_lf", ...)
            // do not exist here.  These selectors are load-bearing.
            Gaus_wide_filter: 'Gaus_wide' + sfx,
            ROI_tight_lf_filter: 'ROI_tight_lf' + sfx,
            ROI_tighter_lf_filter: 'ROI_tighter_lf' + sfx,
            ROI_loose_lf_filter: 'ROI_loose_lf' + sfx,
            Wiener_wide_filters: ['Wiener_wide_U' + sfx,
                                  'Wiener_wide_V' + sfx,
                                  'Wiener_wide_W' + sfx],
            Wire_filters: ['Wire_ind' + sfx, 'Wire_ind' + sfx, 'Wire_col' + sfx],

            r_th_factor: 3.0,
            r_fake_signal_low_th: 375,
            r_fake_signal_high_th: 750,
            r_fake_signal_low_th_ind_factor: 1.0,
            r_fake_signal_high_th_ind_factor: 1.0,
            r_th_peak: 3.0,
            r_sep_peak: 6.0,
            r_low_peak_sep_threshold_pre: 1200,

            // frame tags
            wiener_tag: 'wiener%d' % n,
            decon_charge_tag: 'decon_charge%d' % n,
            gauss_tag: 'gauss%d' % n,
            tight_lf_tag: 'tight_lf%d' % n,
            loose_lf_tag: 'loose_lf%d' % n,
            cleanup_roi_tag: 'cleanup_roi%d' % n,
            break_roi_loop1_tag: 'break_roi_1st%d' % n,
            break_roi_loop2_tag: 'break_roi_2nd%d' % n,
            shrink_roi_tag: 'shrink_roi%d' % n,
            extend_roi_tag: 'extend_roi%d' % n,
            rawdecon_tag: '',
            sparse: false,
            save_negtive_charge: false,

            use_multi_plane_protection: true,
            use_roi_debug_mode: true,
            use_roi_refinement: true,
            mp3_roi_tag: 'mp3_roi%d' % n,
            mp2_roi_tag: 'mp2_roi%d' % n,
            mp_tick_resolution: 10,

            isWrapped: false,
            // PDVD has no plane-swap: natural [U, V, W] layout (contrast PDHD
            // APA0's [0,2,1] bad-APA trick).
            plane2layer: [0, 1, 2],

            Wiener_tight_filters: ['Wiener_tight_U' + sfx,
                                   'Wiener_tight_V' + sfx,
                                   'Wiener_tight_W' + sfx],
        },
    }, nin=1, nout=1, uses=[anode, tpc.fr, tpc.er] + spfilt)
