// PDHD OSP

local pg = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// You won't find anything in this file explicitly mentioned in the config
// because OSP hard-wires the names.
local spfilt = import 'sp-filters.jsonnet';

local detectors = import "detectors.jsonnet";
local detname = "pdhd";
local det = detectors[detname];  // me

function(tpc)

    local resolution = tpc.adc.resolution;
    local fullscale = tpc.adc.fullscale[1] - tpc.adc.fullscale[0];
    local ADC_mV_ratio = ((1 << resolution) - 1 ) / fullscale;
    local anode = tpc.anode;
    pg.pnode({
        type: 'OmnibusSigProc',
        name: tpc.name,
        data: {
            anode: wc.tn(tpc.anode),
            dft: "FftwDFT",
            
            do_not_mp_protect_traditional: true, 
            field_response: wc.tn(tpc.fr),
            filter_responses_tn: [ ], // FIXME: this needs to be special for "bad APA1"
            elecresponse: wc.tn(tpc.er),
            ftoffset: 0.0,
            ctoffset: 1.0*wc.microsecond,
            per_chan_resp: "",
            postgain: 1.0,  // default 1.2
            
            ADC_mV: 1.0 / tpc.adc.lsb_voltage,
            lroi_rebin: 6,
            lroi_th_factor: 3.5,		
            lroi_th_factor1: 0.7,		
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

            tight_lf_tag: 'tight_lf%d' % anode.data.ident,
            loose_lf_tag: 'loose_lf%d' % anode.data.ident,
            cleanup_roi_tag: 'cleanup_roi%d' % anode.data.ident,
            // break_roi_loop1_tag: 'break_roi_1st%d' % anode.data.ident,
            // break_roi_loop2_tag: 'break_roi_2nd%d' % anode.data.ident,
            // shrink_roi_tag: 'shrink_roi%d' % anode.data.ident,
            shrink_roi_tag: "",
            // extend_roi_tag: 'extend_roi%d' % anode.data.ident,
            sparse: false,
            save_negtive_charge: false,

            use_multi_plane_protection: true,
            use_roi_debug_mode: true,   
            use_roi_refinement: true,   
            mp3_roi_tag: 'mp3_roi%d' % anode.data.ident,
            mp2_roi_tag: 'mp2_roi%d' % anode.data.ident,
            mp_tick_resolution: 10,
            // mp_th1: if anode.data.ident==0 then 200 else 1000,
            // mp_th2: if anode.data.ident==0 then 100 else 500,
            
            isWrapped: false,
            // process_planes: if anode.data.ident==0 then [0, 1] else [0, 1, 2],
            plane2layer: [0,1,2], // fixme: bad apa1 has special

            // fixme: "bad apa1" has special
            Wiener_tight_filters: ["Wiener_tight_U", "Wiener_tight_V", "Wiener_tight_W"],

        }
    }, nin=1, nout=1, uses=[anode, tpc.fr, tpc.er] + spfilt)

