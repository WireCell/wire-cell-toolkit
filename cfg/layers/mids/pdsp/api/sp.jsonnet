// The .sp sub-api for signal processing for PDSP.

// Note, we only add these as "uses" because OmnibusSigProc sigproc
// *hard codes* their type/names.....
local spfilt = import "sp-filters.jsonnet";

local low = import "../../../low.jsonnet";
local wc = low.wc;

local resps = import "resps.jsonnet";

// Allow an optional argument "sparse" as this is really an end-user
// decision.  Higher layers may expose this option to the TLA.
function(services, params, options={}) function(anode)
    local pars = std.mergePatch(params, std.get(options, "params", {}));
    local opts = {sparse:true} + options;
    local ident = low.util.idents(anode);
    local resolution = pars.digi.resolution;
    local fullscale = pars.digi.fullscale[1] - pars.digi.fullscale[0];
    local ADC_mV_ratio = ((1 << resolution) - 1 ) / fullscale;

    local res = resps(params).sp;

    low.pg.pnode({
        type: 'OmnibusSigProc',
        name: ident,
        data: {

            /**  
            *  Optimized SP parameters (May 2019)
            *  Associated tuning in sp-filters.jsonnet
            */
            anode: wc.tn(anode),
            dft: wc.tn(services.dft),
            per_chan_resp: "",
            field_response: wc.tn(res.fr),
            elecresponse: wc.tn(res.er),
            ftoffset: 0.0, // default 0.0
            ctoffset: 1.0*wc.microsecond, // default -8.0
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
            wiener_tag: 'wiener' + ident,
            decon_charge_tag: 'decon_charge' + ident,
            gauss_tag: 'gauss' + ident,

            use_roi_debug_mode: false,
            tight_lf_tag: 'tight_lf' + ident,
            loose_lf_tag: 'loose_lf' + ident,
            cleanup_roi_tag: 'cleanup_roi' + ident,
            break_roi_loop1_tag: 'break_roi_1st' + ident,
            break_roi_loop2_tag: 'break_roi_2nd' + ident,
            shrink_roi_tag: 'shrink_roi' + ident,
            extend_roi_tag: 'extend_roi' + ident,

            use_multi_plane_protection: false,
            mp3_roi_tag: 'mp3_roi' + ident,
            mp2_roi_tag: 'mp2_roi' + ident,
            
            isWrapped: false,
            sparse : opts.sparse,
        },
    }, nin=1, nout=1, uses=[anode, services.dft, res.fr, res.er] + spfilt)
