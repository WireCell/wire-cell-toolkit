// microboone signal processing.
//
// This file MUST work for any MB setup.  Any variance should into a
// parameter in variants/*.jsonnet.

local spfilt = import "details/sp-filters.jsonnet";

local low = import "../../../low.jsonnet";
local wc = low.wc;
local pg = low.pg;

local frs = import "frs.jsonnet";


function(services, params, options={}) function(anode)
    local fr = frs(params).nf;

    local cer = params.ductor.binning {
        type: "ColdElecResponse",
        data: params.elec,
    };

    local perchanresp = {
        type: "PerChannelResponse",
        data: {
            filename: params.nf.chresp_file,
        }
    };

    local sigproc_perchan = pg.pnode({
        type: "OmnibusSigProc",
        data: {
            // This class has a HUGE set of parameters.  See
            // OmnibusSigProc.h for the list.  For here, for now, we
            // mostly just defer to the hard coded configuration values.
            // They can be selectively overriddent.  This class also hard
            // codes a slew of SP filter component names which MUST
            // correctly match what is provided in sp-filters.jsonnet.
            anode: wc.tn(anode),
            dft: wc.tn(services.dft),
            field_response: wc.tn(fr),
            elecresponse: wc.tn(cer),
            postgain: 1,  // default 1.2
                          // see elec postgain in params.jsonnet
            per_chan_resp: wc.tn(perchanresp),
            fft_flag: 0,   // 1 is faster but higher memory, 0 is slightly slower but lower memory
        }
    }, nin=1,nout=1, uses=[anode, services.dft, fr, cer, perchanresp] + spfilt);

    local sigproc_uniform = pg.pnode({
        type: "OmnibusSigProc",
        data: {
            anode: wc.tn(anode),
            dft: wc.tn(services.dft),
            field_response: wc.tn(fr),
            elecresponse: wc.tn(cer),
            postgain: 1,  // default 1.2
                          // see elec postgain in params.jsonnet
            per_chan_resp: "",
            shaping: params.elec.shaping,
            fft_flag: 0,    // 1 is faster but higher memory, 0 is slightly slower but lower memory
            // r_fake_signal_low_th: 300,
            // r_fake_signal_high_th: 600,
        }
    }, nin=1,nout=1,uses=[anode, services.dft, fr, cer] + spfilt);

    // ch-by-ch response correction in SP turn off by setting null input
    local sigproc = if std.type(params.nf.chresp_file)=='null'
                    then sigproc_uniform
                    else sigproc_perchan;

    // The L1 SP is a follow-on frame filter but it only takes channels
    // that are considered to be subject to the shorted wire region fields
    // and it also needs both raw and nominal SP waveforms.  A fiddly sub
    // graph is needed to route everything properly.

    local rawsplit = pg.pnode({
        type: "FrameSplitter",
        name: "rawsplitter"
    }, nin=1, nout=2);

    local sigsplit = pg.pnode({
        type: "FrameSplitter",
        name: "sigsplitter"
    }, nin=1, nout=2);

    local chsel = pg.pnode({
        type: "ChannelSelector",
        data: {
            // channels that will get L1SP applied
            channels: std.range(3566,4305),

            // can pass on only the tags of traces that are actually needed.
            tags: ["raw","gauss"]
        }
    }, nin=1, nout=1);
    local l1spfilter = pg.pnode({
        type: "L1SPFilter",
        data: {
            dft: wc.tn(services.dft),
            fields: wc.tn(fr),
            filter: [0.000305453, 0.000978027, 0.00277049, 0.00694322, 0.0153945,
                     0.0301973, 0.0524048, 0.0804588, 0.109289, 0.131334, 0.139629,
                     0.131334, 0.109289, 0.0804588, 0.0524048, 0.0301973, 0.0153945,
                     0.00694322, 0.00277049, 0.000978027, 0.000305453],
            raw_ROI_th_nsigma: 4.2,
            raw_ROI_th_adclimit:  9,
            overall_time_offset : 0,
            collect_time_offset : 3.0,
            roi_pad: 3,
            raw_pad: 15,
            adc_l1_threshold: 6,
            adc_sum_threshold: 160,
            adc_sum_rescaling: 90,
            adc_ratio_threshold: 0.2,
            adc_sum_rescaling_limit : 50,
            l1_seg_length : 120,
            l1_scaling_factor : 500,
            l1_lambda : 5,
            l1_epsilon : 0.05,
            l1_niteration : 100000,
            l1_decon_limit : 100,
            l1_resp_scale : 0.5,
            l1_col_scale : 1.15,
            l1_ind_scale : 0.5,
            peak_threshold : 1000,
            mean_threshold : 500,
            adctag: "raw",                             // trace tag of raw data
            sigtag: "gauss",                           // trace tag of input signal
            outtag: "l1sp",                            // trace tag for output signal
        }
    }, nin=1, nout=1, uses=[services.dft, fr]);

    // merge the split output from NF ("raw" tag) and just the "gauss"
    // from normal SP for input to L1SP
    local rawsigmerge = pg.pnode({
        type: "FrameMerger",
        name: "rawsigmerge",
        data: {
            rule: "replace",

            // note: the first two need to match the order of what data is
            // fed to ports 0 and 1 of this component in the pgraph below!
            mergemap: [
                ["raw","raw","raw"],
                ["gauss","gauss","gauss"],
            ],
        }
    }, nin=2, nout=1);

    // merge out of L1 ("l1sp" tag) with previously split output from
    // regular SP ("wiener" and "gauss" tags).  The "raw" tag is unlikely
    // to be found.
    local l1merge = pg.pnode({
        type: "FrameMerger",
        name: "l1merge",
        data: {
            rule: "replace",

            // note: the first two need to match the order of what data is
            // fed to ports 0 and 1 of this component in the pgraph below!
            mergemap: [
                ["raw","raw","raw"],
                ["l1sp","gauss","gauss"],
                ["l1sp","wiener","wiener"],
            ],
        }
    }, nin=2, nout=1);

    pg.intern([rawsplit], [l1merge], [sigproc, sigsplit, chsel, l1spfilter, rawsigmerge, l1merge],
             edges=[
                 pg.edge(rawsplit, sigproc),
                 pg.edge(sigproc, sigsplit),
                 pg.edge(sigsplit, rawsigmerge),
                 pg.edge(sigsplit, l1merge, 1, 1),

                 pg.edge(rawsplit, rawsigmerge, 1, 1),
                 pg.edge(rawsigmerge, chsel),
                 pg.edge(chsel, l1spfilter),
                 pg.edge(l1spfilter, l1merge),
             ],
             name="L1SP")
