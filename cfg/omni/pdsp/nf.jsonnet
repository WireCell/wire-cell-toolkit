// Return a factory function for a noisefilter node.
// This is assumed to be used in the context of an omni object.

// CAVEAT: this does not try to include all the detail in
// pgrapher/experiment/pdsp!

local wc = import 'wirecell.jsonnet';

{
    noisefilter :: function(anode, binning=$.binning.sp)
        local chndb_data = {
            /// Electronics groups for coherent noise
            groups: [std.range(n * 2560 + u * 40, n * 2560 + (u + 1) * 40 - 1) for u in std.range(0, 19)]
                    + [std.range(n * 2560 + 800 + v * 40, n * 2560 + 800 + (v + 1) * 40 - 1) for v in std.range(0, 19)]
                    + [std.range(n * 2560 + 1600 + w * 48, n * 2560 + 1600 + (w + 1) * 48 - 1) for w in std.range(0, 19)],

            // In reality there is a long list of bad channels.
            bad: [],

            local n = anode.data.ident,

            // last match wins
            channel_info: [ {
                // first entry is default

                channels: {first: n * 2560, last: (n + 1) * 2560 - 1},
                nominal_baseline: 2048.0,  // adc count
                gain_correction: 1.0,  // unitless
                response_offset: 0.0,  // ticks?
                pad_window_front: 10,  // ticks?
                pad_window_back: 10,  // ticks?
                decon_limit: 0.02,
                decon_limit1: 0.09,
                adc_limit: 15,
                roi_min_max_ratio: 0.8, // default 0.8
                min_rms_cut: 1.0,  // units???
                max_rms_cut: 30.0,  // units???

                // parameter used to make "rcrc" spectrum
                rcrc: 1.1 * wc.millisecond, // 1.1 for collection, 3.3 for induction
                rc_layers: 1, // default 2

                // parameters used to make "config" spectrum
                reconfig: {},

                // list to make "noise" spectrum mask
                freqmasks: [],
                
                // field response waveform to make "response" spectrum.
                response: {},

            }, {
                channels: { wpid: wc.WirePlaneId(wc.Ulayer) },
                response: { wpid: wc.WirePlaneId(wc.Ulayer) },
                response_offset: 120, // offset of the negative peak
                pad_window_front: 20,
                decon_limit: 0.02,
                decon_limit1: 0.07,
                roi_min_max_ratio: 3.0,
            }, {
                channels: { wpid: wc.WirePlaneId(wc.Vlayer) },
                response: { wpid: wc.WirePlaneId(wc.Vlayer) },
                response_offset: 124,
                decon_limit: 0.01,
                decon_limit1: 0.08,
                roi_min_max_ratio: 1.5,
            }, {
                channels: { wpid: wc.WirePlaneId(wc.Wlayer) },
                response: { wpid: wc.WirePlaneId(wc.Wlayer) },
                nominal_baseline: 400.0,
                decon_limit: 0.05,
                decon_limit1: 0.08,
            }],                 // channel_info
        };
        local onf_data = {
            /// fixme: fill these in 
            channel_filters: [],
            grouped_filters: [],
            channel_status_filters: [],
        };
        super.noisefilter(anode, binning, chndb_data, onf_data),


}




