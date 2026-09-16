// Overide defaults for specific channels.  If an info is
// mentioned for a particular channel in multiple objects in this
// list then last mention wins.

local wc = import "wirecell.jsonnet";

local handmade = import "chndb-handmade.jsonnet";

// Call with the anode ident number
function(ident, nsamples)
[

    // First entry provides default channel info across ALL
    // channels.  Subsequent entries override a subset of channels
    // with a subset of these entries.  There's no reason to
    // repeat values found here in subsequent entries unless you
    // wish to change them.
    {
        channels: std.range(ident * 2560, (ident + 1) * 2560 - 1),
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

    },

    {
        //channels: { wpid: wc.WirePlaneId(wc.Ulayer) },
        channels: std.range(ident * 2560, ident * 2560 + 800- 1),
        freqmasks: [
            { value: 1.0, lobin: 0, hibin: nsamples - 1 },
            { value: 0.0, lobin: 169, hibin: 173 },
            { value: 0.0, lobin: 513, hibin: 516 },
        ],
        /// this will use an average calculated from the anode
        // response: { wpid: wc.WirePlaneId(wc.Ulayer) },
        /// this uses hard-coded waveform.
        response: { waveform: handmade.u_resp, waveformid: wc.Ulayer },
        response_offset: 120, // offset of the negative peak
        pad_window_front: 20,
        decon_limit: 0.02,
        decon_limit1: 0.07,
        roi_min_max_ratio: 3.0,
    },
    {
        //channels: { wpid: wc.WirePlaneId(wc.Vlayer) },
        channels: std.range(ident * 2560 + 800, ident * 2560 + 1600- 1),
        freqmasks: [
            { value: 1.0, lobin: 0, hibin: nsamples - 1 },
            { value: 0.0, lobin: 169, hibin: 173 },
            { value: 0.0, lobin: 513, hibin: 516 },
        ],
        /// this will use an average calculated from the anode
        // response: { wpid: wc.WirePlaneId(wc.Vlayer) },
        /// this uses hard-coded waveform.
        response: { waveform: handmade.v_resp, waveformid: wc.Vlayer },
        response_offset: 124,
        decon_limit: 0.01,
        decon_limit1: 0.08,
        roi_min_max_ratio: 1.5,
    },

    
    {
        //channels: { wpid: wc.WirePlaneId(wc.Wlayer) },
        channels: std.range(ident * 2560 + 1600, ident * 2560 + 2560- 1),
        nominal_baseline: 400.0,
        decon_limit: 0.05,
        decon_limit1: 0.08,
        // freqmasks: freqbinner.freqmasks(harmonic_freqs, 5.0*wc.kilohertz),
    },

]
    
