// Base channel noise DB object configuration for microboone
// This does not include any run dependent RMS cuts.
// See chndb.jsonnet

local wc = import "wirecell.jsonnet";
local handmade = import "chndb-resp.jsonnet";

function(params, anode, field, rms_cuts=[])
{
    anode: wc.tn(anode),
    field_response: wc.tn(field),

    tick: params.daq.tick,

    // This sets the number of frequency-domain bins used in the noise
    // filtering.  It is expected that time-domain waveforms have the
    // same number of samples.
    nsamples: params.nf.nsamples,

    // For MicroBooNE, channel groups is a 2D list.  Each element is
    // one group of channels which should be considered together for
    // coherent noise filtering.
    groups: [std.range(g*48, (g+1)*48-1) for g in std.range(0,171)],

    // Externally determined "bad" channels.
    bad: [446] + std.range(7136, 7199) + std.range(7201,7214) + std.range(7216,7263),

    // Overide defaults for specific channels.  If an info is
    // mentioned for a particular channel in multiple objects in this
    // list then last mention wins.
    channel_info: [             

        // First entry provides default channel info across ALL
        // channels.  Subsequent entries override a subset of channels
        // with a subset of these entries.  There's no reason to
        // repeat values found here in subsequent entries unless you
        // wish to change them.
        {
            channels: std.range(0, 2400 + 2400 + 3456 - 1),
            nominal_baseline: 2048.0,  // adc count
            gain_correction: 1.0,     // unitless
            response_offset: 0.0,      // ticks?
            pad_window_front: 10,     // ticks?
            pad_window_back: 10,      // ticks?
	    decon_limit: 0.02,
	    decon_lf_cutoff: 0.08,
	    decon_limit1: 0.09,
	    adc_limit: 15,
	    protection_factor: 5,
	    min_adc_limit: 50,
            min_rms_cut: 1.0,         // units???
            max_rms_cut: 5.0,         // units???

            // parameter used to make "rcrc" spectrum
            rcrc: 1.0*wc.millisecond,

            // parameters used to make "config" spectrum
            reconfig : {},

            // list to make "noise" spectrum mask
            freqmasks: [],

            // field response waveform to make "response" spectrum.  
            response: {},

        },

        {
            channels: {wpid: wc.WirePlaneId(wc.Ulayer)},
            freqmasks: [
                { value: 1.0, lobin: 0, hibin: $.nsamples-1 },
                { value: 0.0, lobin: 169, hibin: 173 },
                { value: 0.0, lobin: 513, hibin: 516 },
            ],
            /// this will use an average calculated from the anode
            // response: { wpid: wc.WirePlaneId(wc.Ulayer) },
            /// this uses hard-coded waveform.
            response: { waveform: handmade.u_resp, waveformid: wc.Ulayer },
            response_offset: 79,
            pad_window_front: 20,
	    decon_limit: 0.02,
	    decon_limit1: 0.09,
        },

        {
            channels: {wpid: wc.WirePlaneId(wc.Vlayer)},
            freqmasks: [
                { value: 1.0, lobin: 0, hibin: $.nsamples-1 },
                { value: 0.0, lobin: 169, hibin: 173 },
                { value: 0.0, lobin: 513, hibin: 516 },
            ],
            /// this will use an average calculated from the anode
            // response: { wpid: wc.WirePlaneId(wc.Vlayer) },
            /// this uses hard-coded waveform.
            response: { waveform: handmade.v_resp, waveformid: wc.Vlayer },
            response_offset: 82,
	    decon_limit: 0.025, 
	    decon_lf_cutoff: 0.06,
	    decon_limit1: 0.08,
	},

        {
            channels: {wpid: wc.WirePlaneId(wc.Wlayer)},
            nominal_baseline: 400.0,
	    decon_limit: 0.05,
	    decon_limit1: 0.08,
        },

        {                       // special channel
            channels: 2240,
            freqmasks: [
                { value: 1.0, lobin: 0, hibin: $.nsamples-1 },
                { value: 0.0, lobin: 169, hibin: 173 },
                { value: 0.0, lobin: 513, hibin: 516 },
                { value: 0.0, lobin:  17, hibin:  19 },
            ],
        },

        {                       // these are before hardware fix 
            channels: params.misconfigured.channels,
            reconfig: {
                from: {gain:  params.misconfigured.gain,
                       shaping: params.misconfigured.shaping},
                to:   {gain: params.elec.gain,
                       shaping: params.elec.shaping},
            }
        },
    ] + rms_cuts,
}
    
