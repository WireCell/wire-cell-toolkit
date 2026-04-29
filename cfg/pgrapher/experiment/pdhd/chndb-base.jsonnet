// Base channel noise DB object configuration for microboone
// This does not include any run dependent RMS cuts.
// See chndb.jsonnet

local handmade = import 'chndb-resp.jsonnet';
local wc = import 'wirecell.jsonnet';

// TODO (follow-up): decon_limit, decon_limit1, adc_limit, min_adc_limit and
// roi_min_max_ratio were tuned against the old SBND-copy kernel (peak ~±56 ADC).
// The new PDHD kernel has peak ~±206/247 ADC (~4× larger), so these thresholds
// are not yet re-optimised for PDHD.  Re-tune empirically once NF runs with
// the new kernel.

function(params, anode, field, n, rms_cuts=[], use_freqmask=true, coh_group_shift=3)
  // ADC-domain thresholds below are tuned for FE amplifier gain = 14 mV/fC.
  // For other gains, scale linearly with params.elec.gain.
  local gain_scale = params.elec.gain / (14.0 * wc.mV / wc.fC);
  // chndb-resp.jsonnet stores the FR⊗ER kernel at reference gain=14 mV/fC.
  // Scale element-wise so the kernel tracks the runtime FE gain.
  local scale_resp(arr) = std.map(function(x) x * gain_scale, arr);
  // Frequency-mask toggle threaded through wct-nf-sp.jsonnet's use_freqmask
  // TLA.  When false, all per-channel freqmasks below collapse to [], making
  // the C++ consumer in PDHD::OneChannelNoise a no-op for every channel.
  local freqmask_enabled = use_freqmask;
  {
    anode: wc.tn(anode),
    field_response: wc.tn(field),

    tick: 0.5*wc.us,  // NF always sees 500 ns frames (resampled from 512 ns on data path)

    // This sets the number of frequency-domain bins used in the noise
    // filtering.  It is not necessarily true that the time-domain
    // waveforms have the same number of ticks.  This must be non-zero.
    nsamples: params.nf.nsamples,

    // coh_group_shift: cyclic offset (in offline channels) applied to U and V
    // group boundaries.  Set to 0 to recover the original (pre-fix) definition.
    // Default 3 corrects the +3-channel FEMB-edge misassignment identified in
    // the 027409-evt0-apa0 coherent-noise audit (U/V only; W is unchanged).
    local shift = coh_group_shift,
    local u_group(u) = std.map(function(j) n * 2560 + std.mod(40 * u + shift + j, 800),
                               std.range(0, 39)),
    local v_group(v) = std.map(function(j) n * 2560 + 800 + std.mod(40 * v + shift + j, 800),
                               std.range(0, 39)),
    local w_group(w) = std.range(n * 2560 + 1600 + 48 * w, n * 2560 + 1600 + 48 * (w + 1) - 1),
    groups: [u_group(u) for u in std.range(0, 19)]
            + [v_group(v) for v in std.range(0, 19)]
            + [w_group(w) for w in std.range(0, 19)],


    // Externally determined "bad" channels.
    bad: [2297, 5379, 5472, 5556, 5607, 5608, 5920, 5921, 6072, 7099, 7288, 7679, 2580, 2940, 3347, 3758, 3805, 3866, 4722, 9956, 9986, 9987, 9988, 7876, 9120, 9125, 9126, 9127, 9306, 9307, 9309, 9310, 9534],

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
        channels: std.range(n * 2560, (n + 1) * 2560 - 1),
        nominal_baseline: 2048.0,  // adc count
        gain_correction: 1.0,  // unitless
        response_offset: 0.0,  // ticks?
        pad_window_front: 10,  // ticks?
        pad_window_back: 10,  // ticks?
        decon_limit: 0.02 * gain_scale,
        decon_limit1: 0.09 * gain_scale,
        adc_limit: 60 * gain_scale, // 15,
        min_adc_limit: 200 * gain_scale, // 50,
        roi_min_max_ratio: 0.8, // default 0.8
        min_rms_cut: 10.0 * gain_scale,  // ADC at 14 mV/fC
        max_rms_cut: 30.0 * gain_scale,  // ADC at 14 mV/fC

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
	channels: std.range(n * 2560, n * 2560 + 800- 1),
	// Previous content was U-plane notches at bins [169,173] (~57 kHz)
	// and [513,516] (~171 kHz).  Left empty pending re-analysis.  When
	// ready, populate with wc.freqmasks_phys([...freqs in wc units...], delta)
	// — bins are resolved at runtime from the live frame size and the
	// conjugate-mirror bins are applied automatically.
	freqmasks: [],
        /// this will use an average calculated from the anode
        // response: { wpid: wc.WirePlaneId(wc.Ulayer) },
        /// this uses hard-coded waveform.
        response: { waveform: scale_resp(handmade.u_resp), waveformid: wc.Ulayer },
        response_offset: 127, // argmin of PDHD FR⊗ER kernel (was 120, SBND copy)
        pad_window_front: 20,
        decon_limit: 0.02 * gain_scale,
        decon_limit1: 0.07 * gain_scale,
        roi_min_max_ratio: 3.0,
      },

      {
        //channels: { wpid: wc.WirePlaneId(wc.Vlayer) },
	channels: std.range(n * 2560 + 800, n * 2560 + 1600- 1),
        // Same situation as the U-plane entry above: pending re-analysis.
        // When ready, use wc.freqmasks_phys([...freqs in wc units...], delta)
        // — bins resolved at runtime, conjugate-mirror applied automatically.
        freqmasks: [],
        /// this will use an average calculated from the anode
        // response: { wpid: wc.WirePlaneId(wc.Vlayer) },
        /// this uses hard-coded waveform.
        response: { waveform: scale_resp(handmade.v_resp), waveformid: wc.Vlayer },
        response_offset: 132, // argmin of PDHD FR⊗ER kernel (was 124, SBND copy)
        decon_limit: 0.01 * gain_scale,
        decon_limit1: 0.08 * gain_scale,
        roi_min_max_ratio: 1.5,
      },

      // local freqbinner = wc.freqbinner(params.daq.tick, params.nf.nsamples);
      // local harmonic_freqs = [f*wc.kilohertz for f in
      //   // [51.5, 102.8, 154.2, 205.5, 256.8, 308.2, 359.2, 410.5, 461.8, 513.2, 564.5, 615.8]
      //   [51.5, 77.2, 102.8, 128.5, 154.2, 180.0, 205.5, 231.5, 256.8, 282.8, 308.2, 334.0, 359.2, 385.5, 410.5, 461.8, 513.2, 564.5, 615.8, 625.0]
      // ];
      
      {
        //channels: { wpid: wc.WirePlaneId(wc.Wlayer) },
	channels: std.range(n * 2560 + 1600, n * 2560 + 2560- 1),
        nominal_baseline: 400.0,
        decon_limit: 0.05 * gain_scale,
        decon_limit1: 0.08 * gain_scale,
        // freqmasks: freqbinner.freqmasks(harmonic_freqs, 5.0*wc.kilohertz),
      },

    ] + rms_cuts,
  }
