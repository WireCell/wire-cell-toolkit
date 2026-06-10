// Base channel noise DB object configuration for microboone
// This does not include any run dependent RMS cuts.
// See chndb.jsonnet

local handmade = import 'chndb-resp.jsonnet';
local coh_groups = import 'pdhd-coh-groups.jsonnet';
local femb_params = import 'pdhd-femb-negpulse-groups.jsonnet';
local wc = import 'wirecell.jsonnet';

// TODO (follow-up): decon_limit, decon_limit1, adc_limit, min_adc_limit and
// roi_min_max_ratio were tuned against the old SBND-copy kernel (peak ~±56 ADC).
// The new PDHD kernel has peak ~±206/247 ADC (~4× larger), so these thresholds
// are not yet re-optimised for PDHD.  Re-tune empirically once NF runs with
// the new kernel.

function(params, anode, field, n, rms_cuts=[], use_freqmask=true)
  // ADC-domain thresholds (adc_limit, min/max_rms_cut) are tuned at 14 mV/fC
  // and scale with gain_scale for other gains.  Deconvolved-domain thresholds
  // (decon_limit, decon_limit1) operate on the gain-normalised output and do
  // not scale with FE gain.
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

    // Coherent-noise groups: 60 per anode (20 U-FEMB + 20 V-FEMB + 20 W-FEMB),
    // derived directly from PD2HDChannelMap_WIBEth_visiblewires_v1.txt
    // (DUNE/duneprototypes commit c8f43809).  The per-APA +/-3 wire rotation
    // is already baked into the channel map (the June 30 2025 sign fix on
    // the cathode-crossing correction), so no `coh_group_shift` parameter
    // is needed here.  See pdhd-coh-groups.jsonnet.
    groups: coh_groups.groups[n],

    femb_negpulse_groups: femb_params.groups,

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
        pad_window_front: 20,  // ticks?
        pad_window_back: 20,  // ticks?
        decon_limit: 0.02,
        decon_limit1: 0.09,
        adc_limit: 30 * gain_scale, // 15,
        min_adc_limit: 100 * gain_scale, // 50,
        roi_min_max_ratio: 0.8, // default 0.8
        min_rms_cut: 10.0 * gain_scale,  // ADC at 14 mV/fC
        // min_rms_cut: 8.0 * gain_scale,  // ADC at 14 mV/fC
        max_rms_cut: 50.0 * gain_scale,  // ADC at 14 mV/fC
        // max_rms_cut: 35.0 * gain_scale,  // ADC at 14 mV/fC

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
        decon_limit: 0.01,
        decon_limit1: 0.07,
        roi_min_max_ratio: 0.8,
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
        // APA 0 (n==0) V plane is hardware-faulty and behaves as a collection
        // plane.  Drop the FR⊗ER kernel and zero the response_offset there so
        // PDHD::SignalProtection's deconvolution gate falls through (the gate
        // requires respec.size()>0, respec[0]!=(1,0), and res_offset!=0).
        response: if n == 0 then {} else { waveform: scale_resp(handmade.v_resp), waveformid: wc.Vlayer },
        response_offset: if n == 0 then 0 else 132, // argmin of PDHD FR⊗ER kernel (was 124, SBND copy)
        decon_limit: 0.01,
        decon_limit1: 0.07,
        roi_min_max_ratio: 0.8,
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
        decon_limit: 0.05,
        decon_limit1: 0.08,
        // freqmasks: freqbinner.freqmasks(harmonic_freqs, 5.0*wc.kilohertz),
      },

      // W-plane harmonic noise diagnosed from run 027409 evt 0 (anodes 1 & 3).
      // Two interleaved combs are present:
      //   Comb A: f0 = 18.17 kHz, harmonics h2..h18 (36.34..327.06 kHz)
      //   Comb B: f0 = 28.97 kHz, harmonics h1..h13 (28.97..376.6 kHz)
      // Both are notched on the same channel set.  Comb B was identified
      // post-NF on ch4160/5108/5118/5119/9280/9759/9760/9761 etc. as a
      // residual 29-kHz-spaced peak structure not covered by comb A; a
      // bin-spacing analysis across all listed channels showed best-fit
      // f0 = 28.97 kHz with harmonics 1-12 confirmed at SNR>=4σ.
      // Affected channels: anode1 {4160, 4626-4639, 5107-5119} and
      //                    anode3 {9280-9295, 9758-9780}.
      // Only active when freqmask_enabled=true (use_freqmask TLA).
      // {
      //   channels: if n == 1 then [4160] + std.range(4626, 4639) + std.range(5107, 5119)
      //             else if n == 3 then std.range(9280, 9296) + std.range(9758, 9781)
      //             else [],
      //   freqmasks: if freqmask_enabled && (n == 1 || n == 3) then
      //     wc.freqmasks_phys(
      //       [k * 18.17 * wc.kilohertz for k in std.range(2, 18)] +
      //       [k * 28.97 * wc.kilohertz for k in std.range(1, 13)],
      //       1.0*wc.kilohertz)
      //     else [],
      // },

      // ch10000-10047 (anode 3, W-plane): separate ~20-22 kHz doublet
      // interference absent in all other groups.  Re-measured peaks at
      // 20.7 and 22.0 kHz (previous config used 21.3/22.7 which missed
      // both peaks under ±0.5 kHz notches).
      // {
      //   channels: if n == 3 then std.range(10000, 10048) else [],
      //   freqmasks: if freqmask_enabled && n == 3 then
      //     wc.freqmasks_phys(
      //       [ 20.7*wc.kilohertz,
      //         22.0*wc.kilohertz ],
      //       0.5*wc.kilohertz)
      //     else [],
      // },

      // APA2 (n==2) U/V plane harmonic combs diagnosed from run 027409 evt 0.
      // Two distinct combs, each shared between a U-plane group and a V-plane
      // group (same combs on both planes -> electronics-level pickup):
      //   Comb A: f0 = 16.43 kHz, h2..h12 (32.86..197.16 kHz).  h3 (~49 kHz)
      //           and h5 (~82 kHz) dominate.
      //   Comb B: f0 = 25.57 kHz, h1..h3  (25.57, 51.14, 76.71 kHz).  h3
      //           strongest.
      // Carriers:
      //   Comb A: U ch5202-5203, V ch6243-6282
      //   Comb B: U ch5560-5603, V ch6641-6719
      // G3 (V 5961-5965) had a wide bump 50-90 kHz that is not a discrete
      // comb (peak walks 59->69 kHz channel-to-channel) and is left to the
      // dynamic noisy-channel mask.
      // {
      //   channels: if n == 2 then std.range(5202, 5203) + std.range(6243, 6282)
      //             else [],
      //   freqmasks: if freqmask_enabled && n == 2 then
      //     wc.freqmasks_phys(
      //       [k * 16.43 * wc.kilohertz for k in std.range(2, 12)],
      //       1.0*wc.kilohertz)
      //     else [],
      // },
      // {
      //   channels: if n == 2 then std.range(5560, 5603) + std.range(6641, 6719)
      //             else [],
      //   freqmasks: if freqmask_enabled && n == 2 then
      //     wc.freqmasks_phys(
      //       [k * 25.57 * wc.kilohertz for k in std.range(1, 3)],
      //       1.0*wc.kilohertz)
      //     else [],
      // },

    ] + rms_cuts,
  }
