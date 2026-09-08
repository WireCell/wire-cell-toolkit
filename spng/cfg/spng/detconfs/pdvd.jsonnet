// SPNG ProtoDUNE-VD (pdvd) detector configuration.
//
// Modeled on spng/detconfs/pdhd.jsonnet.  Oracle for all values:
//   reference/dunereco/dunereco/DUNEWireCell/protodunevd/  (params/sp-filters/sp)
//
// PDVD has TWO drift volumes that share one CRP wire geometry but differ in
// electronics response, ADC scaling and noise:
//
//   anode ident 0..3  -> "bottom" drift: analytic ColdElecResponse, 7.8 mV/fC;
//                        ADC baselines [1003.4,1003.4,507.7] mV, fullscale [0.2,1.6] V.
//   anode ident 4..7  -> "top"    drift: measured JsonElecResponse (file), pg 1.36;
//                        ADC baselines [1,1,1] V, fullscale [0,2] V.
//
// This heterogeneity does NOT require a structural change: api.tpc() already
// accepts per-TPC adc/er/fr/filters, and every SPNG job selects a single
// det.tpcs[i] (one TPC = one drift half) via its tpcid.  So each of the 8 TPCs
// simply carries the adc/er for its half.  Choose tpcid in 0..3 for bottom,
// 4..7 for top.
//
// VD wire layout: each anode has two "faces" (two Y-halves of one CRP) with the
// same drift direction ("jumpered", connection=1 for the U/V induction views;
// the W collection view is per-face, connection=0).  There is NO plane swap and
// the anode is NOT wrapped (contrast PDHD).
//
// CAVEATS (documented, need validation against a running FrameToTdm dump):
//  - view_group face order/signs for the jumpered U/V views are a best guess.
//  - decon_roll is derived from response_plane/drift_speed (not independently
//    tuned as PDHD's 129 was).
//  - top-drift electronics response is represented as a JsonElecResponse; SPNG
//    ResponseKernel support for file-based ER at runtime is unverified.
//  - splat smear values are placeholders (no pdvd morse tuning in the oracle).

local wc = import "wirecell.jsonnet";
local api = import "../detector.jsonnet";
local detectors = import "detectors.jsonnet";
local detname = "pdvd";
local det = detectors[detname];  // me

local osp = import "pdvd/osp.jsonnet";

// PDVD works at 500ns internally (raw electronics 512ns are resampled).
local adc_tick = 500 * wc.ns;

// LAr (from protodunevd simparams oracle).
local lar = api.lar(DT = 8.8 * wc.cm2 / wc.s,
                    DL = 4.0 * wc.cm2 / wc.s,
                    lifetime = 1000.0 * wc.ms,
                    drift_speed = 1.473 * wc.mm / wc.us);

local readout_time = 3 * wc.ms;
local tick0_time = -250 * wc.us;
// VD response plane is 18.1 cm from the collection wires (vs PDHD 10 cm).
local response_plane = 18.1 * wc.cm;
local response_time_offset = response_plane / lar.drift_speed;
local response_duration = readout_time + response_time_offset;
local response_start_time = tick0_time - response_time_offset;

// How much to roll the response deconvolution: ~ the ADC tick where FR*ER goes
// to zero, i.e. the response-plane transit in ticks.  FIXME: not independently
// tuned (PDHD hard-codes 129); validate against the FR*ER kernel.
local decon_roll = std.floor(response_time_offset / adc_tick);

local ductor = api.ductor(adc_tick, response_duration, response_start_time);

// No pdvd morse tuning exists; placeholders (only used by depo-splat sim, not by
// the SP / checker path).  FIXME: run wirecell-gen morse-* for real values.
local splat = api.splat(smear_long = [2.7, 2.7, 2.7],
                        smear_tran = [0.75, 0.75, 0.15]);

// Per-half ADC.  resolution/tick shared; baselines+fullscale differ.
local adc_bot = api.adc(tick=adc_tick, resolution=14,
                        baselines=[1003.4*wc.millivolt, 1003.4*wc.millivolt, 507.7*wc.millivolt],
                        fullscale=[0.2*wc.volt, 1.6*wc.volt],
                        readout_duration=readout_time);
local adc_top = api.adc(tick=adc_tick, resolution=14,
                        baselines=[1.0*wc.volt, 1.0*wc.volt, 1.0*wc.volt],
                        fullscale=[0.0*wc.volt, 2.0*wc.volt],
                        readout_duration=readout_time);
// readout_nticks is identical for both halves (same readout_duration).
local readout_nticks = adc_bot.readout_nticks;

// One field response, shared by both halves (oracle uses the same FR file for
// the bottom and top pir entries).
local fr = api.fields_from_name(detname);

// Per-half electronics response.
//  bottom: analytic ColdElecResponse (gain 7.8 mV/fC, shaping 2.2 us).
//  top:    measured JSON shape.  api.elec_response_file() drops the filename, so
//          build the JsonElecResponse directly to keep it faithful.
local er_bot(binning) = api.elec_response(gain=7.8*wc.mV/wc.fC,
                                          shaping=2.2*wc.us,
                                          postgain=1.0,
                                          binning=binning);
local er_top(binning) = {
    type: "JsonElecResponse",
    name: "",
    data: {
        filename: "dunevd-coldbox-elecresp-top-psnorm_400.json.bz2",
        postgain: 1.36,
    } + api.binning_to_time(binning),
};
// SPNG ER binning covers only where ER is nonzero; SIM ER binning spans readout.
local er_spng_binning = api.binning(100, adc_tick);
local er_sim_binning = api.binning(readout_nticks, adc_tick);
local er_spng(ident) = if ident < 4 then er_bot(er_spng_binning) else er_top(er_spng_binning);
local er_sim(ident)  = if ident < 4 then er_bot(er_sim_binning)  else er_top(er_sim_binning);

local wires_obj = api.wires_from_name(detname);

// VD anode faces: two aligned (same-drift) faces per anode.  Both faces of an
// anode share identical anode/response/cathode X planes (vertical-drift CRP).
local apa_cpa = 341.55 * wc.cm;
local cpa_thick = 50.8 * wc.mm;
local apa_w2w = 85.725 * wc.mm;
local apa_g2g = 114.3 * wc.mm;
local apa_plane = 0.5 * apa_g2g;               // grid-wire plane offset
local res_plane = 0.5 * apa_w2w + response_plane;
local cpa_plane = apa_cpa - 0.5 * cpa_thick;
local vd_faces(ident) =
    local sign = if ident > 3 then 1 else -1;  // top:+1 (ident 4..7), bottom:-1
    local centerline = sign * apa_cpa;
    local face = {
        anode:    centerline - sign * apa_plane,
        response: centerline - sign * res_plane,
        cathode:  centerline - sign * cpa_plane,
    };
    [face, face];   // both faces active and identical (VD)

local anodes = [
    api.anode(anode_ident, wires_obj, vd_faces(anode_ident))
    for anode_ident in [0, 1, 2, 3, 4, 5, 6, 7]];

// View groups: U,V are jumpered across the two faces (connection=1); W is a
// single collection plane per face (connection=0).  FIXME: the face order and
// per-face order signs are a best guess pending a FrameToTdm channel-order dump.
local view_groups = [
    api.view_group(0, 1, [0, 1], [1, 1]),  // jumpered U
    api.view_group(1, 1, [0, 1], [1, 1]),  // jumpered V
    api.view_group(2, 0, [0]),             // W face 0
    api.view_group(2, 0, [1]),             // W face 1
];
// Ordering after view groups are concatenated back into 3 per-view tensors.
local view_wpids = [
    api.view_group(0, 1, [0, 1], [1, 1]),
    api.view_group(1, 1, [0, 1], [1, 1]),
    api.view_group(2, 0, [0, 1]),          // W0 + W1 concatenation
];

local pirs(anode) = [
    api.plane_impact_response("", plane,
                              tick=adc_tick,
                              nticks=readout_nticks,
                              fr=fr, er=er_sim(anode.data.ident),
                              rcs=[])
    for plane in [0, 1, 2]];

local noise = api.noise(empirical = api.empirical_noise(det.noise));

// --- SPNG-native time/channel filters (values from sp-filters oracle, _b==_t) ---

// Gaussian smoothing (final signal filter).
local gauss_filter = api.filter_axis([
    api.filter_function(scale=0.12 * wc.megahertz, power=2, kind="lowpass"),
]);

// Wiener "tight" per view (U,V,W).
local hf_tight = [
    api.filter_function(scale=0.148788  * wc.megahertz, power=3.76194),
    api.filter_function(scale=0.1596568 * wc.megahertz, power=4.36125),
    api.filter_function(scale=0.13623   * wc.megahertz, power=3.35324),
];
// Wiener "wide" per view (kept for completeness / parity with PDHD).
local hf_wide = [
    api.filter_function(scale=0.186765 * wc.megahertz, power=5.05429),
    api.filter_function(scale=0.1936   * wc.megahertz, power=5.77422),
    api.filter_function(scale=0.175722 * wc.megahertz, power=4.37928),
];
local lf_loose   = api.filter_function(0.003 * wc.megahertz, kind="highpass");
local lf_tight   = api.filter_function(0.014 * wc.megahertz, kind="highpass");
local lf_tighter = api.filter_function(0.060 * wc.megahertz, kind="highpass");

local wiener_filters = [
    api.filter_axis([hf_tight[0], lf_tighter]),
    api.filter_axis([hf_tight[1], lf_tighter]),
    api.filter_axis([hf_tight[2]]),
];
local dnnroi_filters = [
    api.filter_axis([hf_tight[0], lf_loose]),
    api.filter_axis([hf_tight[1], lf_loose]),
    api.filter_axis([hf_tight[2]]),
];

// Channel (wire) filters used in decon: Wire_ind for U,V; Wire_col for W.
local channel_filters = [
    api.filter_axis([api.filter_function(scale=1.0 / wc.sqrtpi * 5.0)],
                    period=0.5, ignore_baseline=false),
    api.filter_axis([api.filter_function(scale=1.0 / wc.sqrtpi * 5.0)],
                    period=0.5, ignore_baseline=false),
    api.filter_axis([api.filter_function(scale=1.0 / wc.sqrtpi * 10.0)],
                    period=0.5, ignore_baseline=false),
];

local filters = [
    api.view_filters(time_filters=api.time_filters(gauss=gauss_filter,
                                                   wiener=wiener_filters[i],
                                                   dnnroi=dnnroi_filters[i],
                                                   options={
                                                       roll: decon_roll,
                                                       crop: readout_nticks,
                                                   }),
                     channel_filters=api.channel_filters(channel_filters[i]))
    for i in [0, 1, 2]];

// OSP defaults 3/5 sigma for ind/col plus nominal 1.0.
local cvt_ind = api.crossview_threshold(rms_nsigma=3.0, nominal=1);
local cvt_col = api.crossview_threshold(rms_nsigma=5.0, nominal=1);
local cvts = api.crossview_thresholds(cvt_ind, cvt_ind, cvt_col);

// 8 TPCs, each carrying the adc/er for its drift half.
local tpcs = [
    local ident = anode.data.ident;
    local isbot = ident < 4;
    api.tpc(anode, lar=lar, ductor=ductor, splat=splat,
            adc=if isbot then adc_bot else adc_top,
            fr=fr, er=er_spng(ident),
            pirs=pirs(anode), noise=noise,
            view_groups=view_groups,
            view_wpids=view_wpids,
            filters=filters,
            crossview_thresholds=cvts,
            osp_subgraphs=osp)
    for anode in anodes];

api.detector("pdvd", tpcs)
