local wc = import "wirecell.jsonnet";
local api = import "../detector.jsonnet";
local detectors = import "detectors.jsonnet";
local detname = "pdhd";
local det = detectors[detname];  // me

local osp = import "pdhd/osp.jsonnet";

// Caution, PDHD hardware tick is 512ns but we use 500ns internally.  This
// assumes any real DAQ data is resampled.  Changing to 512ns REQUIRES a change
// to the FR among other things.
local adc_tick = 500*wc.ns;

// These were taken from dunereco a0029f0fd0ec8821429abe568338da84fbc06057
local lar = api.lar(DT = 8.8 * wc.cm2 / wc.s,
                    DL = 4.0 * wc.cm2 / wc.s,
                    lifetime = 35.0 * wc.ms,
                    drift_speed = 1.60563* wc.mm / wc.us);


local readout_time = 3*wc.ms;
local tick0_time = -250*wc.us;
//local tick0_time = 0*wc.us;
local response_plane = 10*wc.cm; // relative to collection wires
local response_time_offset = response_plane / lar.drift_speed;
local response_duration = readout_time + response_time_offset;

local response_start_time = tick0_time - response_time_offset;

// How much to roll the response deconvolution. This essentially the ADC tick
// where FR*ER goes to zero.
local decon_roll = 128;

local ductor = api.ductor(adc_tick, response_duration, response_start_time);
local splat = api.splat(
    smear_long= [
        2.691862363980221,
        2.6750200122535057,
        2.7137567141154055
    ],
    smear_tran = [
        0.7377218875719689,
        0.7157764520393882,
        0.13980698710556544
    ]);

local adc = api.adc(tick=adc_tick,
                    resolution=14, 
                    baselines=[1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt],
                    fullscale=[0.2*wc.volt, 1.6*wc.volt],
                    readout_duration=3.0*wc.ms);
local fr = api.fields_from_name(detname);
// SPNG ER duration covers where ER is nonzero
local er_spng = api.elec_response(gain=14.0*wc.mV/wc.fC,
                                  shaping=2.2*wc.us,
                                  binning=api.binning(100, adc_tick));
// SIM ER duration is ridiculously chosen to be the entire readout  which makes no sense.
local er_sim = api.elec_response(gain=14.0*wc.mV/wc.fC,
                                 shaping=2.2*wc.us,
                                 binning=api.binning(adc.readout_nticks, adc_tick));

// rcs...

local wires_obj = api.wires_from_name(detname);

local anodes = [
    api.anode(anode_ident, wires_obj, api.hd_like_faces(anode_ident))
    for anode_ident in [0,1,2,3]];

// Check, eg, FrameToTdm logs to assure correct channel/wire ordering.
local view_groups = [
    // The wrapped U
    api.view_group(0, 2, [1, 0], [-1, 1]),
    // The wrapped V
    api.view_group(1, 2, [1, 0], [-1, 1]),
    // One W face
    api.view_group(2, 0, [0,]),
    // Other W face
    api.view_group(2, 0, [1,]),
];

local pirs(anode) = [
    api.plane_impact_response("", plane,
                              tick=adc.tick,
                              nticks=adc.readout_nticks,
                              fr=fr, er=er_sim,
                              rcs=[])
    for plane in [0,1,2]];

local noise = api.noise(empirical = api.empirical_noise(det.noise));

// Same for all views?
local gauss_filter = api.filter_axis([
    api.filter_function(scale=0.12 * wc.megahertz, power=2, kind="lowpass"),
]);

// by view
local hf_tight = [
    // PDSP Wiener_tight_{U,V,W}
    api.filter_function(scale=0.1487880 * wc.megahertz, power=3.76194),
    api.filter_function(scale=0.1596568 * wc.megahertz, power=4.36125),
    api.filter_function(scale=0.1362300 * wc.megahertz, power=3.35324),
];
local hf_wide = [
    // PDSP's Wiener_wide_{U,V,W}
    api.filter_function(scale=0.186765 * wc.megahertz, power=5.05429),
    api.filter_function(scale=0.193600 * wc.megahertz, power=5.77422),
    api.filter_function(scale=0.175722 * wc.megahertz, power=4.37928),
];
local lf_loose =   api.filter_function(0.002 * wc.megahertz, kind="highpass");
local lf_tight =   api.filter_function(0.016 * wc.megahertz, kind="highpass");
local lf_tighter = api.filter_function(0.080 * wc.megahertz, kind="highpass");

local wiener_filters = [
    api.filter_axis([hf_tight[0], lf_tighter]),
    api.filter_axis([hf_tight[1], lf_tighter]),
    api.filter_axis([hf_tight[2]])
];

// same for all views?
local dnnroi_filters = [
    api.filter_axis([hf_tight[0], lf_loose]),
    api.filter_axis([hf_tight[1], lf_loose]),
    api.filter_axis([hf_tight[2]])
];
    
local channel_filters = [
    api.filter_axis([api.filter_function(scale=1.0 / wc.sqrtpi * 0.75)],
                    period=1.0, ignore_baseline=false),
    api.filter_axis([api.filter_function(scale=1.0 / wc.sqrtpi * 0.75)],
                     period=1.0, ignore_baseline=false),
    api.filter_axis([api.filter_function(scale=1.0 / wc.sqrtpi * 10.0)],
                    period=1.0, ignore_baseline=false)
];

local filters = [
    api.view_filters(time_filters=api.time_filters(gauss=gauss_filter,
                                                   wiener=wiener_filters[i],
                                                   dnnroi=dnnroi_filters[i],
                                                   options={
                                                       roll: decon_roll,
                                                       crop: adc.readout_nticks,
                                                   }),
                     channel_filters=api.channel_filters(channel_filters[i]))
    for i in [0,1,2]];

// OSP defaults are 3/5 sigma for ind/col plus a nominal 1.0 (no units) 
local cvt_ind = api.crossview_threshold(rms_nsigma=3.0, nominal=1);
local cvt_col = api.crossview_threshold(rms_nsigma=5.0, nominal=1);
local cvts = api.crossview_thresholds(cvt_ind, cvt_ind, cvt_col);

// All TPCs are identical except for their anodes
local tpcs = [
    api.tpc(anode, lar=lar, ductor=ductor, splat=splat, adc=adc, fr=fr, er=er_spng,
            pirs=pirs(anode), noise=noise,
            view_groups=view_groups, filters=filters,
            crossview_thresholds=cvts,
            osp_subgraph=osp)
    for anode in anodes ];

api.detector("pdhd", tpcs)

