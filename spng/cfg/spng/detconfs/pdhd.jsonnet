local wc = import "wirecell.jsonnet";
local api = import "../detector.jsonnet";
local detectors = import "detectors.jsonnet";

local detname = "pdhd";
local det = detectors[detname];  // me

// Caution, PDHD hardware tick is 512ns but we use 500ns internally.  This
// assumes any real DAQ data is resampled.  Changing to 512ns REQUIRES a change
// to the FR among other things.
local adc_tick = 500*wc.ns;

local adc = api.adc(tick=adc_tick,
                    resolution=14, 
                    baselines=[1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt],
                    fullscale=[0.2*wc.volt, 1.6*wc.volt]);
local fr = api.fields_from_name(detname);
local er = api.elec_response(gain=14.0*wc.mV/wc.fC,
                             shaping=2.2*wc.us,
                             binning=api.binning(100, adc_tick));

local wires_obj = api.wires_from_name(detname);

local anodes = [
    api.anode(anode_ident, wires_obj, api.hd_like_faces(anode_ident))
    for anode_ident in [0,1,2,3]];

local gauss_filter = api.filter_config(scale=0.12 * wc.megahertz);
local gauss_filters = [
    gauss_filter,
    gauss_filter,
    gauss_filter,
];
local wiener_filters = [
    api.filter_config(scale=0.1487880 * wc.megahertz, power=3.76194),
    api.filter_config(scale=0.1596568 * wc.megahertz, power=4.36125),
    api.filter_config(scale=0.1362300 * wc.megahertz, power=3.35324),
];

local looself = api.filter_config(scale= 0.002 * wc.megahertz);
local dnnroi_filters = [
    looself,
    looself,
    looself,
];
    
local channel_filters = [
    api.filter_config(scale=1.0 / wc.sqrtpi * 0.75, period=1.0, ignore_baseline=false),
    api.filter_config(scale=1.0 / wc.sqrtpi * 0.75, period=1.0, ignore_baseline=false),
    api.filter_config(scale=1.0 / wc.sqrtpi * 10.0, period=1.0, ignore_baseline=false),
];

local filters = [
    api.view_filters(time_filters=api.time_filters(gauss_filters[i], wiener_filters[i], dnnroi_filters[i]),
                     channel_filters=api.channel_filters(channel_filters[i]))
    for i in [0,1,2]];

// All TPCs are identical except for their anodes
local tpcs = [ api.tpc(anode) for anode in anodes ];

local nominal_tpc = api.tpc(adc=adc, fr=fr, er=er, connections=[2,2,0], filters=filters, faces=[0,1]);

api.detector(api.override(nominal_tpc, tpcs))

