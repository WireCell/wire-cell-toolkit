local wc = import "wirecell.jsonnet";
local api = import "../detector.jsonnet";
local detectors = import "detectors.jsonnet";

local detname = "pdhd";
local det = detectors[detname];  // me

local adc = api.adc(resolution=14, 
                    baselines=[1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt],
                    fullscale=[0.2*wc.volt, 1.6*wc.volt]);
local fr = api.fields_from_detector(detname);
local er = api.elec_response(gain=14.0*wc.mV/wc.fC,
                             shaping=2.2*wc.us); // caution, uses default 500ns tick

local wires_obj = api.wires_from_detector(detname);

local anodes = [
    api.anode(anode_ident, wires_obj, api.hd_like_faces(anode_ident))
    for anode_ident in [0,1,2,3]];

local view_layers = [
    api.view_layer(0, connect=2),
    api.view_layer(1, connect=2),
    api.view_layer(2, connect=0),
];

local tpcs = [
    api.tpc(anode, adc, view_layers, fr, er)
    for anode in anodes ];

api.detector(tpcs)

