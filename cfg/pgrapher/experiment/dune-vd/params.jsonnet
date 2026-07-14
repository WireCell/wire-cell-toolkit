// DUNE-VD specific parameters.  This file inerets from the
// generic set of parameters and overrides things specific to PDSP.

local wc = import "wirecell.jsonnet";
local base = import "pgrapher/common/params.jsonnet";

function(params) base {
    // This section will be overwritten in simparams.jsonnet
    det : {

        // The current DUNE-VD goemetry has only one CRP composed by 36
        // independent CRM with side 
        // CRP is on y-z while drift is on x 
        // Only one CRP is defined in this geometry 
        // CRMs are oneside anodes     

        response_plane: params.response_plane,

        local upper_crp_x = 325.00*wc.cm, //300.507*wc.cm,
        local upper_resp_x = upper_crp_x-self.response_plane,
        local cathode_x = -325.00*wc.cm,
        local ncrm = if std.objectHas(params, 'ncrm') then params.ncrm else 36,

        local upper_face = { 
            anode:    upper_crp_x, 
            response: upper_resp_x, 
            cathode:  cathode_x},
       
        volumes: [
            {
                wires: n,       // anode number
                name: "crm%d"%n,
                faces: [ upper_face, upper_face ],
            } for n in std.range(0, ncrm-1)], // std.range is inclusive, i.e. [0, crm-1],
    },

    daq: super.daq {

        tick: 0.5*wc.us, // check this in the TDR, LArSoft
        nticks: params.nticks, //9375, // 1.6 mm/us per 0.5 us assuming 6000 mm drift leght. 

        //readout_time: self.tick*self.nticks,
        //nreadouts: 1,
        //start_time: 0.0*wc.s,
        //stop_time: self.start_time + self.nreadouts*self.readout_time,
        //first_frame_number: 0,
    },

    adc: super.adc {
        
        // ADC resolution (bits): restores the default removed from the
        // shared base params in commit 41e02736 (pre-41e02736 inherited 12).
        resolution: 12,

        // Set 0 for now
        //baselines: [0*wc.millivolt, 0*wc.millivolt, 0*wc.millivolt],
        //resolution: 12,
        //fullscale: [0.2*wc.volt, 1.6*wc.volt],

        // Copied from pdsp. induction plane: 2350 ADC, collection plane: 900 ADC
        // baselines: [1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt],
        // fullscale: [0.2*wc.volt, 1.6*wc.volt],

    },

    // Electronics response synced to protodunevd top-drift elec.  This
    // single-CRP geometry (one upper CRP) is treated as the top drift.
    elec: super.elec {
        // gain is a jsonnet-required placeholder: JsonElecResponse loads
        // its response from file.
        gain: 14*wc.mV/wc.fC,
        type: "JsonElecResponse",
        filename: "dunevd-coldbox-elecresp-top-psnorm_400.json.bz2",
        postgain: 1.36, // 11mV/fC, 1.94 -> 14mV/fC
    },

    sim: super.sim {

        // For running in LArSoft, the simulation must be in fixed time mode. 
        fixed: true,

    },

    overall_short_padding: 0.2*wc.ms,
    sys_status: false,
    sys_resp: {
        start: 0.0 * wc.us,
        magnitude: 1.0,
        time_smear: 1.0 * wc.us,
    },

    files: {

        // Standard wire geometry with 2 wire planes and third dummy induction
        wires: "dunevd10kt-1x6x6-3view30deg-wires-v1.json.bz2",

        // Field response synced to protodunevd (CRP, 18.1 cm response plane).
        fields: [
            "protodunevd_FR_imbalance3p_260501.json.bz2",
        ],

        // Electronics-noise spectra synced to protodunevd top-drift spectra.
        noise: "pdvd-top-noise-spectra-v3.json.bz2",

        chresp: null,

    },

}
