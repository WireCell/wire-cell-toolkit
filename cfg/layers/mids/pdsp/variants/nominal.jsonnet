// This is the file cfg/layers/mids/pdsp/variant/nominal.jsonnet.
// If this file is found at any other path, I will delete it without notice.

// This is the BASE configuration for *pdsp* variant parameter objects.  It
// yields an object of a schema expected throughout cfg/layers/mids/pdsp/.

// Do NOT change this file.

// If the nominal parameter do not suit your needs, follow these steps:
//
// 1. Pick a single, short name that describes your needs.
//
// 2. Start an EMPTY file in cfg/layers/mids/pdsp/variant/<name>.jsonnet
//
// 3. Import the nominal.jsonnet object
//
// 4. Yield the nominal object with any values overridden.
//
// 5. Add this new object to the variants object in mids.jsonnet.

local wc = import "wirecell.jsonnet";
local base_variants = import "../../base/variants.jsonnet";

base_variants.nominal {

    lar : super.lar {
        // Electron lifetime
        lifetime : 10.4*wc.ms,
        // Electron drift speed, assumes a certain applied E-field
        drift_speed : 1.565*wc.mm/wc.us, // at 500 V/cm
    },

    geometry_data : import "geometry.jsonnet",

    binning : {
        tick: 0.5*wc.us,
        nticks: 6000,
    },

    // The electronics response parameters.
    elec : {
        // The FE amplifier gain in units of Voltage/Charge.
        gain : 14.0*wc.mV/wc.fC,

        // The shaping (aka peaking) time of the amplifier shaper.
        shaping : 2.2*wc.us,

        // A realtive gain applied after shaping and before ADC.
        // 
        // Don't use this to fix wrong sign depos.  If you neeed to
        // fix sign of eg larsoft depos its input component has a
        // "scale" parameter which can be given a -1.  Or better, fix
        // the C++ that sets the wrong sign to begin with.
        // 
        // Pulser calibration: 41.649 ADC*tick/1ke
        // theoretical elec resp (14mV/fC): 36.6475 ADC*tick/1ke
        postgain: 1.1365,
    },

    // The "RC" response
    rc : {
        width: 1.1*wc.ms,
    },


    // Digitization model parameters.  
    digi : {
        // A relative gain applied just prior to digitization.  This
        // is not FE gain, see elec for that.
        gain: 1.0,

        // Voltage baselines added to any input voltage signal listed
        // in a per plan (U,V,W) array.
        //
        // Per tdr, chapter 2
        // induction plane: 2350 ADC, collection plane: 900 ADC
        baselines: [1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt],

        // The resolution (bits) of the ADC
        resolution: 12,

        // The voltage range as [min,max] of the ADC, eg min voltage
        // counts 0 ADC, max counts 2^resolution-1.
        //
        // The tdr says, "The ADC ASIC has an input
        // buffer with offset compensation to match the output of the
        // FE ASIC.  The input buffer first samples the input signal
        // (with a range of 0.2 V to 1.6 V)..."
        fullscale: [0.2*wc.volt, 1.6*wc.volt],
    },

    // The noise filter parameter pack
    nf : {
        // In principle, may use a different field response.
        field_file : $.ductor.field_file,

        binning: $.binning,

        // The name of the ChannelNoiseDb
        chndb: "actual",
        filters: {
            channel: ["single"], // sticky, single, gaincalib
            grouped: [],         // grouped,
            status: [],
        },
    },


}
    
