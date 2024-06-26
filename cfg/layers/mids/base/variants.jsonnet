// The base/variants.jsonnet provides base for "nominal" and "actual" parameters.
//
// To define a detector variant param object through INHERITANCE do like:
//
//   // In <detector>/variants.jsonnet
//   local base = import "../../base/variants.jsonnet";
//   local nominal = import "variants/nominal.jsonnet";
//   local actual = import "variants/actual.jsonnet";
//   {
//       nominal: base.nominal + nominal,
//       actual: self.nominal + actual,
//   }
//
// Or to use replacement semantics, change the body to:
//
//   {
//       nominal: std.mergePatch(base.nominal, nominal),
//       actual: std.mergePatch(self.nominal, actual),
//   }
//
// That is, <detector>/variants/nominal.jsonnet overrides some parts of
// base.nominal and then further <detector>/variants/actual.jsonnet overrides
// some parts of the detector's nominal.

local low = import "../../low.jsonnet";
local wc = low.wc;
local detectors = import "detectors.jsonnet";

local nominal = {

    // Note, when constructed via high.params() these may be overridden before
    // the resulting params object is passed to the detector mid api.
    detector_name: "base",
    variant_name: "nominal",
    structure_name: "nominal",
    
    // The high level object summarizing the contents of wire-cell-data
    // corresponding to a canonical detector name.
    detector_data: detectors[self.detector_name],

    // The nominal base uses the same sample time domain binning everywhere.
    // It provides the default binning in more local contexts (such as sim/sp).
    binning : {
        tick: 500*wc.ns,
        nticks: 4096,       // 2ms is enough for anyone!
    },

    // Nominal LAr properties.  Almost certainly a detector should override
    // some or all of these.
    lar : {
        // Longitudinal diffusion constant
        DL : 4.0 * wc.cm2/wc.s,
        // Transverse diffusion constant
        DT : 8.8 * wc.cm2/wc.s,
        // Electron lifetime
        lifetime : 10*wc.ms,
        // Electron drift speed, assumes a certain applied E-field
        drift_speed : 1.6*wc.mm/wc.us, // approx at 500 V/cm
        // LAr density
        density: 1.389*wc.g/wc.centimeter3,
        // Decay rate per mass for natural Ar39.
        ar39activity: 1*wc.Bq/wc.kg,
    },

    
    // This substructure describes the detector geometry.  Almost certainly the
    // geometry.volumes must be supplied by a specific detector.
    //
    // It's main chunk is the "volumes" substructure which must be consistent
    // with the location of the wires in the "wires file".  This can be
    // challenging to get correct by hand.  It can partly be derived from the
    // wires file but that is too slow to do every load into WCT so it is
    // recommended to generate it once like:
    //
    //   $ wirecell-util wires-volumes <detector_name> | jsonnetfmt - > geometry.jsonnet
    //
    // The "anode", "response" and "cathode" planes can be given on that command
    // line but here we override it to match other values from the config.
    geometry_data: import "variants/geometry.jsonnet",
    geometry : std.mergePatch($.geometry_data, {
        xplanes: {
            // Distances between collection plane and a logical plane.
            danode: 0.0,
            dresponse: 10*wc.cm,
            dcathode: nominal.binning.nticks * nominal.binning.tick / nominal.lar.drift_speed,
        }
    }),


    // The electronics response parameters.  Almost certainly, the inheriting
    // detector should override this section.
    elec : {
        // The type of electronics
        type: "cold",

        // The FE amplifier gain in units of Voltage/Charge.
        gain : 14.0*wc.mV/wc.fC,

        // The shaping (aka peaking) time of the amplifier shaper.
        shaping : 2.0*wc.us,

        // A realtive gain applied after shaping and before ADC.
        // 
        // Don't use this to fix wrong sign depos.  If you neeed to
        // fix sign of eg larsoft depos its input component has a
        // "scale" parameter which can be given a -1.  Or better, fix
        // the C++ that sets the wrong sign to begin with.
        postgain: 1.0,
    },

    // The "RC" response.  The detector should override and some may require
    // DFP structural changes to implement RCRC
    rc : {
        width: 1.0*wc.ms,
    },

    // The ductor transforms drifted depos to currents
    ductor: {

        // The (single) field file.  
        field_file: $.detector_data.field,

        // Note, an extension to support "actual" detector may require multiple
        // fields, as in the case of uboone.  The above variable should remain
        // singular.  Access multiple files through another variable such as:
        // field_files: $.detector_data.fields, // (plural)

        // Nominally, we assume simulation time starts at t=0.  However, any
        // ionization deposited between "response" and "anode" planes are
        // "backed up" in space and time.  In order to not trim them out we must
        // start the simulation earlier so that a depo at t=0 at the "anode"
        // plane will have its response in tick >= 0.
        start_time : -($.geometry.xplanes.dresponse - $.geometry.xplanes.danode)/$.lar.drift_speed,

        // Most simply, no custom binning.  
        tbin: 0,
        binning: $.binning,
        readout_time : self.binning.nticks * self.binning.tick,
    },

    // Simulating noise
    noise : {
        model: {
            spectra_file: $.detector_data.noise,

            // These are frequency space binning which are not necessarily
            // same as some time binning - but here they are.
            period: $.binning.tick,     // 1/frequency
            nsamples: $.binning.nticks, // 
            // Optimize binning 
            wire_length_scale: 1.0*wc.cm,
        },
        // Optimize use of randoms
        replacement_percentage: 0.02,
    },

    // Digitization model parameters.  
    digi : {
        // A relative gain applied just prior to digitization.  This
        // is not FE gain, see elec for that.
        gain: 1.0,

        // Voltage baselines added to any input voltage signal listed
        // in a per plan (U,V,W) array.
        baselines: [1*wc.volt, 1*wc.volt, 0.5*wc.volt],

        // The resolution (bits) of the ADC
        resolution: 12,

        // The voltage range as [min,max] of the ADC, eg min voltage
        // counts 0 ADC, max counts 2^resolution-1.
        fullscale: [0*wc.volt, 2*wc.volt],
    },

    // Configure the "splat" (DepoFluxSplat) is an approximation to the
    // combination of ductor+sigproc.
    splat : {
        sparse: true,
        tick: $.binning.tick,
        window_start: 0,
        window_duration: self.tick * $.binning.nticks,
        reference_time: 0.0,
    },

    // The noise filter parameter pack
    nf : {
        // In principle, may use a different field response.
        field_file : $.ductor.field_file,

        binning: $.binning,

        // The name of the ChannelNoiseDb
        chndb: "perfect",
        filters: {
            channel: ["single"], // sticky, single, gaincalib
            grouped: [],         // grouped,
            status: [],
        },
    },

    // Signal processing parameter pack
    sp: {
        // In principle, may use a different field response.
        field_file : $.ductor.field_file,
        binning : $.binning,
    },
};

{
    nominal: nominal,
    actual: nominal,
} + low.ssss(nominal)
