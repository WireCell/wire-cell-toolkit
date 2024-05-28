// Nominal is in name only.  See actual for more realistic.

local wc = import "wirecell.jsonnet";
local base_variants = import "../../base/variants.jsonnet";

base_variants.nominal {

    lar : {
        // Longitudinal diffusion constant
        DL : 6.4 * wc.cm2/wc.s,
        // Transverse diffusion constant
        DT : 9.8 * wc.cm2/wc.s,
        // Electron lifetime
        lifetime : 10.4*wc.ms,
        // given cathode is at 2.548 m, this speed leads to a consistent maximum
        // drifting time with that in data -- 2.560 m / 1.101 mm/us
        drift_speed : 1.098*wc.mm/wc.us,
    },

    // $ wirecell-util wires-volumes uboone \
    //   | jsonnetfmt -n 4 --max-blank-lines 1 - \
    //   > cfg/layers/mids/uboone/geometry.jsonnet
    geometry_data : import "geometry.jsonnet",
    geometry : $.geometry_data + {
        xplanes: {
            // Distances between collection plane and a logical plane.
            danode: 0.0,
            dresponse: 10*wc.cm,
            // Location of cathode measured from collection is based on a dump
            // of ubcore/v06_83_00/gdml/microboonev11.gdml by Matt Toups
            dcathode: 2.5480*wc.m,
        }
    },

    // The nominal time binning for data produced by the detector.
    binning : {
        tick: 0.5*wc.us,
        nticks: 8192,
    },

    elec : {   
        gain : 14.0*wc.mV/wc.fC,
        shaping : 2.2*wc.us,
        postgain: 1.2,
    },

    rc : {
        width: 1.0*wc.ms,
    },

    ductor: super.ductor {
        // These tickle the API
        transform: "standard",  // see overlay
        universe: "single",     // no super-position of different FRs

        // MB uses one or three field files.  The standard variable is
        // "field_file" (singular) but to simplify structure we define the full
        // set even for nominal.  For nominal ("single" universe) the first in
        // the list is used.
        field_files: $.detector_data.fields,

    },

    splat : {
        sparse: true,
        tick: $.ductor.binning.tick,
        window_start: $.ductor.start_time,
        window_duration: $.ductor.readout_time,
    },

    // Simulating noise
    noise : {
        model: {
            spectra_file: "microboone-noise-spectra-v2.json.bz2",
            // These are frequency space binning which are not necessarily
            // same as some time binning - but here they are.
            period: $.binning.tick,     // 1/frequency
            // The MicroBooNE noise filtering spectra runs with a
            // different number of bins in frequency space than the
            // nominal number of ticks in time.
            nsamples: 8192,
            // Optimize binning 
            wire_length_scale: 1.0*wc.cm,
        },
        // Optimize use of randoms
        replacement_percentage: 0.02,
    },

    digi : {
        // There are post-FE amplifiers.  Should this be
        // elec.postgain? (could be, but it's fine here.  There are 2
        // relative gains available).
        gain: 1.0,

        // The resolution (bits) of the ADC
        resolution: 12,

        // These values are pre-amp (ASIC) pedestals. Two additional AC coupling
        // would change the baseline and it also depends on the bias voltage from regulator.
        // +/-10 ADC? variation is expected.
        //baselines: [900*wc.millivolt,900*wc.millivolt,200*wc.millivolt],
        // Actuall MicroBooNE baselines are 2046 ADC for induction planes and 473 for collection.
        baselines: [999.3*wc.millivolt,999.3*wc.millivolt,231.02*wc.millivolt],

        fullscale: [0*wc.volt, 2.0*wc.volt],
    },

    nf: {
        field_file: $.detector_data.field, // "ub-10-half.json.bz2",
        binning: $.binning,
        chresp_file: $.detector_data.chresp,

        // The MicroBooNE noise filtering spectra runs with a
        // different number of bins in frequency space than the
        // nominal number of ticks (binning.nticks) in time.
        nsamples: 9592,

        // to change these, make new variants
        // wct, wcls, wcls_multi
        chndb_type: "wct",
        // dynamic, before, after, perfect
        chndb_epoch: "perfect",
    },        

    sp: {
        field_file: $.detector_data.field, // "ub-10-half.json.bz2",
    },        

    // FIXME: excise this from nominal.
    // Describe misconfigured channels.  Used by NF's chndb and sim's
    // static channel status.  Nominally, everything is perfect.  See
    // other variants for non-perfect content.
    misconfigured: {
        // The list of misconfigured channel IDs 
        channels: [],           
        // Misconfigured parameters, only relevant if there are
        // channels.  Nominal case is there is no misconfig so set
        // same as nominal config.
        gain: $.elec.gain,
        shaping: $.elec.shaping,
    },

    // Imaging paramter pack
    img : {
        charge_error_file: $.detector_data.qerr,
        span: 4,          // Number of ticks to collect into one time slice span
        binning : $.binning,
    },

}


