// This sort of follows
// uboonedata/WireCellData/pgrapher/experiment/uboone/*.jsonnet
// v08_00_00.  However, it assumes a "pure" WCT implementation.  See
// the "overlay" variant for what is needed to configure the official
// MicroBooNE "overlay simulation".

local wc = import "wirecell.jsonnet";

{
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

    geometry: {
        wires_file: "microboone-celltree-wires-v2.1.json.bz2",
        drifts: [{
            wires: 0,
            name: "uboone",
            faces: [
                {
                    // drop any depos w/in this plane.  The exact
                    // choice represents some trade off in
                    // approximations.
                    anode: 0.0*wc.mm, 
                    // plane, arbitrary choice.  Microboone wires
                    // put collection plane at absolute x=-6mm,
                    // response.plane_dx is measured relative to
                    // collection plane wires.
                    response: $.ductor.response_plane - 6*wc.mm,
                    // Location of cathode measured from
                    // collection is based on a dump of
                    // ubcore/v06_83_00/gdml/microboonev11.gdml by
                    // Matt Toups
                    cathode: 2.5480*wc.m,
                },
                null
            ],
        }]
    },

    // The nominal time binning for data produced by the detector.
    binning : {
        tick: 0.5*wc.us,
        // Real detector makes 9595.
        nticks: 9595,
    },

    elec : {   
        gain : 14.0*wc.mV/wc.fC,
        shaping : 2.2*wc.us,
        postgain: 1.2,
    },

    rc : {
        width: 1.0*wc.ms,
    },

    // These parameters only make sense for running WCT simulation on
    // microboone in larsoft.  The "trigger" section is not
    // "standard".  This section is just a set up for use below in
    // "sim".  There is no trigger, per se, in the simulation but
    // rather a contract between the generators of energy depositions
    // (ie, LarG4) and the drift and induction simulation (WCT).  For
    // details of this contract see:
    // https://microboone-docdb.fnal.gov/cgi-bin/private/ShowDocument?docid=12290
    trigger : {

        // A hardware trigger occurs at some "absolute" time but near
        // 0.0 for every "event".  It is measured in "simulation" time
        // which is the same clock used for expressing energy
        // deposition times.  The following values are from table 3 of
        // DocDB 12290.
        hardware_times: {
            none: 0.0,
            
            // BNB hardware trigger time.  Note interactions
            // associated with BNB neutrinos should all be produced
            // starting in the beam gate which begins at 3125ns and is
            // 1600ns in duration.
            bnb : -31.25*wc.ns,

            // Like above but for NUMI.  It's gate is 9600ns long starting
            // at 4687.5ns.
            numi : -54.8675*wc.ns,
            ext : -414.0625*wc.ns,
            mucs: -405.25*wc.ns,
        },
        hw_time: self.hardware_times["bnb"],

        // Measured relative to the hardware trigger time above is a
        // time offset to the time that the first tick of the readout
        // should sample.  This is apparently fixed for all hardware
        // trigger types (?).
        time_offset: -1.6*wc.ms,

        time: self.hw_time + self.time_offset,
    },

    ductor: {

        transform: "standard",  // see overlay

        // We trifurcate the detector to deal with different
        // problematic wire regions.
        field_files: ["ub-10-half.json.bz2",
                      "ub-10-uv-ground-tuned-half.json.bz2",
                      "ub-10-vy-ground-tuned-half.json.bz2"],

        response_plane : 10*wc.cm,


        local tzero = 0*wc.us,
        drift_dt : self.response_plane / $.lar.drift_speed,
        tbin : wc.roundToInt(self.drift_dt / self.binning.tick),
        binning: {
            tick : $.binning.tick,

            nticks : $.ductor.tbin + $.binning.nticks,
        },
        start_time : tzero - self.drift_dt + $.trigger.time,
        readout_time : self.binning.nticks * self.binning.tick,
    },

    // nominal: "sys_status" is FALSE and there is no "sys_resp".
    
    reframer : {
        tbin: (81*wc.us)/($.binning.tick),
        nticks: $.binning.nticks,
        toffset: $.ductor.drift_dt - 81*wc.us,
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
            nsamples: 9592,
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
        field_file: "ub-10-half.json.bz2",
        binning: $.binning,
        chresp_file: "microboone-channel-responses-v1.json.bz2",

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

    sp: {
        field_file: "ub-10-half.json.bz2",
    },        

    // Imaging paramter pack
    img : {
        charge_error_file: "microboone-charge-error.json.bz2",

        // Number of ticks to collect into one time slice span
        span: 4,

        binning : $.binning,

        // fixme: remove old tiling strategy "perfect" and add slicing strategies
    },

}


