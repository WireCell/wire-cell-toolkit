// MicroBooNE-specific parameters.  This file produces a parameters
// object, inheriting from common.params.  

local wc = import "wirecell.jsonnet";
local base = import "pgrapher/common/params.jsonnet";

base {
    lar: super.lar {
        DL : 6.4 * wc.cm2/wc.s,
        DT : 9.8 * wc.cm2/wc.s,
        // given cathode is at 2.548 m, this speed leads to a consistent maximum
        // drifting time with that in data -- 2.560 m / 1.101 mm/us
        drift_speed : 1.098*wc.mm/wc.us, 
    },

    daq: super.daq {
        // The real detector DAQ produces an odd number of ticks.
        // Note, this is different from nsamples defined below!
        // It's different still from nticks used in simulation.
        nticks: 9595,
    },        
    adc: super.adc {
        // There are post-FE amplifiers.  Should this be
        // elec.postgain? (could be, but it's fine here.  There are 2
        // relative gains available).
        gain: 1.0,

        // ADC resolution (bits) -- the MicroBooNE TPC digitises to 12 bits.
        // Previously inherited from the common base params, whose default
        // was removed, so it is set explicitly here.
        resolution: 12,

        // fixme: need double checking
        baselines: [999.3*wc.millivolt,999.3*wc.millivolt,231.02*wc.millivolt],

        // fixme: need double checking
        fullscale: [0*wc.volt, 2.0*wc.volt],
    },

    elec : super.elec {   
        gain : 14.0*wc.mV/wc.fC,
        shaping : 2.2*wc.us,
        postgain: 1.2,

    },

    det : {

        volumes: [
            // Note:
            // $ wirecell-util wires-info microboone-celltree-wires-v2.1.json.bz2
            // anode:0 face:0 X=[-6.00,0.00]mm Y=[-1155.30,1174.70]mm Z=[0.35,10369.60]mm
            //      0: x=-0.00mm dx=6.0000mm
            //      1: x=-3.00mm dx=3.0000mm
            //      2: x=-6.00mm dx=0.0000mm
            // $ wirecell-sigproc response-info ub-10-half.json.bz2 
            // origin:10.00 cm, period:0.10 us, tstart:0.00 us, speed:1.11 mm/us, axis:(1.00,0.00,0.00)
            //      plane:0, location:6.0000mm, pitch:3.0000mm
            //      plane:1, location:3.0000mm, pitch:3.0000mm
            //      plane:2, location:0.0000mm, pitch:3.0000mm

            {
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
                        response: $.elec.fields.start_dx - 6*wc.mm,
                        // Location of cathode measured from
                        // collection is based on a dump of
                        // ubcore/v06_83_00/gdml/microboonev11.gdml by
                        // Matt Toups
                        cathode: 2.5480*wc.m,
                    },
                    
                    null
                ],
            },
        ],

        // The active volume bounding box defined as the location of
        // two extreme corners.  See above note for where these magic
        // numbers come from.
        bounds : {
            tail: wc.point(  -6.0, -1155.30,    0.35, wc.mm),
            head: wc.point(2560.4,  1174.7, 10369.65, wc.mm),
        },

        
    },
    nf: super.nf {

        // The MicroBooNE noise filtering spectra runs with a
        // different number of bins in frequency space than the
        // nominal number of ticks in time.
        nsamples: 9592,
        
        // Add some extra parameters that only make sense for other
        // structures built in uboone/.

        // Also, add this wart which separates "before" and "after"
        // hardware noise fix epochs.  This is used when WCT is run
        // from WCLS.
        run12boundary: 7000,  

        // These are used in chndb and also static channel status db for simulation.
        misconfigured: {
            channels: std.range(2016,2095) + std.range(2192,2303) + std.range(2352,2399),
            gain: 4.7*wc.mV/wc.fC,
            shaping: 1.1*wc.us,
        },

    },

    // Don't bother adding sim overrides here, instead see and use simparams.jsonnet.
    // sim: super.sim {}

    files : {
        wires:"microboone-celltree-wires-v2.1.json.bz2",
        fields:["ub-10-half.json.bz2",
                "ub-10-uv-ground-tuned-half.json.bz2",
                "ub-10-vy-ground-tuned-half.json.bz2"],
        noise: "microboone-noise-spectra-v2.json.bz2",
        chresp: "microboone-channel-responses-v1.json.bz2",
    },

}


