// This is the file cfg/layers/mids/pdsp/variant/simple.jsonnet.
// If this file is found at any other path, I will delete it without notice.

// This files defines a simple VARIANT for pdsp which is NOT MEANT TO BE
// CORRECT.  Do not "improve" it.  It is meant to provide simple parameters to
// make debugging easier.

// If you require a new variant, read and follow the comments at the top of
// "nominal.jsonnet".

local wc = import "wirecell.jsonnet";

local detectors = import "detectors.jsonnet";
local mydet = detectors.pdsp;
local nominal = import "nominal.jsonnet";

local simple = {
    // Define the LAr properties.  In principle, this is NOT subject
    // to variance though in principle it could be (eg, top/bottom of
    // DUNE-VD may have different LAr temps?).  We include this here
    // to keep data out of the the mid code and because this info is
    // used in a few places.
    lar : {
        // Longitudinal diffusion constant
        DL : 4.0 * wc.cm2/wc.s,
        // Transverse diffusion constant
        DT : 8.8 * wc.cm2/wc.s,
        // Electron lifetime
        lifetime : 10.4*wc.ms,
        // Electron drift speed, assumes a certain applied E-field
        // See https://github.com/WireCell/wire-cell-toolkit/issues/266
        // We do not pick a "simple" value (1.6) in order that we match
        // what is *likely* in the FR file. 
        drift_speed : 1.565*wc.mm/wc.us, // at 500 V/cm
        // LAr density
        density: 1.389*wc.g/wc.centimeter3,
        // Decay rate per mass for natural Ar39.
        ar39activity: 1*wc.Bq/wc.kg,
    },

    // Extract all the location for planes used in geometry
    planes : {
        apa_cpa: 3.63075*wc.m, // LArSoft
        // apa_cpa = 3.637*wc.m, // DocDB 203
        cpa_thick: 3.175*wc.mm, // 1/8", from Bo Yu (BNL) and confirmed with LArSoft
        // cpa_thick = 50.8*wc.mm, // DocDB 203
        apa_w2w: 85.725*wc.mm, // between W planes, DocDB 203
        apa_g2g: 114.3*wc.mm, 
        // measured from grid to "anode" plane along drift direction
        offset: 4.76*wc.mm
    },

    // This function returns the minimal geometry needed by WCT.  In
    // the transverse plane the geometry is defined by the "wires
    // file" and in the drift direction by the "drifts" structure
    // which defines a trio of planes: "anode", "response" and
    // "cathode" which are perpendicular to drift.  Each are assumed
    // approximately equipotential along their plane.  If your
    // detector fails this, then use a "variant params" pack for each
    // detector region which does meet this requirement.
    // 
    geometry : {

        // The exhaustive list of the location of every single wire
        // (or strip) segment.
        // wires_file: "protodune-wires-larsoft-v4.json.bz2",
        wires_file: mydet.wires,

        // Comments on how to chose the "anode" plane location: The
        // "anode" cut off plane, here measured from APA centerline,
        // determines how close to the wires do we consider any depo.
        // Anything closer will simply be discarded, else it will
        // either be drifted or "backed up" to the response plane.
        // This is somewhat arbitrary choice.  The "backing up" assume
        // parallel drift points and so some small error is made when
        // applied to ionization produce very close to the wires.  To
        // reduce this small error the "anode" plane can be moved
        // closer to the "response" plane but at the cost of not
        // simulate the volume left beyond the "anode" plane.

        // local apa_plane = 0.5*apa_g2g, // put anode plane at grid plane
        local apa_plane = 0.5*$.planes.apa_g2g - $.planes.offset, // put anode plane at first induction plane

        // Response plane MUST be located so that the wire plane
        // locations used in the FR calculation matches the wire plane
        // locations given to WCT.
        local res_plane = 0.5*$.planes.apa_w2w + $.ductor.response_plane,

        // The cathode plane is like the anode cut off plane.  Any
        // depo not between the two is dropped prior to drifting.  It
        // is demarks the face of the CPA.
        local cpa_plane = $.planes.apa_cpa - 0.5*$.planes.cpa_thick,

        // The drifts are then defined in terms of these above
        // numbers.  You can use "wirecell-util wires-info" or
        // "wirecell-util wires-volumes" or others to understand the
        // mapping of anode number to the 6 locations in global X and
        // Z coordinates.  For Larsoft wires the numbering is column
        // major starting at small X and Z so the centerline is
        // -/+/-/+/-/+.  Also important is that the faces are listed
        // "front" first.  Front is the one with the more positive X
        // coordinates and if we want to ignore a face it is made
        // null.
        drifts : [{
            local sign = 2*(n%2)-1,
            local centerline = sign*$.planes.apa_cpa,
            wires: n,       // anode number
            name: "apa%d"%n,
            faces:
                // top, front face is against cryo wall
                if sign > 0
                then [
                    {
                        anode: centerline + apa_plane,
                        response: centerline + res_plane,
                        cathode: centerline + cpa_plane, 
                    },
                    {
                        anode: centerline - apa_plane,
                        response: centerline - res_plane,
                        cathode: centerline - cpa_plane, 
                    }
                ]
                // bottom, back face is against cryo wall
                else [
                    {
                        anode: centerline + apa_plane,
                        response: centerline + res_plane,
                        cathode: centerline + cpa_plane, 
                    },
                    {
                        anode: centerline - apa_plane,
                        response: centerline - res_plane,
                        cathode: centerline - cpa_plane, 
                    }

                ],
        } for n in std.range(0,5)],
    },

    // The nominal time binning for data produced by the detector.
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

    // The ductor transforms drifted depos to currents
    ductor: {

        # field_file: "dune-garfield-1d565.json.bz2",
        field_file: mydet.fields,

        // The distance from the anode centerline to where the field
        // response calculation begins drifting.  Take care that field
        // response calculations are not always done with the
        // thickness of the anode in mind and may report their
        // "response plane" measured relative to a number of other
        // monuments.  This value is almost certainly required in the
        // "geometry" section.
        response_plane : 10*wc.cm,

        // The time bin where the readout should be considered to
        // actually start given the pre-signal simulated.
        tbin : 0, 
        binning: {
            // Ductor ticks at same speed as ADC
            tick : $.binning.tick,
            // Over the somewhat enlarged domain
            nticks : $.ductor.tbin + $.binning.nticks,
        },
        start_time : -self.response_plane / $.lar.drift_speed,
        readout_time : self.binning.nticks * self.binning.tick,
    },

    // A "splat" (DepoFluxSplat) is an approximation to the combination of
    // ductor+sigproc
    splat : {
        sparse: true,
        tick: $.ductor.binning.tick,
        window_start: $.ductor.start_time,
        window_duration: $.ductor.readout_time,
        reference_time: 0.0,
    },

    // Simulating noise
    noise : {
        model: {
            #spectra_file: "protodune-noise-spectra-v1.json.bz2",
            spectra_file: mydet.noise,

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

    // Signal processing parameter pack
    sp: {
        // In principle, may use a different field response.
        field_file : $.ductor.field_file,
        binning : $.binning,
    },

    dnnroi:  {
        // this is a total fudge factor
        output_scale : 1.0,

        // unet-l23-cosmic500-e50.ts was made for uboone....
        model_file : "unet-l23-cosmic500-e50.ts",
    },
        

    // Imaging paramter pack
    img : {
        // For now we use MicroBooNE's
        #"charge_error_file": "microboone-charge-error.json.bz2",
        "charge_error_file": mydet.qerr,

        // Number of ticks to collect into one time slice span
        span: 4,

        binning : $.binning,

        // The tiling strategy refers to how to deal with substantial
        // dead channels.  The "perfect" strategy simply performs so
        // called "3-plane" tiling.  The "mixed" strategy also
        // includes so called "2-plane" tiling, producing more and
        // more ambiguous blobs but will not have "gaps" due to dead
        // channels.
        tiling_strategy: "perfect",
    },

};
    
local override = nominal + {
    ductor: super.ductor + {
        tbin : 0, 
        binning: {
            tick : $.binning.tick,
            nticks : $.ductor.tbin + $.binning.nticks,
        },
        start_time : -self.response_plane / $.lar.drift_speed,
        readout_time : self.binning.nticks * self.binning.tick,
    }
};

// simple
override

