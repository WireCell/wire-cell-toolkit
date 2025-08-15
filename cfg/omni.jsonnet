/** This is the top-level omni.jsonnet defining the omni object API.

Concrete omni objects may inherit from this object and override attributes and
methods.  The implementer of an omni object MUST read the comments on the object
below to understand what MUST be provided and what MAY or SHOULD be provided.

See omnimap.jsonnet for how concrete omni objects are aggregated.

The omni object API itself only depends on top-level cfg/*.jsonnet files and the
detecgtors.jsonnet file from wire-cell-data.  In particular independent, the
omni object API is independent from any files under pgrapher/.  However,
experiment-specific omni objects utilize whatever they wish.

*/

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

// All detectors should be registered in this index that is provided in
// wire-cell-data.
local detectors = import "detectors.jsonnet";

/// This then defines the omni API:
{
    /// Basic ATTRIBUTES ///
    
    /// MUST supply a unique name to identify the detector variant.
    name : "base",

    /// MAY set. This gives the "detector name" for detectors.jsonnet lookup.  
    detname: self.name,

    /// MAY set a detectors.jsonnet object
    det: detectors[self.detname],


    /// TOOL attributes ///

    /// MAY override to set a the IRandom seed.
    random_seeds: [0,1,2,3,4],
    /// MAY override to set a different kind of IRandom
    random: {
        type: "Random",
        data: {
            generator: "default",
            seeds: $.random_seeds
        }
    },

    /// MAY override to set the IDFT implementation
    dft : {
        type: "FftwDFT",
    },


    /// DETECTOR attributes ///

    /// SHOULD override to define LAr properties.
    lar : {
        // Longitudinal diffusion constant
        DL :  7.2 * wc.cm2/wc.s,
        // Transverse diffusion constant
        DT : 12.0 * wc.cm2/wc.s,
        // Electron lifetime
        lifetime : 8*wc.ms,
        // Electron drift speed, assumes a certain applied E-field
        drift_speed : 1.6*wc.mm/wc.us, // at 500 V/cm
        // LAr density
        density: 1.389*wc.g/wc.centimeter3,
        // Decay rate per mass for natural Ar39.
        ar39activity: 1*wc.Bq/wc.kg,
    },

    /// SHOULD override to provide adc parameters.
    adc : {
        // A relative gain applied just prior to digitization.  This
        // is not FE gain, see elec for that.
        gain: 1.0,

        // Voltage baselines added to any input voltage signal listed
        // in a per plan (U,V,W) array.
        baselines: [900*wc.millivolt,900*wc.millivolt,200*wc.millivolt],

        // The resolution (bits) of the ADC
        resolution: 12,

        // The voltage range as [min,max] of the ADC, eg min voltage
        // counts 0 ADC, max counts 2^resolution-1.
        fullscale: [0*wc.volt, 2.0*wc.volt],
    },        

    /// SHOULD define the drift volumes.  These just show the syntax.
    volumes : [
        {
            wires: 0,           // poorly named.  This is an anode number.
            name: "onesided",
            faces: [ {anode: 0, response: 10*wc.cm, cathode: 2*wc.m},
                     null ],
        },
        {
            wires: 1,
            name: "twosided",
            faces: [ {anode: 0, response: +10*wc.cm, cathode: 2*wc.m},
                     {anode: 0, response: -10*wc.cm, cathode: -2*wc.m} ],
        },
        {
            wires: 2,
            name: "anothertwosided",
            faces: [ {anode: 4*wc.m, response: +10*wc.cm + 4*wc.m, cathode:  2*wc.m + 4*wc.m },
                     {anode: 4*wc.m, response: -10*wc.cm + 4*wc.m, cathode: -2*wc.m + 4*wc.m} ],
        },
    ],

    /// MUST supply various "binning" objects.  A "binning" object describes a
    /// 1-D regular grid of points in a 1-D continuous space: (spacing, number,
    /// start).  The "spacing" gives separation between neighboring points,
    /// "number" is the number of points in the grid and "start" gives the 1-D
    /// coordinate of the first point.
    binning: {
        /// Provides a default time binning.
        default_time_binning: { spacing: 500*wc.ns, number: 1024, start: 0.0*wc.ns },
        /// The time binning for full sim
        sim: self.default_time_binning,
        /// The time binning for full sigproc
        sp: self.default_time_binning,
        /// The time binning for splat, the "fast" version of combined sim+sp "truth".
        splat: self.default_time_binning,
    },

    /// DETECTOR model objects and methods ///

    /// MAY override to supply an array of AnodePlane configuration objects.
    /// This is parameterized by volumes and wires file name.
    anodes : [{
        local wirefile = $.det.wires,
        local wires = { type: 'WireSchemaFile', name: wirefile, data: { filename: wirefile }},
        type: 'AnodePlane',
        name: vol.name,
        data: {
            ident: wires,
            wire_schema: wc.tn(wires),
            faces: vol.faces,
        },
        uses: [wires]
    } for vol in $.volumes],

    /// MAY supply a function to return various response objects.  The "kind"
    /// may be at least one of (sim, sp, splat).  Should return an object with
    /// keys (fr, er, and rc).  All three are ARRAYS of appropriate response
    /// config objects.
    responses(kind, binning=$.binning.sim) :: {
        // By default, we ignore kind.
        local tbinning = { start: binning.start, tick: binning.spacing, nticks: binning.number },
        /// The 2D "field responses (usually only one)
        fr: [{ type: 'FieldResponse', name: $.detname, data: { filename: detectors[$.detname].field } }],
        /// The 1D "long" electronics responses
        rc: [{type: "RCResponse", name: $.detname, data: tbinning}],
        /// The 1D short "electronics" responses
        er: [{ type: 'ColdElecResponse', name: $.detname,
            data: tbinning {
                // The shaping (aka peaking) time of the amplifier shaper.
                shaping : 2.0*wc.us,
                // The FE amplifier gain in units of Voltage/Charge.
                gain: 14.0*wc.mV/wc.fC, 
                // A relative gain applied after shaping and before ADC.
                postgain: 1.0,
            }}]
    },


    /// METHODS to make nodes ///

    /// Generally, methods return a PGraph node.  A omni attributes may be used
    /// with overrides passed as optional arguments.  Where multiple primitive
    /// PGraph nodes are needed, they should be wrapped in a pg.pipeline() or
    /// otherwise assembled into a subgraph with pg.intern(). 

    /// MAY supply an array of configuration objects that implement a pipeline
    /// for drifting a IDepoSet in the given anode.  The results must be
    /// uniquely named in and with the context of the anode.
    drift(anode, time_offset=0.0, fluctuate=true) :: pg.pnode({
        type: 'DepoSetDrifter',
        name: anode.name,
        data: { drifter: "Drifter:"+anode.name }
    }, nin=1, nout=1, uses = [
        pg.pnode({
            type: 'Drifter',
            name: anode.name,
            data: $.lar {
                rng: wc.tn($.random),
                xregions: std.set(std.prune(std.flattenArrays([v.faces for v in $.volumes])), std.toString),
                time_offset: time_offset,
                fluctuate: fluctuate,
            },
        }, nin=1, nout=1, uses=[anode, $.random])]),
        
    /// SHOULD override to return a custom "splat" in the context of an anode.
    /// An optional binning and reference time may be passed by the caller.  See
    /// the "morse" job for a way to estimate the extra smearing needed to match
    /// full sim+sigproc.
    splat(anode, reference_time=0.0, binning = $.binning.splat) ::
        local res = $.responses("splat", binning);
        pg.pnode({
            type: 'DepoFluxSplat',
            name: anode.name,
            data: {
                anode: wc.tn(anode),
                sparse: true,
                tick: binning.spacing,
                window_start: binning.start,
                window_duration: self.tick * binning.number,
                reference_time: reference_time,
                field_response: res.fr[0],
                /// extra smearing, see "morse"
                smear_long: 2.0, // units of tick
                smear_tran: 2.0, // in units of pitch
            },
        }, nin=1, nout=1, uses=[anode, res.fr[0]]),

    /// MAY override to return a trio of pirs for the given field response.
    pirs(binning = $.binning.sim, res = $.responses("sim")) :: [ {
            local tbinning = { start: binning.start, tick: binning.spacing, nticks: binning.number },
            type: 'PlaneImpactResponse',
            name: 'PIR%splane%s' % [res.fr[0].name, plane],
            data: tbinning {
                plane: plane,
                field_response: wc.tn(res.fr[0]),
                short_responses: [wc.tn(er) for er in res.er],
                overall_short_padding: 100*wc.us,
                long_responses: [wc.tn(rc) for rc in res.rc],
            },
            uses: res.fr + res.er + res.rc,
        } for plane in [0,1,2] ],

    /// MAY override to return a custom detector simulation.
    signal(anode, binning = $.binning.sim) ::
        local pirs = $.pirs(binning);
        pg.pnode({
            type: 'DepoTransform',
            name: anode.name,
            data: {
                rng: wc.tn($.random),
                dft: wc.tn($.dft),
                anode: wc.tn(anode),
                pirs: [wc.tn(p) for p in pirs],
                fluctuate: true,
                drift_speed: $.lar.drift_speed,
                first_frame_number: 0,
                readout_time: binning.spacing*binning.number, // a time duration
                start_time: binning.start,
                tick: binning.spacing,
                nsigma: 3,
            },
        }, nin=1, nout=1, uses = pirs + [anode, $.random, $.dft]),


    /// MAY override to return a custom noise adder
    noise(anode, binning = $.binning.sim) :: 
        local model = {
            type: 'EmpiricalNoiseModel',
            name: anode.name,
            data: {
                anode: wc.tn(anode),
                dft: wc.tn($.dft),
                chanstat:"",    // must explicitly empty
                spectra_file: $.det.noise,
                nsamples: binning.number,
                period: binning.spacing,
                wire_length_scale: 1.0*wc.cm,
            },
            uses: [anode, $.dft]
        };
        pg.pnode({
            type: 'AddNoise',
            name: anode.name,
            data: {
                rng: wc.tn($.random),
                dft: wc.tn($.dft),
                model: wc.tn(model),
                nsamples: binning.number, // strictly, this is not necessarily the same as sim.binning
                replacement_percentage: 0.02,
            }}, nin=1, nout=1, uses=[$.random, model]),
    
    /// MAY override to return a custom digitize. 
    digitize(anode) :: pg.pnode({
            type: 'Digitizer',
            name: anode.name,
            data: $.adc {
                anode: wc.tn(anode),
            }
        }, nin=1, nout=1, uses=[anode]),

    /// MUST override to return a noise filter for the anode context.  This is
    /// too messy to provide a reasonable default.
    noisefilter(anode) ::
        std.assertEqual("noisefilter" == "implemented"),

    /// MUST override to return signal processing for the anode context.  This
    /// is too messy to provide a reasonable default.
    sigproc(anode) :: 
        std.assertEqual("sigproc" == "implemented"),

}
