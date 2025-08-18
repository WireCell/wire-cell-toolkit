/** This file defines the "omni object API".

The API is defined by the object produced by this file and the comments giving
expectations for each attribute and function.


Goals: 

The goal of the omni object API is to provide a means by which a user can define
arbitrary WCT configurations in a high-level, detector-independent manner.  This
is possible with an API that defines enough "mid-level" attributes and functions
which then have "low-level" detector-dependent implementations.  The high-level
code need be parameterized by a single string giving the omni object "name",
unique across all omni objects.


Files:

Every canonical or variant detector configuration should implement an omni
object and these should be collected under:

  cfg/omni/<detname>/omni.jsonnet

The "<detname>" is the canonical name as used as a key in the file
"detectors.jsonnet" provided by wire-cell-data.

An "omni.jsonnet" file may produce a single omni object or a list of omni
objects.  If convenient, the "low-level" implementation may place additional
files under cfg/omni/<detname>/.  Implementation is allowed to but discouraged
from using any Jsonnet file under cfg/pgrapher/.  This discouragement is to
encourage a break from the complicated forms in that area.  However, it is
recognized that making a clean break may not be practical and so they are
allowed.

An omni object may inherit from the omni object API object and override
attributes and methods by following the directives defined in comments.  The
implementer of an omni object MUST read and follow the comments on the object
below to understand what MUST be provided and what MAY or SHOULD be provided in
order to produce a valid omni object.


Using:

The high-level configuration file is expected to have content similar to the
following:

  local omnimap = import "omni/all.jsonnet";
  function(name)
      local omni = omnimap[name];
      ...

This parameterized the configuration by the omni name (not necessarily the
canonical detector name).  The "omni" object then provides all the primitives
from which a WCT graph may be constructed.

In some advanced cases, more than one omni name may parameterize a
configuration.  For example, a job that runs both OSP and SPNG will need two
omni objects in order to provide unique "sigproc()" functions.


Evolution:

The omni object API is initially defined with an imperfect understanding of all
current possible and future high-level configurations.  As such, the API is
expected to evolve.  As it does, implementations may require extension or
modification.  To facilitate this, it is strongly recommended to place all
omni.jsonnet files in the wire-cell-toolkit and NOT maintain them elsewhere.
When changes are required, no attempt will be made to accommodate externally
defined omni objects.

*/

local wc = import "wirecell.jsonnet";
local pg = import "pgraph.jsonnet";

// All canonical detectors should be registered in the detectors.jsonnet file
// provided by wire-cell-data.
local detectors = import "detectors.jsonnet";

/// This then defines the omni API:
{
    /** Basic ATTRIBUTES describing the omni object */
    
    /// MUST supply a unique name to identify the detector variant.  This is not
    /// necessarily the canonical "detname".
    name : "base",

    /// SHOULD set. This gives the "canonical detector name" and MUST be a key
    /// in the detectors.jsonnet object.
    detname: self.name,

    /// MAY set a custom detectors.jsonnet type object.  Default will do the
    /// lookup for you.
    detdata: detectors[self.detname],


    /** TOOL attributes giving some basic component configuration.  While the
    tool types are needed by low-level configuration, their implementation is
    more an end-user choice. */

    /// MAY override to set the IRandom seed.  Note, theses seeds may be
    /// interpreted differently by different IRandoms (if/when WCT ever gets
    /// more than just Random).
    random_seeds: [0,1,2,3,4],

    /// MAY override to set a different kind of IRandom.
    random: {
        type: "Random",
        data: {
            generator: "default",
            seeds: $.random_seeds
        }
    },

    /// MAY override to set the IDFT implementation.
    dft: {
        type: "FftwDFT",
    },


    /// DETECTOR attributes.  These define API objects that describe the
    /// detector.  All SHOULD be overridden by detector specific omni objects.

    /// SHOULD override to define LAr properties.
    lar: {
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
    adc: {
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

        // Calculate the ADC range in voltage
        valid_range: self.fullscale[1] - self.fullscale[0],

        // The least significant bit voltage.
        lsb_voltage: self.valid_range / (1 << self.resolution),
    },        

    /// SHOULD define the drift volumes.  This example merely shows the schema.
    volumes : [
        {
            wires: 0,           // obsolete and poorly named when it was used.
            name: "onesided",
            faces: [ {anode: 0, response: $.detdata.response_plane, cathode: 2*wc.m},
                     null ],
        },
        {
            wires: 1,
            name: "twosided",
            faces: [ {anode: 0, response: +$.detdata.response_plane, cathode: 2*wc.m},
                     {anode: 0, response: -$.detdata.response_plane, cathode: -2*wc.m} ],
        },
        {
            wires: 2,
            name: "anothertwosided",
            faces: [ {anode: 4*wc.m, response: +$.detdata.response_plane + 4*wc.m, cathode:  2*wc.m + 4*wc.m },
                     {anode: 4*wc.m, response: -$.detdata.response_plane + 4*wc.m, cathode: -2*wc.m + 4*wc.m} ],
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
    anodes : std.mapWithIndex(function(ident, vol) {
        local wirefile = $.detdata.wires,
        local wires = { type: 'WireSchemaFile', name: wirefile, data: { filename: wirefile }},
        type: 'AnodePlane',
        name: vol.name,
        data: {
            ident: ident,
            wire_schema: wc.tn(wires),
            faces: vol.faces,
        },
        uses: [wires]}, $.volumes),

    /// MAY supply a function to return various response objects.  The "kind"
    /// may be at least one of (sim, sp, splat).  Should return an object with
    /// keys (fr, er, and rc).  All three are ARRAYS of appropriate response
    /// config objects.
    responses(kind, binning=$.binning.sim) :: {

        /// Note, here we ignore "kind" which means user gets same regardless of
        /// application.  Real detector config may be more specific.

        // This is used in common between rc and er
        local tbinning = {
            start: binning.start,
            tick: binning.spacing,
            nticks: binning.number
        },

        /// MAY override, but default should be okay for most. 
        /// The 2D "field responses (usually only one)
        fr: [{ type: 'FieldResponse', name: $.detname, data: { filename: detectors[$.detname].field } }],

        /// SHOULD override, different detectors tend to differ.
        /// The 1D "long" electronics responses
        rc: [{type: "RCResponse", name: $.detname, data: tbinning}],

        /// SHOULD override, different detectors tend to differ.
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
                /// Model the extra smearing that occurs due to real SP filters.
                /// See test/test/test-morse-pdsp.bats for example how to
                /// measure this.  These values are not general.
                smear_long: 2.0, // units of tick
                smear_tran: 2.0, // in units of pitch
            },
        }, nin=1, nout=1, uses=[anode, res.fr[0]]),

    /// MAY override to return a trio of pirs for the given field response.
    /// This is likely okay as-is for most detectors.
    pirs(binning = $.binning.sim, res = $.responses("sim")) :: [ {
        local tbinning = {
            start: binning.start,
            tick: binning.spacing,
            nticks: binning.number
        },
        type: 'PlaneImpactResponse',
        name: 'PIR%splane%s' % [res.fr[0].name, plane],
        data: tbinning {
            plane: plane,
            field_response: wc.tn(res.fr[0]),
            short_responses: [wc.tn(er) for er in res.er],
            // Determines the TOTAL size of the convolution between multiple
            // ER's (or size of one ER), and the ADDITIONAL size to the FR size
            // for FR*ER.  HD tends to pick 100us, VD picks 200us.
            overall_short_padding: 200*wc.us,

            long_responses: [wc.tn(rc) for rc in res.rc],
            // Determines the ADDITIONAL size to the time duration of the signal
            // for the convolution with "RC" aka "long" response.
            long_padding: 1.5 * wc.ms,
        },
        uses: res.fr + res.er + res.rc,
    } for plane in [0,1,2] ],

    /// MAY override to return a custom detector simulation.
    /// This is likely okay as-is for most detectors.
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


    /// MAY override to return a custom noise adder.
    /// Detector MUST override if the default incoherent noise is not sufficient.
    noise(anode, binning = $.binning.sim) :: 
        local model = {
            type: 'EmpiricalNoiseModel',
            name: anode.name,
            data: {
                anode: wc.tn(anode),
                dft: wc.tn($.dft),
                chanstat:"",    // must explicitly empty
                spectra_file: $.detdata.noise,
                nsamples: binning.number,
                period: binning.spacing,
                // Optimization by binning wire lengths.
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
    
    /// MAY override to return a custom digitize.  This likely works for most
    /// detectors as-is.
    digitize(anode) :: pg.pnode({
            type: 'Digitizer',
            name: anode.name,
            data: $.adc {
                anode: wc.tn(anode),
            }
        }, nin=1, nout=1, uses=[anode]),

    /// MUST override to return a noise filter for the anode context.
    ///
    /// The noise filter is generally interpreted as something that is inserted
    /// between raw ADC waveforms and signal processing.  It should consume and
    /// produce ADC-level waveforms.  However, it need not re-apply integer
    /// truncation to its output.
    ///
    /// NF is in general pretty messy and the default supplied here is NOT
    /// adequate.
    ///
    /// Note, detectors that require resampling should include a resampler along
    /// with a noise filter in a pipeline..
    ///
    /// A real detector omni object may pass chndb_data and/or onf_data to
    /// override default OmniChannelNoiseDB and OmnibusNoiseFilter .data
    /// respectively.
    noisefilter(anode, binning=$.binning.sp, chndb_data={}, onf_data={}) ::
        local resp = $.responses("sp", binning);
        local chndb = {
            type: 'OmniChannelNoiseDB',
            name: anode.name,
            uses: [anode, resp.fr[0], $.dft],
            data: {
                anode: wc.tn(anode),
                field_response: wc.tn(resp.fr[0]),
                dft: wc.tn($.dft),
                tick: binning.spacing,


                // WARNING: this is generally not correct.  nsamples is the number
                // of frequency bins in which to filter and not necessarily the same
                // as the number of ticks.
                nsamples: binning.number,

                /// Provide list of electronics groups for coherent noise.
                groups: [],

                // In real detectors there tends to be a long list of bad
                // channels.
                bad: [],

                // A list of channel "info" registered by channel IDs.  Here we
                // give empty but real detectors must give something.  
                channel_info: [ ],
            } + chndb_data,
        };
        // our return value
        pg.pnode({
            type: 'OmnibusNoiseFilter',
            name: anode.name,
            data: {

                // Nonzero forces the number of ticks in the waveform
                nticks: 0,
                // Map internal semantic symbols to cmm mask keys.
                maskmap: {sticky: "bad", ledge: "bad", noisy: "bad"},
                // Real detector should give non-empty list of filters.
                channel_filters: [],
                grouped_filters: [],
                channel_status_filters: [],
                noisedb: wc.tn(chndb),
                intraces: 'orig' + std.toString(anode.data.ident),  // frame tag get all traces
                outtraces: 'raw' + std.toString(anode.data.ident),
            } + onf_data,
        }, uses=[chndb, anode], nin=1, nout=1),

    /// MUST override to return signal processing for the anode context.  This
    /// is too messy to provide a reasonable default.  A real detector omni
    /// object may call this with osp_data to overwrite OmnibusSigProc.data.
    sigproc(anode, binning=$.binning.sp, osp_data={}) :: 
        local resp = $.responses("sp", binning);

        // Use a per-channel response if the detector defines a pcr file.
        local pcr = if std.objectHas($.detdata, "chresp")
                    then {
                        obj: {
                            type: "PerChannelResponse",
                            data: { filename: $.detdata.chresp, },
                        },
                        tn: wc.tn(self.obj),
                        uses: [self.obj]
                    }
                    else {
                        obj: null,
                        tn: "",
                        uses: []
                    };

        pg.pnode({
            type: 'OmnibusSigProc',
            name: anode.name,
            data: {
                anode: wc.tn(anode),
                dft: wc.tn($.dft),
                field_response: wc.tn(resp.fr[0]),
                elecresponse: wc.tn(resp.er[0]),
                per_chan_resp: pcr.tn,
                fft_flag: 0,  // 1 is faster but higher memory, 0 is slightly slower but lower memory
                // Note, this parameter wrongly named.  It is ADC count per *WCT voltage unit* which is explicitly NOT mV.
                ADC_mV: 1.0/$.adc.lsb_voltage,

                // frame tags
                wiener_tag: 'wiener%d' % anode.data.ident,
                decon_charge_tag: 'decon_charge%d' % anode.data.ident,
                gauss_tag: 'gauss%d' % anode.data.ident,

                use_roi_debug_mode: false,
                tight_lf_tag: 'tight_lf%d' % anode.data.ident,
                loose_lf_tag: 'loose_lf%d' % anode.data.ident,
                cleanup_roi_tag: 'cleanup_roi%d' % anode.data.ident,
                break_roi_loop1_tag: 'break_roi_1st%d' % anode.data.ident,
                break_roi_loop2_tag: 'break_roi_2nd%d' % anode.data.ident,
                shrink_roi_tag: 'shrink_roi%d' % anode.data.ident,
                extend_roi_tag: 'extend_roi%d' % anode.data.ident,

                use_multi_plane_protection: false,
                mp3_roi_tag: 'mp3_roi%d' % anode.data.ident,
                mp2_roi_tag: 'mp2_roi%d' % anode.data.ident,
                
                isWrapped: false,
            } + osp_data,
        }, nin=1, nout=1, uses=[anode, $.dft, resp.fr[0], resp.er[0]] + pcr.uses),



    /// MUST override to return imaging for the anode context.  FIXME: this may
    /// not be the right abstraction.
    img(anode)::
        std.assertEqual("img" == "implemented"),        

}
