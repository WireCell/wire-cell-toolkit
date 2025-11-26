/// This abstracts how to construct a detector.
local wc = import "wirecell.jsonnet";
local known_detectors = import "detectors.jsonnet";

{
    /// A binning summarizes a regular sampling along one dimension.
    binning(nsteps, step, start=0):: {
        nsteps:nsteps, step:step, start:start
    },
    /// Many components have time binning using specific variable names.  This converts
    binning_to_time(b):: { nticks: b.nsteps, tick: b.step, start: b.start },

    /// Note, do NOT change this even if you target some other value.
    default_adc_period: 500*wc.ns,


    /// Define some liquid argon.
    lar(DL,                     // Longitudinal diffusion constant, O(10 * wc.cm2/wc.s)
        DT,                     // Transverse diffusion constant, O(10 * wc.cm2/wc.s)
        lifetime,               // Electron lifetime
        drift_speed=1.6*wc.mm/wc.us, // Electron drift speed, assumes a certain applied E-field
        density=1.389*wc.g/wc.centimeter3, // LAr density
        ar39activity=1*wc.Bq/wc.kg):: // Decay rate per mass for natural Ar39.
        {
            DL : DL,
            DT : DT,
            lifetime : lifetime,
            drift_speed : drift_speed,
            density: density,
            ar39activity: ar39activity,
        },

    /// Describe one "anode" aka AnodePlane.  An AnodePlane is what holds "wires"
    /// (which may be in the shape of strips).  An AnodePlane has one (uboone) or
    /// two (DUNE HD and VD) "faces".  Faces can be opposing (DUNE HD) or aligned
    /// (DUNE VD).  The wires and the faces together imply one or two drift volumes
    ///
    /// @param ident The identifier in the wires file for wires of this anode.
    /// @param wires_object A WireSchemaFile config.
    /// @param faces An array of faces.  This MUST have a null place holder for insensitive faces, order matters.
    anode(ident, wires_object, faces):: {
        type : "AnodePlane",
        name : "a"+std.toString(ident),
        data : {
            // This IDENT must match a set of wires
            ident : ident,
            // The wire schema file
            wire_schema: wc.tn(wires_object),
            faces : faces,
        },
        uses: [wires_object],
    },

    /// Provide parameters governing response simulation (DepoTransform).
    ductor(tick, readout_duration, start_time=0):: {
        tick: tick, readout_duration: readout_duration, start_time: start_time
    },

    /// Describe an ADC
    adc(tick=$.default_adc_period, // sampling period
        resolution=14,          // number of bits 
        gain=1.0,               // unitless pregain VIN->VADC
        /// voltage baseline biases
        baselines=[1*wc.volt,1*wc.volt,1*wc.volt],
        /// Full scale of VADC
        fullscale=[0*wc.volt,2*wc.volt],
        /// Duration of one ADC waveform
        readout_duration=1*wc.ms)::
        {
            tick: tick,
            resolution: resolution, // how many bits in the ADC
            gain: gain,           // a relative, unitless gain
            baselines: baselines, // per-plane baseline volutabes
            fullscale: fullscale, // ADC fullscale voltage range

            readout_duration: readout_duration,
            // It would probably be better to use std.round but I don't have it
            // in my current version of Jsonnet.  This function came in >to
            // Jsonnet 0.21.0, April 2023.
            readout_nticks: std.floor(readout_duration/tick),

            max_count: 1 << self.resolution,

            
            valid_range: self.fullscale[1] - self.fullscale[0],

            /// How much voltage the ADC sees for each count, ie, after relative gain.
            vadc_per_count: self.valid_range / (self.max_count),
            /// Line-level Voltage of least-significant bit, ie before relative gain.
            lsb_voltage: self.vadc_per_count/self.gain,

        },

    /// Describe the longitudinal drift dimension of one face by giving absolute X
    /// locations of three planes.
    sensitive_face(ident, anode, response, cathode):: {
        anode:anode,           // where we consider drifting to stop
        response:response,     // boundary between constant drift and field response
        cathode:cathode,  // the furthest distance from anode to consider ionization
    },

    /// An anode face may have no wires or be considered insensitive (eg "back" of an APA).
    /// It is represented as a null value in the anode.faces
    innactive_face():: null,

    /// Make a faces array for DUNE hd-like detectors.  This is offered without
    /// much discussion.  See various "params.jsonnet" in wct/cfg/pgrapher for
    /// understanding.
    hd_like_faces(anode_ident,  // determines direction.
                  apa_cpa = 3.5734*wc.m,
                  cpa_thick = 3.175*wc.mm,
                  apa_w2w = 85.87*wc.mm,
                  plane_gap = 4.76*wc.mm)::
        local apa_g2g = apa_w2w + 6*plane_gap;
        local apa_plane = 0.5*apa_g2g - plane_gap;
        // relative to collection wires.  WARNING this must be consistent with the FR!
        local response_plane = 10*wc.cm;
        local res_plane = 0.5*apa_w2w + response_plane;
        local cpa_plane = apa_cpa - 0.5*cpa_thick;
        local sign = 2*(anode_ident%2)-1;
        local centerline = sign*apa_cpa;
        if sign > 0
        then [
            // top, front face is against cryo wall
            null,
            {
                anode: centerline - apa_plane,
                response: centerline - res_plane,
                cathode: centerline - cpa_plane, 
            }
        ]
        else [
            // bottom, back face is against cryo wall
            {
                anode: centerline + apa_plane,
                response: centerline + res_plane,
                cathode: centerline + cpa_plane, 
            },
            null
        ],

    /// The wire endpoints and channel map are typically provided by wires file.
    /// This is large and shared by multiple anodes so we use a "service" type
    /// component to hold and share the data.  Normally there is only one wires file
    /// per job but if you have multiple, give each a unique name.
    wires_from_file(filename, name=""):: {
        type: "WireSchemaFile",
        name: name,
        data: { filename: filename }
    },
    
    /// Get a wires file registered for the named detector
    wires_from_name(detname):: $.wires_from_file(known_detectors[detname].wires),
    
    /// The field response object
    fields_from_file(filename):: {
        type: "FieldResponse",
        name: filename,
        data: { filename: filename }
    },

    /// Get a field response registered by detector name.  In general, a list of
    /// fields can be registered and may be referenced by index.  The convention is
    /// not strong but usually index=0 is "nominal".
    fields_from_name(detname, index=0):: $.fields_from_file(known_detectors[detname].fields[index]),

    default_er_binning: $.binning(100, $.default_adc_period),

    /// Define a cold electronics response
    ///
    /// @param kind: The kind of electronics response model: "cold" or "warm". 
    /// @param gain The amplifier converting gain in units of [voltage/charge]
    /// @param shaping The amplifier shaping time in units of [time]
    /// @param postgain A relative, unitless gain
    /// @param tick A nominal sample period in units of [time].  This almost always should match ADC tick.
    /// @param nticks A nominal number of samples that at least covers where ER is non-zero.
    /// @param start A start time (probably always should be left at 0) in units of [time].
    ///
    /// Caution: The "nticks" is not well defined.  Older OSP code wants this to
    /// be the readout while newer SPNG more correctly wants this to be just
    /// large enough to cover where ER is nonzero.
    ///
    elec_response(gain, shaping, binning=$.default_er_binning, kind='cold', postgain=1.0, name=""):: {
        type: if kind == 'cold'
              then 'ColdElecResponse'
              else 'WarmElecResponse',
        name: name,
        data: {
            gain: gain, shaping: shaping, postgain: postgain,
        } + $.binning_to_time(binning)
    },

    /// An arbitrary electronics response from JSON file.  See elec_response() for common args.
    elec_response_file(filename, binning=$.default_er_binning, postgain=1.0, name=""):: {
        type: 'JsonElecResponse',
        name: name,
        data: {
            postgain: postgain, 
        } + $.binning_to_time(binning)
    },
        
    // fixme: add rc
    
    // Create a plane impact response (PIR) for one plane.  Any null values will
    // be supplied by tpc.  The context_name can be explicitly passed as "" if a
    // PIR for a given plane is the same for all anodes.
    plane_impact_response(context_name, plane,
                          tick,   // usually adc.tick
                          nticks, // usually adc.readout_nticks
                          fr, er, // as used in sim, may differ than osp and spng
                          rcs=[], // one or two "RC" responses.  Padding is
                          // typically not customized, but it needs to be big
                          // enough to avoid DFT cyclic artifacts.
                          long_padding=1.5*wc.ms,
                          overall_short_padding=200*wc.us,
                         )::
        {
            type: "PlaneImpactResponse",
            name: context_name + "p" + std.toString(plane),
            data: {
                plane : plane,
                field_response : wc.tn(fr),

                start: 0,
                tick: tick,
                nticks: nticks,

                overall_short_padding: overall_short_padding,
                short_responses: [wc.tn(er)],

                long_padding : long_padding,
                long_responses: [wc.tn(rc) for rc in rcs],
            },
            uses: [fr, er] + rcs

        },
    /// Each TPC gets one PIR per plane collected into an array.
    plane_impact_responses(u, v, w):: [u, v, w],

    /// A detector's empirical noise model parameters not derived from tpc or control.
    empirical_noise(spectra_file, wire_length_scale=10, replacement_percentage=0.02, chanstat=""):: {
        spectra_file: spectra_file,
        wire_length_scale: wire_length_scale,
        replacement_percentage: replacement_percentage,
        chanstat: chanstat,
    },
    /// Bundle the types of noise model parameters
    noise(empirical):: {
        empirical: empirical,
    },


    /// Define a filter configuration.  Kind can be "lowpass" (aka hf) or
    /// "highpass" (aka lf).  The kind maps to an analytical function of the
    /// other parameters.
    filter_config(scale,               // 
                  power=2.0,           // exponential power
                  kind="lowpass",
                  period=$.default_adc_period, // sampling period
                  // if true, zero-frequency is zeroed
                  ignore_baseline=true)::
        {
            scale: scale,
            power: power,
            kind: kind,
            period: period,
            ignore_baseline: ignore_baseline
        },

    /// Various filters in the time dimension are defined with semantic labels.
    /// These each should be defined with filter_config().  The filter set are 
    /// specified on a per-view basis.
    time_filters(gauss, wiener, dnnroi):: {
        /// The filter for final output signals from SPNG.  This should preserve
        /// signal.
        gauss: gauss,
        // The filter for input to threshold-based intial ROI construction.  It
        // should maximize signal/noise.
        wiener: wiener,         
        // The filter used to provide the non-MP2/MP3 image input to DNNROI. 
        dnnroi: dnnroi,
    },
    
    /// There is just one type of channel filter  per view which is used in decon.
    channel_filters(decon):: {
        decon: decon
    },

    /// Collect the time and channel filters for one view.
    view_filters(time_filters, channel_filters):: {
        time: time_filters,
        channel: channel_filters
    },

    /// Collect all the per-view filters in view index order
    tpc_filters(u, v, w):: [ u, v, w ],

    /// See Threshold.h for what these mean.
    crossview_threshold(nominal=0.0, rms_nsigma=0.0, rms_axis=-1, rms_max_value=0, binary=true):: {
        nominal: nominal,
        rms_nsigma: rms_nsigma,
        rms_axis: rms_axis,
        rms_max_value: rms_max_value,
        binary: binary,
    },
    /// Per-view list
    crossview_thresholds(cvt_u, cvt_v, cvt_w):: [cvt_u, cvt_v, cvt_w],

    /// A view group collects the wire planes in a common layer / view index
    /// in one face or two faces that are connected in an unspecified anode.
    ///
    /// The term "view" is introduced here to distinguish from (wire) "plane" (a
    /// physical unit) or "layer" (a logical/relative location).  A "view" is
    /// the combination of both.  When a view has two faces, the face ident
    /// array gives an ordering.
    ///
    /// @param connection A number describing face connectivity:
    ///        0: one face not connected to the other side
    ///        1: two faces connected on one edge (jumper)
    ///        2: two faces connected on both edges (wrapped)
    /// @param view_index The view/plane index. 0:U, 1:V 2:W
    /// @param face_idents An array of one or two face ID numbers.
    view_group(connection, view_index, face_idents):: {
        name: 'group_v' +  std.toString(view_index) +
              'f' + std.join('f',std.map(std.toString, face_idents)) +
              'c' + std.toString(connection),
        view_index: view_index,
        view_layer: wc.wpid_index_to_layer(view_index),
        face_idents: face_idents,
        connection: connection,

        // Get the wpids for this view group in a given anode.
        wpids(anode_ident):: [wc.WirePlaneId(self.view_layer, fid, anode_ident)
                              for fid in self.face_idents],
    },

    /// Helper to resolve case of 1 vs 2 face groups
    view_groups_one_layer(connection, view_index, face_idents)::
        if connection > 0
        then [$.view_group(connection, view_index, face_idents)]
        else [$.view_group(connection, view_index, [fid]) for fid in face_idents],

    /// Partition the set of faces and layers by connections and face_idents
    view_groups(connections,     // per-plane connection number: 0:none, 1:jumper, 2:wrapped
                face_idents)::   // which face(s), and provide an ordering
        wc.flatten([ $.view_groups_one_layer(connections[view_index], view_index, face_idents)
          for view_index in [0,1,2]]),

    /// Name for TPC of given ident.
    tpc_name(tpc_ident):: "tpc" + std.toString(tpc_ident),

    /// Define a TPC API.
    ///
    /// This rolls up all params required to cover all the higher-level
    /// construction.
    ///
    /// In general, each TPC may have unique parameters though typically they
    /// will all be identical.  
    ///
    tpc(anode=null,
        lar=null,
        ductor=null,
        adc=null, fr=null, er=null, rcs=[],
        pirs=null, noise=null,
        connections=null, filters=null, faces=null,
        crossview_thresholds=null)::
        {
            // the AnodePlane config object
            anode: anode,
            // The liquid argone
            lar: lar,
            // Induction parameters
            ductor: ductor,
            // ADC related parameters
            adc:adc,
            // response object configs, rcs may be an array of "RC" long responses.
            fr:fr, er:er, rcs:rcs,
            // how same-layer views on connected across faces
            connections: connections,
            // The decon filters config
            filters: filters,
            // Ordered enumeration of face idents
            faces: faces,

            // plane impact responses
            pirs: pirs,

            // noise models and their params
            noise: noise,

            // Configuration for Threshold for input to CrossViews
            crossview_thresholds: crossview_thresholds,

            // derived info below

            // Note, these can not be accessed while null!
            ident: anode.data.ident,
            name: $.tpc_name(self.ident),

            view_groups: $.view_groups(connections, faces),
    },


    /// Construct a detector as a named collection of tpcs.
    ///
    detector(name, tpcs, lar=null):: {
        name: name,
        tpcs: tpcs,
        tpc: {[t.name]:t for t in tpcs},
        lar: if std.type(lar) == "null"
             then tpcs[0].lar
             else lar,
    },

    /// Return a detector with a subset of TPCs selected by ID numbers.
    subset(det, tpc_idents)::
        if std.length(tpc_idents) == 0
        then det
        else det { tpcs: [det.tpc[$.tpc_name(ident)] for ident in tpc_idents] },


    // /// Some components accept various control parameters.  This bundles them.
    // control(device="cpu", verbosity=0):: {
    //     verbosity: verbosity, // 0:silent 1:one-per-exec 2:many-per-exec
    //     device: device,   // "cpu", "cuda", "gpu0", "gpu1", ...
    // },

    test: {
        adc: $.adc(),
    }
}

