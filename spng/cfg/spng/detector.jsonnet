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


    /// Describe an ADC
    adc(resolution=14, gain=1.0, baselines=[1*wc.volt,1*wc.volt,1*wc.volt], fullscale=[0*wc.volt,2*wc.volt]):: {
        resolution: resolution, // how many bits in the ADC
        gain: gain,           // a relative, unitless gain
        baselines: baselines, // per-plane baseline volutabes
        fullscale: fullscale, // ADC fullscale voltage range

        max_count: std.pow(2,self.resolution),
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
    wires_from_detector(detname):: $.wires_from_file(known_detectors[detname].wires),
    
    /// The field response object
    fields_from_file(filename):: {
        type: "FieldResponse",
        name: filename,
        data: { filename: filename }
    },

    /// Get a field response registered by detector name.  In general, a list of
    /// fields can be registered and may be referenced by index.  The convention is
    /// not strong but usually index=0 is "nominal".
    fields_from_detector(detname, index=0):: $.fields_from_file(known_detectors[detname].fields[index]),

    default_er_binning: $.binning(100,500*wc.ns),

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
        name : std.toString(ident),
        data : {
            // This IDENT must match a set of wires
            ident : ident,
            // The wire schema file
            wire_schema: wc.tn(wires_object),
            faces : faces,
        },
        uses: [wires_object],
    },

    
    /// An anode face is composed of 3 layers wire planes labeled "u", "v" and
    /// "w".  When an anode has two faces, some number of the 3 layers on one
    /// face may be electrically connected to wires on the other face.  Eg, DUNE
    /// VD have 2 faces side-by-side with 1 face edge of U and V connected
    /// ("jumpered" connect=1) while DUNE HD has 2 opposing faces and both face
    /// edges connected ("wrapped" connect=2) while both have independent
    /// collection faces (connect=0).  The view_layer object defines these
    /// associations.
    ///
    /// @param view_index: 0 for U, 1 for V and 2 for W layers.
    /// @param connect: 0 for independent faces, 1 for one-side jumper, 2 for two-side wrapped.
    /// @param face_idents: the order of face IDENTS to make when connect>0. 
    ///
    /// Given face_idents=[f1,f2] then the channels may be arranged in
    /// wire-attachment-number (WAN) of each face in order [f1,f2] such that the
    /// two channels at the boundary between the two blocks of channels "f1" and
    /// "f2" have neighboring wires.  This is most important for connect=1.  For
    /// connect=2 both [f1,f2] and [f2,f1] meet the condition as the wires array
    /// is cyclic.  For connect=2 case, it merely sets an order convention.
    ///
    /// The face ident numbers are as given in the wires file.  
    view_layer(view_index, connect, face_idents=[0,1]) :: {
        view_index: view_index, connect: connect, face_idents: face_idents
    },

    /// This groups configuration related to an anode and its drift volume(s).
    ///
    /// This is not directly mapped to any one WCT component but holds context
    /// used by many.
    ///
    /// @param anode An anode config object.
    /// @param adc An adc() object.
    /// @param view_layers An array of view_layers objects, one for u, v and w views.
    /// @param fr A field response config object
    /// @param er An electronics response config object
    /// @param rc An array of zero or more long-time response objects.
    tpc(anode, adc, view_layers, fr, er, rc=[], name=""):: {
        name: name,
        anode: anode, adc:adc,
        view_layers: view_layers,
        fr:fr, er:er, rc:rc,
        // fixme: add readout binning, various time offsets.
    },

    /// Construct a detector as a named collection of tpcs.
    ///
    /// This is what each detectors/<name>.jsonnet should export.
    detector(tpcs, name=""):: {
        name: name,
        tpcs: tpcs,
    },

    test: {
        adc: $.adc(),
    }
}

