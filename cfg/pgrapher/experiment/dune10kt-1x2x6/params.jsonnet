local wc = import "wirecell.jsonnet";
local base = import "pgrapher/dune/params.jsonnet";

function(params) base {
    // This section will be overwritten in simparams.jsonnet
    det : {

        // The "faces" is consumed by, at least, the Drifter and
        // AnodePlane.  The "wires" number is used to set
        // AnodePlane.ident used to lookup the anode in WireSchema.
        // It corresponds to the anode number.

        // between center lines
        // local apa_cpa = 3.637*wc.m, // DocDB 203
        local apa_cpa = 3.63075*wc.m, // LArSoft
        // local cpa_thick = 50.8*wc.mm, // DocDB 203
        local cpa_thick = 3.175*wc.mm, // 1/8", from Bo Yu (BNL) and confirmed with LArSoft

        // X positions from dune10kt gdml:
        // U: 39.5355 mm, V: 34.7755 mm, W: 30.0155 mm
        local apa_w2w = 60.031*wc.mm,
        local plane_gap = 4.76*wc.mm,
        local apa_g2g = apa_w2w + 6*plane_gap, // 88.591*wc.mm,

        // The "anode" cut off plane, here measured from APA
        // centerline, determines how close to the wires do we
        // consider any depo.  Anything closer will simply be
        // discarded, else it will either be drifted or "backed up" to
        // the response plane.  This is somewhat arbitrary choice.
        // Placing it w/in the response plane means any depos that are
        // "backed up" won't have proper field response.  But, the
        // tighter this is made, the less volume is simulated.
        // local apa_plane = 0.5*apa_g2g, // pick it to be at the grid wires
        local apa_plane = 0.5*apa_g2g - plane_gap, // pick it to be at the first induction wires

        // The "response" plane is where the field response functions
        // start.  Garfield calcualtions start somewhere relative to
        // something, here's where that is made concrete.  This MUST
        // match what field response functions also used.
        response_plane: 10*wc.cm, // relative to collection wires
        local res_plane = 0.5*apa_w2w + self.response_plane,

        // The cathode plane is like the anode cut off plane.  Any
        // depo not between the two is dropped prior to drifting.
        local cpa_plane = apa_cpa - 0.5*cpa_thick,


        // The volumes are then defined in terms of these above
        // numbers.  You can use "wirecell-util wires-info" or
        // "wirecell-util wires-volumes" or others to understand the
        // mapping of anode number to the 6 locations in X and Z.  For
        // Larsoft wires the numbering is column major starting at
        // small X and Z so the centerline is -/+/-/+/-/+.  Also
        // important is that the faces are listed "front" first.
        // Front is the one with the more positive X coordinates and
        // if we want to ignore a face it is made null.
        volumes: [
            {
                local sign = 2*(n%2)-1,
                local centerline = sign*apa_cpa,
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
            } for n in std.range(0,11)], // twelve anodes

        // This describes some rough, overall bounding box.  It's not
        // directly needed but can be useful on the Jsonnet side, for
        // example when defining some simple kinematics.  It is
        // represented by a ray going from extreme corners of a
        // rectangular solid.  Again "wirecell-util wires-info" helps
        // to choose something.
        bounds : {
          tail: wc.point(-363.376, -600.019, -0.87625, wc.cm),
          head: wc.point( 363.376,  600.019, 1393.46, wc.cm),
        }
    },

    daq: super.daq {
        nticks: 6000,
    },

    adc: super.adc {
        // per tdr, chapter 2
        // induction plane: 2350 ADC, collection plane: 900 ADC
        baselines: [1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt],

        // check this.  The tdr says, "The ADC ASIC has an input
        // buffer with offset compensation to match the output of the
        // FE ASIC.  The input buffer first samples the input signal
        // (with a range of 0.2 V to 1.6 V)..."
        fullscale: [0.2*wc.volt, 1.6*wc.volt],

    },

    // This sets a relative gain at the input to the ADC.  Note, if
    // you are looking to fix SimDepoSource, you are in the wrong
    // place.  See the "scale" parameter of wcls.input.depos() defined
    // in pgrapher/common/ui/wcls/nodes.jsonnet.
    // also, see later overwriting in simparams.jsonnet
    elec: super.elec {
      postgain: 1.1365, // pulser calibration: 41.649 ADC*tick/1ke
                       // theoretical elec resp (14mV/fC): 36.6475 ADC*tick/1ke
      shaping: 2.2 * wc.us,
    },

    sim: super.sim {

        // For running in LArSoft, the simulation must be in fixed time mode. 
        fixed: true,

        // The "absolute" time (ie, in G4 time) that the lower edge of
        // of final readout tick #0 should correspond to.  This is a
        // "fixed" notion.
        local tick0_time = if std.objectHas(params, 'G4RefTime') then params.G4RefTime else 0,

        // Open the ductor's gate a bit early.
        local response_time_offset = $.det.response_plane / $.lar.drift_speed,
        local response_nticks = wc.roundToInt(response_time_offset / $.daq.tick),

        ductor : {
            nticks: $.daq.nticks + response_nticks,
            readout_time: self.nticks * $.daq.tick,
            start_time: tick0_time - response_time_offset,
        },

        // To counter the enlarged duration of the ductor, a Reframer
        // chops off the little early, extra time.  Note, tags depend on how 
        reframer: {
            tbin: response_nticks,
            nticks: $.daq.nticks,
        }
        
    },

    files: {
        wires: "dune10kt-1x2x6-wires-larsoft-v1.json.bz2", // "protodune-wires-larsoft-v4.json.bz2",

        fields: [
            // "garfield-1d-3planes-21wires-6impacts-dune-v1.json.bz2",
            // "garfield-1d-boundary-path-rev-dune.json.bz2",
            // "dune-garfield-1d565.json.bz2",
            "dune-garfield-1d60563.json.bz2",
        ],

        // fixme: this is for microboone and probably bogus for
        // protodune because (at least) the span of MB wire lengths do
        // not cover pdsp's.
        noise: "protodune-noise-spectra-v1.json.bz2",


        chresp: null,
    },

}

