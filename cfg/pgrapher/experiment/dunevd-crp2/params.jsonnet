// ProtoDUNE-SP specific parameters.  This file inerets from the
// generic set of parameters and overrides things specific to PDSP.

local wc = import "wirecell.jsonnet";
local base = import "pgrapher/common/params.jsonnet";

base {
    // This section will be overwritten in simparams.jsonnet
    det : {

        // See:  wirecell-util wire-volumes protodune-wires-larsoft-v3.json.bz2
        // to help with defining these parameters.

        // between center lines
        local apa_cpa = 313.03*wc.cm,
        local cpa_thick = 50.8*wc.mm,
        local apa_w2w = 85.725*wc.mm,
        local plane_gap = 4.76*wc.mm,
        local apa_g2g = 114.3*wc.mm, 

        local apa_plane = 0.5*apa_g2g, // pick it to be at the grid wires

        // The "response" plane is where the field response functions
        // start.  Garfield calcualtions start somewhere relative to
        // something, here's where that is made concrete.  This MUST
        // match what field response functions also used.
        response_plane: 18.92*wc.cm, // relative to collection wires
        local res_plane = 0.5*apa_w2w + self.response_plane,

        // The cathode plane is like the anode cut off plane.  Any
        // depo not between the two is dropped prior to drifting.
        local cpa_plane = apa_cpa - 0.5*cpa_thick,

        volumes: [
            {
                local world = 100,
                local split = s*10, // 1: left, 2: right
                local anode = a, // physical anode number
                wires: world + split + anode,
                name: "anode%d"%(world + split + anode),

                faces: [
                    {
                        anode:     15.07*wc.cm,
                        response:  15.07*wc.cm - 18.92*wc.cm,
                        cathode:   -300*wc.cm,
                    },
                    null
                ],
            } for a in std.range(0,1) for s in std.range(1,2)
        ],

        // This describes some rough, overall bounding box.  It's not
        // directly needed but can be useful on the Jsonnet side, for
        // example when defining some simple kinematics.  It is
        // represented by a ray going from extreme corners of a
        // rectangular solid.  Again "wirecell-util wires-info" helps
        // to choose something.
        bounds : {
            tail: wc.point( 0.0, -150,   0.0, wc.cm),
            head: wc.point(23.0,  150, 300.0, wc.cm),
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

    elec: super.elec {
      type: "JsonElecResponse",
      filename: "dunevd-coldbox-elecresp-top-psnorm_400.json.bz2",
      postgain: 1.0,
    },

    sim: super.sim {

        // For running in LArSoft, the simulation must be in fixed time mode. 
        fixed: true,

        // The "absolute" time (ie, in G4 time) that the lower edge of
        // of final readout tick #0 should correspond to.  This is a
        // "fixed" notion.
        local tick0_time = -250*wc.us,

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
        wires: "dunevdcrp2-wires-larsoft-v1.json.bz2",

        fields: [
            "dunevd-resp-isoc3views-18d92.json.bz2",
        ],

        noise: "protodune-noise-spectra-v1.json.bz2",


        chresp: null,
    },

}

