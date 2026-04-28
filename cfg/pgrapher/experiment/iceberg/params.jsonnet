local base = import 'pgrapher/dune/params.jsonnet';
local wc = import 'wirecell.jsonnet';

function(params) base {
  det: {

    local apa_cpa = 360.299 * wc.mm, // between center lines, from Shekar Mishra
    local cpa_thick = 3.175 * wc.mm, // 1/8 inch
    local apa_w2w = 3.00155 * 2 * wc.cm,
    local plane_gap = 4.76 * wc.mm,
    local apa_u2u = apa_w2w + 4*plane_gap,

    local apa_plane = 0.5 * apa_u2u,  // pick it to be at the fist induction wires

    // The "response" plane is where the field response functions
    // start.  Garfield calcualtions start somewhere relative to
    // something, here's where that is made concrete.  This MUST
    // match what field response functions also used.
    response_plane: 10 * wc.cm,  // relative to collection wires
    local res_plane = 0.5 * apa_w2w + self.response_plane,

    // The cathode plane is like the anode cut off plane.  Any
    // depo not between the two is dropped prior to drifting.
    local cpa_plane = apa_cpa - 0.5 * cpa_thick,


    volumes: [
      {
        wires: 0,  // anode number
        name: 'apa0',
        // face 0 receives +x depos
        faces:
          [
            {
              anode:    apa_plane,
              response: res_plane,
              cathode:  cpa_plane,
            },
            {
              anode:    -apa_plane,
              response: -res_plane,
              cathode:  -cpa_plane,
            }
          ],
      }
    ],

    // This describes some rough, overall bounding box.  It's not
    // directly needed but can be useful on the Jsonnet side, for
    // example when defining some simple kinematics.  It is
    // represented by a ray going from extreme corners of a
    // rectangular solid.  Again "wirecell-util wires-info" helps
    // to choose something.
    bounds: {
      tail: wc.point(-30,  0.0,   0.0, wc.cm),
      head: wc.point( 30, 93.8, 115.3, wc.cm),
    },
  },

  daq: super.daq {
    nticks: 8256,
  },

  adc: super.adc {
    baselines: [900 * wc.millivolt, 900 * wc.millivolt, 200 * wc.millivolt],

    // The resolution (bits) of the ADC
    resolution: 14,

    // The voltage range as [min,max] of the ADC, eg min voltage
    // counts 0 ADC, max counts 2^resolution-1.
    fullscale: [0.0 * wc.volt, 1.4 * wc.volt],
  },

  elec: super.elec {
    postgain: 1.0,
    shaping: 2.0 * wc.us,
  },


  sim: super.sim {

    // For running in LArSoft, the simulation must be in fixed time mode.
    fixed: true,

    // The "absolute" time (ie, in G4 time) that the lower edge of
    // of final readout tick #0 should correspond to.  This is a
    // "fixed" notion.
    tick0_time: if std.objectHas(params, 'G4RefTime') then params.G4RefTime else 0,

    // Open the ductor's gate a bit early.
    local response_time_offset = $.det.response_plane / $.lar.drift_speed,
    local response_nticks = wc.roundToInt(response_time_offset / $.daq.tick),

    ductor: {
      nticks: $.daq.nticks + response_nticks,
      readout_time: self.nticks * $.daq.tick,
      start_time: $.sim.tick0_time - response_time_offset,
    },

    // To counter the enlarged duration of the ductor, a Reframer
    // chops off the little early, extra time.  Note, tags depend on how
    reframer: {
      tbin: response_nticks,
      nticks: $.daq.nticks,
    },

  },

  files: {
    wires: 'iceberg-wires-larsoft-v2.json.bz2',

    fields: [
      'dune-garfield-1d565.json.bz2',
    ],

    wiregroups: "iceberg_group_to_channel_map_incoh.json.bz2",
    noise: "iceberg_noise_model_incoh.json.bz2",
    // noise: 'protodune-noise-spectra-v1.json.bz2',


    chresp: null,
  },

}
