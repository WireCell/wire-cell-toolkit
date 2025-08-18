// Here we override params.jsonnet to provide simulation-specific params.

local base = import 'params.jsonnet';
local wc = import 'wirecell.jsonnet';

base {

  // redefine the detector volumes with the cryostat side included
  det : {

      // The "faces" is consumed by, at least, the Drifter and
      // AnodePlane.  The "wires" number is used to set
      // AnodePlane.ident used to lookup the anode in WireSchema.
      // It corresponds to the anode number.

      // Numbers determined from protodunehd_v6 gdml
      local apa_cpa = 3.5734*wc.m,
      local cpa_thick = 3.175*wc.mm, // 1/8", from Bo Yu (BNL) and confirmed with LArSoft
      local apa_w2w = 85.87*wc.mm,
      local plane_gap = 4.76*wc.mm,
      local apa_g2g = apa_w2w + 6*plane_gap,

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
          } for n in std.range(0,3)],

      // This describes some rough, overall bounding box.  It's not
      // directly needed but can be useful on the Jsonnet side, for
      // example when defining some simple kinematics.  It is
      // represented by a ray going from extreme corners of a
      // rectangular solid.  Again "wirecell-util wires-info" helps
      // to choose something.
      bounds : {
          tail: wc.point(-4.0, 0.0, 0.0, wc.m),
          head: wc.point(+4.0, 6.1, 7.0, wc.m),
      }
  },

  daq: super.daq {
    nticks: 6000,
  },

  sim: super.sim {

    // For running in LArSoft, the simulation must be in fixed time mode.
    fixed: true,
    continuous: false,
    fluctuate: true,

  },

  sys_status: false,
  sys_resp: {
    // overall_short_padding should take into account this offset "start".
    start: -10 * wc.us,
    magnitude: 1.0,
    time_smear: 1.0 * wc.us,
  },

  rc_resp: {
    width: 1.1*wc.ms,
    rc_layers: 0,
  }
}
