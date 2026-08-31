// ProtoDUNE-SP specific parameters.  This file inerets from the
// generic set of parameters and overrides things specific to PDSP.

local wc = import "wirecell.jsonnet";
local base = import "pgrapher/dune/params.jsonnet";

// Bottom-drift electronics-noise file, selected by the front-end gain.  The
// cold electronics has four gain settings -- 4.7 / 7.8 / 14 / 25 mV/fC -- and
// noise-spectra files currently exist only for 7.8 and 14.  Any other gain
// (a setting with no file, or a value that is not a valid setting) aborts the
// configuration rather than silently using a wrong-gain spectrum.  See
// pdvd/nf_plot/electronics_gain_and_noise.md.
local pdvd_bottom_noise(gain) =
    local g = gain / (wc.mV / wc.fC);
    if std.abs(g - 7.8) < 0.05 then "pdvd-bottom-noise-spectra-7d8mVfC-v1.json.bz2"
    else if std.abs(g - 14.0) < 0.05 then "pdvd-bottom-noise-spectra-14mVfC-v1.json.bz2"
    else error ("PDVD bottom noise: no spectra file for elec.gain = " + g
                + " mV/fC.  Valid cold-electronics gain settings are"
                + " 4.7/7.8/14/25 mV/fC; spectra files exist only for 7.8 and 14.");

base {
    // This section will be overwritten in simparams.jsonnet
    det : {

        // See:  wirecell-util wire-volumes protodunevd-wires-larsoft-v6.json.bz2
        // to help with defining these parameters.  (v6 = default wire file, see
        // files.wires below and pdvd/docs/qlmatch/pdvd-crp-anode-plane-geometry.md)

        // between center lines
        local apa_cpa = 341.55*wc.cm,
        // GDML CathodeBlock (protodunevd_v4_refactored.gdml); corrected 2026-07
        // from the legacy 50.8 mm (= DocDB 203 / ProtoDUNE-SP nominal, which no
        // PDVD GDML uses).  Drift-facing surface -> |x| = 3.0 cm, drift distance
        // cpa_plane = 338.55 cm ~ GDML CRMActive 338.5.  See
        // pdvd/docs/pdvd-tpc-geometry-fiducial.md.
        local cpa_thick = 60.0*wc.mm,
        local apa_w2w = 85.725*wc.mm,
        local plane_gap = 4.76*wc.mm,
        local apa_g2g = 114.3*wc.mm, 

        // FV / sensitive-volume anode edge.  2026-07-13: moved to the PHYSICAL
        // CRP shield plane, the drift-facing boundary of the active LAr.  The
        // confirmed CRP stack (W collection fixed) is W -3.2mm-> V -10mm-> U
        // -3.2mm-> Shield, so the shield sits 16.4 mm below W (toward the drift
        // volume); the 16.4 mm W..shield is PCB stack, not active drift region.
        // anode = +-(341.55 - 1.64) = +-339.91 cm.  Supersedes the 2026-07-12
        // U-plane choice (2 x 0.2 mm = +-341.51 cm), which used the v5 wire
        // file's LArSoft 0.2 mm-step convention.  Drives DetectorVolumes::
        // inner_bounds -> QLMatching FV/u-coordinate + boundary flags and the
        // inner_bounds-based clustering.  NOT byte-identical.  See
        // pdvd/docs/qlmatch/pdvd-crp-anode-plane-geometry.md.
        local apa_plane = 16.4*wc.mm, // W -> shield plane (3.2 + 10 + 3.2 mm)

        // The "response" plane is where the field response functions
        // start.  Garfield calcualtions start somewhere relative to
        // something, here's where that is made concrete.  This MUST
        // match what field response functions also used.
        response_plane: 18.1*wc.cm, // relative to collection wires
                                    // synced to params.files.fields
        local res_plane = 0.5*apa_w2w + self.response_plane,

        // The cathode plane is like the anode cut off plane.  Any
        // depo not between the two is dropped prior to drifting.
        local cpa_plane = apa_cpa - 0.5*cpa_thick,

        volumes: [
            {
                // local world = 100,
                // local split = s*10, // 1: left, 2: right
                local anode = a, // physical anode number
                // wires: world + split + anode,
                // name: "anode%d"%(world + split + anode),
                wires: anode,
                name: "anode%d" %a,

                local sign = if a>3 then 1 else -1,
                local centerline = sign * apa_cpa,
                faces:
                // top drift volume
                if sign > 0
                then [
                    {
                        anode: centerline - apa_plane,
                        response: centerline - res_plane,
                        cathode: centerline - cpa_plane, 
                    },

                    {
                        anode: centerline - apa_plane,
                        response: centerline - res_plane,
                        cathode: centerline - cpa_plane, 
                    }
                ]
                // bottom drift volume
                else [
                    {
                        anode: centerline + apa_plane,
                        response: centerline + res_plane,
                        cathode: centerline + cpa_plane, 
                    },

                    {
                        anode: centerline + apa_plane,
                        response: centerline + res_plane,
                        cathode: centerline + cpa_plane, 
                    }
                ],
            } for a in std.range(0,7)
        ],

        // This describes some rough, overall bounding box.  It's not
        // directly needed but can be useful on the Jsonnet side, for
        // example when defining some simple kinematics.  It is
        // represented by a ray going from extreme corners of a
        // rectangular solid.  Again "wirecell-util wires-info" helps
        // to choose something.
        bounds : {
            tail: wc.point(-3.15, -3.42, 0, wc.m),
            head: wc.point(3.13, 3.42, 3.04, wc.m),
        }
    },

    lar: super.lar {
        // Calibrated from PDVD data: anode->cathode crossing tracks' reconstructed
        // drift x-span (both drift volumes, 142 events) vs the collection-plane ->
        // cathode-surface distance gives v = v_reco * D / S.  The cathode surface
        // moved from 2.54 -> 3.0 cm (cpa_thick 50.8 -> 60 mm, GDML), so D:
        // 339.01 -> 338.55 cm and v rescales 1.57 -> 1.568 (v proportional to D).
        // Consistent with the Walkowiak nominal at the PDVD field.  See
        // pdvd/docs/clus-workflow.md (drift-velocity calibration).
        drift_speed: 1.568 * wc.mm / wc.us,  // was 1.57 (D=339.01); 1.6 orig default
    },

    daq: super.daq {
        nticks: 6000,
    },

    // Real data frames are 8000 ticks after the 512->500 ns Resampler
    // (7813 raw ticks * 512/500 -> 8000).  OmnibusNoiseFilter now pushes
    // the actual frame size into OmniChannelNoiseDB on first frame, so this
    // override is no longer load-bearing — but it is kept as a documented
    // default and to pre-size freqbinner() calls in chndb-base.jsonnet.
    nf: super.nf {
        nsamples: 8000,
    },

    adc: super.adc {

        resolution: 14,

        // reuse ProtoDUNE SP values for bottom drift
        baselines: [1003.4*wc.millivolt,1003.4*wc.millivolt,507.7*wc.millivolt],
        fullscale: [0.2*wc.volt, 1.6*wc.volt],
        // will rewrite top drift values in sim/digitizer and sigproc
        // FIXME: need a more elegant way
        // top drift: ~1 volt baselines, 0-2volt full scale
    },

    // This sets a relative gain at the input to the ADC.  Note, if
    // you are looking to fix SimDepoSource, you are in the wrong
    // place.  See the "scale" parameter of wcls.input.depos() defined
    // in pgrapher/common/ui/wcls/nodes.jsonnet.
    // also, see later overwriting in simparams.jsonnet
    elecs: [
      super.elec { // bottom drifter
        postgain: 1.0,
        shaping: 2.2 * wc.us,
        gain:7.8*wc.mV/wc.fC,
      },
      super.elec { // top
        type: "JsonElecResponse",
        filename: "dunevd-coldbox-elecresp-top-psnorm_400.json.bz2",
        postgain: 1.36, // 11mV/fC, 1.94 -> 14mV/fC
      },
    ],
    elec: $.elecs[0], // nominal 

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
        // v6 = v5 with the U/V induction planes moved to their PHYSICAL CRP
        // spacing (W kept fixed at |x|=341.55 cm): V at W-3.2 mm (341.23),
        // U at W-13.2 mm (340.23), vs the v5 LArSoft 0.2 mm-step convention.
        // See pdvd/docs/qlmatch/pdvd-crp-anode-plane-geometry.md.  NOT
        // byte-identical to v5 reco (sigproc/imaging induction-plane geometry).
        wires: "protodunevd-wires-larsoft-v6.json.bz2",
        strip_length: "PDVD_strip_length.json.bz2",

        fields: [
            // "dunevdcrp2-FR-fixcoll-adjustind.json.bz2",
            // "dunevdcrp2-FR-fixcoll-adjustind.json.bz2", // repeat for top drifter
            // "protodunevd_FR_3view_speed1d55.json.bz2", // remember to sync response plane position above
            // "protodunevd_FR_3view_speed1d55.json.bz2",
            "protodunevd_FR_imbalance3p_260501.json.bz2",
            "protodunevd_FR_imbalance3p_260501.json.bz2",
        ],

        // Electronics-noise spectra.  Bottom: selected by front-end gain
        // (`pdvd_bottom_noise`, above -- set the gain in $.elecs[0]); the
        // 7.8 mV/fC file is data-retuned from run039324, see
        // noise_spectrum_comparison.md.  Top: single gain setting.
        noises: [
            pdvd_bottom_noise($.elec.gain),
            "pdvd-top-noise-spectra-v3.json.bz2",
        ],

        // chresp: "protodunevd-params-channel-responses-v0.json.bz2",
        chresp: null,
    },

}

