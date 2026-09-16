// cfg/layers/mids/pdvd/variants/nominal.jsonnet
//
// BASE ("nominal") variant parameter object for *pdvd* (ProtoDUNE Vertical
// Drift).  Yields the schema expected throughout cfg/layers/mids/pdvd/.
//
// SCOPE: the layers/omnijob workflow (and test/scripts/spdir) is a single-APA
// job -- it uses anodes()[0].  For PDVD that is anode ident 0, a BOTTOM-drift
// CRP.  So this nominal describes the bottom drift: analytic ColdElecResponse at
// 7.8 mV/fC and the bottom ADC scaling.  The top drift (ident 4..7, a measured
// JsonElecResponse and [0,2] V ADC) is out of scope here; it is represented
// fully in the SPNG detconf (spng/detconfs/pdvd.jsonnet).
//
// Values follow the dunereco protodunevd oracle.

local wc = import "wirecell.jsonnet";
local base_variants = import "../../base/variants.jsonnet";

base_variants.nominal {

    lar: super.lar {
        // PDVD diffusion (base defaults DL 4.0 / DT 8.8 already match).
        lifetime: 1000.0 * wc.ms,
        drift_speed: 1.473 * wc.mm / wc.us,
    },

    geometry_data: import "geometry.jsonnet",

    // The VD field response begins 18.1 cm from the collection plane (vs the
    // base/HD 10 cm).  Override the response-plane distance the base derived.
    geometry+: {
        xplanes+: {
            dresponse: 18.1 * wc.cm,
        },
    },

    binning: {
        tick: 0.5 * wc.us,
        nticks: 6000,
    },

    // Bottom-drift electronics response.
    elec: {
        type: "cold",
        gain: 7.8 * wc.mV / wc.fC,
        shaping: 2.2 * wc.us,
        postgain: 1.0,
    },

    rc: {
        width: 1.1 * wc.ms,
    },

    // Bottom-drift digitization (14-bit ADC).
    digi: {
        gain: 1.0,
        baselines: [1003.4 * wc.millivolt, 1003.4 * wc.millivolt, 507.7 * wc.millivolt],
        resolution: 14,
        fullscale: [0.2 * wc.volt, 1.6 * wc.volt],
    },

}
