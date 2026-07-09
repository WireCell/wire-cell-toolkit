// PDVD charge-light (Q/L) matching constants -- SKELETON, NOT PRODUCTION.
//
// PDVD has no Q/L matching chain wired yet (clus.jsonnet runs without it);
// this file collects the PDVD-decided QLMatching light-model constants so the
// future chain (opflash_source -> flash_attach -> matching, as in
// pdhd/qlmatching.jsonnet) can import them.  Everything light-model-related
// below is derived and validated in pdvd/photlib/ (wcp-porting repo) and
// documented in pdvd/docs/pdvd-photon-model.md:
//
//  - light_model 'library': per-point visibilities interpolated from a 10 cm
//    grid sampled off the official PDFastSimANN v5 computable graph
//    (protodune_vd_v5_128nm_tf2.6) -- the as-built geometry that matches the
//    raw-data opdet positions exactly.  Library channel == flash-chain OpDet.
//  - semimodel_file: the semi-analytical FIT to the same ANN (fallback /
//    cross-check; cathode XAs fit to ~13%, PMTs ~40%, membrane XAs need the
//    unported lateral-cosine branch -- see the JSON's _comment).
//  - The library mode needs no same-TPC x-sign gate: the ANN encodes cathode
//    opacity and the double-sided cathode XAs.
//
// Still to calibrate before production: offset_us (light<->charge time base),
// QtoL / absolute PE scale, per-PD gains, and the DAPHNE<->module assignment
// check for Arapuca mirror pairs (jjo map vs PDS_Mapping_v04152025 -- see
// pdvd-photon-model.md "open items").

local wc = import 'wirecell.jsonnet';

{
    nchan: 40,

    // Dead in data: 24/27/28/34.  Ar-blind (official eff_Ar=0, no/quenched
    // WLS at 128nm): 16 (membrane, no PTP), 29/39 (PEN+Q), 32 (no coating).
    // Keep the Ar-blind ones masked for pure-Ar running; unmask 16/29/39 for
    // Xe-doped running with the 175nm library.
    ch_mask: [16, 24, 27, 28, 29, 32, 34, 39],

    // Official per-OpDet detection efficiencies (PDVD_PDS_Mapping_v04152025):
    // XAs 0.03, TPB-coated PMTs 0.12, PEN PMTs 0.036, Ar-blind 0.  The
    // absolute PE scale on data still needs a crosser calibration on top
    // (PDHD-style); these set the relative PD-type weighting.
    VUVEfficiency: [
        0.03, 0.03, 0.03, 0.03,                  // 0-3   membrane XA
        0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,  // 4-11 cathode XA
        0.12, 0.036, 0.12, 0.12,                 // 12-15 TCO PMTs
        0.0, 0.03,                               // 16-17 membrane XA (16 no-WLS)
        0.12, 0.036, 0.12, 0.12,                 // 18-21 TCO PMTs
        0.03, 0.03,                              // 22-23 membrane XA
        0.036, 0.036, 0.036, 0.036, 0.036, 0.0,  // 24-29 bottom PMTs (29 PEN+Q)
        0.036, 0.036, 0.0, 0.036, 0.036, 0.036,  // 30-35 (32 no coating)
        0.036, 0.036, 0.036, 0.0,                // 36-39 (39 PEN+Q)
    ],
    VISEfficiency: std.makeArray(40, function(i) 0.0),

    match_data: {
        nchan: $.nchan,
        ch_mask: $.ch_mask,
        // both flat XAs (type 0) and PMTs (type 1) are active PDs
        active_opdet_types: [0, 1],
        doReflectedLight: false,   // library vis is total photon arrival
        QtoL: 1.0,                 // placeholder until PE-scale calibration
        VUVEfficiency: $.VUVEfficiency,
        VISEfficiency: $.VISEfficiency,

        light_model: 'library',
        photon_library_file: 'pdvd/photodet/pdvd-photlib-vis-v5-128nm.json',
        semimodel_file: 'pdvd/photodet/semi-analytical-pdvd.json',
    },
}
