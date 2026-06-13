// PDHD charge-light (Q/L) matching helper.
//
// PDHD counterpart of cfg/pgrapher/experiment/sbnd/qlmatching.jsonnet.  Builds the
// matching graph nodes for ProtoDUNE-HD; the matching-only detector constants live
// here (not in params.jsonnet).
//
// Per drift-side group:  TensorFileSource -> FlashTensorToOpticalPCs -> QLMatching
//   - opflash_source   reads the opflash tensor archive (opflash_pdhd-wct.tar.gz).
//   - flash_attach     (Aux fan-in, 2->1) expands that matrix into the canonical
//                      flash / light / flashlight point clouds on the cluster root.
//   - matching         QLMatching reads those flash PCs and writes a per-cluster
//                      matched-flash scalar (so Cluster::get_flash() works).
//
// PDHD differs from SBND in three ways that matter here:
//   - 160 OpDets, all flat X-ARAPUCAs (active_opdet_types=[0], NOT the SBND PMT [1]);
//   - no reflected/VIS light (doReflectedLight=false; see pdhd/photodet README in
//     wire-cell-data), so the semi-analytical model needs only the VUV tables;
//   - central cathode at x=0 (the C++ default cathode_x), set in semi-analytical-pdhd.json.
//
// VUVEfficiency / QtoL are placeholders (uniform) sufficient to RUN the matching;
// the absolute visibility->PE scale still needs calibration against PDHD data/MC.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// trigger_offset: per-event readout-vs-trigger offset (opflash metadata offset_us,
// ~250us) folded into the matching geometry so the (offset-free, time_offset=0)
// charge x lands on the same trigger time base as the flash times.  Default 0 =>
// bit-identical (e.g. SBND, which bakes the offset into x_raw at imaging time).
function(params, trigger_offset=0 * wc.us) {
    // --- PDHD matching constants (matching-only) ---
    local nchan = 160,
    local ch_mask = [],   // no dead OpDets masked yet

    // visibility->PE efficiency.  Uniform placeholder for all 160 X-ARAPUCA windows;
    // VIS unused (reflected light off) but kept the right length for the predictor.
    local vuv_eff = 0.03,
    local VUVEfficiency = std.makeArray(nchan, function(i) vuv_eff),
    local VISEfficiency = std.makeArray(nchan, function(i) 0.0),

    // Opflash archive reader (the single per-event PDHD opflash file).  `inname`
    // is the full path to opflash_pdhd-wct.tar.gz (caller prefixes the input dir).
    opflash_source(tag, inname):: g.pnode({
        type: 'TensorFileSource',
        name: 'opflash_src_%s' % tag,
        data: {
            inname: inname,
            prefix: 'opflash_',
        },
    }, nin=0, nout=1),

    // Opflash matrix -> canonical flash/light/flashlight PCs (2->1 fan-in).
    // port 0 = cluster pctree, port 1 = opflash matrix.
    flash_attach(tag):: g.pnode({
        type: 'FlashTensorToOpticalPCs',
        name: 'flash_attach_%s' % tag,
        data: {
            nchan: nchan,
            // PDHD opflash has no per-frame CAF time offset -> use the raw flash time.
            correct_flash_time: false,
        },
    }, nin=2, nout=1),

    // Charge-light matching for one drift-side group.  `anodes` are the group's
    // APAs; anodes[0] is the representative (drift geometry / per-TPC OpDet mask --
    // every APA in a PDHD drift group shares one drift volume / TPC side, and its
    // ident, 0 or 1, selects the correct cathode side).  All group anodes are
    // registered on the Grouping (grouping_anodes) so blobs from every APA resolve.
    // `dv` is the group's DetectorVolumes node.  `calib_dump` (default '') is the
    // hand-scan calibration JSON path: when non-empty QLMatching writes the per-event
    // candidate-bundle universe there for the pdhd/ql_scan viewer; '' => off, no dump,
    // production output bit-identical.
    matching(anodes, dv, tag, face, calib_dump=''):: g.pnode({
        type: 'QLMatching',
        name: 'matching_%s' % tag,
        data: {
            anode: wc.tn(anodes[0]),
            grouping_anodes: [wc.tn(a) for a in anodes],
            tpc_face: face,   // imaging face of this drift side (0 = -x / APA0,2; 1 = +x / APA1,3)
            detector_volumes: wc.tn(dv),
            data: if std.objectHas(params, 'reality') && params.reality == 'sim' then false else true,
            QtoL: 1.0,
            doReflectedLight: false,
            drift_speed: params.lar.drift_speed,
            trigger_offset: trigger_offset,
            nchan: nchan,
            ch_mask: ch_mask,
            flash_minPE: 50,
            active_opdet_types: [0],   // X-ARAPUCA (flat), not the SBND PMT default [1]
            semimodel_file: 'pdhd/photodet/semi-analytical-pdhd.json',
            VUVEfficiency: VUVEfficiency,
            VISEfficiency: VISEfficiency,
            calib_dump: calib_dump,
        },
    }, nin=1, nout=1, uses=[dv] + anodes),
}
