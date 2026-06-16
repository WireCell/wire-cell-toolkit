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
// readout_window_ticks: post-resample frame length used by the window-truncation
// flag.  PDHD's window (~5999) is far longer than the C++ default (SBND's 3427),
// so run_clus_evt.sh reads the real value from the SP frame and passes it in; the
// 6000 default is a sane PDHD fallback if it is not supplied.
function(params, trigger_offset=0 * wc.us, readout_window_ticks=6000) {
    // --- PDHD matching constants (matching-only) ---
    local nchan = 160,
    // Static optical dead-channel mask (OpChannel == OpDet, 0..159).
    //   - 3                      noisy (LArSoft IgnoreChannels)
    //   - 86,87,97,107,116,117   dead  (LArSoft v1 IgnoreChannels; confirmed 0 PE / 115 evts)
    //   - 120..159               DAPHNE full-stream readout, skipped by the snippet decoder
    //                            => 0 PE in every event; never enter the WCT opflash.
    // Run-to-run dead channels (e.g. the x<0 side, instrumented only in run 27980) are
    // caught per-event by auto_mask below, not statically.  See
    // pdhd/docs/pds-opchannel-opdet-mapping.md and pdhd-light-raw-data.md.
    local ch_mask = [3, 86, 87, 97, 107, 116, 117] + std.range(120, 159),

    // visibility->PE efficiency.  Uniform across all 160 X-ARAPUCA windows; VIS
    // unused (reflected light off) but kept the right length for the predictor.
    // 0.023 = first data calibration (run 27305 crosser anchors, lambda=100 cm):
    // direct-PMT scale ~0.77 x the old 0.03 placeholder.  Provisional (pinned mainly
    // by one bright crosser); see pdhd/docs/ql-light-normalization-study.md.
    local vuv_eff = 0.023,
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
            // Assemble the round-1/2 LASSO matrices sparse (block-sparse P/PF): on the
            // bright outlier (run 29107 evt 1015, ~440 flashes) the dense path's P/PT
            // spike and dense Gram build dominate QLMatching's time and memory. Sparse
            // products change the FP accumulation order, so the diagnostic LASSO strength
            // can move at the ULP level while every match assignment is unchanged (verified
            // byte-identical mabc output across all 30 evts of run 29107). C++ default OFF
            // (dense, byte-identical to history). See match/docs/qlmatching-perf-evt1015-pdhd.md.
            sparse_lasso: true,
            drift_speed: params.lar.drift_speed,
            trigger_offset: trigger_offset,
            nchan: nchan,
            ch_mask: ch_mask,
            // Per-event dynamic dead-channel auto-mask: drops a channel that never fires
            // this event while its nearest live neighbours do -- catches run-dependent
            // dead channels (x<0 side cabling varies by run) absent from the static
            // ch_mask above.  Thresholds tuned to PDHD: a live channel reaches >=30 PE
            // event-max (p1=31.6) while a dead one sits at 0, so pe_low=10 cleanly
            // separates them (cleared by 99.83% of live channels).  C++ default OFF.
            auto_mask: true,
            auto_mask_pe_low: 10,      // event-max PE below this => dead candidate
            auto_mask_pe_bright: 50,   // neighbour-median PE meaning "light present" (per-flash p75=55)
            auto_mask_neighbors: 4,
            auto_mask_min_contrast: 1,
            auto_mask_min_flash: 3,
            flash_minPE: 50,
            // No candidate-flash time cut: include EVERY flash from the event.
            // flash->get_time() is the RAW (trigger-relative) flash time.  The
            // C++/SBND default +-1.5 ms clipped the LATE half of PDHD's ~3 ms
            // readout; an earlier fix widened the window to the full readout,
            // but a finite window is still a cut.  Set bounds wider than any
            // conceivable PDHD readout (+-1 s) so the window can never exclude a
            // flash -- all flashes reach matching and are ranked there.  (On the
            // current data this is bit-identical to the full-readout window: 0 of
            // 707 flashes across the 23-event sample fell outside it.)  The
            // flash_minPE = 50 floor still applies.  PDHD-only; SBND keeps its
            // own qlmatching.jsonnet.
            flash_mintime: -1 * wc.s,
            flash_maxtime: 1 * wc.s,
            // Real PDHD readout window (post-resample SP frame length, ~5999),
            // supplied by run_clus_evt.sh from the SP frame.  Without it the C++
            // default (SBND's 3427, below the PDHD cathode at ~4498) would falsely
            // flag every mid-drift cluster as window-truncated.
            readout_window_ticks: readout_window_ticks,
            active_opdet_types: [0],   // X-ARAPUCA (flat), not the SBND PMT default [1]
            semimodel_file: 'pdhd/photodet/semi-analytical-pdhd.json',
            VUVEfficiency: VUVEfficiency,
            VISEfficiency: VISEfficiency,
            calib_dump: calib_dump,
        },
    }, nin=1, nout=1, uses=[dv] + anodes),
}
