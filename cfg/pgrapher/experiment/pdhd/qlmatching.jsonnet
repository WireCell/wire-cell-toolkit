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
    // Channels 120..159 are the DAPHNE full-stream PDs (the entire -x / side-0
    // z<250 half).  They USED to be dropped by the snippet decoder (0 PE) and were
    // statically masked here; the full-160 re-extraction (all-PD light reco) now
    // decodes them with real light (run 29107: ~38/40 fire, ~70k PE/evt), so they
    // are NO LONGER statically masked -- masking them discarded half of side 0.
    // Runs whose light reco lacks the full stream leave 120..159 at 0 PE; those are
    // caught per-event by auto_mask below (as are the run-dependent x<0 dead
    // channels), not statically.  See pdhd/docs/pds-opchannel-opdet-mapping.md and
    // pdhd-light-raw-data.md.
    local ch_mask = [3, 86, 87, 97, 107, 116, 117],

    // visibility->PE efficiency.  Uniform across all 160 X-ARAPUCA windows; VIS
    // unused (reflected light off) but kept the right length for the predictor.
    // 0.0145 = data calibration on run 29107 (28 clean two-boundary anode-cathode
    // crossers, lambda=300 cm): median direct-PMT scale 0.63 x the 0.023 prior.
    // Supersedes the 27305 tuning (0.023, lambda=100) -- that run had only ~5-7
    // anchors and over-concentrated the model.  See ql-light-normalization-study.md.
    // 0.01254 = run-29107 evt-983 hand-scan label retune (self-consistent with the
    // APA0 measured_pe_scale below): the +x side-1 anchor over-predicted by ~16%
    // (meas/pred 0.865) on the sizable+low-ks matches, so vuv_eff = 0.0145 x 0.865.
    // See ql-light-normalization-study.md / ql_light_calib/fit_labels.py.
    local vuv_eff = 0.01254,
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
    // Shared QLMatching `data` block -- everything independent of the per-node
    // anode/face wiring.  matching() (one drift side per node) and matching_joint()
    // (both drift sides in one node) each splice this in with their own anode/face keys.
    local match_data(dv, calib_dump) = {
            detector_volumes: wc.tn(dv),
            data: if std.objectHas(params, 'reality') && params.reality == 'sim' then false else true,
            QtoL: 1.0,
            doReflectedLight: false,
            // Per-channel MEASURED-PE gain correction (length nchan; 1.0 = identity,
            // C++ default empty = byte-identical for SBND/ICARUS). The -x full-data-
            // stream half ("APA0", optical ch 120-159) under-reports PE; scale its
            // MEASUREMENT up by 1.57 so it matches the (recalibrated) prediction.
            // 1.57 = g x median(pred_old/meas|APA0) = 0.865 x 1.814, self-consistent
            // with the vuv_eff retune above (run-29107 evt-983 sizable+low-ks labels;
            // the brighter ~2.3 tail is high-ks saturation, excluded). APA2 (ch80-119)
            // shows the same elevation but is left to the common model per the
            // APA0-only scope. ql_light_calib/fit_labels.py.
            measured_pe_scale: std.makeArray(nchan, function(i) if i >= 120 then 1.57 else 1.0),
            // Per-channel light-error model sigma = (PE<knee) ? floor : frac*PE.
            // frac 0.3 -> 0.44 from the evt-983 label residuals of the corrected model
            // (intrinsic per-PMT scatter on side1+APA0; floor/knee keep C++ defaults).
            pe_err_frac: 0.44,
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

            // --- Bundle prefilters (PDHD default ON; C++ defaults OFF, so any config
            // that does not set these is bit-identical) ----------------------------
            // (a) Containment: drop a (flash, cluster) bundle whose cluster is not
            // contained in the TPC drift box once the flash T0 x-offset is applied.
            // The box is the TWO-APA union (compute_geometry unions inner_bounds over
            // grouping_anodes above) -- the SAME box used for the at_x_boundary /
            // two_boundary flag walk -- so on each -x/+x drift side it spans both APAs
            // (z ~ 0..4.6 m, verified on run 29107: geometry z-span 462 cm). It also
            // prunes non-contained junk candidates from the pre-LASSO pool.
            require_containment: true,

            // Opaque-cathode mismatched-candidate filter. With an opaque cathode a flash
            // lit on one drift side cannot illuminate a cluster on the other side, so a
            // (cluster, opposite-side flash) bundle is non-physical -- EXCEPT for a
            // genuine cathode-crosser: at the flash's T0 the cluster's end reaches the
            // cathode (at_x_boundary), so its far half can be the source of the lit
            // side's light (whose own-side flash may be missing/dark). Keep same-side and
            // cross-side cathode-crossers; drop every other cross-side bundle. (Brightness
            // is irrelevant -- a bright crosser is as valid as a dim one; a mid-drift
            // cluster merely contained at some opposite flash's T0 is a coincidence.)
            // Cleans the pre-LASSO pool AND the hand-scan candidate tables (same test in
            // dump_calib). C++ default OFF (bit-identical for SBND/ICARUS).
            cross_side_filter: true,

            // cathode_ext1 / cathode_ext2: the cathode-end window [u_cathode+ext2,
            // u_cathode+ext1).  ext1 is the containment edge (how far PAST the cathode a
            // cluster's trimmed end may sit and still count as contained / at_x_boundary;
            // = high_x_cut_ext1, also the PE-inclusion edge, QLMatching.cxx:1082,2474).
            // ext2 is the FLAG-ONLY lower edge (how far SHORT of the cathode an end may
            // sit and still flag at_x_boundary, QLMatching.cxx:2506) -- it touches neither
            // containment nor PE inclusion.
            //
            // PDHD's drift residual is degenerate t0/velocity/SCE at ~+-2 cm
            // (project_cathode_crossing_offset): genuine cathode crossers scatter ~+-1.75 cm
            // around the cathode and that spread is irreducible by any single velocity.  We
            // set drift_speed=1.576 so the most-overshooting crosser sits just INSIDE the
            // cathode (run 29107 evt 983: clus62 +0.84, clus44 +0.33 cm) rather than past it,
            // which lets ext1 revert from the earlier 2.5 cm to 1.5 cm -- ~the C++ default
            // 1.2, just +0.3 cm of noise cushion (clus62 has 0.66 cm margin to the edge).
            // The undershoot residual then lands entirely on ext2: at 1.576 the worst
            // undershooter is clus96 at -2.61 cm, so ext2 widens from the C++ default -2.0
            // to -3.0 cm (0.4 cm margin) to keep flagging it at_x_boundary.  ext2 is the
            // benign lever -- widening it only labels more near-cathode ends as boundary
            // crossers (exempting them from reject_overpred); it cannot drop a bundle.
            cathode_ext1: 1.5 * wc.cm,
            cathode_ext2: -3.0 * wc.cm,

            // Robust containment endpoint.  compute_endpoint_flags' gap-based trims only
            // fire for an ISOLATED deep straggle (separated from the body by a >0.75 cm
            // gap); a thin off-axis OVERCLUSTERING tail that stays within 0.75 cm of the
            // body is never trimmed, so a few stray points dragging the drift endpoint
            // past the cathode falsely fail containment and silently drop that candidate
            // (run 29107 evt 983: cluster Bee-89/ident-55, an 11-point / 0.3% off-axis
            // tail reaching 2.3 cm past the cathode, killed every cathode-time candidate).
            // This gap-independent pass snaps each endpoint back inside the in-edge iff
            // the outer material is SPARSE -- its point mass < max(frac*cluster_points,
            // count).  Point mass (not blob count) is the sparsity measure, so a dense
            // genuine track-end / at_x_boundary crosser (hundreds of points => the cluster
            // is grossly mis-t0'd, not overclustered) exceeds the allowance and is left to
            // fail containment.  Verified on evt 983: rescues ident-55, and NO cluster
            // loses at_x_boundary (the calibration crossers keep their flags + strengths).
            // 1% OR 15 points: the 15-pt floor catches small clusters (where 1% < 1 pt),
            // the 1% scales for large ones.  C++ default OFF (robust_endpoint_trim=false)
            // => byte-identical for every other config (verified vs baseline, 0 bundles
            // differ).  See docs/qlmatching-chain.md.
            robust_endpoint_trim: true,
            robust_endpoint_frac: 0.01,
            robust_endpoint_count: 15,
            // Charge-weighted straggle judge (OR with the point-count one above).  Overclustered
            // satellites can carry many points yet a tiny charge fraction (evt 991 ident-31:
            // 91 pts past the anode, >5x the 15-pt count floor, but ~0.1% of cluster charge and
            // ~50x lower per-point q than the real track body), so the point-count judge misses
            // them and require_containment vetoes an otherwise-anode-piercing candidate.  Trim
            // outer material whose charge < 1% of cluster charge; a dense genuine track-end
            // exceeds it and is left to fail.  C++ default 0.0 => disabled (byte-identical).
            robust_endpoint_charge_frac: 0.005,
            // Absolute per-point charge-density ceiling (size-independent): only trim
            // charge-sparse outer material (diffuse overclustered satellites, ~150 q/pt),
            // never a charge-dense real track tip into the boundary (~2500-8800 q/pt).  This
            // is what stops the trim from shaving a genuine track tail to force containment,
            // independent of cluster size.  C++ default 0 => disabled.
            robust_endpoint_charge_abs: 1500,

            // (b) Over-prediction reject: before the chi2 fit, drop a bundle whose
            // predicted light hugely exceeds the measured light over the masked PMT set:
            //   reject if  sum(pred)/sum(meas) > overpred_total_ratio
            //          or  pred/meas at the brightest predicted PMT > overpred_maxch_ratio
            // Boundary/window-truncated bundles are EXEMPT (incl. every at_x_boundary
            // two-boundary light anchor), so this never touches the calibration set.
            // Both cuts are PRE-LASSO (they prune candidates, triggering a global
            // re-solve), so their effect is NOT the per-bundle cull count -- on run
            // 29107 they swap ~44% of the (noise-dominated) matched flashes while keeping
            // the count flat (1399->1428) and IMPROVING quality (the removed matches had
            // median best-KS ~0.87 = junk; stable matches held at paired median dKS=0;
            // strict ks<0.2 two-boundary anchors grew 3->14). See docs/qlmatching-chain.md.
            // GT-tightened ceilings.  The run-29107 evt-983 hand scan now provides PDHD
            // ground truth: of its 31 labeled matches, 24 are boundary/window-truncated
            // (overpred-EXEMPT) and 7 are subject to this cut; their worst values are
            // R_total=1.90 and R_max=5.04 (mc 1000020, a real match that over-predicts one
            // PMT 5x).  The optical-model retune (vuv_eff 0.0145->0.01254 + APA0 measured
            // scale) also tightened genuine over-prediction, so the earlier deliberately-
            // LOOSE provisional ceilings (5.0, 25.0 -- sized above the rough model's range)
            // are no longer needed.  (3.0, 10.0) keeps every GT match with ~1.6x / 2.0x
            // margin over the worst GT case, while culling the egregious tail (auto-selected
            // winners run out to R_max~18).  R_max is NOT dropped below ~8: GT mc 1000020
            // legitimately over-predicts one channel 5x, and this is single-event GT, so 2x
            // headroom is kept.  Verified on evt 983: all 31 GT matches survive.  SBND keeps
            // its own (2.9/4.3, 10 GT events).  See ql_light_calib/containment_overpred_check.py.
            reject_overpred: true,
            overpred_total_ratio: 3.0,
            overpred_maxch_ratio: 10.0,

            // --- Flag fit-advantage (SBND method; default OFF in C++, byte-identical
            // for any config that omits them) ----------------------------------------
            // (a) LASSO L1 down-weight: a boundary/near-PD/window-truncated bundle (incl.
            // the xTPC cathode-crossers) has its measured light UNDER-reported (charge
            // drifts past the cathode / out the readout window), so a raw PE-mismatch L1
            // penalty wrongly shrinks it below the strength cutoff. Multiplying its L1
            // weight by lasso_boundary_weight (0.2 = C++/SBND default) lets it survive.
            lasso_flag_weight: true,
            lasso_boundary_weight: 0.2,
            // (b) Per-bundle chi2 relaxation (prototype FlashTPCBundle.cxx:480): widens the
            // chi2 denominator for a near-PD channel with a big measured excess over
            // prediction, and drops the single worst-chi2 channel if it is a dead PD
            // (pe==0, pred>0). The excess-widening thresholds (chi2_pmt_excess/ratio/inflate,
            // C++ defaults 350/1.3/0.5) are SBND-PMT-scaled and largely inert at PDHD
            // ARAPUCA PE levels (no channel reaches a 350-PE excess), so the ACTIVE effect
            // here is the dead-PD worst-channel drop (detector-agnostic); the excess term is
            // left for a later PDHD retune. Verified on evt 983: all 31 hand-scan matches
            // survive both levers. See docs/qlmatching-chain.md.
            chi2_relax: true,

            active_opdet_types: [0],   // X-ARAPUCA (flat), not the SBND PMT default [1]
            semimodel_file: 'pdhd/photodet/semi-analytical-pdhd.json',
            VUVEfficiency: VUVEfficiency,
            VISEfficiency: VISEfficiency,
            calib_dump: calib_dump,
    },

    // One-drift-side matcher (historical PDHD path; kept for back-compat).  `anodes`
    // are the group's APAs; anodes[0] is the representative drift geometry / OpDet-mask
    // side; `face` its imaging face.
    matching(anodes, dv, tag, face, calib_dump=''):: g.pnode({
        type: 'QLMatching',
        name: 'matching_%s' % tag,
        data: {
            anode: wc.tn(anodes[0]),
            grouping_anodes: [wc.tn(a) for a in anodes],
            tpc_face: face,   // imaging face of this drift side (0 = -x / APA0,2; 1 = +x / APA1,3)
        } + match_data(dv, calib_dump),
    }, nin=1, nout=1, uses=[dv] + anodes),

    // JOINT both-sides matcher (SBND-style): both drift volumes enter ONE node so the
    // cross-cathode (xTPC) consistency pass can pair a cathode-crosser's two halves (one
    // per side).  `sides` = per-side anode lists ([[APA0,APA2],[APA1,APA3]]); each side's
    // representative (sides[i][0]) is input port i and faces[i] its imaging face.  Both
    // grouping_anodes and tpc_faces are PER-INPUT, so each run's active-volume box and
    // wire face come from its OWN side.  Emits one merged tree + one combined calib JSON.
    matching_joint(sides, dv, faces, calib_dump=''):: g.pnode({
        type: 'QLMatching',
        name: 'matching_joint',
        data: {
            anodes: [wc.tn(side[0]) for side in sides],
            grouping_anodes: [[wc.tn(a) for a in side] for side in sides],
            tpc_faces: faces,
            // Merge per-side optical-flash display PCs into the joint grouping (mirrors
            // clus_all_tpc PointTreeMerging.root_pcs_to_merge) -> one flash list.
            root_pcs_to_merge: ['opflash'],
            // Cross-cathode cathode-crossing consistency flag pass (needs both sides in
            // one node): pairs coincident at_x_boundary / window-truncated halves across
            // the central cathode (drift side 0 vs side 1).  CONSUMED by matching, not
            // observation-only: cull_cross_tpc flags the coincident crosser pairs (run
            // 29107 evt 983: 349 coincident pairs) and the cull_inconsistent ladder then
            // keeps a cluster's scenario-1 xTPC crosser over its rivals (same code as
            // SBND).  Geometry cuts SBND-seeded; the flag fit-advantage (lasso_flag_weight
            // + chi2_relax above) is now also enabled, matching SBND's method.
            xtpc_flag: true,
            xtpc_dmax: 5 * wc.cm,
            xtpc_dmax2: 300 * wc.cm,
            xtpc_angle_max: 20,
            xtpc_hough_radius: 15 * wc.cm,
        } + match_data(dv, calib_dump),
    }, nin=std.length(sides), nout=1, uses=[dv] + [a for side in sides for a in side]),
}
