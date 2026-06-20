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
    // 0.01281 = run-29107 4-event hand-scan label retune (evts 983/991/999/1007, 25
    // flag-clean sizable+low-ks anchors; flag_PMT/xboundary/wtrunc dropped). The +x/
    // apa1 INTEGRAL meas/pred g=1.022 at lambda=300 (concentration-matched; the spread
    // is N90-correct so the integral is the clean norm handle), vuv_eff = 0.01254 x
    // 1.022. Supersedes the evt-983-only 0.01254. (KS rides to the grid edge on these
    // non-crosser anchors and cannot refit lambda; N90+integral keep lambda=300.)
    // See ql-light-normalization-study.md / ql_light_calib/fit_labels_multi.py.
    local vuv_eff = 0.01281,
    local VUVEfficiency = std.makeArray(nchan, function(i) vuv_eff),
    local VISEfficiency = std.makeArray(nchan, function(i) 0.0),

    // --- Per-channel (per-PD) MEASURED-PE gain calibration ---
    // The optical SP chain deconvolves all 160 PDs with just TWO SPE templates
    // (FBK / HPK, one shape per type); residual per-channel SiPM gain is not removed
    // and biases the bundle chi2/KS.  measured_pe_scale (Opflash::init) absorbs it.
    // Fit from the cleanest run-29107 hand-scan matches (54 low-ks, no flag_PMT/
    // flag_wtrunc; ql_light_calib/fit_perchannel_scale.py).  Granularity = grouped
    // block x type base + individual breakout only for well-sampled, tight-scatter
    // outliers (the prior per-channel SPE-shape study overfit -- pdhd-spe-template-
    // tuning.md; grouping generalises).  fbk_ch = template_index 0 in
    // pdhd-spe-templates.json; the rest are HPK.
    local fbk_ch = [4, 14, 24, 34, 40, 42, 45, 46, 47, 49, 50, 52, 55, 56, 57, 59,
                    60, 62, 65, 66, 67, 69, 70, 72, 75, 76, 77, 79, 84, 85, 86, 87,
                    94, 95, 96, 97, 104, 105, 106, 107, 114, 115, 116, 117, 120, 121,
                    124, 125, 127, 129, 130, 131, 134, 135, 137, 139, 140, 141, 144,
                    145, 147, 149, 150, 151, 154, 155, 157, 159],
    // Block x type base = current_block_scale x light-weighted median(pred/meas) over
    // each group's clean anchors.  The dominant effect is per-TYPE: FBK reads ~12-21%
    // low everywhere (the documented FBK tail over-subtraction) -> FBK scaled up.  The
    // old uniform APA0 (ch120-159) 1.14 splits cleanly into FBK 1.20 / HPK 1.00 -- the
    // block average conflated an FBK-only defect with HPK that needs no boost.  The +x
    // (ch0-79) bases are renormalised (meas-weighted mean held at 1.0, k=1.128) so the
    // +x integral -- hence vuv_eff -- is unchanged; APA2 (ch80-119) ~1.0; -x bases
    // subsume/refine the old 1.14/1.0.
    local pd_scale = function(i)
        local fbk = std.member(fbk_ch, i);
        if i >= 120 then (if fbk then 1.20 else 1.00)            // APA0 full-stream
        else if i >= 80 then (if fbk then 0.98 else 0.96)       // APA2
        else (if fbk then 1.12 else 0.98),                      // +x (norm-preserving)
    // Individual overrides: channels >3 standard-errors AND >0.20 off their group base,
    // well sampled (N>=12) and tight (MAD<0.25) -- a stable per-channel gain, not anchor
    // scatter.  Mostly high-gain HPK PDs reading ~1.5-2x too high (scale DOWN); ch23 is
    // a low-gain HPK (UP).  The 3-5x factors the raw fit wanted on ch40/50/60/70 are
    // statistical artifacts (N<=5, MAD up to 12) -> left at group default, NOT scaled.
    local perch_override = {
        '20': 0.71, '23': 1.47, '25': 0.59, '30': 0.73, '35': 0.68, '48': 0.70,
        '88': 0.60, '98': 0.56, '108': 0.59, '118': 0.53, '139': 0.67, '149': 0.73,
    },

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
            // C++ default empty = byte-identical for SBND/ICARUS).  Built from the
            // block x type base (pd_scale) with sparse per-channel overrides above.
            // Supersedes the uniform-1.14-APA0 retune (which conflated the FBK defect
            // with HPK).  See ql_light_calib/fit_perchannel_scale.py + the per-channel
            // calibration section of ql-light-normalization-study.md.
            measured_pe_scale: std.makeArray(nchan, function(i)
                if std.objectHas(perch_override, std.toString(i))
                then perch_override[std.toString(i)]
                else pd_scale(i)),
            // Per-channel light-error model sigma = (PE<knee) ? floor : frac*PE.
            // frac = 0.40: the high-PE method-of-moments value (E[(pred-meas)^2]=meas+frac^2*pred^2
            // on the bright tail of the run-29107 hand-accepted matches). We do NOT calibrate the
            // bundle chi2/ndf to ~1 (the earlier frac=0.60 did that -- median chi2/ndf 1.66->1.13 --
            // but chi2/ndf~1 is not required here, and the larger frac is a looser matching error
            // that over-tolerates mismatch; the matcher is KS-led, so the intrinsic per-PMT 0.40 is
            // the right scale).  floor/knee + the low-PE inflation below are unchanged.
            // ql_light_calib/fit_perchannel_scale.py + validate_chain.py.
            pe_err_frac: 0.40,
            // Compute the bundle-chi2 per-PMT error from the PREDICTED pe (not the
            // measured-based flash error). Required by the low-PE inflation below, and
            // on its own already cures the catastrophic "predicted light, measured ~0"
            // penalty the measured-based branch gives (perr=floor there -> chi2 ~ pred^2).
            pe_err_on_pred: true,
            // Efficiency-aware low-PE error inflation: PDs detect zero far more often than
            // Poisson when little light is predicted (run-29107 hand scans: ~50% of pred
            // 2-5 PE channels measure 0). Grow the relative error as pred falls so these
            // are tolerated: rel = frac + (lowpe_frac-frac)*exp(-pred/lowpe_knee). Tuned so
            // a typical measured-zero channel contributes ~1 to chi2 (median 5.2 -> 1.0)
            // while pred>=50 is untouched. ql_light_calib/fit_lowpe.py on evts 983/991/999/1007.
            pe_err_lowpe_frac: 1.55,
            pe_err_lowpe_knee: 5.5,
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

            // Perf only (matching bit-identical): apply the cross_side_filter drop BEFORE
            // the per-point visibility loop instead of after it, skipping the
            // SemiAnalyticalModel evaluation of cross-side bundles that get discarded
            // anyway.  After the A+B+B2 LASSO work the vis_loop is ~90% of QLMatching on
            // busy events; on evt 1015 this cuts vis_loop ~38% and the process wall ~28%
            // with the surviving candidate set (pre_bundles) unchanged.  C++ default OFF
            // (legacy post-loop drop, byte-for-byte the old path).  See
            // match/docs/qlmatching-perf-evt1015-pdhd.md sec 9-10.
            crossside_skip_vis: true,

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
            // Gap-aware (detachment) anode trim.  Some overclustered junk is moderately
            // dense (~2000 q/pt, e.g. evt 1007 uid-126: 7 gap-separated groups marching
            // from u=-111cm up to the body) and so overlaps the genuine track-end band
            // the charge_abs gate protects -- the density judge cannot trim it.  But it
            // is DETACHED from the contiguous body by large drift gaps (a real anode-
            // piercing track end is continuous, slices ~0.08cm apart).  Trim sub-anode
            // material that contains a gap > 3cm and carries < 1% of cluster charge.
            // C++ default 0 => disabled (byte-identical).
            robust_endpoint_gap: 3.0 * wc.cm,
            robust_endpoint_gap_charge_frac: 0.01,

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
            // chi2_relax excess-widening threshold, re-scaled from the SBND/MicroBooNE 350 PE
            // to a PDHD ARAPUCA level.  Path-1 widens the chi2 denom of a close_to_PMT channel
            // whose MEASURED pe exceeds prediction by > chi2_pmt_excess PE AND pe > ratio*pred
            // (un-modeled direct light / saturation the semi-analytical model under-predicts).
            // PDHD's lower light yield means genuine excess is smaller in absolute PE than
            // SBND's: on the run-29107 close_to_PMT bundles the qualifying per-channel excess
            // has median ~22 PE, p90 ~165, with a saturation tail to ~10k; 100 PE targets that
            // genuine-excess tail (ratio/inflate kept at the faithful C++ 1.3/0.5).  Measured
            // benign for selection (the ladder is KS-led, c2n ceilings 35 loose) -- a faithful
            // PDHD value, not a fit.  ql_light_calib/validate_chain.py.
            chi2_pmt_excess: 100.0,

            // --- High-consistency ladder (SBND method; C++ default OFF => single-branch,
            // byte-identical for any config that omits it) ----------------------------
            // flag_high_consistent gates the PRE-LASSO cull_inconsistent purity cull: when a
            // cluster has a consistent bundle (ladder-flagged OR xtpc), its NON-consistent
            // rivals are dropped before the fit (an xtpc scenario-1 bundle takes absolute
            // priority).  With the ladder OFF the single-branch (highconsist_ks_max/min_ndf)
            // is used.  The multi-branch ladder (TimingTPCBundle::examine_bundle) is KS-LED:
            // a true match has low KS, a wrong one high KS; chi2/ndf only fences the tail.
            // Branches (OR): B1 clean (ks<hc_clean_ks, c2n<hc_clean_c2); B2 good; B3 needs
            // flag_two_boundary (auto-computed when the ladder is on); B4
            // at_x_boundary|close_to_PMT|window_truncated, chi2 RELAXED (missing charge).
            // KS ceilings are kept at the SBND values (KS is the purity lever -- on the
            // run-29107 4-event hand scan ks<0.10 is 88% pure: 57 GT vs 8 junk); the chi2/ndf
            // ceilings are RAISED to PDHD scale.  PDHD's semi-analytical light model is rougher
            // than SBND's, so a clean-KS GT match runs chi2/ndf up to ~32 (vs SBND's ~1-2) while
            // the few low-KS junk matches sit at chi2/ndf >= 40.  The SBND chi2 ceilings (6/4/8)
            // amputated clean-KS PDHD GT matches (e.g. evt991 clu11 ks=0.033 c2n=6.4, evt1007
            // clu87 ks=0.051 c2n=29) -- they fell non-consistent and were culled when their
            // cluster had another ladder-passing bundle, losing 6 clean GT winners.  c2n<35 on
            // the low-KS branches keeps every clean-KS GT (max ~32) and fences the low-KS junk
            // tail (c2n>=40).  High-KS boundary matches (missing charge, ks 0.15-0.40) cannot be
            // separated from junk by KS, so the ladder does NOT flag them consistent -- they win
            // via the LASSO + lasso_flag_weight boundary down-weight as before (cull case-3
            // keep-all leaves them untouched unless a rival is consistent).  GUARDRAIL-tuned,
            // not refit; ql_light_calib/validate_chain.py.  The 12 xtpc-consistent GT winners
            // are protected by the scenario-1/consistent priority, so the cull cannot drop them.
            // Apply the ch_mask (+ saturation mask) to the KS shape metric too, not just
            // chi2/LASSO -- functional parity with SBND (bundle_mask_ks:true; C++ default
            // false). A dead PD reading 0 where the model predicts light should not be charged
            // to a match's KS.  PDHD's 7 dead channels (3/86/87/97/107/116/117) carry ~0
            // predicted light for the current hand-scan matches, so this is near-neutral on the
            // 4-event FOM (correctness/SBND-alignment, not a tuning gain).  Ladder ks ceilings
            // held (no GT match's KS crosses a ceiling under the mask).
            bundle_mask_ks: true,
            highconsist_ladder: true,
            highconsist_ks_max: 0.06,    // single-branch fallback (ladder-off); unused here
            highconsist_min_ndf: 3,      // B1/B2/B3 ndf floor
            hc_clean_ks: 0.06, hc_clean_c2: 35.0,   // SBND ks 0.06; c2n 6->35 (PDHD chi2 scale)
            hc_good_ks:  0.10, hc_good_c2:  35.0,   // SBND ks 0.09->0.10 (ks<0.10 = pure band); c2n 4->35
            hc_tb_ks:    0.10, hc_tb_c2:    35.0,   // two_boundary; c2n 8->35
            hc_miss_ks:  0.08, hc_miss_c2:  60.0,   // boundary/PMT/truncated: tight ks, relaxed chi2 (SBND)
            hc_miss_min_ndf: 5,

            // --- Empty-flash light rescue (SBND method; C++ default OFF => byte-identical) ---
            // After the LASSO, a flash left EMPTY (no bundle above the strength cutoff) adopts
            // its best light-quality candidate from the pre-fit snapshot if
            // metric = ks*(chi2/ndf)^rescue_exponent  (x rescue_boundary_weight per at_x_boundary
            // then close_to_PMT)  <  rescue_metric_max.  One-flash-per-cluster is enforced by
            // reassignment only when the empty flash is a strictly better light match.
            // Unlike SBND (where hand-scan misses were timing-degenerate, ~1/5 light-recoverable),
            // PDHD has many clean GT matches sitting on LASSO-emptied flashes: on the run-29107
            // 4-event scan 16 empty flashes have a GT-accept best candidate at metric < 0.16, with
            // a clear gap before the ks=1.0 cross-side no-flash crossers (metric > 0.8 -- those are
            // xTPC-annotation cases that must NOT be light-rescued).  The lowest non-GT candidate
            // sits at 0.057 (off-scan, neutral) and NO human-rejected match appears at low metric.
            // rescue_metric_max = 0.20 captures the clean GT recoveries, excludes the crossers, and
            // admits only neutral off-scan adoptions; exponent/boundary_weight kept at SBND 0.8/0.8.
            // ql_light_calib/validate_chain.py.
            empty_rescue: true,
            rescue_metric_max: 0.20,
            rescue_exponent: 0.8,
            rescue_boundary_weight: 0.8,

            // --- Cluster-centric rescue (§J; C++ default OFF => byte-identical) ---
            // The empty-flash rescue above only fills flashes left WHOLLY empty. Big
            // charge clusters with an excellent candidate (ks ~0.03-0.15, chi2/ndf ~1-2,
            // pred/meas ~0.8-1.2) but driven to strength 0 by the LASSO L1 sparsity (a
            // rival already explains that flash) stay UNMATCHED even when the flash is
            // non-empty -- the prototype's many-clusters-per-flash case the strength cut
            // + best-per-cluster pipeline drops. This pass adopts the best ACCEPTED
            // candidate for each still-unmatched cluster, attaching it even onto an
            // already-non-empty flash (multiple clusters per flash is physical and
            // GT-endorsed -- the hand-scan labels list several cluster_idents per flash).
            // Acceptance bar (PE-scale-aware, unlike §I's ks-only metric): ks<ks_max AND
            // chi2/ndf<chi2ndf_max AND ratio_lo < pred/meas < ratio_hi. Tuned on the
            // run-29107 4-event hand scan (evts 983/991/999/1007): of the 22 GT matches
            // the LASSO left unmatched, the end-to-end pass recovers 15 to their EXACT
            // hand-scan flash (16 of 22 strands now matched; 1 lands on a non-GT but
            // non-rejected flash), re-introduces 0 scanner-REJECTED matches, and leaves 6
            // (low-quality candidates, or candidates cull_inconsistent removed before the
            // rescue snapshot). Loosening the bar costs purity (a rejected_auto pair
            // reappears) with no recall gain. calib-dump vs labels-evt*.json comparison.
            cluster_rescue: true,
            cluster_rescue_ks_max: 0.20,
            cluster_rescue_chi2ndf_max: 8.0,
            cluster_rescue_ratio_lo: 0.4,
            cluster_rescue_ratio_hi: 2.5,

            // Both drift sides read the SAME global opflash archive (all-PD light
            // reco -> one opflash_pdhd-wct.tar.gz), so each per-side node would dump
            // the full flash list and the Bee op display would double every flash.
            // Key the Bee-op flash gid by the flash's physical side so the two nodes
            // emit ONE gid per physical flash (duplicate collapses; cross-side xTPC
            // matches still resolve).  SBND/PDVD per-TPC flashes are unaffected.
            opflash_phys_gid: true,

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
