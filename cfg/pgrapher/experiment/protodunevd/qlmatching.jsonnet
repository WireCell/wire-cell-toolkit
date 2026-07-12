// PDVD charge-light (Q/L) matching helper.
//
// PDVD counterpart of cfg/pgrapher/experiment/pdhd/qlmatching.jsonnet.  Builds
// the matching graph nodes for ProtoDUNE-VD; the matching-only detector
// constants live here (not in params.jsonnet).  Light-model provenance:
// pdvd/photlib/ (wcp-porting repo) + pdvd/docs/pdvd-photon-model.md; chain
// design + knob rationale: pdvd/docs/pdvd-qlmatching.md.
//
// PDVD differs from PDHD in the ways that matter here:
//   - 40 OpDets of MIXED type: 8 cathode X-ARAPUCAs (x=0, double-sided),
//     8 membrane-wall XAs (y=+-417.6 cm), 8 z-wall PMTs, 16 bottom PMTs
//     (x=-336.5 cm, behind the BOTTOM anode; the top CRP has no PDs).
//     active_opdet_types=[0,1] (both XA and PMT).
//   - ONE all-PD flash (no per-side split; the cathode XAs see both drift
//     volumes) => shared_flash joint LASSO + opdet_all_volumes (all 40 PDs on
//     both drift-side runs) + NO cross-side/xTPC machinery.
//   - vd_surface_flags: flag_close_to_PMT fires for activity near a PD-bearing
//     surface (bottom anode -> bottom PMTs; +-y walls -> that wall's membrane
//     XAs) and the chi2 relaxation applies only to those channels; the cathode
//     end sets the inert flag_at_cathode (treatment designed later).
//   - light model = 'library' (v5 PDFastSimANN sampled on a 10 cm grid); the
//     fitted semi-analytical JSON is the fallback ('semi'), where the membrane
//     XAs must be masked (the port's cosine=|dx|/d is wrong for y-normal PDs).
//
// Channel roles (flash-chain OpDet order == v5 ANN order; from
// pdvd/photlib/pdvd-photlib-chanmap.json, three-way gated against the v5 GDML,
// the official PDVD_PDS_Mapping_v04152025 and the raw-data opdet positions):
//   0-3    membrane XA, top volume    (0,2 at y=+417.6; 1,3 at y=-417.6)
//   4-11   cathode XA (x~0, double-sided)
//   12,13  membrane XA, bottom volume (12 y+; 13 y-, NO WLS -> Ar-blind)
//   14-17  z-wall PMTs, bottom volume (14,15 z=+409; 16,17 z~-100; 14,16,17
//          TPB 0.12; 15 PEN 0.036)
//   18,19  membrane XA, bottom volume (18 y+; 19 y-)
//   20-23  z-wall PMTs (20,22,23 TPB 0.12; 21 PEN 0.036)
//   24-39  bottom PMTs (x=-336.5; PEN 0.036 except 29/39 PEN+Q and 32 uncoated
//          -> Ar-blind; 24/27/28/34 dead in data)
//
// QtoL = 0.070 from the beam-flash gold-pair calibration (see the knob comment);
// trigger_offsets = [bottom, top] per-CRATE light<->charge offsets (the PDVD
// analogue of PDHD's opflash offset_us; per event from the rawwf trigoff tree
// via the archive metadata offset_bot_us/offset_top_us -- the BDE/TDE charge
// windows open up to ~32 us apart).  Any per-run residual constant
// (data/ql_trigger_offset.txt, normally 0) is folded into BOTH values by
// run_clus_evt.sh before they get here.

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';

// require_containment / flash_minPE are function parameters so the driver can
// loosen them for the trigger-offset diagnostic runs (containment is meaningless
// until the time base is calibrated) without editing this file.
function(params, trigger_offset=0 * wc.us, readout_window_ticks=10000,
         light_model='library', require_containment=true, flash_minPE=25,
         trigger_offsets=null) {
    // Per-input [bottom, top] offsets; null => scalar trigger_offset for both
    // (the C++ per-input array, when set, REPLACES the scalar).
    local trigoffs = if trigger_offsets == null
                     then [trigger_offset, trigger_offset]
                     else trigger_offsets,
    local nchan = 40,

    // Static optical dead-channel mask:
    //   dead in data (not in the DAPHNE readout): 24, 27, 28, 34
    //   32 (uncoated PMT): official eff_Ar = eff_Xe = 0 -- masked, though the
    //     data show it responding at its PEN peers' level (ablib_gold.py
    //     Ar-blind closure 1.85x peers); revisit if an official 175 nm
    //     efficiency for it appears.
    // The Ar-blind-only channels 13 (membrane XA, no PTP) and 29/39 (PEN+Q
    // PMTs) are LIVE under the Xe/175 nm model (unmasked 2026-07-11 with the
    // library switch; they respond in data at their peers' level -- see
    // pdvd-questions-dune.md sec 3).
    local ch_mask_base = [24, 27, 28, 32, 34],
    // Semi-analytical mode: additionally mask the LIVE membrane XAs -- the WCT
    // port fixes cosine=|dx|/d (orientation-0), which is wrong/divergent for
    // the y-normal wall XAs (see pdvd-photon-model.md sec 6).  13 (now live
    // under Xe) is a membrane XA too.
    local ch_mask = if light_model == 'semi'
        then ch_mask_base + [0, 1, 2, 3, 12, 13, 18, 19]
        else ch_mask_base,

    // Official per-OpDet detection efficiencies, Xe/175 nm column eff_Xe
    // (PDVD_PDS_Mapping_v04152025, read per channel from
    // pdvd-photlib-chanmap.json; identical to eff_Ar except 13/29/39 which
    // are Ar-blind-only): XA (PTP) 0.03, TPB-coated PMT 0.12, PEN PMT 0.036,
    // uncoated 0.  A newer official map exists (PDVD_PDS_Mapping_v09162025,
    // dunecore) -- re-check these values against it when convenient.  These
    // set the relative PD-type weighting; the absolute scale rides on QtoL
    // (data calibration pending).  Masked channels keep their nominal value
    // (inert).
    local VUVEfficiency = [
        0.03, 0.03, 0.03, 0.03,                    // 0-3   membrane XA (top)
        0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,  // 4-11 cathode XA
        0.03, 0.03,                                // 12,13 membrane XA (13 Xe-only)
        0.12, 0.036, 0.12, 0.12,                   // 14-17 z-wall PMTs (15 PEN)
        0.03, 0.03,                                // 18,19 membrane XA
        0.12, 0.036, 0.12, 0.12,                   // 20-23 z-wall PMTs (21 PEN)
        0.036, 0.036, 0.036, 0.036, 0.036, 0.036,  // 24-29 bottom PMTs (29 PEN+Q Xe-live)
        0.036, 0.036, 0.0, 0.036, 0.036, 0.036,    // 30-35 (32 uncoated, eff_Xe=0 too)
        0.036, 0.036, 0.036, 0.036,                // 36-39 (39 PEN+Q Xe-live)
    ],
    local VISEfficiency = std.makeArray(nchan, function(i) 0.0),

    // PD-surface channel sets for vd_surface_flags (see channel roles above).
    local wall_ylo_channels = [1, 3, 13, 19],    // y=-417.6 membrane XAs
    local wall_yhi_channels = [0, 2, 12, 18],    // y=+417.6 membrane XAs
    local bottom_pmt_channels = std.range(24, 39),

    // Opflash archive reader (the single per-event PDVD opflash file).  `inname`
    // is the full path to opflash_pdvd-wct.tar.gz (caller prefixes the input dir).
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
            // PDVD opflash has no per-frame CAF time offset -> raw flash time.
            correct_flash_time: false,
        },
    }, nin=2, nout=1),

    // Shared QLMatching `data` block -- everything independent of the per-node
    // anode/face wiring.
    local match_data(dv, calib_dump) = {
            detector_volumes: wc.tn(dv),
            data: if std.objectHas(params, 'reality') && params.reality == 'sim' then false else true,
            // Data PE scale, anchored on GROUND-TRUTH pairs: the 80 beam-
            // trigger flashes (position known to ~1 us from the trigoff tc_us,
            // see check_trigger_flash.py) paired with their beam cluster
            // (largest predicted-light bundle on that flash).  Under the
            // Xe/175 nm library + eff_Xe + the Xe-live mask, geometry-fix
            // dumps (tpc_extra_faces union, full +-336.4cm Y): gold-pair
            // median measured/predicted = 0.070, [16,84]% = [0.024,0.176]
            // (80 pairs, ql_light_calib/fit_qtol_gold.py; superseded values
            // 0.082 pre-geometry-fix Xe, 0.11 128nm-era -- fixing the Y
            // truncation raised predicted light, so QtoL moved down).  The
            // library+official-eff model OVER-predicts ~14x as one global
            // normalization (units/recombination/SPE-scale product).
            // DO NOT refit this from auto-selected "good-KS" bundles: with a
            // mis-scaled QtoL the LASSO is amplitude-inert and the selection
            // is dominated by accidentals ~40x brighter than prediction (that
            // route gave QtoL~40, off by ~350x -- see pdvd-qlmatching.md).
            QtoL: 0.070,
            doReflectedLight: false,   // library vis is total photon arrival
            nchan: nchan,
            ch_mask: ch_mask,
            active_opdet_types: [0, 1],   // XAs AND PMTs
            VUVEfficiency: VUVEfficiency,
            VISEfficiency: VISEfficiency,

            // Visibility backend (see header).  In library mode the semimodel
            // is still loaded for the OpDet table.
            light_model: light_model,
            // Xe/175 nm as of 2026-07-11: the runs are Xe-doped by every
            // data test (Ar-blind channels live at peer level, 175 nm wins
            // the gold-pair shape A/B 56/80, mixture scan monotonic --
            // pdvd-questions-dune.md sec 3 + ql_light_calib/ablib_gold.py).
            photon_library_file: 'pdvd/photodet/pdvd-photlib-vis-v5-175nm.json',
            semimodel_file: 'pdvd/photodet/semi-analytical-pdvd.json',

            // --- PDVD single-flash mode (all C++-default-OFF knobs) ---
            // One all-PD flash covers both drift volumes: the two drift-side
            // runs share each physical flash and ONE joint LASSO explains it
            // with clusters from both sides (a per-side fit would let each side
            // independently absorb the full PE -- biased for cathode-crossers).
            shared_flash: true,
            // Keep all 40 PDs on both runs (cathode XAs at x=0 are double-sided;
            // the legacy per-TPC split by cathode x would halve the flash).
            opdet_all_volumes: true,
            // PD-surface endpoint flags: bottom-anode proximity -> bottom PMTs;
            // +-y wall proximity -> that wall's membrane XAs; chi2 relaxation
            // restricted to the triggering surface's channels; the cathode end
            // sets the inert flag_at_cathode.  The z-wall PMTs (~1% of flash PE)
            // are deliberately NOT flagged in round 1 (add their walls later if
            // the hand scan shows a need).
            vd_surface_flags: true,
            pd_wall_cushion: 10 * wc.cm,
            pd_wall_channels_ylo: wall_ylo_channels,
            pd_wall_channels_yhi: wall_yhi_channels,
            anode_pd_channels: [bottom_pmt_channels, []],  // [bottom volume, top volume]

            // Dead-PD self-check: per-event dynamic auto-mask on top of the
            // static ch_mask.  Same-type neighbour pool (XA vs PMT efficiencies
            // differ ~4x, so a cross-type (Y,Z) neighbour is not a meaningful
            // brightness reference).  pe_bright/min_contrast tuned on the 120
            // DIAG=2 dumps by emulating the C++ gate with ch24 dropped from the
            // static mask: bottom-PMT neighbour medians almost never reach the
            // earlier bright=20 (ch24 caught in 55% of events); at 10/3 the
            // known-dead ch24 is caught in 91/120 events with healthy channels
            // masked in <6%.  Two genuinely DIM channels, ch16 and ch33
            // (event-max PE median ~5.6, ~50x below their peers), get
            // per-event masked when quiet -- the safe direction; their status
            // is flagged in pdvd-questions-dune.md (per-run bad channels).
            auto_mask: true,
            auto_mask_same_type: true,
            auto_mask_pe_low: 5,
            auto_mask_pe_bright: 10,
            auto_mask_neighbors: 3,
            auto_mask_min_contrast: 3,
            auto_mask_min_flash: 3,

            // Flash admission: no time cut (flash times run 0..7.5 ms on the
            // light base; the SBND default +-1.5 ms would clip the tail) and
            // the light-chain-study PE knee (~20-30) as the floor.
            flash_minPE: flash_minPE,
            flash_mintime: -1 * wc.s,
            flash_maxtime: 1 * wc.s,
            // Real PDVD readout window (post-resample SP frame length, 10000
            // ticks x 0.5 us = 5 ms), supplied by run_clus_evt.sh from the SP
            // frame.  Without it the C++ default (SBND's 3427) would falsely
            // flag every mid-drift cluster as window-truncated.
            readout_window_ticks: readout_window_ticks,

            // Containment prefilter (drop bundles whose cluster leaves the
            // drift box at the flash's T0).  Meaningless until the light<->
            // charge time base is calibrated -- the driver disables it for the
            // trigger-offset diagnostic runs.
            require_containment: require_containment,

            // Generic SBND/PDHD-proven levers (all C++-default-OFF):
            // sparse LASSO (perf; the joint system is bigger than a per-side one),
            sparse_lasso: true,
            // flag-aware L1 down-weight (a boundary/near-PD bundle's measured
            // light is an underestimate; don't let the LASSO shrink it away),
            lasso_flag_weight: true,
            lasso_boundary_weight: 0.2,
            // mask the KS shape metric on the same channels the chi2/LASSO drop,
            bundle_mask_ks: true,
            // per-bundle chi2 relaxation: with vd_surface_flags the excess
            // widening applies only to the near-surface PD channels; the
            // dead-PD worst-channel drop is detector-agnostic.  chi2_pmt_excess
            // starts at the PDHD ARAPUCA scale (retune after the hand scan).
            chi2_relax: true,
            chi2_pmt_excess: 100.0,
            // Per-channel light-error model, PDHD-proven structure: error on
            // the PREDICTED pe (cures the "predicted light, measured ~0"
            // catastrophe; the measured-based branch gives chi2/ndf ~ 10^3-10^4
            // on the gold pairs) with efficiency-aware low-PE inflation
            // rel = frac + (lowpe_frac-frac)*exp(-pred/lowpe_knee).  Fractions
            // widened vs PDHD's 0.40/1.55/5.5 because the PDVD per-channel
            // measured-PE calibration is still raw (within the cathode-XA
            // group alone the gold-pair per-channel scale spans x13 -- the
            // known "no 1-PE peak" cathode SPE question): on the QtoL=0.11
            // gold pairs this model gives median chi2/ndf 9.7, 71% under the
            // hc ladder ceiling 35.  Retune (and fit measured_pe_scale) after
            // the hand scan provides multi-topology ground truth.
            pe_err_on_pred: true,
            pe_err_floor: 2.0,
            pe_err_frac: 0.60,
            pe_err_lowpe_frac: 2.0,
            pe_err_lowpe_knee: 10.0,
            // KS-led high-consistency ladder (purity cull before the fit).  KS
            // ceilings = SBND/PDHD; chi2/ndf ceilings at the loose PDHD scale --
            // with the uncalibrated PDVD PE scale the chi2 is inflated, so the
            // ladder simply fires less (safe direction).
            highconsist_ladder: true,
            highconsist_ks_max: 0.06,
            highconsist_min_ndf: 3,
            hc_clean_ks: 0.06, hc_clean_c2: 35.0,
            hc_good_ks:  0.10, hc_good_c2:  35.0,
            hc_tb_ks:    0.10, hc_tb_c2:    35.0,
            hc_miss_ks:  0.08, hc_miss_c2:  60.0,
            hc_miss_min_ndf: 5,

            // --- Cathode-crosser (xTPC) machinery, ENABLED 2026-07-11 ---
            // Works under shared_flash: cull_cross_tpc pairs candidate bundles
            // across the two drift sides by flash-TIME coincidence (identical
            // times when the flash is shared) and splits sides by anode_x vs
            // the cathode.  Cuts validated on the 17 hand/finder-confirmed
            // crosser pairs of run 039252 evts 298567/298581/298595 (wcp
            // pdvd/ql_display/docs/ql-cathode-crosser-recipe.md).
            xtpc_flag: true,                // C++ default false
            // PDVD-specific scenario-1 distance ceiling: each volume's active
            // edge sits ~3 cm from x=0 (cathode_x = +-3 cm per side) and the
            // top/bottom trigger crates skew by up to ~30 us (~5 cm), so the
            // true pairs meet at 10-22 cm -- NOT the PDHD/SBND 5 cm.
            xtpc_dmax: 25 * wc.cm,          // C++ default 5 cm
            xtpc_dmax2: 300 * wc.cm,        // C++ default (scenario-2 ceiling)
            xtpc_angle_max: 20,             // C++ default (scenario-2 angles, deg)
            xtpc_hough_radius: 15 * wc.cm,  // C++ default (local-direction radius)
            // Bind each direction-confirmed scenario-1 pair to ONE coincident
            // flash; pinned bundles are exempt from the round-1/2 strength
            // prunes (fit_round{1,2}_shared already honor the exemption).
            xtpc_joint_pin: true,           // C++ default false
            xtpc_pin_angle: 20,             // deg, min(local vhough, global PCA)
            // Widened at-cathode flag window (C++ default -2 cm, PDHD -3 cm):
            // a clean untruncated crosser half whose endpoint stops 2-7 cm
            // short of the cathode (e.g. the evt298567 gid83 hand pick, whose
            // halves carry neither at_x_boundary nor window_truncated with the
            // default window) must acquire at_x_boundary to enter the xtpc
            // candidate pool (admission = at_x_boundary || window_truncated).
            // Side effect (intended): those bundles also join the
            // lasso_flag_weight boundary down-weight group, like every other
            // boundary bundle (~2% of evt298567 bundles flip).
            cathode_ext2: -12 * wc.cm,

            // DELIBERATELY OFF for round 1 (C++ defaults):
            //  - reject_overpred: the gold-pair pred/meas scatter is still
            //    ~x3 either way around QtoL and the per-channel PE scale is
            //    uncalibrated; enable with data-tuned ceilings after the hand
            //    scan (SBND 2.9/4.3, PDHD 3.0/10 for reference).
            //  - measured_pe_scale: beam-gold per-channel fit exists (see
            //    pdvd-qlmatching.md) but is single-topology (bright beam
            //    showers; possible saturation bias) -- refit on hand-scan GT.
            //  - empty_rescue / cluster_rescue: per-run "empty flash" concepts,
            //    not shared_flash-aware (the C++ skips them with a warning).
            //  - cross_side_filter / crossside_skip_vis / opflash_phys_gid:
            //    per-side flash concepts; PDVD has one flash.  (xtpc_* is NOT
            //    in this list any more -- enabled above; the joint LASSO alone
            //    proved insufficient to keep both crosser halves together.)
            //  - robust_endpoint_trim / pmt_nonlinearity: PDHD/SBND-tuned;
            //    revisit after the hand scan.

            drift_speed: params.lar.drift_speed,
            trigger_offset: trigger_offset,
            trigger_offsets: trigoffs,
            calib_dump: calib_dump,
    },

    // JOINT both-sides matcher: the two drift volumes enter ONE node so the
    // shared-flash joint LASSO sees both sides' candidate bundles.  `sides` =
    // per-side anode lists ([[anode0..3],[anode4..7]]); each side's
    // representative (sides[i][0]) is input port i and faces[i] its imaging
    // face (PDVD: both 0). Each CRP anode's two faces share x-bounds but split
    // the Y range in half (adjacent, disjoint) rather than duplicating it like
    // PDHD/SBND -- tpc_extra_faces unions in face 1 so the active-volume box
    // (and the light-prediction inclusion gate) covers the full Y extent
    // instead of truncating to face 0's half.
    matching_joint(sides, dv, faces, calib_dump=''):: g.pnode({
        type: 'QLMatching',
        name: 'matching_joint',
        data: {
            anodes: [wc.tn(side[0]) for side in sides],
            grouping_anodes: [[wc.tn(a) for a in side] for side in sides],
            tpc_faces: faces,
            tpc_extra_faces: [1 for side in sides],
            // Merge the (single, input-0) optical-flash display PC into the
            // joint grouping.
            root_pcs_to_merge: ['opflash'],
        } + match_data(dv, calib_dump),
    }, nin=std.length(sides), nout=1, uses=[dv] + [a for side in sides for a in side]),
}
