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
// QtoL = 0.094 from the beam-flash gold-pair calibration (see the knob comment);
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
         trigger_offsets=null, drift_speed=null, drift_speeds=null,
         cathode_ext1=null, anode_ext1_margin=null, use_saturation_flag=false,
         saturation_mask_fit=true, chi2_sat_inflate=null,
         use_coverage_flag=false, coverage_min=1.0,
         coverage_mask_fit=true, pe_err_nodata=null,
         robust_endpoint_trim=false, robust_endpoint_frac=null,
         robust_endpoint_count=null, robust_endpoint_charge_frac=null,
         robust_endpoint_charge_abs=null, robust_endpoint_gap=null,
         robust_endpoint_gap_charge_frac=null,
         robust_endpoint_walk_to_floor=false,
         robust_endpoint_gap_cathode=false,
         xtpc_cathode_tol=null, xtpc_cathode_qfrac=null,
         mask_wall_xa=false, wall_flags=true,
         flash_sel_cathode=false, flash_sel_minPE=null,
         flash_sel_min_fired=null, flash_sel_fired_pe=null,
         reject_overpred=false, overpred_total_ratio=null,
         overpred_maxch_ratio=null,
         // Per-PD-family PE-error overrides (C++ pe_err_family_* knob,
         // families = [cathode XAs 4-11, PMTs]).  All null (default) => keys
         // omitted => byte-identical global error model above.  A single
         // non-null value activates the family arrays; null members become
         // -1 = "use the global value" in C++.
         pe_err_cath_floor=null, pe_err_cath_frac=null,
         pe_err_cath_lowpe_frac=null, pe_err_cath_lowpe_knee=null,
         pe_err_pmt_floor=null, pe_err_pmt_frac=null,
         pe_err_pmt_lowpe_frac=null, pe_err_pmt_lowpe_knee=null,
         // xtpc / selection quality gates (C++ defaults legacy-inert; doc 19).
         // null / false => keys omitted => byte-identical.
         xtpc_pin_min_strength=null,
         xtpc_sc1_light_gate=false, xtpc_sc1_ks_max=null, xtpc_sc1_c2n_max=null,
         xtpc_cathode_ks_max=null,
         postcull_unflagged=false, postcull_ks_max=null, postcull_c2n_max=null,
         // Sweepable ladder ceilings (doc 19 phase 4).  null => the literal
         // operating values below (compiled JSON unchanged).
         hc_clean_ks=null, hc_clean_c2=null, hc_good_ks=null, hc_good_c2=null,
         hc_tb_ks=null, hc_tb_c2=null, hc_miss_ks=null, hc_miss_c2=null,
         hc_miss_min_ndf=null,
         // Sweepable LASSO regularization (doc 19 phase 4).  null => key
         // omitted => the C++ defaults (lambda 0.1, delta_charge 0.01,
         // delta_light 0.025, delta_shape 0.01, bkg_weight 0.5,
         // strength_cutoff 0.05, lasso_boundary_weight 0.2).
         lasso_lambda=null, delta_charge=null, delta_light=null,
         delta_shape=null, bkg_weight=null, strength_cutoff=null,
         lasso_boundary_weight=null) {
    // Per-input [bottom, top] offsets; null => scalar trigger_offset for both
    // (the C++ per-input array, when set, REPLACES the scalar).
    local trigoffs = if trigger_offsets == null
                     then [trigger_offset, trigger_offset]
                     else trigger_offsets,
    local nchan = 40,

    // Static optical dead-channel mask (Ar/128 nm model as of 2026-07-13):
    //   dead in data (not in the DAPHNE readout): 24, 27, 28, 34
    //   Ar-blind (official eff_Ar = 0; no/quenched WLS at 128 nm):
    //     13 (membrane XA, no PTP), 29/39 (PEN+Q PMTs), 32 (uncoated PMT).
    //   persistently faulty/dim -- ~8-50x below their same-type peers in EVERY
    //     run examined (039252/039253/039349), so a stable hardware state, not a
    //     per-event fluctuation; previously only caught by the per-event
    //     auto_mask below.  Masked here so their biased light never enters the
    //     LASSO/chi2/KS: 2 (top membrane XA), 16/17 (z-wall PMTs), 33 (bottom
    //     PMT).  Data justification: pdvd/docs/qlmatch/pdvd-pd-functionality-run39252.md
    //     (wcp-porting repo).  ch0 is left LIVE (dimmer than its bright top-wall
    //     peers but meas/pred ~0.83, i.e. functional); the other live PMTs are
    //     also kept -- they respond in proportion to the library and carry the
    //     only real -x light away from the walls.
    // These are the channels that do NOT contribute usable Ar 128 nm light.  (For
    // Xe-doped running unmask 13/29/39 and switch to the 175 nm library --
    // that Xe verdict was adopted 2026-07-11 in 0adb15fa and REVERTED here;
    // see pdvd-questions-dune.md sec 3.)
    local ch_mask_base = [2, 13, 16, 17, 24, 27, 28, 29, 32, 33, 34, 39],
    // Semi-analytical mode: additionally mask the LIVE membrane XAs -- the WCT
    // port fixes cosine=|dx|/d (orientation-0), which is wrong/divergent for
    // the y-normal wall XAs (see pdvd-photon-model.md sec 6).  13 is already
    // in the base mask.
    //
    // mask_wall_xa (library mode): mask the SAME live membrane/wall XAs on
    // reliability grounds -- the evt298567 hand scan shows they are BIMODAL
    // (zero response in 9/15 matches with >=20 PE predicted, x2.5 when they
    // do respond; no scale factor fixes it) and dropping them improves the
    // match/reject separation on every metric (AUC(ks) 0.751->0.795, good-rung
    // matches 43%->59%; pdvd/docs/qlmatch/17_pdvd-pd-family-reliability-
    // evt298567.md).  mask_ks (bundle_mask_ks below) is already true, so the
    // masked channels leave the chi2, LASSO AND the KS shape test together.
    // Default false => byte-identical pre-study config.
    local wall_xa_live = [0, 1, 3, 12, 18, 19],   // 2, 13 already in ch_mask_base
    local ch_mask = if light_model == 'semi'
        then ch_mask_base + [0, 1, 2, 3, 12, 18, 19]
        else if mask_wall_xa then ch_mask_base + wall_xa_live
        else ch_mask_base,
    // Cathode X-ARAPUCAs: the one PD family with full-stream readout (100%
    // flash coverage), proportional response, and per-channel medians within
    // +-30% of each other -- the ruler family (doc 17).  Used by the
    // flash_sel_* admission and the overpred_channels scope below.
    local cathode_channels = std.range(4, 11),

    // Per-PD-family PE-error override plumbing (scan-tuning doc 19).  Hoisted
    // to function scope: computed field names below cannot see object-locals.
    // Family order [cathode XAs (4-11), PMTs (14-17, 20-39)]; nn(null) = -1 =
    // "use the global value" in C++.
    local pe_err_pmt_channels = std.range(14, 17) + std.range(20, 39),
    local pe_err_nn(v) = if v == null then -1 else v,
    local pe_err_fam_on =
        pe_err_cath_floor != null || pe_err_cath_frac != null ||
        pe_err_cath_lowpe_frac != null || pe_err_cath_lowpe_knee != null ||
        pe_err_pmt_floor != null || pe_err_pmt_frac != null ||
        pe_err_pmt_lowpe_frac != null || pe_err_pmt_lowpe_knee != null,

    // Official per-OpDet detection efficiencies, Ar/128 nm column eff_Ar
    // (PDVD_PDS_Mapping_v04152025, read per channel from
    // pdvd-photlib-chanmap.json): XA (PTP) 0.03, TPB-coated PMT 0.12, PEN PMT
    // 0.036, Ar-blind 0 (13/29/32/39).  These set the relative PD-type
    // weighting.  Masked channels keep their nominal value (inert).
    //
    // Per-type DATA scale factors (owner-adopted 2026-07-14, Ar/128 nm kept):
    // median Sum(meas)/Sum(pred) per PD type over the 192 strict geometric
    // cathode-crosser anchors of the saturation-fixed 120-event `_satrep`
    // reprocess (veto OFF + flag + twoside repair; railed channels excluded
    // from BOTH sums -- the veto had biased the cathode ratio x5 low):
    //   cathode XA  10.116 [7.62, 24.4]  n=191
    //   membrane XA  1.655 [0.27,  8.1]  n=113
    //   PMTs         0.352 [0.009, 2.0]  n=163 (TPB and PEN share the PMT
    //     factor; the official 0.12/0.036 relative weighting is kept)
    // fit: pdvd/ql_light_calib/fit_qtol_crossers.py, doc
    // pdvd/docs/qlmatch/pdvd-qtol-recalibration.md.  The factors multiply the
    // official effs so the per-type median meas/pred is 1 at UNCHANGED
    // QtoL = 0.094; they are effective weights, not physical efficiencies
    // (the ~x4 far-channel rise of the 128 nm visibility shape remains).
    local eff_scale_cathode = 10.116,
    local eff_scale_membrane = 1.655,
    local eff_scale_pmt = 0.352,
    local eff_xa_cath = 0.03 * eff_scale_cathode,
    local eff_xa_memb = 0.03 * eff_scale_membrane,
    local eff_pmt_tpb = 0.12 * eff_scale_pmt,
    local eff_pmt_pen = 0.036 * eff_scale_pmt,
    local VUVEfficiency = [
        eff_xa_memb, eff_xa_memb, eff_xa_memb, eff_xa_memb,  // 0-3   membrane XA (top)
        eff_xa_cath, eff_xa_cath, eff_xa_cath, eff_xa_cath,  // 4-7   cathode XA
        eff_xa_cath, eff_xa_cath, eff_xa_cath, eff_xa_cath,  // 8-11  cathode XA
        eff_xa_memb, 0.0,                                    // 12,13 membrane XA (13 no-WLS)
        eff_pmt_tpb, eff_pmt_pen, eff_pmt_tpb, eff_pmt_tpb,  // 14-17 z-wall PMTs (15 PEN)
        eff_xa_memb, eff_xa_memb,                            // 18,19 membrane XA
        eff_pmt_tpb, eff_pmt_pen, eff_pmt_tpb, eff_pmt_tpb,  // 20-23 z-wall PMTs (21 PEN)
        eff_pmt_pen, eff_pmt_pen, eff_pmt_pen,               // 24-26 bottom PMTs
        eff_pmt_pen, eff_pmt_pen, 0.0,                       // 27-29 (29 PEN+Q)
        eff_pmt_pen, eff_pmt_pen, 0.0,                       // 30-32 (32 uncoated)
        eff_pmt_pen, eff_pmt_pen, eff_pmt_pen,               // 33-35
        eff_pmt_pen, eff_pmt_pen, eff_pmt_pen, 0.0,          // 36-39 (39 PEN+Q)
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
            // (largest predicted-light bundle on that flash).
            // Ar/128 nm value ~= 0.094 (2026-07-13).  ESTIMATE, not a fresh
            // gold-pair refit: the jjo beam-flash timestamp table
            // (data/jjo_triglight_offsets.txt) needed by fit_qtol_gold.py /
            // ablib_gold.py has been removed and the raw files changed schema
            // (triglight, no charge_bde_us), so the 80-pair anchor cannot be
            // rebuilt without reconstructing t_expect.  0.094 = the recorded
            // pre-Y-truncation gold-pair Ar value 0.11 scaled by the Xe
            // truncation factor 0.070/0.082 (0.11*0.070/0.082=0.0939); the
            // library pred-ratio S175/S128 on the current-geometry beam-like
            // bundles corroborates (brightest-flash-per-event median 1.29 ->
            // 0.070*1.29=0.090, trending into 0.094 as the population is made
            // more beam/cathode-like).  A full 128 nm gold-pair refit is a
            // flagged follow-up once the trigger table is restored.  The
            // library+official-eff model OVER-predicts ~10x as one global
            // normalization (units/recombination/SPE-scale product).
            // DO NOT refit this from auto-selected "good-KS" bundles: with a
            // mis-scaled QtoL the LASSO is amplitude-inert and the selection
            // is dominated by accidentals ~40x brighter than prediction (that
            // route gave QtoL~40, off by ~350x -- see pdvd-qlmatching.md).
            // 2026-07-14: the crosser-anchor recalibration on the
            // saturation-fixed dumps keeps this value; the per-type residuals
            // are absorbed by the VUVEfficiency scale factors above.
            QtoL: 0.094,
            doReflectedLight: false,   // library vis is total photon arrival
            nchan: nchan,
            ch_mask: ch_mask,
            active_opdet_types: [0, 1],   // XAs AND PMTs
            VUVEfficiency: VUVEfficiency,
            VISEfficiency: VISEfficiency,

            // Visibility backend (see header).  In library mode the semimodel
            // is still loaded for the OpDet table.
            light_model: light_model,
            // Ar/128 nm as of 2026-07-13: reverted from the Xe/175 nm verdict
            // (0adb15fa, 2026-07-11) -- the runs are treated as pure-Argon
            // scintillation for this pass (see pdvd-questions-dune.md sec 3).
            photon_library_file: 'pdvd/photodet/pdvd-photlib-vis-v5-128nm.json',
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
            // wall_flags=false empties BOTH wall channel lists, which is the C++
            // "wall inactive" idiom (QLMatching.h pd_wall_channels doc): the +-y
            // wall proximity then sets no flag_close_to_PMT and adds no relax
            // channels.  Adopted with mask_wall_xa: the wall XAs the relaxation
            // was protecting are masked out of the fit anyway, and the wall flag
            // only exempted junk from the overpred prefilter / widened the chi2
            // on wall-PMT-free bundles for nothing.  The bottom-anode proximity
            // flag (anode_pd_channels) is NOT touched.  Default true =>
            // byte-identical pre-study config.
            pd_wall_channels_ylo: if wall_flags then wall_ylo_channels else [],
            pd_wall_channels_yhi: if wall_flags then wall_yhi_channels else [],
            anode_pd_channels: [bottom_pmt_channels, []],  // [bottom volume, top volume]

            // Dead-PD self-check: per-event dynamic auto-mask on top of the
            // static ch_mask.  Same-type neighbour pool (XA vs PMT efficiencies
            // differ ~4x, so a cross-type (Y,Z) neighbour is not a meaningful
            // brightness reference).  pe_bright/min_contrast tuned on the 120
            // DIAG=2 dumps by emulating the C++ gate with ch24 dropped from the
            // static mask: bottom-PMT neighbour medians almost never reach the
            // earlier bright=20 (ch24 caught in 55% of events); at 10/3 the
            // known-dead ch24 is caught in 91/120 events with healthy channels
            // masked in <6%.  The persistently DIM channels (ch2/16/17/33) are
            // now in the static ch_mask_base above (confirmed stable across
            // 039252/253/349, pdvd-pd-functionality-run39252.md); auto_mask
            // remains the per-event safety net for anything transient/new.
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
            // Cathode-scoped flash admission ON TOP of flash_minPE: a flash must
            // show real light on the cathode XAs (sum >= flash_sel_minPE PE and
            // >= flash_sel_min_fired channels at >= flash_sel_fired_pe PE) to
            // enter matching -- the all-PD floor alone admits flashes whose PE
            // lives on the bimodal wall XAs / self-trigger-silent PMTs.
            // Operating point 5 PE / 2 fired @ 1 PE keeps 289/297 (97.3%) of the
            // keep-round confirmed crosser+boundary candle flashes (the 8 lost
            // are dim wall-hugging tracks, <=160 PE total, two with literally
            // zero cathode signal) and cuts ~7% of the previously admitted
            // flashes (pdvd/docs/qlmatch/18_*.md).  C++ defaults: empty channel
            // list / 0 / 0 / 1.0 => keys omitted when off => byte-identical
            // pre-study config.
            [if flash_sel_cathode then 'flash_sel_channels']: cathode_channels,
            [if flash_sel_cathode && flash_sel_minPE != null then 'flash_sel_minPE']:
                flash_sel_minPE,
            [if flash_sel_cathode && flash_sel_min_fired != null then 'flash_sel_min_fired']:
                flash_sel_min_fired,
            [if flash_sel_cathode && flash_sel_fired_pe != null then 'flash_sel_fired_pe']:
                flash_sel_fired_pe,
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

            // Light-pattern over-prediction prefilter, scoped to the CATHODE XAs
            // (overpred_channels): a candidate bundle predicting far more cathode
            // light than the flash measured cannot be the match.  The scope
            // matters -- judged on all PDs, a self-trigger-silent PMT / bimodal
            // wall XA reads 0 and fakes over-prediction on good bundles; the
            // cathode XAs are always read out.  Ceilings from the 18-event
            // keep-round scan (pdvd/docs/qlmatch/18_*.md): R_total <= 15 /
            // R_max <= 50 keeps every evt298567 full-scan match (max 1.43/4.97)
            // and culls ~16% of all non-exempt candidate bundles; the only
            // confirmed candles above them are 3 wrong-amplitude picks
            // (meas ~ 0, pred in the thousands).  Boundary/truncated bundles
            // stay exempt (C++).  C++ defaults false / empty / 1e9 => keys
            // omitted when off => byte-identical pre-study config.
            [if reject_overpred then 'reject_overpred']: true,
            [if reject_overpred then 'overpred_channels']: cathode_channels,
            [if reject_overpred && overpred_total_ratio != null then 'overpred_total_ratio']:
                overpred_total_ratio,
            [if reject_overpred && overpred_maxch_ratio != null then 'overpred_maxch_ratio']:
                overpred_maxch_ratio,

            // Generic SBND/PDHD-proven levers (all C++-default-OFF):
            // sparse LASSO (perf; the joint system is bigger than a per-side one),
            sparse_lasso: true,
            // flag-aware L1 down-weight (a boundary/near-PD bundle's measured
            // light is an underestimate; don't let the LASSO shrink it away),
            lasso_flag_weight: true,
            lasso_boundary_weight: 0.2,
            // mask the KS shape metric on the same channels the chi2/LASSO drop,
            bundle_mask_ks: true,
            // Per-flash DAPHNE-rail channel masking (needs a light archive made
            // with OpHitFinder flag_saturation; see pdvd/docs/qlmatch/
            // pdvd-saturation-recovery.md).  C++ default false.  Key omitted
            // when off => byte-identical pre-fix config.
            [if use_saturation_flag then 'use_saturation_flag']: true,
            // ...and whether a railed channel is DROPPED from the chi2/KS or
            // kept there at its clipped PE.  Keeping it is right: the clipped
            // value is a LOWER BOUND on the true light (11_pdvd-saturation-
            // recovery.md §2.1 `clip` = a tight deterministic underestimate;
            // and the production chain repairs it via OpDecon
            // saturation_repair), so it carries real information -- dropping
            // it also makes a WRONG prediction free on exactly the bright
            // cathode channels that anchor Q/L matching.  The LASSO rows stay
            // zeroed regardless (a lower bound must not pull the fitted charge
            // down).  C++ default true = drop; key omitted then =>
            // byte-identical.
            [if !saturation_mask_fit then 'saturation_mask_fit']: false,
            // Extra chi2 width on a railed channel: denom += (pe*inflate)^2,
            // the same form as chi2_pmt_inflate but gated on the rail flag
            // alone.  C++ default 0 = none; key omitted when null.
            [if chi2_sat_inflate != null then 'chi2_sat_inflate']: chi2_sat_inflate,
            // Per-flash readout-coverage tracking: a membrane-XA / PMT channel
            // is a 16.4-us self-trigger snippet stream (duty ~5-30%); with no
            // snippet over a flash's window it reads measured = 0 (needs a
            // light archive made with OpHitFinder emit_coverage; see
            // pdvd/docs/qlmatch/14_pdvd-lightpattern-sp-investigation.md).
            // C++ defaults false/1.0.  Keys omitted when off => byte-identical
            // pre-fix config.
            [if use_coverage_flag then 'use_coverage_flag']: true,
            [if use_coverage_flag && coverage_min != 1.0 then 'coverage_min']: coverage_min,
            // ...and whether an uncovered channel is DROPPED from the fit or
            // kept at that 0.  Keeping it is right: DAPHNE has no dead time
            // (min inter-snippet gap 0.336 us) and the self-trigger threshold
            // is ~1 PE, so the silence measures "this channel saw < ~1 PE"
            // rather than nothing at all -- dropping it lets a prediction
            // over-shoot there for free (doc 14 section 12).  The channel is
            // still labelled "nodata" in the calib dump / ql_scan either way.
            // C++ default true = drop; key omitted then => byte-identical.
            [if !coverage_mask_fit then 'coverage_mask_fit']: false,
            // PE_err floor (PE) for kept no-data channels: their measurement is
            // the upper limit "< self-trigger threshold", not "0 +- pe_err_floor".
            // PDVD does NOT need this -- pe_err_floor below is already 2.0 with
            // pe_err_knee at the C++ default 1.0, so a measured 0 carries +-2 PE,
            // about twice the ~1 PE threshold.  Present for detectors whose
            // pe_err_floor is tighter (the C++ default is 0.3, at which a 5 PE
            // over-prediction on a no-data channel would cost chi2 ~278).
            // C++ default -1 = disabled; key omitted when null.
            [if pe_err_nodata != null then 'pe_err_nodata']: pe_err_nodata,
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
            // Per-PD-family override of the model above (scan-tuning doc 19).
            // C++ default: empty arrays.  Keys omitted when all family params
            // are null => byte-identical pre-knob config.  Family order:
            // [cathode XAs (4-11), PMTs (14-17, 20-39)]; -1 => global value.
            [if pe_err_fam_on then 'pe_err_family_channels']:
                [cathode_channels, pe_err_pmt_channels],
            [if pe_err_fam_on then 'pe_err_family_floor']:
                [pe_err_nn(pe_err_cath_floor), pe_err_nn(pe_err_pmt_floor)],
            [if pe_err_fam_on then 'pe_err_family_frac']:
                [pe_err_nn(pe_err_cath_frac), pe_err_nn(pe_err_pmt_frac)],
            [if pe_err_fam_on then 'pe_err_family_lowpe_frac']:
                [pe_err_nn(pe_err_cath_lowpe_frac), pe_err_nn(pe_err_pmt_lowpe_frac)],
            [if pe_err_fam_on then 'pe_err_family_lowpe_knee']:
                [pe_err_nn(pe_err_cath_lowpe_knee), pe_err_nn(pe_err_pmt_lowpe_knee)],
            // KS-led high-consistency ladder (purity cull before the fit).  KS
            // ceilings = SBND/PDHD; chi2/ndf ceilings at the loose PDHD scale --
            // with the uncalibrated PDVD PE scale the chi2 is inflated, so the
            // ladder simply fires less (safe direction).
            highconsist_ladder: true,
            highconsist_ks_max: 0.06,
            highconsist_min_ndf: 3,
            // Ceilings sweepable via the null-default function args (doc 19
            // phase 4); null (default) compiles to exactly these literals.
            hc_clean_ks: if hc_clean_ks == null then 0.06 else hc_clean_ks,
            hc_clean_c2: if hc_clean_c2 == null then 35.0 else hc_clean_c2,
            hc_good_ks:  if hc_good_ks == null then 0.10 else hc_good_ks,
            hc_good_c2:  if hc_good_c2 == null then 35.0 else hc_good_c2,
            hc_tb_ks:    if hc_tb_ks == null then 0.10 else hc_tb_ks,
            hc_tb_c2:    if hc_tb_c2 == null then 35.0 else hc_tb_c2,
            hc_miss_ks:  if hc_miss_ks == null then 0.08 else hc_miss_ks,
            hc_miss_c2:  if hc_miss_c2 == null then 60.0 else hc_miss_c2,
            hc_miss_min_ndf: if hc_miss_min_ndf == null then 5 else hc_miss_min_ndf,
            // Sweepable LASSO regularization (doc 19 phase 4); keys omitted
            // when null => C++ defaults => byte-identical pre-knob config.
            [if lasso_lambda != null then 'lasso_lambda']: lasso_lambda,
            [if delta_charge != null then 'delta_charge']: delta_charge,
            [if delta_light != null then 'delta_light']: delta_light,
            [if delta_shape != null then 'delta_shape']: delta_shape,
            [if bkg_weight != null then 'bkg_weight']: bkg_weight,
            [if strength_cutoff != null then 'strength_cutoff']: strength_cutoff,
            [if lasso_boundary_weight != null then 'lasso_boundary_weight']:
                lasso_boundary_weight,

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
            // At-cathode flag window inner edge (C++ default -2 cm, PDHD -3 cm).
            // Previously widened to -12 cm so that clean untruncated crosser
            // halves stopping 2-7 cm short of the cathode (e.g. the evt298567
            // gid83 hand pick) still acquire at_x_boundary and enter the xtpc
            // candidate pool (admission = at_x_boundary || window_truncated).
            // That window over-flagged ordinary non-crosser cluster ends near
            // the cathode (evt298567 bot fl62 c3/c75 shown atCATH/xbound though
            // several cm short), so it is narrowed back to the C++ default -2 cm.
            // Trade-off: a genuine crosser half whose end stops >2 cm short of
            // the cathode no longer acquires at_x_boundary via this window and
            // must rely on window_truncated for xtpc admission -- revisit if the
            // gid83-style short-half pool needs recovering.
            cathode_ext2: -2 * wc.cm,
            // xtpc CATHODE RESCUE (C++ default 0 = OFF; keys omitted when null =>
            // byte-identical pre-knob config). A crosser pair whose halves meet a
            // few cm off the nominal cathode fails BOTH xtpc admission gates at
            // once: the overshooting half fails containment, the short half misses
            // the at_cathode window above (the -12 cm cathode_ext2 experiment tried
            // to fix the short half but over-flagged non-crossers BECAUSE it went
            // through at_x_boundary; the rescue widens a dedicated xtpc-candidacy
            // flag instead, leaving at_x_boundary and its ladder/LASSO consumers
            // legacy). tol = length past the containment gate / below the
            // at_cathode window tolerated FOR XTPC CANDIDACY ONLY; qfrac = charge
            // fraction discardable as overclustered junk when measuring the
            // overshoot. A provisionally kept uncontained half is purged before the
            // fit unless cull_cross_tpc confirms it (scenario 1, one half
            // contained). PDVD 039252 evt298567 top-97+bot-139 vs flash 96
            // (pdvd/docs/qlmatch/16_pdvd-clus97-crosser-evt298567.md §10).
            [if xtpc_cathode_tol != null then 'xtpc_cathode_tol']: xtpc_cathode_tol,
            [if xtpc_cathode_qfrac != null then 'xtpc_cathode_qfrac']: xtpc_cathode_qfrac,

            // xtpc / selection quality gates (scan-tuning doc 19).  All C++
            // defaults are legacy-inert; every key omitted when its arg is
            // null/false => byte-identical pre-knob config.
            // Pin strength floor: pinned bundle loses the strength-cutoff
            // exemption when its LASSO solution <= floor (C++ default 0 = off).
            [if xtpc_pin_min_strength != null then 'xtpc_pin_min_strength']:
                xtpc_pin_min_strength,
            // Scenario-1/xtpc-consistent light gate (C++ default false;
            // ks/c2n ceilings default 0.3/50 in C++).
            [if xtpc_sc1_light_gate then 'xtpc_sc1_light_gate']: true,
            [if xtpc_sc1_light_gate && xtpc_sc1_ks_max != null then 'xtpc_sc1_ks_max']:
                xtpc_sc1_ks_max,
            [if xtpc_sc1_light_gate && xtpc_sc1_c2n_max != null then 'xtpc_sc1_c2n_max']:
                xtpc_sc1_c2n_max,
            // Cathode-rescue ks ceiling (C++ default 0 = off).
            [if xtpc_cathode_ks_max != null then 'xtpc_cathode_ks_max']:
                xtpc_cathode_ks_max,
            // Post-fit cull of unflagged low-quality selections (C++ default
            // false; ks/c2n ceilings default 0.30/20 in C++).
            [if postcull_unflagged then 'postcull_unflagged']: true,
            [if postcull_unflagged && postcull_ks_max != null then 'postcull_ks_max']:
                postcull_ks_max,
            [if postcull_unflagged && postcull_c2n_max != null then 'postcull_c2n_max']:
                postcull_c2n_max,

            // Cathode-side containment tolerance (how far past the cathode a
            // cluster's drift end may reach and still count as contained;
            // cut is last_u < u_cathode + cathode_ext1).  C++ default +1.2 cm
            // (QLMatching.h:238).  Key omitted when the arg is null =>
            // byte-identical pre-knob config; the driver passes a value only
            // for the anode-pull / cushion study (run 039252 evt298567).
            [if cathode_ext1 != null then 'cathode_ext1']: cathode_ext1,

            // Anode-side containment slack: extra tolerance BELOW anode_ext1,
            // subtracted to form the anode floor for BOTH the containment gate
            // (first_u > anode_ext1 - margin) and the anode flag-window inner edge
            // (at_x_boundary / close_to_PMT).  The C++ ties the two deliberately --
            // a cluster admitted by the slack must still acquire the boundary flags.
            // C++ default 1.0 cm (QLMatching.h, prototype ProtectOverClustering.cxx
            // :296 hard-codes it).  Key omitted when the arg is null =>
            // byte-identical pre-knob config; the driver passes 2.0 cm for PDVD
            // production (floor -2 -> -4 cm), so a full-gap cathode-crosser whose
            // anode end pokes a few cm past the anode under the correct flash stays
            // a candidate instead of being prefiltered onto a neighbouring flash
            // (run 039252 evt298567 top:22 missed by 0.501 cm; the 2 cm anode pull
            // is what pushes it out -- see pdvd-cathode-containment-flash-demotion).
            [if anode_ext1_margin != null then 'anode_ext1_margin']: anode_ext1_margin,

            // Robust-endpoint trim: snap a drift endpoint back inside the detector
            // edge when the outer material is overclustered junk rather than a real
            // track end, before the containment gate reads the endpoints.  The
            // master switch is robust_endpoint_trim; the rest are its judges, and
            // NONE of them do anything while it is false (QLMatching.cxx:3699 gates
            // the whole block).  C++ default false => key suppressed when off =>
            // byte-identical pre-knob config.
            //
            // PDVD leaves this OFF pending an A/B + a census of affected clusters.
            // The driver turns it on for the evt298567 clus 34 demo: 8 points
            // carrying 0.08% of that cluster's charge sit 28.4 cm off its anode end
            // and stretch its drift extent to 371.7 cm -- past the 336.9 cm drift
            // depth -- so no flash T0 can contain it and every candidate bundle is
            // dropped before the fit.  Only the GAP judge can trim them: they are
            // detached, not diffuse, so the density judge (charge_abs) protects them
            // as it would a genuine track tip.
            // See pdvd/docs/qlmatch/15_pdvd-clus34-unmatched-evt298567.md.
            // PDHD's proven operating point is pdhd/qlmatching.jsonnet:306-332.
            [if robust_endpoint_trim then 'robust_endpoint_trim']: true,
            // Point-count judge: trim when outer points < max(frac*npoints, count).
            // C++ defaults 0.05 / 0.0.
            [if robust_endpoint_frac != null then 'robust_endpoint_frac']: robust_endpoint_frac,
            [if robust_endpoint_count != null then 'robust_endpoint_count']: robust_endpoint_count,
            // Charge-fraction judge, OR'd with the point-count one.  C++ default 0 =>
            // disabled (the point-count judge stands alone).
            [if robust_endpoint_charge_frac != null then 'robust_endpoint_charge_frac']:
                robust_endpoint_charge_frac,
            // Per-point charge-DENSITY ceiling guarding the charge-fraction path: only
            // diffuse material is trimmed, never a charge-dense real track tip.
            // C++ default 0 => disabled.
            [if robust_endpoint_charge_abs != null then 'robust_endpoint_charge_abs']:
                robust_endpoint_charge_abs,
            // Gap (detachment) judge: trim outer material split from the body by a gap
            // wider than robust_endpoint_gap that carries at most
            // robust_endpoint_gap_charge_frac of the cluster charge.  Keyed on
            // detachment, not density, so it reaches junk the density judge protects.
            // Applies to the ANODE end only unless robust_endpoint_gap_cathode is set.
            // C++ defaults 0 / 0 => disabled.
            [if robust_endpoint_gap != null then 'robust_endpoint_gap']: robust_endpoint_gap,
            [if robust_endpoint_gap_charge_frac != null then 'robust_endpoint_gap_charge_frac']:
                robust_endpoint_gap_charge_frac,
            // Mirror that gap judge at the CATHODE end, reusing the same two thresholds.
            // The knobs above were added for anode-end cases, but nothing in their
            // justification is direction-specific, and the upstream cause is symmetric:
            // clustering_isolated merges any "small" cluster into the nearest big one
            // within a hardcoded 80 cm (no angle/direction/gap test), at whichever end
            // it lies.  cf. PDHD evt 1007 uid-126 (anode) vs PDVD evt298567 apa-4
            // ident 97 (cathode): same pathology, only the first had a judge.
            // C++ default false => key omitted => byte-identical pre-knob config, and
            // PDHD (which sets the two keys above) is unaffected until it opts in.
            // See pdvd/docs/qlmatch/16_pdvd-clus97-crosser-evt298567.md.
            [if robust_endpoint_gap_cathode then 'robust_endpoint_gap_cathode']: true,
            // Where the anode-end walk stops calling material "outside".  Default
            // anode_in (= anode_ext1, -2 cm), which is anode_ext1_margin cm ABOVE the
            // floor the containment gate uses (first_u > anode_ext1 - margin, -4 cm
            // for PDVD) -- so a cluster whose body starts in that 2 cm dead band has
            // the walk march past its detached junk INTO the body, and the judges then
            // correctly refuse to trim real track charge.  Whether the rescue works
            // ends up depending on cluster SIZE, which is not physics: evt298567 apa-0
            // ident 34 (2649 pts) is rescued, apa-4 ident 1 (1375 pts) is not, on the
            // same 20 swallowed points.  When true the walk breaks at that same floor.
            // C++ default false => key omitted => byte-identical pre-knob config.
            // See pdvd/docs/qlmatch/15_pdvd-clus34-unmatched-evt298567.md.
            [if robust_endpoint_walk_to_floor then 'robust_endpoint_walk_to_floor']: true,

            // DELIBERATELY OFF for round 1 (C++ defaults):
            //  (reject_overpred is NOT in this list any more -- the hand scan
            //  provided the data-tuned cathode-scoped ceilings above.)
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

            // Scalar override (sec 8.12 recalibration): null => the common
            // params.lar.drift_speed => byte-identical pre-knob config.  Also
            // the value the calib dump exports as d["drift_speed"].
            drift_speed: if drift_speed == null then params.lar.drift_speed
                         else drift_speed,
            // Per-input [bottom, top] drift speeds (the two PDVD volumes can
            // carry different calibrated values, sec 8.12 of
            // pdvd-anode-time-consistency.md). C++ default [] => the scalar
            // drift_speed for every input. Key omitted when null =>
            // byte-identical pre-knob compiled config.
            [if drift_speeds != null then 'drift_speeds']: drift_speeds,
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
