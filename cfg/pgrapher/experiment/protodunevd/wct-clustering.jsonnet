local g = import "pgraph.jsonnet";
local f = import "pgrapher/common/funcs.jsonnet";
local wc = import "wirecell.jsonnet";

local tools_maker = import 'pgrapher/common/tools.jsonnet';

local reality = 'data';
local data_params = import 'pgrapher/experiment/protodunevd/params.jsonnet';
local simu_params = import 'pgrapher/experiment/protodunevd/simparams.jsonnet';
local params = if reality == 'data' then data_params else simu_params;
local tools_all = tools_maker(params);

// Top-level function: parameters overridable via --tla-str / --tla-code
function(
    input = ".",
    // Indices into tools_all.anodes to process; default = all
    anode_indices = std.range(0, std.length(tools_all.anodes) - 1),
    // Directory for output mabc-*.zip and bee data files ('' means current directory)
    output_dir = '',
    run = 1,
    subrun = 1,
    event = 1,
    // Stepped-sampler fallback: blobs the stepped grid leaves point-less get
    // one point at the blob center.  Default off -> bit-identical output.
    stepped_center_fallback = false,
    // Wires file for this job's anodes (doc qlmatch/29: cluster July's v6 imaging with v6 geometry).  '' => the
    // params.jsonnet file (v7-uvwfit) => compiled config byte-identical.  Same idiom as wct-nf-sp-dnnroi.jsonnet.
    // run_clus_evt.sh: PDVD_CLUS_WIRES.
    wires_file = '',
    // doc pdvd/28 round 2: the SBND Q/L-tail fast flavors, off by default.
    // ql_po_fast => ProtectOverclustering busy_num_threshold=200 (per-apa
    // stage); ql_dg_fast => Deghost on the 'ctpc_fast' graph flavor (per-apa
    // and per-group stages).  false => keys omitted in clus.jsonnet =>
    // byte-identical compiled config (proof: doc pdvd/28 sec 12).
    ql_po_fast = false,
    ql_dg_fast = false,
    // doc 31 round 3: resolve induction charge by channel IDENT on a wrapped
    // strip's continuation (PDVD 1568 U/V wires).  false => key omitted =>
    // byte-identical.  See pdvd/docs/nf_sp_img_clus/31_*.md.
    // PDVD PRODUCTION ON, owner flip 2026-09-03 (doc 31 round 6).
    wrapped_channel_charge = true,
    // Event-T0 / readout-tick0 compensation for the live BlobSampler drift-x
    // conversion.  PDVD has no per-event T0, so default 0 (no preset T0).
    // -250us would place trigger-time activity at its true x; see
    // pgrapher/experiment/protodunevd/clus.jsonnet.
    time_offset = 0,
    // Accept any-volume containment in the switch_scope T0Correction filter
    // (out-of-time clusters keep participating in all-APA merge passes).
    relax_containment_filter = true,
    // Interpose charge-light (Q/L) matching between the per-drift-group
    // clustering (stage 3) and the final all-TPC clustering (stage 4).  Default
    // false => the historical no-matching chain, bit-identical.  PDVD runs ONE
    // all-PD flash list (opflash_input) against both drift volumes via the
    // shared-flash joint QLMatching.  See pgrapher/experiment/protodunevd/
    // qlmatching.jsonnet and pdvd/docs/pdvd-qlmatching.md.
    do_qlmatch = false,
    // Full path to the per-event opflash tensor archive (the light chain's
    // work/<RUN>_light<EVENTNO>/opflash_pdvd-wct.tar.gz — NOT in the charge
    // input dir, unlike PDHD).  Required when do_qlmatch.
    opflash_input = '',
    // Hand-scan calibration dump: when true (and do_qlmatch) QLMatching writes
    // work/<output_dir>/calib-evt<event>.json for the pdvd/ql_scan viewer.
    calib = false,
    // Dump the optical "op" bee instance (measured flash PE + Q/L predicted PE
    // per matched cluster) from the all-TPC MABC.  Needs do_qlmatch.
    save_opflash = false,
    // doc pdvd/39 round 2: record what clustering_isolated merged, so the PR job
    // can undo it (pr.jsonnet unmerge_assoc).  Writes the
    // isolated/assoc_cluster_id/assoc_cluster_main perblob arrays; cluster
    // membership is untouched.  NOT free of output change: the saved pctree
    // gains three arrays (its clustering-global Bee member must stay identical --
    // that is the gate).  false => both keys omitted => byte-identical.
    clus_save_assoc_id = false,
    // Per-event light-vs-charge time-base offsets in MICROSECONDS, PER CRATE
    // (opflash metadata offset_bot_us/offset_top_us + the per-run residual;
    // see run_clus_evt.sh).  The BDE (bottom volume, anodes 0-3) and TDE (top
    // volume, anodes 4-7) charge windows open up to ~32 us apart, hence two
    // values.  Applied DOWNSTREAM (matching geometry + T0Correction x_t0cor).
    // Defaults 0.
    trigger_offset_bot_us = 0,
    trigger_offset_top_us = 0,
    // Post-resample readout-window length (ticks) for the Q/L window-truncation
    // flag.  run_clus_evt.sh reads the real value from the SP frame (10000).
    readout_window_ticks = 10000,
    // Q/L visibility backend: 'library' (v5 ANN grid, default) or 'semi'
    // (fitted semi-analytical fallback; membrane XAs auto-masked).
    light_model = 'library',
    // Q/L knob overrides for the trigger-offset diagnostic runs (see
    // pdvd/ql_light_calib/): containment is meaningless until the time base is
    // calibrated, and a higher PE floor bounds the un-pruned bundle count.
    ql_require_containment = true,
    ql_flash_minpe = 25,
    // Per-side calibrated drift speeds in MM/US (bottom volume = anodes 0-3,
    // top = anodes 4-7), sec 8.12 of pdvd-anode-time-consistency.md.  null =>
    // the common legacy 1.568 everywhere -- compiled config byte-identical.
    // When set, BOTH clustering (BlobSampler x_raw + dvm x_t0cor) and the
    // joint QLMatching (per-input drift_speeds) get the same per-side values.
    drift_speed_bot_mmus = null,
    drift_speed_top_mmus = null,
    // Q/L cathode-side containment tolerance in CM (how far past the cathode a
    // drift end may reach and still count contained; C++ default +1.2 cm).
    // null => key omitted => C++ default => compiled config byte-identical.
    // Set only for the anode-pull / cushion study (run 039252 evt298567).
    ql_cathode_ext1_cm = null,
    // Q/L anode-side containment slack in CM: extra tolerance below anode_ext1,
    // forming the anode floor (anode_ext1 - margin) for BOTH the containment gate
    // and the at_x_boundary / close_to_PMT flag window.  C++ default 1.0 cm
    // (floor -3 cm); null => key omitted => compiled config byte-identical.
    // run_clus_evt.sh passes 2.0 cm (floor -4 cm) for PDVD production.
    ql_anode_margin_cm = null,
    // Per-flash DAPHNE-rail channel masking in QLMatching.  Needs a light
    // archive made with flag_saturation (run_light_evt.sh
    // PDVD_FLAG_SATURATION=1).  C++ default false; key suppressed when off =>
    // compiled config byte-identical.  docs/qlmatch/pdvd-saturation-recovery.md.
    ql_use_saturation_flag = false,
    // Whether a rail-flagged channel is also DROPPED from the chi2/KS (the
    // 2026-07-14..16 operating point) or kept there at its clipped PE.
    // Keeping it is right -- the clipped PE is a LOWER BOUND on the true
    // light, not a missing measurement, and the production chain additionally
    // repairs it (PDVD_SAT_REPAIR=1); dropping it makes a wrong prediction
    // free on the bright cathode channels that anchor Q/L matching
    // (docs/qlmatch/11_pdvd-saturation-recovery.md section 6).  The LASSO rows
    // stay zeroed either way.  C++ default true = drop; key suppressed then =>
    // compiled config byte-identical.
    ql_saturation_mask_fit = true,
    // Extra chi2 width on a railed channel: denom += (pe*inflate)^2 (the
    // chi2_pmt_inflate form, gated on the rail flag alone).  null => key
    // suppressed => C++ default 0 = none.  PDVD's pe_err_on_pred/pe_err_frac
    // 0.60 already puts a ~60%-of-predicted sigma on these channels (doc 11
    // section 6), so this is a spare lever, not a requirement.
    ql_chi2_sat_inflate = null,
    // Per-flash readout-coverage tracking in QLMatching.  Needs a light
    // archive made with emit_coverage (run_light_evt.sh
    // PDVD_EMIT_COVERAGE=1).  C++ default false; key suppressed when off =>
    // compiled config byte-identical.
    // docs/qlmatch/14_pdvd-lightpattern-sp-investigation.md.
    ql_use_coverage_flag = false,
    // Whether a readout-uncovered channel is also DROPPED from the fit (the
    // 2026-07-14..16 operating point) or kept at its measured 0.  Keeping it
    // is right -- DAPHNE has no dead time (min inter-snippet gap 0.336 us)
    // and the self-trigger threshold is ~1 PE, so the silence measures
    // "< ~1 PE" rather than nothing (doc 14 section 12;
    // docs/qlmatch/scripts/analyze_coverage_deadtime.py).  The channel is
    // labelled "nodata" in the dump / ql_scan either way.  C++ default true
    // = drop; key suppressed then => compiled config byte-identical.
    ql_coverage_mask_fit = true,
    // PE_err floor (PE) for kept no-data channels.  PDVD does not need it:
    // qlmatching.jsonnet already sets pe_err_floor 2.0 with pe_err_knee at
    // the C++ default 1.0, so a measured 0 carries +-2 PE (~2x the ~1 PE
    // threshold).  null => key suppressed => C++ default -1 = disabled.
    ql_pe_err_nodata = null,
    // Robust-endpoint trim: snap a drift endpoint back inside the detector edge
    // when the outer material is overclustered junk rather than a real track end,
    // before the containment gate reads the endpoints.  ql_robust_trim is the
    // MASTER switch -- the judges below are inert while it is false.  OFF => key
    // suppressed => compiled config byte-identical.  PDVD has not A/B'd or censused
    // this; run_clus_evt.sh turns it on only for the evt298567 clus 34 demo (see
    // docs/qlmatch/15_pdvd-clus34-unmatched-evt298567.md).  Pure pass-through: null
    // => qlmatching omits the key => C++ default.
    ql_robust_trim = false,
    ql_robust_frac = null,
    ql_robust_count = null,
    ql_robust_charge_frac = null,
    ql_robust_charge_abs = null,
    // Gap (detachment) judge, anode end only, in CM.
    ql_robust_gap_cm = null,
    ql_robust_gap_charge_frac = null,
    // Break the anode-end walk at the floor the containment gate actually uses
    // (anode_ext1 - anode_ext1_margin) instead of anode_ext1, so the walk stops
    // counting material the gate would have accepted as "outside".  Fixes the case
    // where a body starting in that dead band drags real track charge into the
    // judges' pile and the trim then correctly refuses (evt298567 apa-4 ident 1 vs
    // flash 120).  Needs ql_robust_trim to do anything.  OFF => key suppressed =>
    // compiled config byte-identical.
    ql_robust_walk_to_floor = false,
    ql_robust_gap_cathode = false,
    // xtpc cathode rescue (doc 16 §10): tolerance in CM past the containment gate /
    // below the at_cathode window granted FOR XTPC CANDIDACY ONLY; a provisionally
    // kept uncontained crosser half is purged pre-fit unless cull_cross_tpc confirms
    // it (scenario 1, one half contained).  qfrac = charge fraction discardable as
    // overclustered junk when measuring the overshoot.  null => keys suppressed =>
    // compiled config byte-identical (C++ default 0 = OFF).  NOT censused; the demo
    // is run_clus_evt.sh PDVD_QL_XTPC_CATHODE_TOL_CM / _QFRAC (evt298567 clus 97).
    ql_xtpc_cathode_tol_cm = null,
    ql_xtpc_cathode_qfrac = null,
    // Cathode-XA-anchored operating point (docs/qlmatch/18_*.md; toolkit knobs in
    // cfg/.../protodunevd/qlmatching.jsonnet).  All four default to the LEGACY
    // behaviour here => compiled config byte-identical when the runner does not
    // pass them; run_clus_evt.sh turns them ON for PDVD production.
    // Mask the live membrane/wall XAs (0,1,3,12,18,19) out of the chi2/LASSO/KS
    // (bimodal response, doc 17).
    ql_mask_wall_xa = false,
    // +-y wall proximity flag (vd_surface_flags wall component): false empties the
    // wall channel lists => no wall flag_close_to_PMT / relax channels.  With the
    // wall XAs masked the relaxation protected nothing and only exempted junk from
    // the overpred prefilter.
    ql_wall_flags = true,
    // Cathode-scoped flash admission on top of ql_flash_minpe: sum(PE over cathode
    // XAs) >= minpe AND >= min_fired cathode channels at >= fired_pe PE.
    ql_flash_sel_cathode = false,
    ql_flash_sel_minpe = null,
    ql_flash_sel_min_fired = null,
    ql_flash_sel_fired_pe = null,
    // Cathode-scoped over-prediction prefilter (R_total / R_max ceilings).
    ql_reject_overpred = false,
    ql_overpred_total_ratio = null,
    ql_overpred_maxch_ratio = null,
    // Per-PD-family PE-error overrides (scan-tuning doc 19; C++
    // pe_err_family_* knob, families [cathode XAs 4-11, PMTs]).  All null
    // (default) => keys suppressed => byte-identical global error model.
    ql_peerr_cath_floor = null,
    ql_peerr_cath_frac = null,
    ql_peerr_cath_lowpe_frac = null,
    ql_peerr_cath_lowpe_knee = null,
    ql_peerr_pmt_floor = null,
    ql_peerr_pmt_frac = null,
    ql_peerr_pmt_lowpe_frac = null,
    ql_peerr_pmt_lowpe_knee = null,
    // Optional third pe_err family: the live membrane/wall XAs
    // (0,1,3,12,18,19), wall-XA inclusion study (docs/qlmatch/25_*.md sec 8;
    // meaningful with ql_mask_wall_xa=false).  All null => two-family (or
    // no-family) arrays unchanged => byte-identical.
    ql_peerr_wall_floor = null,
    ql_peerr_wall_frac = null,
    ql_peerr_wall_lowpe_frac = null,
    ql_peerr_wall_lowpe_knee = null,
    // Per-channel multiplier on the MEASURED PE (length-40 array; Opflash
    // read).  null => key suppressed => byte-identical (C++ empty = identity).
    ql_measured_pe_scale = null,
    // QLMatching QtoL, forwarded to qlmatching.jsonnet `qtol` (toolkit default
    // 0.094; the key is always emitted, so 0.094 here => compiled config
    // byte-identical).  run_clus_evt.sh PDVD_QTOL overrides.  Doc pdvd/100 sec 8.
    ql_qtol = 0.094,
    // xtpc / selection quality gates (scan-tuning doc 19); null/false =>
    // keys suppressed => byte-identical legacy behaviour.
    ql_xtpc_pin_min_strength = null,
    ql_xtpc_sc1_light_gate = false,
    ql_xtpc_sc1_ks_max = null,
    ql_xtpc_sc1_c2n_max = null,
    ql_xtpc_cathode_ks_max = null,
    ql_xtpc_pin_confirms_rescue = false,
    ql_cathode_rescue_solo_ks = null,
    ql_cathode_rescue_solo_c2n = null,
    ql_sc1_evict_ks_margin = null,
    ql_postcull_unflagged = false,
    ql_postcull_ks_max = null,
    ql_postcull_c2n_max = null,
    // Rescue blind-spot fix (doc 23 phase 1a): early postcull pass before the
    // rescues. C++ default false; key omitted when off => byte-identical.
    ql_postcull_before_rescue = false,
    // Window-truncated overprediction cull (doc 23 phase 2); false/null =>
    // keys suppressed => byte-identical.
    ql_postcull_wtrunc = false,
    ql_postcull_wtrunc_ratio_hi = null,
    ql_postcull_wtrunc_sat_frac = null,
    // xtpc-pin overprediction cull (doc 23 phase 2); false/null => keys
    // suppressed => byte-identical.
    ql_postcull_pin = false,
    ql_postcull_pin_ratio_hi = null,
    // Sweepable ladder ceilings + LASSO regularization (doc 19 phase 4);
    // null => the operating literals / C++ defaults, compiled JSON unchanged.
    ql_hc_clean_ks = null, ql_hc_clean_c2 = null,
    ql_hc_good_ks = null, ql_hc_good_c2 = null,
    ql_hc_tb_ks = null, ql_hc_tb_c2 = null,
    ql_hc_miss_ks = null, ql_hc_miss_c2 = null,
    ql_hc_miss_min_ndf = null,
    ql_lasso_lambda = null, ql_delta_charge = null, ql_delta_light = null,
    ql_delta_shape = null, ql_bkg_weight = null, ql_strength_cutoff = null,
    ql_lasso_boundary_weight = null,
    // Shared-flash-aware rescues (doc 19 phase 5); false/null => keys
    // suppressed => byte-identical.
    ql_empty_rescue_shared = false,
    ql_rescue_metric_max = null,
    ql_cluster_rescue_shared = false,
    ql_cluster_rescue_ks_max = null,
    ql_cluster_rescue_c2n_max = null,
    ql_cluster_rescue_ratio_lo = null,
    ql_cluster_rescue_ratio_hi = null,
    ql_cluster_rescue_precull = false,
    ql_cluster_rescue_precull_additive = false,
    // Relaxed second-chance rescue tier (doc 21); false/null => keys
    // suppressed => byte-identical.  min_len in CM (converted to native
    // units here).
    ql_cluster_rescue_relaxed = false,
    ql_cluster_rescue_relaxed_ks_max = null,
    ql_cluster_rescue_relaxed_c2n_max = null,
    ql_cluster_rescue_relaxed_ratio_lo = null,
    ql_cluster_rescue_relaxed_ratio_hi = null,
    ql_cluster_rescue_relaxed_minlen_cm = null,
    // Saturation-aware rescue ratio-high extension (doc 23 phase 1b);
    // false/null => keys suppressed => byte-identical.
    ql_cluster_rescue_sat_relax = false,
    ql_cluster_rescue_sat_frac_min = null,
    ql_cluster_rescue_sat_ratio_mult = null,
    // Post-QLMatching cathode-crossing STITCH tip-touch relaxation (stage-4
    // ClusteringCathodeConnect, cfg/.../protodunevd/clus.jsonnet).  Ports the
    // PDHD tip-touch branch to PDVD: when a genuine crosser's two halves reach
    // the cathode and nearly touch (~1cm gap), cc_tip_touch_cut (CM) drops the
    // uninformative cc_pca connection-alignment term and cc_tip_touch_angle_cut
    // (DEG) accepts on the local charge-weighted Hough even when a curved half
    // inflates the global PCA above angle_cut.  null => keys suppressed => C++
    // defaults tip_touch_cut=0 / tip_touch_angle_cut=angle_cut (OFF) =>
    // compiled config byte-identical.  run_clus_evt.sh sets the PDVD operating
    // point.  Only meaningful with do_qlmatch (needs the matched cluster_t0).
    // _cm is converted to internal units below; _angle_cut is in DEGREES.
    cc_tip_touch_cut_cm = null,
    cc_tip_touch_angle_cut = null,
    // Post-QLMatching cathode-crossing STITCH distance gates (stage-4
    // ClusteringCathodeConnect).  PDVD's cathode is ~6 cm thick, so a genuine
    // crosser's two halves stop at their own cathode face (~+-3 cm) and meet
    // 6-15 cm apart -- past the SBND-tuned dis_cut=5 / drift_cut=8 /
    // cathode_x_cut=5.  Data-driven from the 151 QL-confirmed xtpc_pin pairs
    // (28 evts): dis p95=15, dX p95=12, tip|x| p95=6.6.  null => clus.jsonnet
    // keeps 5/8/5 => byte-identical.  Set in CM by run_clus_evt.sh census point.
    cc_cathode_x_cut_cm = null,
    cc_drift_cut_cm = null,
    cc_dis_cut_cm = null,
    // 6cm-cathode crosser cc_pca relaxation (degrees): loose connection-vector
    // bound for drift-gated crossers.  null => clus.jsonnet keeps C++ default 0
    // (relaxation OFF) => byte-identical.  Census operating point ~65 deg.
    cc_crosser_conn_relax = null,
    // 6cm-cathode crosser tt_pca bound (degrees): admit bent crossers. null =>
    // clus.jsonnet keeps C++ default 0 (bound stays angle_cut) => byte-identical.
    // Census operating point ~15 deg.
    cc_crosser_pca_angle = null,
    // 6cm-cathode near-cathode closest-approach retry (cm). null => clus.jsonnet keeps
    // C++ default 0 (retry OFF) => byte-identical. Census operating point ~10 cm.
    cc_cathode_band_dis = null,
    // Emit a second pre-pipeline raw-imaging Bee set split by drift side
    // (img-side-bot / img-side-top) alongside img-global, for per-side video
    // frames.  false => bee_points_sets unchanged in clus.jsonnet => byte-identical.
    bee_img_per_side = false,
    // Persist the post-Q/L point-cloud tree for the PR job (doc pdvd/25 M1).
    // '' => the inert dump_mode sink in clus.jsonnet is unchanged => byte-identical.
    // Runner: run_clus_evt.sh -save-pctree => work/<RUN6>_<EVT>/pctree-evt<ID>.tar.gz.
    pctree_outname = '',
)

local params_w = if wires_file == '' then params else params { files+: { wires: wires_file } };
local tools = if wires_file == '' then tools_all else tools_maker(params_w);
local anodes = [tools.anodes[i] for i in anode_indices];

local cluster_source(fname) = g.pnode({
    type: "ClusterFileSource",
    name: fname,
    data: {
        inname: fname,
        anodes: [wc.tn(a) for a in anodes],
    }
}, nin=0, nout=1, uses=anodes);

local trigger_offset_bot = trigger_offset_bot_us * wc.us;
local trigger_offset_top = trigger_offset_top_us * wc.us;
// Per-side drift speeds (null => legacy common 1.568 inside clus.jsonnet /
// scalar drift_speed inside QLMatching; keys/args omitted => byte-identical).
local drift_speed_bot = if drift_speed_bot_mmus == null then null
                        else drift_speed_bot_mmus * wc.mm / wc.us;
local drift_speed_top = if drift_speed_top_mmus == null then null
                        else drift_speed_top_mmus * wc.mm / wc.us;
local clus = import 'pgrapher/experiment/protodunevd/clus.jsonnet';
local clus_maker = clus(output_dir=output_dir, runNo=run, subRunNo=subrun, eventNo=event, stepped_center_fallback=stepped_center_fallback,
                        wrapped_channel_charge=wrapped_channel_charge,
                        time_offset=time_offset, relax_containment_filter=relax_containment_filter,
                        trigger_offset=trigger_offset_bot, trigger_offset_top=trigger_offset_top,
                        drift_speed_b=drift_speed_bot, drift_speed_t=drift_speed_top);

// Drift-side groups: anodes 0-3 (bottom drift) and anodes 4-7 (top drift).
// With a subset anode_indices only non-empty groups are built and the final
// merge multiplicity shrinks to match.
local group_defs = [
    { name: "group0123", anodes: [a for a in anodes if a.data.ident < 4] },
    { name: "group4567", anodes: [a for a in anodes if a.data.ident >= 4] },
];
local groups = [gd for gd in group_defs if std.length(gd.anodes) > 0];
local ngroups = std.length(groups);

// One drift-side pipe: per-anode sources -> per_apa (stages 1+2) -> per_group (stage 3).
local group_pipe(gd) =
    local n = std.length(gd.anodes);
    local actives = [cluster_source("%s/clusters-apa-anode%d-ms-active.tar.gz"%[input, a.data.ident]) for a in gd.anodes];
    local maskeds = [cluster_source("%s/clusters-apa-anode%d-ms-masked.tar.gz"%[input, a.data.ident]) for a in gd.anodes];
    local apa_pipes = [clus_maker.per_apa(gd.anodes[i], dump=false, po_fast=ql_po_fast, dg_fast=ql_dg_fast) for i in std.range(0, n - 1)];
    local pg = clus_maker.per_group(gd.anodes, gd.name, dump=false, dg_fast=ql_dg_fast,
                                    save_assoc_id=clus_save_assoc_id);
    g.intern(
        innodes = actives + maskeds,
        centernodes = apa_pipes,
        outnodes = [pg],
        edges =
            [g.edge(actives[i], apa_pipes[i], 0, 0) for i in std.range(0, n - 1)] +
            [g.edge(maskeds[i], apa_pipes[i], 0, 1) for i in std.range(0, n - 1)] +
            [g.edge(apa_pipes[i], pg, 0, i) for i in std.range(0, n - 1)]
    );

local group_pipes = [group_pipe(gd) for gd in groups];
// With Q/L matching on, ONE joint QLMatching node emits a single pre-merged tree, so
// the all-TPC stage skips its input PointTreeMerging (premerged).  Without matching
// the two per-side clustering outputs still fan into the ngroups-way merge.
// null => byte-identical (key suppressed in clus.jsonnet); else cm -> internal.
local cc_tip_touch_cut = if cc_tip_touch_cut_cm == null then null
                         else cc_tip_touch_cut_cm * wc.cm;
// Distance gates: null => the clus.jsonnet defaults (5/8/5 cm) => byte-identical.
local cc_cathode_x_cut = if cc_cathode_x_cut_cm == null then 5*wc.cm else cc_cathode_x_cut_cm * wc.cm;
local cc_drift_cut = if cc_drift_cut_cm == null then 8*wc.cm else cc_drift_cut_cm * wc.cm;
local cc_dis_cut = if cc_dis_cut_cm == null then 5*wc.cm else cc_dis_cut_cm * wc.cm;
local clus_all_tpc = if do_qlmatch
    then clus_maker.all_tpc(anodes, premerged=true, save_opflash=save_opflash,
                            cc_tip_touch_cut=cc_tip_touch_cut, cc_tip_touch_angle_cut=cc_tip_touch_angle_cut,
                            cc_cathode_x_cut=cc_cathode_x_cut, cc_drift_cut=cc_drift_cut, cc_dis_cut=cc_dis_cut,
                            cc_crosser_conn_relax=cc_crosser_conn_relax, cc_crosser_pca_angle=cc_crosser_pca_angle,
                            cc_cathode_band_dis=cc_cathode_band_dis, bee_img_per_side=bee_img_per_side,
                            tensor_outname=pctree_outname, save_assoc_id=clus_save_assoc_id)
    else clus_maker.all_tpc(anodes, ngroups=ngroups,
                            cc_tip_touch_cut=cc_tip_touch_cut, cc_tip_touch_angle_cut=cc_tip_touch_angle_cut,
                            cc_cathode_x_cut=cc_cathode_x_cut, cc_drift_cut=cc_drift_cut, cc_dis_cut=cc_dis_cut,
                            cc_crosser_conn_relax=cc_crosser_conn_relax, cc_crosser_pca_angle=cc_crosser_pca_angle,
                            cc_cathode_band_dis=cc_cathode_band_dis, bee_img_per_side=bee_img_per_side,
                            tensor_outname=pctree_outname, save_assoc_id=clus_save_assoc_id);

// JOINT Q/L matching (shared-flash): both drift sides enter ONE QLMatching node and
// each reads the SAME all-PD opflash archive.  Per side: opflash source ->
// flash_attach (2->1 fan-in: port0 = group cluster tree, port1 = opflash matrix);
// flash_attach[i] feeds the joint node's input port i.  The joint node runs ONE
// LASSO over the shared flashes (clusters from both sides explain each flash) and
// emits one merged tree (-> premerged clus_all_tpc).  PDVD imaging face is 0 for
// both drift sides (the two faces of a CRP share x-bounds).
local qlm = import 'pgrapher/experiment/protodunevd/qlmatching.jsonnet';
local qlm_maker = qlm(params_w, trigger_offset_bot, readout_window_ticks, light_model,
                      ql_require_containment, ql_flash_minpe,
                      trigger_offsets=[trigger_offset_bot, trigger_offset_top],
                      // Single recalibrated velocity => scalar override (the
                      // calib dump's d["drift_speed"] then carries it for the
                      // viewers/scripts); genuinely split values => per-input
                      // array [bottom, top] (readers use per-geometry speeds).
                      drift_speed=if drift_speed_bot != null
                                     && drift_speed_bot == drift_speed_top
                                  then drift_speed_bot else null,
                      drift_speeds=if drift_speed_bot == null
                                      || drift_speed_top == null
                                      || drift_speed_bot == drift_speed_top
                                   then null
                                   else [drift_speed_bot, drift_speed_top],
                      // null => qlmatching omits the key => C++ default +1.2 cm.
                      cathode_ext1=if ql_cathode_ext1_cm == null then null
                                   else ql_cathode_ext1_cm * wc.cm,
                      // null => qlmatching omits the key => C++ default 1.0 cm.
                      anode_ext1_margin=if ql_anode_margin_cm == null then null
                                        else ql_anode_margin_cm * wc.cm,
                      use_saturation_flag=ql_use_saturation_flag,
                      saturation_mask_fit=ql_saturation_mask_fit,
                      chi2_sat_inflate=ql_chi2_sat_inflate,
                      use_coverage_flag=ql_use_coverage_flag,
                      coverage_mask_fit=ql_coverage_mask_fit,
                      pe_err_nodata=ql_pe_err_nodata,
                      // Master switch off => qlmatching suppresses the key and the
                      // whole trim block stays inert whatever the judges say.
                      robust_endpoint_trim=ql_robust_trim,
                      robust_endpoint_frac=ql_robust_frac,
                      robust_endpoint_count=ql_robust_count,
                      robust_endpoint_charge_frac=ql_robust_charge_frac,
                      robust_endpoint_charge_abs=ql_robust_charge_abs,
                      robust_endpoint_gap=if ql_robust_gap_cm == null then null
                                          else ql_robust_gap_cm * wc.cm,
                      robust_endpoint_gap_charge_frac=ql_robust_gap_charge_frac,
                      robust_endpoint_walk_to_floor=ql_robust_walk_to_floor,
                      robust_endpoint_gap_cathode=ql_robust_gap_cathode,
                      xtpc_cathode_tol=if ql_xtpc_cathode_tol_cm == null then null
                                       else ql_xtpc_cathode_tol_cm * wc.cm,
                      xtpc_cathode_qfrac=ql_xtpc_cathode_qfrac,
                      // Cathode-XA operating point (docs/qlmatch/18_*.md); all
                      // legacy-default here, turned on by run_clus_evt.sh.
                      mask_wall_xa=ql_mask_wall_xa,
                      wall_flags=ql_wall_flags,
                      flash_sel_cathode=ql_flash_sel_cathode,
                      flash_sel_minPE=ql_flash_sel_minpe,
                      flash_sel_min_fired=ql_flash_sel_min_fired,
                      flash_sel_fired_pe=ql_flash_sel_fired_pe,
                      reject_overpred=ql_reject_overpred,
                      overpred_total_ratio=ql_overpred_total_ratio,
                      overpred_maxch_ratio=ql_overpred_maxch_ratio,
                      // Per-PD-family PE-error overrides (doc 19); all null
                      // => qlmatching suppresses the family keys.
                      pe_err_cath_floor=ql_peerr_cath_floor,
                      pe_err_cath_frac=ql_peerr_cath_frac,
                      pe_err_cath_lowpe_frac=ql_peerr_cath_lowpe_frac,
                      pe_err_cath_lowpe_knee=ql_peerr_cath_lowpe_knee,
                      pe_err_pmt_floor=ql_peerr_pmt_floor,
                      pe_err_pmt_frac=ql_peerr_pmt_frac,
                      pe_err_pmt_lowpe_frac=ql_peerr_pmt_lowpe_frac,
                      pe_err_pmt_lowpe_knee=ql_peerr_pmt_lowpe_knee,
                      pe_err_wall_floor=ql_peerr_wall_floor,
                      pe_err_wall_frac=ql_peerr_wall_frac,
                      pe_err_wall_lowpe_frac=ql_peerr_wall_lowpe_frac,
                      pe_err_wall_lowpe_knee=ql_peerr_wall_lowpe_knee,
                      measured_pe_scale=ql_measured_pe_scale,
                      qtol=ql_qtol,
                      // xtpc / selection quality gates (doc 19).
                      xtpc_pin_min_strength=ql_xtpc_pin_min_strength,
                      xtpc_sc1_light_gate=ql_xtpc_sc1_light_gate,
                      xtpc_sc1_ks_max=ql_xtpc_sc1_ks_max,
                      xtpc_sc1_c2n_max=ql_xtpc_sc1_c2n_max,
                      xtpc_cathode_ks_max=ql_xtpc_cathode_ks_max,
                      xtpc_pin_confirms_rescue=ql_xtpc_pin_confirms_rescue,
                      cathode_rescue_solo_ks=ql_cathode_rescue_solo_ks,
                      cathode_rescue_solo_c2n=ql_cathode_rescue_solo_c2n,
                      sc1_evict_ks_margin=ql_sc1_evict_ks_margin,
                      postcull_unflagged=ql_postcull_unflagged,
                      postcull_ks_max=ql_postcull_ks_max,
                      postcull_c2n_max=ql_postcull_c2n_max,
                      postcull_before_rescue=ql_postcull_before_rescue,
                      postcull_wtrunc_overpred=ql_postcull_wtrunc,
                      postcull_wtrunc_ratio_hi=ql_postcull_wtrunc_ratio_hi,
                      postcull_wtrunc_sat_frac=ql_postcull_wtrunc_sat_frac,
                      postcull_pin_overpred=ql_postcull_pin,
                      postcull_pin_ratio_hi=ql_postcull_pin_ratio_hi,
                      // Sweepable ladder + LASSO regularization (doc 19).
                      hc_clean_ks=ql_hc_clean_ks, hc_clean_c2=ql_hc_clean_c2,
                      hc_good_ks=ql_hc_good_ks, hc_good_c2=ql_hc_good_c2,
                      hc_tb_ks=ql_hc_tb_ks, hc_tb_c2=ql_hc_tb_c2,
                      hc_miss_ks=ql_hc_miss_ks, hc_miss_c2=ql_hc_miss_c2,
                      hc_miss_min_ndf=ql_hc_miss_min_ndf,
                      lasso_lambda=ql_lasso_lambda,
                      delta_charge=ql_delta_charge,
                      delta_light=ql_delta_light,
                      delta_shape=ql_delta_shape,
                      bkg_weight=ql_bkg_weight,
                      strength_cutoff=ql_strength_cutoff,
                      lasso_boundary_weight=ql_lasso_boundary_weight,
                      // Shared-flash-aware rescues (doc 19 phase 5).
                      empty_rescue_shared=ql_empty_rescue_shared,
                      rescue_metric_max=ql_rescue_metric_max,
                      cluster_rescue_shared=ql_cluster_rescue_shared,
                      cluster_rescue_ks_max=ql_cluster_rescue_ks_max,
                      cluster_rescue_chi2ndf_max=ql_cluster_rescue_c2n_max,
                      cluster_rescue_ratio_lo=ql_cluster_rescue_ratio_lo,
                      cluster_rescue_ratio_hi=ql_cluster_rescue_ratio_hi,
                      cluster_rescue_precull=ql_cluster_rescue_precull,
                      cluster_rescue_precull_additive=ql_cluster_rescue_precull_additive,
                      // Relaxed second-chance tier (doc 21).
                      cluster_rescue_relaxed=ql_cluster_rescue_relaxed,
                      cluster_rescue_relaxed_ks_max=ql_cluster_rescue_relaxed_ks_max,
                      cluster_rescue_relaxed_chi2ndf_max=ql_cluster_rescue_relaxed_c2n_max,
                      cluster_rescue_relaxed_ratio_lo=ql_cluster_rescue_relaxed_ratio_lo,
                      cluster_rescue_relaxed_ratio_hi=ql_cluster_rescue_relaxed_ratio_hi,
                      cluster_rescue_relaxed_min_length=
                          if ql_cluster_rescue_relaxed_minlen_cm == null then null
                          else ql_cluster_rescue_relaxed_minlen_cm * wc.cm,
                      // Saturation-aware rescue ratio-high extension (doc 23 phase 1b).
                      cluster_rescue_sat_ratio_relax=ql_cluster_rescue_sat_relax,
                      cluster_rescue_sat_frac_min=ql_cluster_rescue_sat_frac_min,
                      cluster_rescue_sat_ratio_mult=ql_cluster_rescue_sat_ratio_mult);
local calib_dump_joint =
    if calib then '%s/calib-evt%s.json' % [output_dir, std.toString(event)]
    else '';

local graph = if do_qlmatch then
    local srcs = [qlm_maker.opflash_source(groups[i].name, opflash_input)
                  for i in std.range(0, ngroups - 1)];
    local atts = [qlm_maker.flash_attach(groups[i].name) for i in std.range(0, ngroups - 1)];
    local dv_all = clus_maker.detector_volumes(anodes);
    local sides  = [gd.anodes for gd in groups];
    local faces  = [0 for gd in groups];
    local mat = qlm_maker.matching_joint(sides, dv_all, faces, calib_dump_joint);
    g.intern(
        innodes=group_pipes,
        centernodes=srcs + atts + [mat],
        outnodes=[clus_all_tpc],
        edges=
            [g.edge(srcs[i], atts[i], 0, 1) for i in std.range(0, ngroups - 1)] +
            [g.edge(group_pipes[i], atts[i], 0, 0) for i in std.range(0, ngroups - 1)] +
            [g.edge(atts[i], mat, 0, i) for i in std.range(0, ngroups - 1)] +
            [g.edge(mat, clus_all_tpc, 0, 0)]
    )
else g.intern(
    innodes=group_pipes,
    outnodes=[clus_all_tpc],
    edges=[g.edge(group_pipes[i], clus_all_tpc, 0, i) for i in std.range(0, ngroups - 1)]
);

local app = {
  type: 'Pgrapher',
  data: {
    edges: g.edges(graph),
  },
};

local cmdline = {
    type: "wire-cell",
    data: {
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellSio", "WireCellSigProc", "WireCellImg", "WireCellRoot", "WireCellTbb", "WireCellClus", "WireCellMatch", "WireCellAux"],
        apps: ["Pgrapher"]
    }
};

[cmdline] + g.uses(graph) + [app]
