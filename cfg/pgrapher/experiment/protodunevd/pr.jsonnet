// PDVD pattern-recognition (PR) builder -- the tail that runs AFTER Q/L matching
// on the persisted point-cloud tree (doc pdvd/25 secs 4-5).  Forked BY DUPLICATION
// from the SBND pr() builder in pgrapher/experiment/sbnd/clus.jsonnet (HEAD
// 0913743f, lines 847-2665); the SBND file is untouched.  The geometry objects
// (DetectorVolumes, PCTransformSet, the per-(anode,face) 'stepped' samplers with
// the per-crate drift speeds, the T0 scope coordinates) are taken from the PDVD
// clustering module so the PR job sees exactly what the clustering built.
//
// PDVD deltas vs SBND, each stated where it lives below:
//   * no unmerge_bundle / unmerge_assoc stages (PDVD never runs examine_bundles
//     nor writes assoc_cluster_id);
//   * 16 retiler samplers (8 two-sided anodes), per-crate drift speed;
//   * pdvd_pr_fv: one BoxFiducial spanning both drift volumes;
//   * pdvd_box_recomb at 0.44 kV/cm with the table parameter set;
//   * nticks is an argument (PDVD readout window 10000, SBND 3427);
//   * Bee detector 'protodunevd'; single-event RSE from the job TLAs.
// Every knob keeps the SBND key-suppression idiom: an unset knob is absent from
// the compiled config and means "C++ default".  The SBND-tuned DEFAULTS of the
// pr() arguments are kept verbatim (they document the SBND operating point);
// the PDVD operating point is set in pdvd/wct-pr-perevt.jsonnet, which passes
// every physics knob explicitly.
//
// Import as
//   local pr_mod = import 'pgrapher/experiment/protodunevd/pr.jsonnet';
//   local pr_maker = pr_mod(output_dir=..., runNo=..., subRunNo=..., eventNo=...,
//                           drift_speed_b=..., drift_speed_t=..., trigger_offset=...,
//                           trigger_offset_top=..., time_offset=...);
//   pr_maker.pr(anodes, pipeline_names=[...], ...)

local wc = import 'wirecell.jsonnet';
local g = import 'pgraph.jsonnet';
local clus = import 'pgrapher/common/clus.jsonnet';
local clus_mod = import 'pgrapher/experiment/protodunevd/clus.jsonnet';

function(output_dir='', runNo=1, subRunNo=1, eventNo=1, stepped_center_fallback=false,
         time_offset=0 * wc.us, relax_containment_filter=true,
         trigger_offset=0 * wc.us, trigger_offset_top=null,
         drift_speed_b=null, drift_speed_t=null) {
    // The PDVD clustering module, configured EXACTLY as the Q/L job configured it
    // (same drift speeds, trigger offsets, time offset): switch_scope re-derives
    // x_t0cor from cluster_t0 through these, so they must match for the M2
    // round-trip gate.
    local clus_maker = clus_mod(output_dir=output_dir, runNo=runNo, subRunNo=subRunNo, eventNo=eventNo,
                                stepped_center_fallback=stepped_center_fallback,
                                time_offset=time_offset, relax_containment_filter=relax_containment_filter,
                                trigger_offset=trigger_offset, trigger_offset_top=trigger_offset_top,
                                drift_speed_b=drift_speed_b, drift_speed_t=drift_speed_t),
    local bee_dir = if output_dir == '' then 'data' else output_dir,
    local evt_out_prefix = bee_dir + '/',

    pr(anodes, dump=true, pipeline_names=[], tensor_outname='',
              // Readout length in ticks for the Magnify/PrDisplay writers (SBND 3427;
              // PDVD production window 10000, run_clus_evt.sh readout_window_ticks).
              nticks=10000,
              // flag_mains_min_length (internal units): PDVD's Q/L never builds a
              // main+associated decomposition, so QLMatching flags no main; the
              // 'flag_mains' stage (ClusteringFlagMatchedMains) flags every
              // flash-matched cluster (matched_flash_gid >= 0 with a real t0)
              // at least this long as a main.  0 = every matched cluster.
              flag_mains_min_length=0,
              // PDVD boundary vetoes for the STM verdict (doc 25 M3).  All C++
              // default OFF; keys omitted when off => byte-identical config.
              // readout_edge_guard: the stop's fitted arrival tick within
              // stm_readout_edge_ticks of [0, nticks) = readout-window truncation.
              // stm_cathode_guard_cm: overrides TaggerCheckSTM guard_cathode_cm (C++
              // 5 cm; PDVD's cathode is a 6 cm slab, so 12 cm) when cathode_guard is on.
              stm_readout_edge_guard=false,
              stm_readout_edge_ticks=60,
              stm_cathode_guard_cm=null,
              // Steiner-terminal per-point charge floor (electrons; CreateSteinerGraph
              // terminal_charge_threshold, C++ default 4000 = prototype).  null =>
              // key omitted => byte-identical.  PDVD needs its own (doc 25 M3).
              steiner_terminal_charge=null,
              // save_in_scope (doc 87): add the per-cluster T_cluster tree to
              // tracking-pr.root -- the in-scope set (switch_scope's scope_filter,
              // the SAME predicate the Bee clustering layer is gated on) plus the
              // per-bundle summary.  It is what lets an arm run with mabc-pr.zip
              // and the pctree suppressed and still produce a full nusel table.
              // C++ default false.  Key omitted when off => byte-identical config.
              save_in_scope=false,
              // pr_bee (doc 87): write pr_evt<ID>/mabc-pr.zip.  false => bee_zip is the
              // empty string, which MultiAlgBlobClustering reads as "write no Bee
              // zip at all" (an empty name used to raise IOError, so it was never
              // a legal value).  DEFAULT TRUE = today's behaviour, byte-identical.
              pr_bee=true,
              // trackfitting_config_file: the SBND TrackFitting parameter JSON.
              // DEFAULT = the canonical in-tree file, resolved through
              // WIRECELL_PATH by TaggerCheckSTM/TaggerCheckNeutrino
              // (Persist::resolve; an absolute path is passed through unchanged).
              // '' selects the C++ preset defaults, which are uBooNE-hard-coded --
              // never right for SBND, which is why this no longer defaults to ''.
              trackfitting_config_file='pgrapher/experiment/protodunevd/pdvd_track_fitting.json',
              particle_dataset=null, extra_uses=[],
              // dl_weights: SCN (DL) neutrino-vertex weights, WIRECELL_PATH-resolved.
              // DEFAULT = the uBooNE-trained net, i.e. the DL vertex is ON for SBND
              // (owner adopted 2026-07-30 on nueCC48 evt 18253/1/172230, where the
              // geometric vertex sat at the far end of a proton track and the DL
              // vertex moved it 9.7 cm onto the true interaction point -- docs/pr/4).
              // '' restores the geometric-only vertex (still what every identity
              // gate must pass, CLAUDE.md M4).  The weights remain uBooNE-TRAINED:
              // SBND retraining is docs/pr/2 gap G3 and stays open.
              // REQUIRES libpython preloaded in the job (LD_PRELOAD=libpython3.11.so.1.0);
              // without it the SCN module import fails and the code falls back to the
              // geometric vertex with a single WARN line -- grep the log for
              // "DL vertex failed" (expect none).
              dl_weights='uboone/scn_vtx/t48k-m16-l5-lr5d-res0.5-CP24.pth',
              // DL re-rank operating point, threaded for A/B runs (doc pr/79).
              // min_accept: rerank-total acceptance threshold for the DL route
              // (TaggerCheckNeutrino); top_k: DL voxel candidates admitted to
              // the reranker.
              // min_accept 4.0 -> 10.0 adopted 2026-08-15 (owner, doc pr/79):
              // live A/B on the 473-label hand-scan measured +36/473 correct
              // (322 -> 358; 51 fixed / 15 regressed, every regression's
              // baseline DL total in [4,10)).  Pass 4.0 for the pre-flip arm.
              dl_vtx_min_accept_score=10.0,
              dl_vtx_top_k=5,
              // doc pr/105: dl_vtx_rerank=false selects the legacy single-argmax
              // DL branch (top voxel snapped, dl_vtx_cut gate, else the
              // traditional vertex) -- strategy 3.1 of the vertex-selection
              // comparison.  Default true = the pinned production value, so the
              // compiled JSON is byte-identical when unset.
              dl_vtx_rerank=true,
              // beam_window: internal-unit [low, high] on the matched flash time
              // (cluster_t0).  DEFAULT = the SBND BNB gate after the
              // frame_apply_at_caf correction, which is what production passes
              // (run_nusel_evt.sh BEAM_WINDOW="0.2,2.2").  [0,0] disables it and
              // makes beam_window_only inert.
              beam_window=[0.2 * wc.us, 2.2 * wc.us],
              // --- TGM/FC knobs.  Every default below is the SBND PRODUCTION
              // operating point (run_full1k_nusel.sh's flag set: -chord -rescue
              // -rescue-chord -fvz 5 -fvzi 3 -main-pair-real -fvx 2.5 -fvy 3,
              // plus NUCAND/CHORD_MODE defaults), adopted as the module defaults
              // 2026-07-27 (owner; sbnd_xin/docs/64_cfg-sync.md).  They used to
              // sit at the pre-adoption values here and were reached only via
              // runner flags, so cfg/ advertised a configuration nobody ran.
              // Each C++ knob still defaults off/legacy, so other detectors are
              // unaffected; pass the old value here for a pre-adoption A/B:
              //   tgm_neutrino_candidate false, tgm_chord_charge false,
              //   tgm_chord_mode 'chord', tgm_component_extremes false,
              //   tgm_component_rescue false, tgm_rescue_chord false,
              //   tgm_main_pair false, tgm_main_pair_mode 'path',
              //   tgm_fv_zmax_margin 3, tgm_fv_zmax_margin_interior 0,
              //   tgm_fv_x_margin 2, tgm_fv_y_margin 2.5. ---
              tgm_neutrino_candidate=true, tgm_chord_charge=true,
              tgm_chord_mode='path', tgm_component_extremes=true,
              tgm_component_rescue=true, tgm_rescue_chord=true,
              tgm_main_pair=true, tgm_main_pair_mode='real',
              tgm_fv_zmax_margin=5, tgm_fv_zmax_margin_interior=3,
              tgm_fv_x_margin=2.5, tgm_fv_y_margin=3,
              save_stm_fit=false, unmerge_bundle_mode='real',
              // doc pr/34 §10 particle-flow (Bee mc tree) port-fidelity knobs.
              // C++ defaults false; keys omitted when off => byte-identical
              // pre-knob config.  Display-only stage: mc.json is the artifact.
              pf_track_main_cluster_only=false,
              pf_track_bridged_clusters=false,  // doc pr/40 round 9 B2; C++ default false. Key omitted when off => byte-identical.
              pf_shower_vertex_barrier=false,
              pf_shower_parent_precedence=false,
              pf_pi0_node_per_id=false,
              pf_pdg_name_prototype_fallback=false,
              // doc pr/38 Round 4: graph-faithful parentage for barrier-orphaned
              // PF track segments.  C++ default false; key omitted when off =>
              // byte-identical pre-knob config.  Display-only (mc.json).
              pf_orphan_track_parentage=false,
              // doc pr/65 round 3: audit-only orphan net -- log unclaimed
              // segments, fabricate no PF-root node.  C++ default false; key
              // omitted when off => byte-identical.  Display-only (mc.json).
              pf_orphan_audit_only=false,
              // doc pr/84 round 2 (F1/F2): vertex-touching pseudo-parent
              // suppression + remote-gap anchor.  C++ defaults false/3cm/8cm;
              // keys omitted when off/null => byte-identical pre-knob config.
              // Display-only (mc.json).  Scalar params are in cm.
              pf_direct_when_touching=false,
              pf_touch_max=null,
              // doc 77 round 1 (2026-08-24): pf_touch_cross_main/_max
              // removed -- zero movers, F1.0 probe failure (pr/84:607/622).
              pf_pseudo_gap_from_main=false,
              // doc pr/84 round 3 (G1): guarantee unique jsTree node ids in
              // mc.json.  C++ default false; key omitted when off =>
              // byte-identical pre-fix config.
              pf_unique_node_ids=false,
              // pf_drop_stray_satellites (doc pr/92): skip showers the kine
              // side's kine_drop_stray_satellites gate dropped from Enu, so
              // the PF tree and Enu describe the same particle set.  Inert
              // unless that kine knob is also on.  C++ default false; key
              // omitted when off => byte-identical pre-pr/92 config.
              pf_drop_stray_satellites=false,
              // pf_orphan_confident_track (doc pr/93 round 4): emit a root PF
              // node for an unclaimed main-cluster segment with a confident
              // non-electron template PID, length > pf_orphan_track_min_cm
              // (null => C++ 50cm), and straight-long (SBND 18255-315167's
              // freed 150.7cm proton).  C++ default false; keys omitted when
              // off => byte-identical.
              pf_orphan_confident_track=false,
              pf_orphan_track_min_cm=null,
              // pf_orphan_guard_freed (doc pr/123 round 2): root PF node
              // for a pass4-track-guard-freed segment (kPass4GuardFreed);
              // C++ default false. Key omitted when off => byte-identical.
              pf_orphan_guard_freed=false,
              // pf_orphan_near_cross_cluster (doc pr/128): PF node for an
              // unclaimed CROSS-CLUSTER track segment whose fit points come
              // within pf_orphan_near_gap_cm (null => C++ 5cm) of the emitted
              // candidate and which is longer than pf_orphan_near_min_len_cm
              // (null => C++ 30cm).  Every orphan pool is same_cluster gated,
              // so this class reaches no output (SBND 18255-137238's muon).
              // C++ default false; keys omitted when off => byte-identical.
              pf_orphan_near_cross_cluster=false,
              pf_orphan_near_gap_cm=null,
              pf_orphan_near_min_len_cm=null,
              // Continuation terms (doc pr/128 §3.1): the touch must be
              // END-to-END and run STRAIGHT ON.  Proximity alone admitted two
              // cosmics on SBND 18255-72786 (+1151 MeV on a 701 MeV
              // candidate).  null => C++ 10cm / 30deg.
              pf_orphan_near_end_tol_cm=null,
              pf_orphan_near_kink_deg=null,
              // pf_conn4_near_candidate (doc pr/128): stop skipping a conn-4
              // shower whose closest approach to the main cluster is within
              // pf_conn4_near_gap_cm (null => C++ 20cm) -- the candidate's own
              // material that pr/74 conn3_unreachable or the pr/123+pr/124
              // prune re-seeds stamped conn-4 (SBND 18255-105074, 377 MeV).
              // Far conn-4 (490 of 514 showers, >=50cm) stays skipped.
              // C++ default false; keys omitted when off => byte-identical.
              pf_conn4_near_candidate=false,
              pf_conn4_near_gap_cm=null,
              // pf_track_owns_loose_vertex (doc pr/93 round 4): a vertex the
              // track BFS reached via a real segment is not claimable by a
              // root shower whose only tie to it is the loose fill_sets view
              // (SBND 18264-69314's 67 MeV daughter stolen from its muon).
              // C++ default false; key omitted when off => byte-identical.
              pf_track_owns_loose_vertex=false,
              // restore_demoted_mains (doc pr/20 Part I P2; C++ default false,
              // key omitted when null => byte-identical pre-knob config): tag a
              // split-off part that was ITSELF a matched bundle main before the
              // flash-time merge with flag_demoted_main.  Requires the Q/L stage
              // to have run with save_bundle_main_provenance, else the visitor
              // warns and flags nothing.  Only the OUTER (flash) un-merge takes
              // it -- unmerge_assoc undoes the isolated grouping, a different
              // question.
              restore_demoted_mains=null,
              // require_provenance (doc pr/23 sec 4.2; C++ default false, key
              // omitted when null => byte-identical pre-knob config): with
              // restore_demoted_mains on, ABORT on a pctree with no wasmain
              // array (a legacy pre-pr/20 tree) instead of warn-and-skip.
              // Legacy-tree runs opt out with false explicitly.
              require_provenance=null,
              // evaluate_demoted_mains (doc pr/20 Part I P3; C++ default false,
              // key omitted when false => byte-identical pre-knob config): let
              // TaggerCheckTGM / STM / FC evaluate a cluster carrying
              // flag_demoted_main.  Inert unless restore_demoted_mains above put
              // the flag there.
              evaluate_demoted_mains=false,
              // tgm_exempt_demoted_main (doc pr/25, SBND evt 320029; C++ default
              // false, key omitted when false => byte-identical pre-knob
              // config): with tgm_main_pair on, TaggerCheckTGM's
              // main_pair_rejects vetoes EVERY demoted-main pair
              // unconditionally (its own real_cluster_main/path-component
              // provenance is all-zero by construction post-carve) -- this
              // exempts flag_demoted_main clusters from that veto so the
              // usual CASE-A/CASE-B boundary geometry actually runs on them.
              // DESIGNED, NOT the SBND default -- changes cosmic verdicts,
              // owner sign-off pending (escalation rule 1).
              tgm_exempt_demoted_main=false,
              // skip_cosmic_companions / cosmic_companion_min_length (doc pr/20
              // Part I P4; C++ defaults false / 0, keys omitted when off =>
              // byte-identical): act on the verdict P3 produced -- drop a TGM/
              // STM-tagged companion from the neutrino's other_clusters, unless
              // it is shorter than the floor (cm).  Inert without P3.
              // nu_fallback_demoted_mains (docs/73 sec 12, round 3): when NO
              // candidate survives the primary loop, consider demoted mains
              // (same gates).  Inert without restore_demoted_mains upstream;
              // pairs with evaluate_demoted_mains (P3) so the candidates carry
              // tagger verdicts.  SBND PRODUCTION ON since 2026-08-17
              // (docs/73 sec 12 owner flip).
              // sp_photon_flag: store the single-photon tagger's verdict in
              // TaggerInfo::photon_flag, as prototype NeutrinoID.cxx:271 does.
              // The port ran singlephoton_tagger() and filled its shw_sp_*
              // features but discarded the verdict (docs/pr/26 sec. 8.2).
              // C++ default false = that gap; key omitted when off =>
              // byte-identical.  Only the uBooNE tagger ntuple's photon_flag
              // branch changes when on -- nothing in the chain reads it.
              // mip_dqdx: SBND MIP dQ/dx scale in e/cm handed to
              // TaggerCheckSTM AND (since docs pr/7-pr/8) to
              // tagger_check_neutrino as the PR chain's flat-template/cal_4mom
              // amplitude (C++ default 50000 = MicroBooNE).  56000 is the
              // SBND value derived in sbnd_xin/docs/48; pass 50000 to isolate
              // the reference-table change from the MIP-scale change in an A/B.
              mip_dqdx=56000,
              // stm_consistent_fv: TaggerCheckSTM's containment gate uses the
              // same fiducial + margins as tagger_check_tgm / tagger_check_fc.
              // Default TRUE = SBND production; pass false to restore the
              // pre-doc-49 FiducialUtils fallback for an A/B.  Details at the
              // tagger_check_stm call below.
              stm_consistent_fv=true,
              // stm_accept_guards: the doc-63 round-1 STM acceptance guards
              // (charge-desert one-objectness veto, spike-not-ramp nu-vertex
              // veto, eval ratio2 normalization cap).  C++ default false.
              // DEFAULT TRUE = SBND production as of doc 63 (owner
              // 2026-07-26); pass false for a pre-campaign A/B (key then
              // omitted => byte-identical pre-knob config).
              stm_accept_guards=true,
              // stm_proton_muon_guard: the doc-63 round-2 muon-consistency
              // guard on detect_proton (C++ default false; key omitted when
              // off => byte-identical).  DEFAULT TRUE = SBND production as
              // of doc 63 (owner 2026-07-26).
              stm_proton_muon_guard=true,
              // stm_cathode_guard: the doc-63 round-3 cathode-truncation veto
              // (C++ default false; key omitted when off => byte-identical).
              // DEFAULT TRUE = SBND production as of doc 63 (owner
              // 2026-07-26).
              stm_cathode_guard=true,
              // stm_anode_dist_fix: the doc-63 round-4a fix of the inverted
              // face selection in check_stm_conditions' dist_to_anode helper
              // (the anode-clipped-TGM catch fired at the SBND cathode).
              // C++ default false; key omitted when off => byte-identical.
              // DEFAULT TRUE = SBND production as of doc 63 (owner
              // 2026-07-26, after validation).
              stm_anode_dist_fix=true,
              // stm_second_track_guard: the doc-63 round-4b/4c second-track
              // vetoes (long straight leftover past the kink; long MIP-like
              // other-track segment).  C++ default false; key omitted when
              // off => byte-identical.  DEFAULT TRUE = SBND production as of
              // doc 63 (owner 2026-07-26, after validation).
              stm_second_track_guard=true,
              // stm_deficit_guard / stm_vertex_kink_guard: the doc-63
              // round-5 stop-region vetoes (charge-deficient end with no
              // Bragg = truncation; sharp end turn into a hot prong =
              // vertex).  C++ defaults false; keys omitted when off =>
              // byte-identical.  DEFAULT TRUE = SBND production as of doc 63
              // (owner 2026-07-26, after validation).
              stm_deficit_guard=true,
              stm_vertex_kink_guard=true,
              // stm_descent_guard: doc-94 round-1 veto on a stop reached
              // travelling UPWARD or near-horizontally.  A cosmic stopping
              // muon arrived from the sky, so it entered a boundary face
              // ABOVE the point where it stopped.  DEFAULT FALSE -- the cut
              // is set from the measured population, and stm_descent_cos_y
              // +1.01 is above the feature's range, so turning the boolean on
              // WITHOUT lowering the cut is a pure probe that changes no
              // verdict.  C++ defaults false/+1.01/10; keys omitted when off
              // => byte-identical.  See sbnd_xin/docs/94.
              stm_descent_guard=false,
              stm_descent_cos_y=1.01,
              stm_descent_min_cm=10.0,
              // stm_vertex_hadron_guard: doc-94 veto on a LONG, HEAVILY
              // IONIZING prong off the fitted main -- a proton from a neutrino
              // vertex, which check_other_tracks (a second-MUON predicate)
              // lets through because protons scatter and every clause there
              // wants straightness > 0.975 or MIP-band charge.
              // DEFAULT TRUE = SBND production as of doc 94 (owner
              // 2026-09-02, after the validation below).  NOT byte-identical.
              //   - recovers 3 of the 4 events the owner adjudicated as
              //     neutrinos: 966-2-22, 304-6-28, 146-60-31 all flip
              //     STM -> nu-candidate.  (707-18-12 was adjudicated a
              //     genuine STM on 2026-09-02 then re-adjudicated a NEUTRINO
              //     the same day, the owner having re-read the truth; it is
              //     the one of the five that no guard recovers -- docs/94.)
              //   - measured A/B over all 3067 SBND data events: 1 bundle
              //     flips of 34,827, and the owner adjudicated that bundle
              //     (64475:23) a NEUTRINO, so the measured cost is zero;
              //   - breaks 0 of the 36 owner-adjudicated correct STMs
              //     (doc 62, scan-d59k/stm-baseline.tsv);
              //   - causal negative control: len bar raised to 30 cm => 0
              //     fires and every verdict returns to baseline.
              // C++ defaults stay false/12/1.5, so every other detector is
              // untouched; set false for the A/B.
              stm_vertex_hadron_guard=true,
              stm_hadron_len_cm=12.0,
              stm_hadron_mip=1.5,
              // stm_entry_rise_guard: doc-94 round-2 veto on a CONTIGUOUS
              // elevated dQ/dx run ANCHORED at the boundary end of the fit
              // that then DECAYS to the body level -- two particles sharing
              // that stretch, one of which left the detector, so the boundary
              // point is an EXIT and the vertex sits a few tens of cm inside.
              // Every other predicate in TaggerCheckSTM reads the STOP end;
              // this is the only one that reads where the fit starts.
              // Round 3 (owner hand-scan, 2026-09-02) adds stm_entry_kink_deg:
              // the fitted path must also TURN somewhere along the muon.  Two
              // particles meeting at a vertex turn; a single muon carrying a
              // delta-ray fluctuation at the boundary does not.  The AND is
              // what works -- the kink is worthless alone (104 of the 246
              // evaluated STM bundles, 42.3%, clear 22 deg).
              // stm_entry_max_cm 30 -> 60 in the same round: the owner
              // adjudicated 350099 (shoulder 48.8 cm) a NEUTRINO, so the
              // "no decay above 30 cm" premise was wrong, and the kink now
              // carries what that bound was there for.
              //
              // DEFAULT TRUE = SBND production as of doc 94 round 3 (owner
              // 2026-09-02, after hand-adjudicating every bundle it moves).
              // NOT byte-identical.  Validation:
              //   - recovers 827-27-4, the one of the owner's five that the
              //     round-1 vertex_hadron_guard could not: STM -> nu-candidate;
              //   - all 3067 SBND data events: 34,825 bundles identical, 2
              //     flipped, 0 one-arm-only, 0 bundles gain a cosmic tag, 414
              //     of 416 STM bundles keep their tag, and BOTH flips
              //     (164466:7, 350099:15) are bundles the owner adjudicated
              //     NEUTRINOS -- measured contamination ZERO;
              //   - right on 10 of the 11 hand-adjudicated bundles; the one
              //     error is a MISS (707-18-12, docs/94 sec 14), not a false
              //     release;
              //   - probe arm with the cut disabled reproduces the baseline on
              //     3067 of 3067 events;
              //   - causal negative control: frac 1.3 -> 5.0 so no elevated
              //     run can form => 0 fires, every verdict back to baseline.
              // C++ defaults stay false/1.3/5/60/70/22, so every other
              // detector is untouched; set false for the A/B.
              stm_entry_rise_guard=true,
              stm_entry_frac=1.3,
              stm_entry_min_cm=5.0,
              stm_entry_max_cm=60.0,
              stm_entry_min_len_cm=70.0,
              stm_entry_kink_deg=22.0,
              // stm_d66_cuts: the doc-66 sec 12 diffusion-margin cut package
              // (Michel-veto res_length floor 6 -> 6.5 cm, detect_proton
              // track_medium gate 1.0 -> 1.05, block-B ks2 entry 0.05 ->
              // 0.055, C1 peak clause 4.3 -> 4.1).  The 4.0/8.8 diffusion
              // revert moved four owner-adjudicated bundles across the
              // prototype constants by hairline margins; these values were
              // chosen from a full-1000-event margin sweep with zero measured
              // collateral (sbnd_xin/docs/66 sec 12.2).  C++ defaults are the
              // prototype constants; false omits the keys => byte-identical
              // pre-package config.  DEFAULT TRUE = SBND production as of doc
              // 66 (owner 2026-07-27).
              stm_d66_cuts=true,
              stm_michel_res_cm=6.5,
              stm_proton_tm_max=1.05,
              stm_proton_b_ks2_max=0.055,
              stm_proton_c_peak_max=4.1,
              // beam_window_only: run the PR tail (steiner + TGM + STM + FC)
              // ONLY on the beam-coincident bundle -- the main clusters whose
              // matched flash time (cluster_t0) falls in beam_window, plus (for
              // steiner) the companions sharing their matched_flash_gid.
              // DEFAULT TRUE = SBND production as of doc 56: the taggers are
              // validated cosmic/containment verdicts for the beam bundle, and
              // the out-of-time bundles they used to be run on are cosmics by
              // construction (~11 bundles/event, 1 in window).  Pass false to
              // restore the evaluate-every-bundle behavior (C++ knobs default
              // off, keys omitted => compiled config byte-identical to the
              // pre-doc-56 one).  Inert if beam_window is empty.
              beam_window_only=true,
              // nu_skip_cosmic: tagger_check_neutrino skips in-window mains
              // already tagged cosmic upstream (flag_TGM, flag_STM, or Q/L
              // lm_flag>0), so neutrino PR runs only on untagged nu candidates.
              // Per-main: a cosmic-tagged longest bundle does not veto a clean
              // runner-up.  C++ default false (key omitted when off =>
              // byte-identical uBooNE config).  DEFAULT TRUE for SBND (owner
              // 2026-07-30, sbnd_xin/docs/pr/3); zero production impact while
              // tagger_check_neutrino is not in pipeline_names.
              // nu_skip_cosmic_bundle: lift that verdict from the main to the
              // whole flash bundle -- if ANY in-window main sharing a
              // matched_flash_gid is cosmic-tagged, no main of that bundle is
              // eligible, so neutrino PR does not run on the bundle at all.
              // Needed because a bundle can hold a second main: SBND evt
              // 18255/52195 gid 6 holds the TGM-tagged 513 cm cathode-crosser
              // AND a 5-point 1.7 cm shard, and under the per-main rule the
              // shard became the nu candidate and pulled the bundle's untagged
              // 400 cm associated muon through a full PR + track fit.
              // C++ default false (key omitted when off => byte-identical
              // uBooNE config).  DEFAULT TRUE for SBND (owner 2026-08-01,
              // sbnd_xin/docs/pr/3 sec. 8).  NOT bit-identical: it removes PR
              // output on cosmic bundles that previously produced it.
              // nu_skip_cosmic_bundle_min_length (cm): design-A guard on the
              // bundle veto (docs/pr/16 sec. 7).  An UNTAGGED in-window main at
              // least this long survives the veto: the TGM/STM taggers examined
              // it and declined to tag it, and the cosmic-tagged bundle-mate
              // stays out of the PR ensemble regardless (companions are
              // associated-only).  Unexamined shards (out-of-scope mains carry
              // no verdict; evt 52195's 1.3 cm shard) stay vetoed.  C++ default
              // 0 = veto every bundle-mate (key omitted => byte-identical).
              // DEFAULT 15 for SBND (owner 2026-08-01): restores evt 10550's
              // 18.5 cm candidate, keeps 13 cm and below vetoed.
              nu_skip_cosmic_bundle_min_length=15,
              // dir_weak_use_score: route the PR chain's direction-weakness
              // reads through segment_is_dir_weak() -- the faithful port of the
              // prototype's ONLY public accessor ProtoSegment::is_dir_weak()
              // (score>0.07/0.13 thresholds; sentinel score 100 = weak).  The
              // original port read the raw flag at all ~83 sites (docs/pr/6).
              // C++ default false (key omitted => byte-identical uBooNE
              // config).  DEFAULT TRUE for SBND (owner 2026-07-30); zero
              // production impact while tagger_check_neutrino is not in
              // pipeline_names.
              dir_weak_use_score=true,
              // mip_dqdx_median: the PR chain's median-dQ/dx threshold scale in
              // e/cm (C++ default 43000 = uBooNE; every x1.2/1.3/1.4/1.75...
              // ratio cut is quoted against it -- docs pr/7 sec 5, pr/8).  The
              // 50k-role flat-template amplitude reuses the existing mip_dqdx
              // arg above (56000).  48000 keeps the uBooNE 43/50 internal ratio
              // (owner 2026-07-30, placeholder pending an SBND median-MIP
              // measurement).  null falls back to the uBooNE default.
              mip_dqdx_median=48000,
              // Proton-template direction vote (doc pr/8): when the muon-vs-flat
              // direction gate abstains both ways, a proton-consistent, clearly
              // asymmetric template fit toward a free end declares the direction
              // (guards G1-G4).  C++ default false (key omitted => byte-identical
              // uBooNE config).  DEFAULT TRUE for SBND (owner 2026-07-30);
              // score/asym thresholds stay at the C++ defaults 0.25/1.3 pending
              // the pr/8 sec 6 calibration.
              proton_dir_vote=true,
              // Endpoint-trim retry (doc pr/9 sec 6 F1): on double abstention,
              // retry the dQ/dx PID once with exactly 1 sample excluded at the
              // hypothesized stopping end -- an ill-defined endpoint dilutes or
              // inflates the tip dQ/dx, which is compared against the template's
              // Bragg maximum (evt 172230: tip -37% vetoed a clean proton;
              // trimmed retry passes at score 0.077).  Dynamic: trims 0 samples
              // when the untrimmed comparison already decides, never more than
              // 1.  Runs before the proton_dir_vote fallback.  C++ default
              // false (key omitted => byte-identical uBooNE config).  DEFAULT
              // TRUE for SBND (owner 2026-07-30).
              endpoint_trim_retry=true,
              // fit_vertex short-segment exclusion (doc pr/9 sec 11 F3c), cm.
              // Segments with wcpt-path length below the cut do not count as
              // vertex-fit legs (they stay in the graph/particle flow): a
              // 0.62 cm drift-blur vertex-activity stub dragged the evt 172230
              // main vertex 2.4 cm via the post-stub refit.  >=3 surviving legs
              // => fit on the survivors; <=2 => skip the fit (the two-leg
              // position was already fit).  null = C++ default 0 = legacy
              // include-all (key omitted => byte-identical uBooNE config).
              // DEFAULT 1.0 cm for SBND (owner 2026-07-30, deliberate
              // prototype divergence).
              fit_vertex_min_seg_length=1.0,
              // cathode_x / cathode_kink_xcut (cm): the cathode kink veto,
              // doc pr/20 Part II B0.  segment_search_kink's four accept
              // criteria are each guarded by para_angle > 7.5-15 deg, a guard
              // against isochronous imaging artifacts that is INVERTED at the
              // cathode -- the crossing is drift-x dominated (para_angle
              // 61-78 deg, wide open) while the ~2 cm transverse mismatch
              // across the gap supplies the turn, so one cosmic leaves the PR
              // graph as two segments.  cathode_kink_xcut suppresses kink
              // ACCEPT candidates within that distance of x = cathode_x; the
              // angle arithmetic itself is untouched.
              // null = C++ default 0 = OFF (key omitted => byte-identical).
              cathode_x=0,
              // cathode_wide_kink_angle (deg) / _skirt / _baseline (cm): the
              // wide-baseline cathode kink ACCEPT, doc pr/47 sec 8 (O1) --
              // the converse of the veto above.  At a cathode-crossing fit
              // index the shipped index-windowed refl_angle statistic is
              // suppressed by the gap/distortion (52085: a genuine 33-38 deg
              // two-prong junction measures only ~23 deg locally), so a fifth
              // accept path fires when the skirt-excluded PCA turn angle
              // across the crossing exceeds the angle cut.  null = C++
              // default 0 = OFF (key omitted => byte-identical).
              // shower_topo_demote_len (cm, doc pr/25 sec 3): demote any
              // kShowerTopology segment longer than this to a track.  Owner
              // hand-scan 2026-08-03: 10/10 long shower-topology segments on
              // a selected nu-candidate main cluster are tracks, none showers.
              // null = C++ default 0 = OFF (key omitted => byte-identical).
              // ---- doc sbnd_xin/docs/pr/30 §11 port-fidelity knobs ------------
              // All five default to the pre-pr/30 behaviour, so the compiled JSON is
              // byte-identical until one is set.  See cfg/pgrapher/common/clus.jsonnet
              // for what each does.  fit_exclusion (P1) and oov_prototype_parity (F2)
              // turn NEW behaviour on (default off); first_seg_local_pca (P2) and
              // other_seg_relaxed_accept (P4) gate behaviour that is ALREADY production,
              // so null = on = legacy and false restores the prototype's narrower
              // version.  doc 77 round 1 (2026-08-24): graph_endpoint_strict (P8)
              // removed -- "must stay OFF" (pr/30 P8; pr/86:450).
              // shower_topo_proto_dir (doc pr/31 sec 11, F2 was P2): true skips the
              // stage-3 segment_determine_shower_direction call, leaving the topology
              // shower with the direction segment_is_shower_topology set -- the
              // prototype's state.  C++ default false = today's path = byte-identical.
              // doc pr/32 sec 11: the four stage-4 vertex-ID port fixes.  C++
              // default false = today's path; keys omitted => byte-identical.
              // The SBND operating point lives in wct-pr-perevt.jsonnet.
              // doc pr/31 sec 12: the sec 10.12 topology/PID/direction port
              // fixes (F5 cont-muon 30cm dir3, F6 empty-window abstain, F3
              // shower-topo reset, F1 preserve-4mom, F4 local median, F7
              // vertex-by-index -- F7 deliberately dormant pending pr/30 F4).
              // C++ default false = today's path; keys omitted =>
              // byte-identical.  The SBND operating point lives in
              // wct-pr-perevt.jsonnet.
              // Steiner TERMINAL filter fidelity (doc pr/29 D1 and D12).  Both
              // OFF here = the historical toolkit behaviour, keys omitted =>
              // byte-identical config.  Turning either on can only ADD Steiner
              // terminals, never remove them, because both restore a way for a
              // terminal to PASS the reference-containment test:
              //   steiner_terminal_wire_tol=1     the prototype's one wire of
              //     slack on both sides of all three planes.  Applies ONLY to
              //     the terminal filter; get_extreme_wcps shares the C++ helper
              //     and correctly keeps the exact test the prototype gives it.
              //   steiner_terminal_adjacent_slice=true   makes the t+-1 slice
              //     fallback actually resolve.  The time_blob_map key is in
              //     TICKS, so the historical step of 1 matches nothing on SBND
              //     (4 ticks per slice) and the whole branch is dead.
              // NOT flipped on by default: that is an unconditional production
              // output change (CLAUDE.md sec.5 rule 1) and is the owner's call.
              steiner_terminal_wire_tol=0,
              steiner_terminal_adjacent_slice=false,
              // Steiner EDGE-WEIGHT charge fidelity (doc pr/29 D2).  OFF here =
              // the historical toolkit behaviour, key omitted => byte-identical.
              //   steiner_edge_charge_forward_dead_mix=true   weights steiner
              //     edges with charges computed under the disable_dead_mix_cell
              //     this chain actually passes (false), as the prototype does;
              //     the toolkit dropped the argument at the call and always
              //     took the true branch of calc_charge_wcp.
              // Unlike the two above this can move a weight in EITHER
              // direction, so it can add or drop tree edges.
              steiner_edge_charge_forward_dead_mix=false,
              // Isochronous first-segment endpoint finding (doc pr/24 round 2,
              // SBND evt 271851): principal-axis endpoints for filled 2-D
              // sheet clusters instead of the wire-footprint boundary metric.
              // false/nulls = C++ defaults (false / 40 cm / 25 cm / 0.35 /
              // 0.02) = OFF (keys omitted => byte-identical).
              // examine_vertices_3 extension-retraction guard (doc pr/24
              // round 5).  false/null = C++ defaults (false / -1.0 cm) = OFF
              // (keys omitted => byte-identical).
              // doc pr/67: log-only trajectory-coverage probe + the
              // counterfactual override for find_proto_vertex's hardcoded
              // main-cluster branch-search round budget.  false/null =
              // C++ defaults = OFF (keys omitted => byte-identical).
              // protect_bundle stage knobs (doc pr/23): the PR-stage
              // overclustering protection (uboone's second graph examination,
              // ProtectOverClustering.cxx).  The stage only acts when
              // 'protect_bundle' is named in pipeline_names -- which the
              // production jobs do by DEFAULT since the doc pr/23 sec 9 flip
              // (owner 2026-08-02).  protect_graph_name: connected_blobs flavor
              // (null => 'relaxed').  protect_cathode_*: the SBND cathode
              // re-join divergence -- INTERNAL units (unlike cathode_kink_xcut
              // one block up, which is cm); nulls => C++ defaults, i.e. the
              // re-join pass disabled (prototype-faithful).
              protect_graph_name=null,
              // C++ default true (doc pr/23 ordering): convicted in-window
              // mains (TGM/STM/lm) do not open their bundle.  null => key
              // omitted => C++ default.
              protect_skip_convicted=null,
              // doc pr/94 round 3: let a convicted main open its bundle so the
              // bundle's unconvicted members are graph-examined.  C++ default
              // false; key omitted when null.
              protect_open_convicted_bundles=null,
              protect_cathode_x=0,
              protect_cathode_rejoin_xcut=5 * wc.cm,
              protect_cathode_rejoin_dyz=4 * wc.cm,
              protect_cathode_rejoin_dis=8 * wc.cm,
              // Direction-agreement fallback for a pair that fails ONLY
              // cathode_rejoin_dyz (doc pr/25, SBND evt 489327): dyz is a
              // frame-aligned transverse-offset bound and rejects an OBLIQUE
              // crosser by construction (the offset across the cathode gap is
              // real track travel, not misalignment).  nulls => C++ default
              // 0 = fallback disabled => byte-identical to the block above.
              // INTERNAL units except protect_cathode_rejoin_angle (degrees).
              protect_cathode_rejoin_perp=null,
              protect_cathode_rejoin_angle=null,
              protect_cathode_rejoin_dir_radius=null,
              protect_cathode_rejoin_dir_npts=null,
              // cosmic_y_* (cm): cosmic_tagger()'s four "reaches the top of the
              // detector" tests, re-anchored from uBooNE's top face (y = +117 cm)
              // to SBND's (y = +200 cm, sbnd_y_top above).  The uBooNE literals
              // 100 / 102 / 80 / 50 cm are 17 / 15 / 37 / 67 cm below its top --
              // an entry-point tolerance for a downward cosmic, which does not
              // scale with detector height -- so SBND keeps the same offsets and
              // gets 183 / 185 / 163 / 133 cm.  At the uBooNE values these cuts
              // are near-meaningless on SBND (100 cm is mid-detector), which is
              // what docs/pr/2 sec 2e(iv) flagged.  DEFAULT ON for SBND (owner
              // 2026-07-30); null on any one of them restores that cut's uBooNE
              // value.  NOT bit-identical: cosmic_tagger verdicts can change.
              // vertex_z_prior_scale (cm): denominator of the upstream-z penalty
              // (z - min_z)/scale that ranks main-vertex candidates against the
              // +0.25-per-track bonuses.  uBooNE's 200 cm spans ~5.2 penalty
              // units over its 1037 cm detector; SBND is 501 cm long, so the
              // same 200 cm halves the prior's dynamic range -- most visible in
              // compare_main_vertices_global, which ranks candidates from
              // DIFFERENT clusters of the beam bundle, where separations do run
              // toward the full detector length.  DEFAULT = the length-scaled
              // 200 x 501/1037 = 96.6 -> 100 cm (owner 2026-07-30).  The
              // alternative reading -- that the prior is a per-cm trade-off
              // against track bonuses and transfers unchanged -- keeps 200;
              // pass null for that (docs/pr/2 sec 2e(iv)).
              // NOT bit-identical: vertex ranking can change.
              // ssm_target_dir / ssm_absorber_dir: the SSM tagger's beam-line
              // reference directions [x,y,z] in the detector frame (docs/pr/2
              // sec 2e(i)).  null = the C++ defaults = the prototype's uBooNE
              // BNB-target / NuMI-absorber vectors, so the compiled config is
              // unchanged.  SBND HAS NO VALUE FOR EITHER YET: the BNB-target
              // direction must be re-derived in the SBND frame, and the
              // NuMI-absorber features have no obvious SBND meaning.  Until
              // then the 8 ssm_*_angle_{target,absorber} features carry uBooNE
              // geometry -- they are only reachable now, not fixed.
              // kine_*: the charge -> kinetic-energy calibration constants of
              // NeutrinoEnergyReco (docs/pr/2 sec 2e(iii)).  null = the C++
              // defaults = the uBooNE-tuned literals they replaced.
              // The three RECOMBINATION survival factors carry SBND values
              // since 2026-07-30 (docs/pr/10 sec 6): the uBooNE empirical
              // 0.7 / 0.5 / 0.35 scaled by the table-integrated survival
              // ratio between the official uBooNE Box at 0.273 kV/cm and the
              // doc-55 free-power SBND fit (C excluded -- degenerate with
              // the fudge factors, which deliberately STAY uBooNE: they
              // absorb gain/lifetime normalization, not field physics).
              // NOT bit-identical: kine-derived BDT features move; the
              // T_kine tree of the pr/3-style tracking-pr.root shows the
              // effect (Enu -12..-14% on nuecc48 172230/235435/444187).
              // null on any one restores that factor's uBooNE value.
              // The plane weights [0.25,0.25,1.0], the 0.04 asymmetry
              // switch and kine_w_value (23.6 eV) still have NO SBND value.
              // doc pr/35 sec 10.2 (F1 = P1+P8): live start-segment PDG at the
              // four fill_kine_tree sites (prototype parity).  C++ default
              // false; key omitted when off => byte-identical pre-knob config.
              // ---- doc sbnd_xin/docs/pr/36 sec 10 tagger-stage knobs -------
              // F1 (= P1): give the match_isFC recompute the SAME fiducial +
              // margins tagger_check_{stm,tgm,fc} use (pdvd_pr_fv +
              // pdvd_pr_fv_margins), mirroring stm_consistent_fv above.
              // false => keys omitted => the historical FiducialUtils
              // fallback, byte-identical pre-knob config.
              neutrino_consistent_fv=false,
              // sbnd_xin/docs/74 G1/G2: give cosmic_tagger()'s containment
              // tests (its inside_fv lambda + the flag-1 vertex test) the
              // SAME pdvd_pr_fv + pdvd_pr_fv_margins as TGM/STM/FC, instead
              // of the FiducialUtils zero-margin sensitive union (which has
              // no wall inset and excludes the CPA slab |x| < 0.45 cm).
              // C++ default false; key omitted when off => byte-identical.
              cosmic_consistent_fv=false,
              // sbnd_xin/docs/75: give the nue/single-photon taggers'
              // containment tests (angular_cut, shower_to_wall,
              // bad_reconstruction_2/_2_sp) the SAME pdvd_pr_fv +
              // pdvd_pr_fv_margins as cosmic_consistent_fv above -- same
              // zero-margin FiducialUtils inconsistency, same fix pattern.
              // C++ default false; key omitted when off => byte-identical.
              nue_sp_consistent_fv=false,
              // F3 (= P2): single-photon SCE correction gate.  Vacuous on
              // SBND today (clus_geom_helper is ''); kept OFF by owner
              // decision 2026-08-04 so a future SBND SCE helper enables it
              // as its own explicit step.  C++ default false.
              // F4 (= P3+P5): graph-index-ordered tagger accumulation sets
              // (M4 house-rule determinism fix).  C++ default false.
              // F5 (= P6): prototype wcpt-identity stem-endpoint rule at the
              // 18 seg_endpoint_near sites.  C++ default false.
              // F6 (= P8): broken_muon_id counts distinct cluster ids.  C++
              // default false.
              // F7 (= P4): neutrino_type verdict bitmask + its T_tagger
              // branch (threaded to BOTH tagger_check_neutrino and
              // tagger_output).  C++ default false.
              neutrino_type_bitmask=false,
              // doc pr/94 Phase 1: per-bundle identity + per-activity
              // cosmic-flag T_tagger branches (tagger_output only, plumbing
              // only -- nothing populates them yet).  C++ default false.
              nu_per_bundle=false,
              nu_per_bundle_min_length=null,  // doc pr/94 Phase 5b; cm; null => C++ default 0 (no floor)
              // doc 80: MCS muon momentum.  ONE argument derives both the
              // TaggerCheckNeutrino computation key (mcs_enable, via the knob
              // bag) and the UbooneTaggerOutputVisitor T_kine branch key
              // (mcs_output), so the two gates can never disagree.  C++
              // default false on both; keys omitted when off =>
              // byte-identical pre-MCS config AND schema.
              mcs_enable=false,
              // doc pr/94 round 3: the selected candidate gets the main-cluster
              // PR treatment for its own pass even when it is a demoted main.
              // C++ default false; key omitted when off.  NOT gated on
              // nu_per_bundle -- the legacy demoted-main fallback needs it too.
              // doc 75: closes the DL-swap flag leak nu_selected_as_main's
              // own guard leaves open (see common/clus.jsonnet comment).
              // C++ default false; key omitted when off.
              // ---- doc sbnd_xin/docs/pr/33 §10 EM-shower-clustering knobs.
              // All C++ default false = keys omitted = byte-identical
              // pre-knob config.
              // F1 (= P1): prototype calculate_num_daughter_tracks callee at
              // the main-vertex proton-skip site / the examine_showers
              // daughter_length site.
              // F2 (= P2): read the PDG off the object the prototype reads
              // (4 sites start-segment; 1 inverted site shower-type; 2 sites
              // exact ==13 muon test).  Parity at the :170 site needs
              // from_start_segment AND exact_muon_test together.
              // F3 (= P3): shared pi0-id allocation stream across the two
              // pi0 finders (prevents pio_id collision in the nue tagger
              // pi0 block and the Bee mc.json grouping).
              // F4 (= P6): is_shower gains the prototype's abs(pdg)==11
              // disjunct at the cluster-center-point site.
              // F5 (= P12): shower_less same-index tie-break by stable
              // shower id (house-rule determinism fix, prototype n/a).
              // doc pr/39: exclude a shower's own start vertex from the
              // end_point farthest-vertex search (prototype map_vtx_segs
              // parity).  Ships OFF pending owner gate review.
              // doc pr/91 round 3: flood-fill frontier = visited, not merely
              // present in the view.  Ships OFF pending owner gate review.
              // doc sbnd_xin/docs/pr/40 -- track (proton/pion/muon)
              // mis-identified as electron.  F1 (persistence), F2/F3 (dQ/dx
              // guards on wholesale track-to-electron conversion).  All
              // default false = legacy = byte-identical.
              // doc sbnd_xin/docs/pr/40 round 2 -- two follow-on defects from the
              // pr/40 round: F1 zero-KE persistence stub, F2 proton-
              // daughter -> pion, F3 reclass_pinfo negative-KE stub.  All
              // default false = legacy = byte-identical.
              // doc sbnd_xin/docs/pr/40 round 4 -- two follow-on defects from
              // round 2/3's F5: F7 clears the shower flags a relabelled pion
              // still carried (it was still being wrapped as a Shower); F8
              // relabels a muon segment at a multi-proton (>=2, charge-
              // confirmed) non-neutrino-vertex hadronic vertex to pion.  Both
              // default false = legacy = byte-identical.
              // doc sbnd_xin/docs/pr/40 round 5 -- muon mis-identified as
              // electron, three independent mechanisms.  F9 narrows F1 so it
              // no longer rescues an undirected electron guess; F10 excludes
              // a long, straight candidate from the main-vertex EM-shower
              // selection; F11 gives segment_is_shower_trajectory the same
              // straightness veto F3 gave segment_is_shower_topology's
              // dQ/dx.  All three default false = legacy = byte-identical.
              // doc sbnd_xin/docs/pr/40 round 6 -- boundary-level fixes the
              // round-5 measurement demanded.  F12 keeps the shower flood-
              // fill from absorbing a confident straight non-electron track;
              // F13 keeps connecting_to_main_vertex from force-setting a
              // proton-daughter pion to electron; F14 widens the Michel
              // stopping-muon rescue past its weak-dir degree-2 limits.
              // All three default false = legacy = byte-identical.
              // doc 77 round 1 (2026-08-24): shower_connect_protected_pion_
              // guard removed -- measured dead, never flipped (pr/40:1459).
              // doc sbnd_xin/docs/pr/74 round 2 -- P1 cascade guard + P2
              // Michel-terminal check.  C++ defaults false/40cm/1.3/40cm.
              // Keys omitted when off/null => byte-identical pre-pr/74.
              // doc sbnd_xin/docs/pr/40 round 9 -- straight-track PID guard
              // family + B2 cross-cluster bridge.  C++ defaults false
              // (scalars 25 deg / 1.8 cm live in C++).  Keys omitted when
              // off/null => byte-identical pre-round-9 config.
              // doc pr/97 D1; C++ default false => legacy indeterminate main_pi read.
              // doc pr/92 -- stray-satellite drop from kine/PF.  false/null =
              // C++ defaults = OFF (20 MeV / 8 cm / 60 deg / 45 deg / 90 cm /
              // 30 cm / 25 deg); keys omitted => byte-identical pre-pr/92.
              // doc pr/84 round 2 (F3): stitch radius in cm; null = C++
              // default 0 = OFF, key omitted => byte-identical.
              // doc pr/84 round 3 (S1): collapse showers that share a start
              // segment.  C++ default false = OFF, key omitted when off =>
              // byte-identical pre-fix config.
              // doc pr/44: a multi-segment long-muon pseudo-shower keeps its
              // muon start segment (prototype parity; the update_particle_type
              // majority vote at the in_main_cluster seeding site is a
              // toolkit-only addition).  C++ default false; key omitted when
              // off => byte-identical pre-fix config.
              // doc pr/40 round 10; false = C++ default = OFF, key omitted =>
              // byte-identical.
              // doc pr/93 round 3 -- C++ defaults false; key-suppressed when off.
              // doc pr/93 round 4 -- C++ defaults false / null = C++ defaults
              // (orphan floor 50cm; sccc 5cm/15deg base + 12cm/7.5deg aligned).
              // Keys suppressed when off => byte-identical.
              // doc pr/99 round 3; C++ defaults; keys suppressed when off.
              // doc pr/101 Enu accounting round; C++ defaults; keys suppressed when off.
              // doc pr/43 round 2 -- C++ defaults false; keys suppressed when off.
              // doc pr/45 -- find_other_segments empty-2D-tree sentinel guard
              // (SBND 18255-56463 isochronous tail).  C++ default false; key
              // suppressed when off => byte-identical.
              // doc pr/46 -- long-muon stub bridge in find_cont_muon_segment
              // (18255-55595 broken muon behind a 2.4 cm vertex stub).  C++
              // default false; key suppressed when off => byte-identical.
              // doc pr/48 -- back-to-back track fixes
              // (18255-51513/56211/57903/59335/57485: nu vertex mid-segment
              // on one unbroken track, dQ/dx rising at BOTH ends, no angular
              // kink).  two_end_break = the two-end residual-range break
              // pass; kink_walk_dqdx_stop / kink_break_protect = the 59335
              // walk-overshoot + EV4-absorption fixes.  C++ defaults false;
              // keys suppressed when off => byte-identical.  teb_* numerics:
              // null = C++ defaults (doc pr/48 sec 9 operating point), inert
              // while two_end_break is off.
              teb_min_len=null,
              teb_min_arm=null,
              teb_min_arm_pts=null,
              teb_stub_max=null,
              teb_accept_range=null,
              teb_rise_r1=null,
              teb_rise_r2=null,
              teb_abs_end_min=null,
              teb_dip_floor=null,
              teb_score_cap_r1=null,
              teb_score_cap_r2=null,
              teb_turn_angle=null,
              teb_turn_baseline=null,
              teb_turn_skirt=null,
              // doc pr/90 round 2: R2 argmax arm-fill guard.  null = C++
              // default 0 = legacy, key suppressed.  doc 77 round 1
              // (2026-08-24): teb_second_max removed -- negative on its own
              // motivating events (pr/90 sec 8.5).
              // doc pr/90 round 4: chain-topology admission (D1), route R3
              // turn/activity (D3), R2 bragg veto (D4).  false/null = C++
              // default = legacy, key suppressed.
              kink_dqdx_hot_ratio=null,
              // doc pr/49 -- cross-cluster projection-ghost deweighting in
              // the trajectory fit's 2D charge association (18255-57441
              // V-plane ghost: an unrelated cluster's charge aliases with
              // the fitted track in one view and detours the fit; live cells
              // outside the fitted cluster's own blob coverage that sit
              // inside a 3D-distant foreign cluster's keep their measurement
              // at reduced weight -- cells covered by nobody, or claimed by
              // a genuinely touching neighbor, keep full weight, so
              // no-ghost events are untouched).  C++ default -1 = off; null
              // omits the key => byte-identical.  >= 0 = on, value =
              // wire/slice tolerance in cells (0 = strict; the 57441
              // contamination is ONE cell away, so >= 1 re-admits it).
              // doc pr/50 -- suspend the pr/49 deweighting during
              // find_proto_vertex (the partition-forming stage; its recursive
              // kink walk is globally sensitive to fit perturbations --
              // 18255-172230 lost its true-kink main vertex to a 2.7 cm
              // neighbor).  All later fitting stages keep the deweighting.
              // C++ default false = pr/49 behavior; false omits the key =>
              // byte-identical.
              // doc pr/50 -- main-vertex kink-consistency snap (172230-class
              // near-vertex robustness; C++ defaults in TaggerCheckNeutrino.h).
              // false/null omit the keys => byte-identical.
              // doc pr/85 -- carry the old vertex's arms through the snap
              // residual below this arc (cm).  C++ default 0 = off; null
              // omits the key => byte-identical.
              // doc pr/104 -- main-vertex junction snap (C++ defaults in
              // TaggerCheckNeutrino.h).  false/null omit the keys => byte-identical.
              // doc pr/51 -- main-vertex graph audit (near-vertex graph-shape
              // repair: dup-corridor merge / charge-less-bridge removal /
              // micro-stub absorb + re-seat / one refit; C++ defaults in
              // TaggerCheckNeutrino.h).  false/null omit the keys =>
              // byte-identical.
              // esva_ignore_empty_2d (docs/73 sec 12, round 3, evt 78242):
              // eliminate_short_vertex_activities case 5 must not read the
              // empty-2D-index sentinel (-1) as "covered" on cathode-crossing
              // clusters.  SBND PRODUCTION ON since 2026-08-17 (docs/73
              // sec 12 owner flip).
              // doc pr/51 round 3 -- op3 satellite-anchor radius (cm).
              // C++ default 0 (main-vertex-only, round 2); null omits the
              // key => byte-identical.
              // doc pr/85 -- op3 interposed-stub absorb at the main-vertex
              // anchor (angle in deg, C++ default 150).  false/null omit
              // the keys => byte-identical.
              // doc pr/86: interposed-splice candidate ceiling, cm (C++
              // default 0 = use mvga_stub).  null omits the key =>
              // byte-identical.
              // doc pr/86 P4: satellite-anchor op3 overlap threshold (C++
              // default 0 = use mvga_dup_frac).  null omits the key =>
              // byte-identical.
              // doc pr/86 P1b: interposed splice at degree-1 main anchors
              // (C++ default false).  false omits the key => byte-identical.
              // doc pr/86 round 2: op3 post-carry straighten reach (cm) and
              // op3.5 junction-collapse radius (cm).  C++ defaults 0 (off).
              // Keys omitted when null => byte-identical.
              // doc pr/83 r3: op1 scope/threshold decouple (cm / fraction;
              // radius -1 = unscoped), post-op3 dup pass, carry cap,
              // abandoned-cluster dup audit.  C++ defaults 0/0/false/0/false
              // (all legacy).  Keys omitted when null/false =>
              // byte-identical.
              // doc 77 round 1 (2026-08-24): mvga_carry_max removed -- not
              // needed, class A cleared 8/8 with it OFF (pr/83 r3 sec 8.5).
              // doc pr/83 r4 -- projective dup collapse; null omits => byte-identical.
              // doc 77 round 1 (2026-08-24): dl_vtx_swap_guard removed --
              // live A/B -36/1014 (pr/89 round 5).
              // doc pr/106 sec 10: exclusion-free charge cloud for the DL vertex net. C++ default false.
              // doc pr/112 sec 11: dual chain. C++ defaults false/"snap"/false/2.0/true/0; null/false => key omitted.
              // doc pr/107: dQ/dx fit keeps every trajectory point (prototype parity). C++ default false.
              // doc 77 round 1 (2026-08-24): dl_vtx_topo_weight/_center
              // (pr/89 Arm C2 rule-1 topology term) removed -- live A/B
              // -8/1014.
              // doc pr/51 round 3 -- apply the traditional-path swap
              // decision instead of discarding it.  C++ default false;
              // false omits the key => byte-identical.
              // doc pr/51 round 4 -- diagnostic-only rough-path probe.
              // C++ default false; false omits the key => byte-identical.
              // doc pr/51 round 5 -- steiner gap penalty (H1 short-cut fix):
              // do_rough_path routes on the support-penalized
              // "steiner_graph_gap" flavor when scale > 0.  C++ defaults:
              // scale 0 (off), dead_alpha 0.25, min_edge 0.5 cm,
              // sample_step 0.3 cm, point_radius 0.2 cm.  null omits the
              // keys => byte-identical.
              // doc pr/73: per-edge DEBUG sentinel for the steiner_graph_gap
              // scan (endpoints, midpoint, w, bad, both vertex charges,
              // deficit).  Log-only diagnostic.  C++ default false; false
              // omits the key => byte-identical compiled config.
              // doc pr/51 round 6 -- weak-charge deficit term on the same
              // gap flavor.  C++ defaults: weak_scale 0 (off), weak_qref
              // 2000 (charge units).  null omits the keys => byte-identical.
              // doc pr/73 round 2 F3a -- do_rough_path route excursion cap,
              // cm.  C++ default -1 = off; null omits the key.
              // doc pr/83 -- oriented break_segment splits (find_vertices, not
              // boost source/target).  C++ default false; key omitted when
              // off => byte-identical pre-fix config.
              // doc pr/51 round 7 -- robust vertex fit (dynamic per-leg
              // direction windows for MyFCN).  C++ defaults: robust false,
              // main_only true, min_len 10, rin_margin 2, rout_frac 0.5,
              // rout_min 9, rout_max 18, angle 20, min_pts 5, min_aniso 3,
              // prior_range 1 (lengths cm, angle deg).  false/null omit the
              // keys => byte-identical.
              // doc pr/54 -- keep well-supported isolated residual segments
              // in find_other_segments (18255-142421 separated EM shower with
              // no fitted trajectory).  C++ defaults: keep false, floors
              // 25 points / 3 cm.  false/null omit the keys => byte-identical.
              // doc pr/102 P1 -- OR-disjuncts on the pr/54 keep: min_nnf
              // (terminal not-faked floor) and len_admit (cm).  C++ defaults
              // 0 = off.  null omits the keys => byte-identical.
              // doc 77 round 1 (2026-08-24): other_seg_uncover_3d (pr/102
              // P2, 3-D uncovered-charge radius) removed -- 23 ADVERSE
              // movers, stays OFF.
              // doc pr/67 round 3 (S2) -- isochronous-snap size gate, cm.
              // C++ default 10.0 = legacy.  null omits the key.
              // doc pr/65 round 3 -- offer graph-unreachable main-cluster
              // segments (kept-isolated pr/54 residuals) to the shower
              // absorbers (reachability-relaxed guards).  C++ default false;
              // false omits the key => byte-identical.
              // doc pr/59 round 2 -- per-cluster orphaned-associate_points
              // rescue.  C++ default false; false omits the key =>
              // byte-identical.
              // doc pr/64 round 7 -- reassign same-cluster association
              // orphans that Stage C of clustering_points_segments would
              // otherwise drop (18259-18625: 12-18 pt blob at PF segment
              // 126042's own fit endpoint, in img charge but absent from
              // shower_track/associate_points).  C++ default false; false
              // omits the key => byte-identical.
              // doc pr/64 round 8 -- clear a merge survivor's
              // associate_points when examine_structure_final_1/_1p/_3
              // deletes a segment that had non-empty associate_points, so
              // pr/59's reassociate_cluster_orphans any_orphan trigger
              // correctly re-fires.  C++ default false; false omits the key
              // => byte-identical.
              // doc pr/72 round 2 -- guard examine_structure_3 against
              // merging a genuine near-vertex track stub into an unrelated
              // shower/track trunk (18255-196649).  C++ default false;
              // false/null omits the key => byte-identical.  Numeric
              // defaults null = component keeps its own C++ default
              // (fitted from a 117-event census).
              // doc pr/45 -- paint muon-typed (+-13) pseudo-showers as track in
              // the Bee shower_track layer + PrDisplayDump (18255-56463: 411 cm
              // muon painted red).  C++ default false; key suppressed when off
              // => byte-identical.
              pseudo_shower_track_paint=false,
              // muon_dqdx_curve [c0, c1, pivot_cm, power]: the muon
              // median-dQ/dx-vs-length envelope used by nine tagger cuts, as a
              // multiple of mip_dqdx_median.  DEFAULT = the SBND fit
              // (2026-07-30, docs/pr/10 sec 4): the SBND stopping-muon table
              // median (0.5 kV/cm, /48000) times the uBooNE empirical/table
              // margin g(L), refit with the prototype's functional form.
              // NOT bit-identical: nine tagger cuts move (looser at short L
              // -- the higher field keeps more of the Bragg rise in dQ/dx).
              // Scales as 1/mip_dqdx_median: re-derive when the 48000
              // placeholder becomes a measurement.  null restores the uBooNE
              // refit [0.8866, 0.9533, 18, 0.4234] (byte-identical pre-knob).
              // use_power_recomb: hand the taggers (STM + neutrino PR) the
              // free-power Modified-Box recombination fitted to SBND stopping
              // tracks (docs/55 sec 7g canonical, PowerBoxRecombination
              // defaults) instead of the plain pdvd_box_recomb above.
              // DEFAULT ON for SBND (owner 2026-07-30, docs/pr/10 sec 7):
              // every model-driven dQ/dx -> dE/dx conversion now uses the
              // measured SBND recombination.  false restores pdvd_box_recomb,
              // which is byte-identical to the pre-pr/10 config AND (since
              // doc 88 re-typed it to PracticalBoxRecombination) computes the
              // same numbers it did then -- upstream e6fb7ef3 changed what
              // Gen::BoxRecombination does with practical-unit parameters.
              use_power_recomb=true,
              // sp_dedx_use_recomb_model: route the single-photon stem dE/dx
              // through the configured recombination model instead of the
              // inline uBooNE-field (0.273 kV/cm) inverse Box.  DEFAULT ON
              // for SBND (owner 2026-07-30) with sp_mean_dedx_cut=2.23, the
              // physical-scale transfer of the legacy 2.3 (docs/pr/10 sec 5).
              // false/null restore the inline formula and the 2.3 literal.
              // dl_vtx_cut: max distance (mm; C++ default 25.0 = 2.5 cm) from
              // the DL SCN prediction to accept a candidate vertex.  Threaded
              // for configurability (docs/pr/2 sec 7.4); null keeps the C++
              // default, which is coupled to the uBooNE-trained net (gap G3).
              // fast_xgb_forest: book the two XGB BDT combiners with
              // TmvaGradForest instead of TMVA::Reader -- same scores, ~4 s
              // and ~0.9 GB less per PR job (sbnd_xin/docs/76 round 2).
              // C++ default false.  Key omitted when off => byte-identical
              // pre-knob config.
              fast_xgb_forest=false,
       // doc 77 round 2: the pattern-recognition knob bag, built once by the
       // job (wct-pr-perevt.jsonnet) from its TLAs and handed to
       // TaggerCheckNeutrino as-is.  An absent key is that knob's C++ default,
       // exactly as the per-knob `+ (if x then {x: ...})` clauses it replaces.
       tcn_knobs={}):: {
        // Only gate when the caller actually supplied a window; beam_window=[0,0]
        // (the arg default, i.e. "no beam window") must not silently drop every
        // cluster's tagger evaluation.
        local beam_gate = beam_window_only && beam_window[1] > beam_window[0],
        local dv = clus_maker.detector_volumes(anodes),
        local pcts = clus_maker.pc_transforms(dv),
        // DetectorVolumes implements IFiducial -- used by MakeFiducialUtils / the
        // taggers' inside_fiducial_volume().  NOT the FV_* metadata below:
        // DetectorVolumes::contained() is contained_by(p).valid(), i.e. the union of
        // the per-face IAnodeFace::sensitive() boxes, which AnodePlane builds as
        // x in [anode_x, cathode_x] and so run to the W plane with no margin
        // (SBND |x| <= 201.45, |y| <= 199.965, z in [0, 501.0]; |x| < 0.45 is a hole).
        // Any tagger given no explicit fiducial therefore tests a volume ~3 cm more
        // permissive at every wall than pdvd_pr_fv + pdvd_pr_fv_margins below --
        // TaggerCheckSTM and TaggerCheckNeutrino's match_isFC still do; see
        // sbnd_xin/docs/49_stm-containment-fv-inconsistency.md.
        local cm_old = clus.clustering_methods(
            prefix='pr', detector_volumes=dv, pc_transforms=pcts, fiducial=dv, coords=clus_maker.scope_coords),
        local cm = clus.clustering_methods(
            prefix='pr', detector_volumes=dv, pc_transforms=pcts, fiducial=dv, coords=clus_maker.t0cor_coords),
        // Box-model recombination at the SBND drift field (uBooNE used 0.273 kV/cm).
        // doc 88: PracticalBoxRecombination, NOT BoxRecombination -- the
        // numbers below are in practical units (kV/cm, (kV/cm)(g/cm^2)/MeV,
        // g/cm^3).  Upstream e6fb7ef3 made Gen::BoxRecombination consistently
        // WCT-unit, so handing it these numbers is wrong by a factor
        // units::cm/units::MeV = 10 in the quenching term.  The Practical
        // class carries the pre-e6fb7ef3 arithmetic verbatim, so this path is
        // bit-identical to what it computed before that merge.
        // PDVD (doc pdvd/25 sec 7c): Modified Box with the SAME parameter set the
        // dQ/dx reference tables were generated with (convert_field.C: alpha 0.93,
        // beta 0.212, rho 1.38) at the PDVD field 0.44 kV/cm implied by the
        // data-calibrated drift velocity (pdvd/stm/pdvd_transport.py).  SBND's
        // A=1.0/B=0.255 is a different parameter set (energy_loss/docs/
        // dqdx_consistency_check.md sec 6); PDVD keeps tables and model consistent.
        local pdvd_box_recomb = {
            type: 'PracticalBoxRecombination',
            name: 'pdvd_box_recomb',
            data: { A: 0.93, B: 0.212, Efield: 0.44, rho: 1.38, Wi: 23.6e-6 },
        },
        // Free-power Modified Box fitted to SBND stopping-track dQ/dx vs residual
        // range (docs/55 sec 7g; canonical parameters in stm_ref_dqdx.json:
        // R = ln(A+u)/u, u = k*(dEdx/2.1)^p, times normalization C).  The data
        // block restates the C++ defaults so the operating point is visible here;
        // selected via use_power_recomb (docs/pr/10).
        // PDVD has NO fitted power-box yet (doc pdvd/25 M7: fit it on the PDVD
        // stopping-muon sample).  The block restates the SBND fit / C++ defaults so
        // use_power_recomb=true is well-defined; PDVD production keeps it false.
        local pdvd_power_recomb = {
            type: 'PowerBoxRecombination',
            name: 'pdvd_power_recomb',
            data: { A: 0.93, k: 0.282371, p: 1.362179, C: 0.855175, pivot: 2.1, Wi: 23.6e-6, dedx_max: 77.0 },
        },
        // The recombination model the taggers actually receive.  With
        // use_power_recomb=false this compiles to exactly the pre-knob config.
        local pdvd_recomb = if use_power_recomb then pdvd_power_recomb else pdvd_box_recomb,
        // TGM fiducial: ONE box spanning BOTH TPCs (the overall FV bounds of
        // dvm above), so a cathode-crossing track is not an "exiter" at x=0.
        // The default fiducial=dv cannot serve here: DetectorVolumes::contained()
        // is the union of per-face sensitive volumes, which excludes the CPA slab
        // (|x| < 0.45 cm).  Margins go in via the tagger's fv_tolerance instead
        // of the box, mirroring the metadata *_margin values.
        // PDVD: ONE box spanning BOTH drift volumes (bottom anodes 0-3 at x<0, top
        // anodes 4-7 at x>0, the 6 cm cathode slab centred on x=0), so a cathode
        // crosser is not an "exiter" at x=0.  x = the shield-plane (anode face)
        // position +-339.91 cm where the sensitive volume starts (protodunevd
        // clus.jsonnet dvm a0f0pA/a4f0pA FV_x = +-3399.1 mm); y/z = the active
        // volume (+-336.4 cm, 0.05..299.25 cm: dvm overall FV_y/z are these inset
        // by 15 cm).  Margins enter through fv_tolerance below, whose defaults
        // (2 / 2.5 / 3 cm) equal dvm.overall's FV_*_margin values.
        local pdvd_pr_fv = {
            type: 'BoxFiducial',
            name: 'pdvd_pr_fv',
            data: {
                bounds: {
                    tail: { x: -339.91 * wc.cm, y: -336.4 * wc.cm, z: 0.05 * wc.cm },
                    head: { x: 339.91 * wc.cm, y: 336.4 * wc.cm, z: 299.25 * wc.cm },
                },
            },
        },
        // Margin vector convention (inside_fv / Clustering_Util fc_check): index 4
        // is applied as contained(z - tv[4]), so with a NEGATIVE entry it insets
        // the DOWNSTREAM (z ~ 500 cm) face; index 5 insets the upstream face.
        // tgm_fv_zmax_margin (cm; default 3 = byte-identical legacy value)
        // parametrizes only the downstream inset -- shared by tagger_check_tgm AND
        // tagger_check_fc below, so "contained" keeps one meaning across both
        // verdicts (docs/27_fc-tgm-consistent-fv.md).
        // tgm_fv_x_margin / tgm_fv_y_margin (cm; defaults 2 / 2.5 = byte-identical
        // legacy values) parametrize the drift-x and vertical-y insets, both faces
        // symmetric.  Shared by tagger_check_tgm AND tagger_check_fc, same as the
        // downstream-z knob, so "contained" keeps one meaning across both verdicts.
        local pdvd_pr_fv_margins = [-tgm_fv_x_margin * wc.cm, -tgm_fv_x_margin * wc.cm, -tgm_fv_y_margin * wc.cm, -tgm_fv_y_margin * wc.cm, -tgm_fv_zmax_margin * wc.cm, -3 * wc.cm],
        // tgm_fv_zmax_margin_interior (cm; default 0 = OFF, key omitted =>
        // byte-identical): when > 0, check_tgm's CASE-A interior-support tests
        // (chord midpoints + waypoint re-check) use THIS downstream-z inset
        // instead of tgm_fv_zmax_margin, i.e. the doc-32 widening becomes
        // endpoint-only.  Rationale (doc 32 caveat, doc 35): a corner clipper
        // running ALONG the downstream wall inside the widened 3->5 cm band
        // (evt287517 cluster 16, evt289805 cluster 9) keeps its midpoint support
        // at the legacy 3 cm interior.  TGM only -- tagger_check_fc containment
        // and the ENDPOINT exit tests keep pdvd_pr_fv_margins unchanged.
        // x/y track the endpoint vector above: the doc-35 endpoint-only widening
        // applies to the downstream-z inset only.
        local pdvd_pr_fv_margins_interior = [-tgm_fv_x_margin * wc.cm, -tgm_fv_x_margin * wc.cm, -tgm_fv_y_margin * wc.cm, -tgm_fv_y_margin * wc.cm, -tgm_fv_zmax_margin_interior * wc.cm, -3 * wc.cm],
        // Retiler for the steiner stage: same 'stepped' samplers that built the 3d
        // PC (PointTreeBuilding), one per (anode, face) -- PDVD anodes are two-sided,
        // so 16 samplers, each with its crate's drift speed (clus.jsonnet live_sampler).
        local improve2 = cm.improve_cluster_2(
            anodes=anodes,
            samplers=[clus.sampler(clus_maker.live_sampler(a, f), apa=a.data.ident, face=f)
                      for a in anodes for f in [0, 1]]),
        // Visitors available to the PR pipeline, by name.  switch_scope re-applies
        // the per-cluster T0 correction on the loaded tree (the corrected scope is
        // runtime state and does not persist through the tarball); it recomputes
        // deterministically from cluster_t0.  fiducialutils MUST precede any
        // tagger (they silently no-op without it).
        local cm_by_name = {
            switch_scope: cm_old.switch_scope(),
            // PDVD (doc 25 M3): every flash-matched cluster becomes a main so the
            // taggers and the per-bundle neutrino PR evaluate it (QLMatching only
            // flags MicroBooNE-style decomposed mains, which PDVD never builds).
            // Placed right after switch_scope.  Absent from pipeline_names =>
            // absent from the compiled config.
            flag_mains: {
                type: 'ClusteringFlagMatchedMains',
                name: 'pr',
                data: {
                    grouping: 'live',
                    require_t0: true,
                    min_length: flag_mains_min_length,
                    skip_flagged: true,
                },
            },
            // PDVD has NO unmerge_bundle / unmerge_assoc stages (doc pdvd/25 sec 4):
            // examine_bundles is not run by the PDVD Q/L chain and its isolated()
            // pass does not write assoc_cluster_id, so there is nothing to undo.
            // PR-stage overclustering protection (doc pr/23): uboone's SECOND
            // graph-examination round (WCPPID::Protect_Over_Clustering ->
            // Examine_graph, run after Q/L matching and before NeutrinoID) --
            // split each beam-bundle cluster at graph-component boundaries so
            // vertex-proximate photon fragments merged by Clustering_neutrino are
            // fit as separate clusters instead of one MST-bridged trajectory
            // (doc pr/22 sec 8, evt 386948).  Position (doc pr/23 ordering
            // decision): AFTER the cosmic taggers (tagger_check_tgm/stm/fc) and
            // BEFORE tagger_check_neutrino, with 'steiner' named a second time
            // right after it so the split clusters' steiner products are rebuilt
            // (the visitor drops the pre-split ones; CreateSteinerGraph
            // replace=false touches nothing else).  This matches uboone -- cosmic
            // verdicts on unsplit clusters (wire-cell-prod-stm.cxx:806-810),
            // Protect_Over_Clustering only in the nue executable
            // (wire-cell-prod-nue.cxx:1322).  The cathode re-join
            // knobs keep it from splitting cathode crossers, an SBND geometry
            // uboone did not have (doc pr/20).  In the production pipeline_names
            // by DEFAULT since doc pr/23 sec 9; dropped from the list => absent
            // from the compiled config => the pre-pr/23 chain.
            protect_bundle: cm.protect_bundle(
                graph_name=if protect_graph_name == null then 'relaxed' else protect_graph_name,
                beam_window_only=beam_gate,
                beam_window_low=beam_window[0],
                beam_window_high=beam_window[1],
                skip_convicted=protect_skip_convicted,
                open_convicted_bundles=protect_open_convicted_bundles,
                cathode_x=protect_cathode_x,
                cathode_rejoin_xcut=protect_cathode_rejoin_xcut,
                cathode_rejoin_dyz=protect_cathode_rejoin_dyz,
                cathode_rejoin_dis=protect_cathode_rejoin_dis,
                cathode_rejoin_perp=protect_cathode_rejoin_perp,
                cathode_rejoin_angle=protect_cathode_rejoin_angle,
                cathode_rejoin_dir_radius=protect_cathode_rejoin_dir_radius,
                cathode_rejoin_dir_npts=protect_cathode_rejoin_dir_npts),
            // SBND has no beam_flash flag (QLMatching sets main/associated_cluster
            // instead) -- process every scope-passing cluster, narrowed to the
            // beam-coincident bundle when beam_window_only is on (the default; see
            // the pr() argument).  Keys omitted when off => byte-identical.
            steiner: cm.steiner(retiler=improve2, perf=true, require_beam_flash=false,
                                beam_window_only=beam_gate,
                                beam_window_low=beam_window[0],
                                beam_window_high=beam_window[1],
                                terminal_wire_tol=steiner_terminal_wire_tol,
                                terminal_adjacent_slice=steiner_terminal_adjacent_slice,
                                edge_charge_forward_dead_mix=steiner_edge_charge_forward_dead_mix)
              + { data+: { [if steiner_terminal_charge != null then 'terminal_charge_threshold']: steiner_terminal_charge } },
            // The doc pr/23 second steiner pass, named right after protect_bundle:
            // replace=false rebuilds ONLY the clusters protect_bundle purged
            // (split retained + fragments).  A replace=true second pass would
            // erase+emplace every in-window cluster's steiner_graph and dangle
            // the GraphAlgorithms the STM fit cached in the tagger stage
            // (bad_alloc from garbage num_vertices, SBND evt 54095).  Distinct
            // component name => distinct instance; in the production
            // pipeline_names by DEFAULT since doc pr/23 sec 9, always paired
            // with protect_bundle.
            steiner_refresh: cm.steiner(name='refresh',
                                retiler=improve2, perf=true, require_beam_flash=false,
                                beam_window_only=beam_gate,
                                beam_window_low=beam_window[0],
                                beam_window_high=beam_window[1],
                                // Same terminal-filter settings as the first pass:
                                // the refresh rebuilds only the clusters
                                // protect_bundle purged, and they must be built the
                                // same way as their peers (doc pr/29 D1, D12).
                                terminal_wire_tol=steiner_terminal_wire_tol,
                                terminal_adjacent_slice=steiner_terminal_adjacent_slice,
                                edge_charge_forward_dead_mix=steiner_edge_charge_forward_dead_mix,
                                replace=false)
              + { data+: { [if steiner_terminal_charge != null then 'terminal_charge_threshold']: steiner_terminal_charge } },
            fiducialutils: cm.fiducialutils(),
            tagger_check_stm: cm.tagger_check_stm(
                evaluate_demoted_mains=evaluate_demoted_mains,
                trackfitting_config_file=trackfitting_config_file,
                particle_dataset=wc.tn(particle_dataset),
                recombination_model=wc.tn(pdvd_recomb),
                require_in_scope=true,
                // save_stm_fit (C++ default false; key omitted when off =>
                // byte-identical): persist per-pass STM fits as cluster PCs +
                // grouping slot "stm" for the Bee stm_fit layer and the
                // stm_magnify ROOT dump below.  Runner flag: -stm-fit.
                save_stm_fit=save_stm_fit,
                // mip_dqdx: SBND MIP dQ/dx scale in e/cm, replacing the inherited
                // MicroBooNE 50000.  Anchored to the muon reference table the same
                // way 50000 was anchored to MicroBooNE's: the uBooNE table plateau
                // (rr = 59.5 cm) is 48879.4 e/cm and 50000/48879.4 = 1.02293; the
                // SBND table regenerated at 0.5 kV/cm plateaus at 54657.7 e/cm and
                // 56000/54657.7 = 1.02456, so the same 2.3-2.5% headroom is kept
                // and 56000 is as round a number as 50000 was.
                // NOT byte-identical -- see sbnd_xin/docs/48.
                mip_dqdx=mip_dqdx,
                // accept_guards (C++ default false; key omitted when off =>
                // byte-identical): doc-63 round-1 acceptance guards.
                accept_guards=stm_accept_guards,
                // proton_muon_guard (C++ default false; key omitted when off =>
                // byte-identical): doc-63 round-2 muon-consistency guard.
                proton_muon_guard=stm_proton_muon_guard,
                // cathode_guard (C++ default false; key omitted when off =>
                // byte-identical): doc-63 round-3 cathode-truncation veto.
                cathode_guard=stm_cathode_guard,
                // anode_dist_fix (C++ default false; key omitted when off =>
                // byte-identical): doc-63 round-4a dist_to_anode face fix.
                anode_dist_fix=stm_anode_dist_fix,
                // second_track_guard (C++ default false; key omitted when off =>
                // byte-identical): doc-63 round-4b/4c second-track vetoes.
                second_track_guard=stm_second_track_guard,
                // deficit_guard / vertex_kink_guard (C++ defaults false; keys
                // omitted when off => byte-identical): doc-63 round-5 vetoes.
                deficit_guard=stm_deficit_guard,
                vertex_kink_guard=stm_vertex_kink_guard,
                // descent_guard (C++ default false; keys omitted when off =>
                // byte-identical): doc-94 round-1 travel-direction veto.
                descent_guard=stm_descent_guard,
                guard_descent_cos_y=(if stm_descent_guard then stm_descent_cos_y else null),
                guard_descent_min_cm=(if stm_descent_guard then stm_descent_min_cm else null),
                // vertex_hadron_guard (C++ default false; keys omitted when off
                // => byte-identical): doc-94 round-1 vertex-hadron veto.
                vertex_hadron_guard=stm_vertex_hadron_guard,
                guard_hadron_len_cm=(if stm_vertex_hadron_guard then stm_hadron_len_cm else null),
                guard_hadron_mip=(if stm_vertex_hadron_guard then stm_hadron_mip else null),
                // entry_rise_guard (C++ default false; keys omitted when off
                // => byte-identical): doc-94 round-2 entry-rise veto.
                entry_rise_guard=stm_entry_rise_guard,
                guard_entry_frac=(if stm_entry_rise_guard then stm_entry_frac else null),
                guard_entry_min_cm=(if stm_entry_rise_guard then stm_entry_min_cm else null),
                guard_entry_max_cm=(if stm_entry_rise_guard then stm_entry_max_cm else null),
                guard_entry_min_len_cm=(if stm_entry_rise_guard then stm_entry_min_len_cm else null),
                guard_entry_kink_deg=(if stm_entry_rise_guard then stm_entry_kink_deg else null),
                // doc-66 sec 12 cut package (C++ defaults = prototype constants;
                // keys omitted when stm_d66_cuts=false => byte-identical).
                michel_res_length_cut=(if stm_d66_cuts then stm_michel_res_cm * wc.cm else null),
                proton_tm_max=(if stm_d66_cuts then stm_proton_tm_max else null),
                proton_b_ks2_max=(if stm_d66_cuts then stm_proton_b_ks2_max else null),
                proton_c_peak_max=(if stm_d66_cuts then stm_proton_c_peak_max else null),
                // stm_consistent_fv (default TRUE = SBND production): give the
                // containment gate the SAME fiducial + margins tagger_check_tgm and
                // tagger_check_fc get, so "contained" means one thing across all
                // three verdicts.  Without it cluster_fc_check falls back to
                // FiducialUtils -> the union of per-face SENSITIVE volumes, which
                // exceeds this box at every wall even before the margins (SBND 0.40
                // cm x / 0.65 y / 0.85 z) and is holed at the CPA slab: on a
                // 30-event sample 96 of the 147 clusters the STM tagger skipped as
                // "fully contained" were called exiters by tagger_check_fc, so no
                // dQ/dx fit was ever attempted for them.  The prototype has no such
                // split -- check_stm and check_tgm are members of ONE ToyFiducial
                // and share its boundaries, inset once by boundary_dis_cut = 3 cm --
                // so this restores parity rather than inventing a volume.
                // NOT byte-identical.  Set false for the A/B (runner: -no-stm-fv).
                // See sbnd_xin/docs/49_stm-containment-fv-inconsistency.md.
                // Only the ENDPOINT tests exist in cluster_fc_check, so the
                // interior_fv_tolerance vector tagger_check_tgm needs has no
                // counterpart here.
                fiducial=(if stm_consistent_fv then wc.tn(pdvd_pr_fv) else null),
                fv_tolerance=(if stm_consistent_fv then pdvd_pr_fv_margins else []),
                // Beam-window gate on the main-cluster loop; see the pr()
                // beam_window_only argument.  Keys omitted when off =>
                // byte-identical.
                beam_window_only=beam_gate,
                beam_window_low=beam_window[0],
                beam_window_high=beam_window[1])
              // PDVD boundary vetoes (pr() args above); the common helper
              // returns a plain component object, so the keys are merged into
              // its data (no pnode wrapper here).  Omitted when off.
              + { data+: {
                    [if stm_readout_edge_guard then 'readout_edge_guard']: true,
                    [if stm_readout_edge_guard then 'guard_readout_edge_ticks']: stm_readout_edge_ticks,
                    [if stm_readout_edge_guard then 'readout_nticks']: nticks,
                    [if stm_cathode_guard_cm != null then 'guard_cathode_cm']: stm_cathode_guard_cm,
                  } },
            // STM-stage Magnify-tracking ROOT dump (doc sbnd_xin/docs/40): reads
            // the stm_fit/stm_pass cluster PCs and the "stm" TrackFitting slot,
            // writes tracking-stm.root (T_rec_charge/T_proj_data/T_bad_ch/Trun)
            // with the two-TPC concatenated-per-plane channel convention.  Only
            // active when named in pipeline_names (needs -stm-fit; the WireCellRoot
            // plugin must be loaded by the job).
            stm_magnify: {
                type: 'PdvdMagnifyTrackingVisitor',
                name: 'pr',
                data: {
                    grouping: 'live',
                    track_fitting_name: 'stm',
                    output_filename: evt_out_prefix + 'tracking-stm.root',
                    runNo: runNo,
                    subRunNo: subRunNo,
                    eventNo: eventNo,
                    anodes: [wc.tn(a) for a in anodes],
                    detector_volumes: wc.tn(dv),
                    dQdx_scale: 0.1,
                    dQdx_offset: -1000.0,
                    // Readout length in ticks; only used to clamp T_bad_ch time
                    // ranges (a dead region with unbounded x converts to +-1e9
                    // ticks and makes the Magnify projection pads unreadable).
                    // SBND 3427; PDVD passes its readout window (pr() nticks arg).
                    nticks: nticks,
                },
            },
            // Through-going-muon tagger (prototype check_tgm port).  Runs on every
            // matched main cluster.  tgm_neutrino_candidate (C++ default false;
            // key omitted when off => byte-identical): enable the ported
            // check_neutrino_candidate veto so in-beam-window bundles may be
            // tagged; when off in-beam bundles are never tagged (conservative v1
            // behavior).  Must run after fiducialutils (dead-region /
            // signal-processing checks) and before tagger_check_stm (which skips
            // TGM-flagged mains).
            // require_in_scope: evaluate only clusters that pass switch_scope's
            // active-volume filter.  Without it the tagger also sees the
            // out-of-volume shards switch_scope splits off (which keep an inherited
            // flag_main_cluster) -- on the SBND 10-event reco1 sample those were 43%
            // of all evaluated mains and tagged TGM at 85% vs 44% for real clusters.
            tagger_check_tgm: cm.tagger_check_tgm(
                evaluate_demoted_mains=evaluate_demoted_mains,
                fiducial=wc.tn(pdvd_pr_fv),
                fv_tolerance=pdvd_pr_fv_margins,
                // Key omitted when the knob is 0 (empty list) => byte-identical.
                interior_fv_tolerance=(if tgm_fv_zmax_margin_interior > 0
                                       then pdvd_pr_fv_margins_interior else []),
                beam_window_low=beam_window[0],
                beam_window_high=beam_window[1],
                // Beam-window gate on the main-cluster loop -- the SAME window the
                // in-beam protection above uses; see the pr() beam_window_only
                // argument.  Key omitted when off => byte-identical.
                beam_window_only=beam_gate,
                require_in_scope=true,
                check_neutrino_candidate=tgm_neutrino_candidate,
                require_chord_charge=tgm_chord_charge,
                // C++ default "chord". Key omitted then => byte-identical.
                chord_charge_mode=tgm_chord_mode,
                component_extremes=tgm_component_extremes,
                // C++ default false. Key omitted when off => byte-identical.
                component_rescue=tgm_component_rescue,
                // C++ default false. Key omitted when off => byte-identical
                // doc-32 rescue behavior.  Rescued-end pairs must also pass the
                // straight-chord test (evt288727 two-cosmic composite).
                rescue_chord_check=tgm_rescue_chord,
                // C++ default false. Key omitted when off => byte-identical.
                // A pair may tag only when an end lies in the main cluster --
                // the TGM verdict follows the main cluster, not a merged-in
                // through-going fragment (evt289343 cluster 9).
                main_component_pairs=tgm_main_pair,
                // C++ default "path" (largest-path-component proxy). Key omitted
                // then => byte-identical.  "real" = per-blob flash-merge
                // provenance (needs a save_real_cluster_id pctree; falls back to
                // the proxy on old tarballs).
                main_component_mode=tgm_main_pair_mode,
                // C++ default false. Key omitted when off => byte-identical.
                // See the pr() arg comment (doc pr/25, SBND evt 320029).
                exempt_demoted_main_pairs=tgm_exempt_demoted_main),
            // Fully-contained tagger.  Independent of TGM/STM: it evaluates every
            // in-scope main cluster and only records a containment verdict (flag
            // "FC"), so it neither vetoes nor is vetoed by them.  Placed LAST in
            // the pipeline on purpose -- cluster_fc_check lazily populates
            // PCA/hough/steiner-boundary caches, so running it after the other
            // taggers leaves their inputs pristine rather than relying on those
            // caches being order-independent.  Coverage is unaffected by position.
            // fiducial/fv_tolerance are the SAME objects tagger_check_tgm gets, so
            // "contained" means one thing across both verdicts.  Without them FC
            // fell back to FiducialUtils (per-face sensitive volumes, no margin),
            // which is both more permissive at every wall than TGM's inset box and
            // holed at the CPA slab -- see docs/27_fc-tgm-consistent-fv.md.
            tagger_check_fc: cm.tagger_check_fc(
                evaluate_demoted_mains=evaluate_demoted_mains,
                fiducial=wc.tn(pdvd_pr_fv),
                fv_tolerance=pdvd_pr_fv_margins,
                require_in_scope=true,
                // Beam-window gate on the main-cluster loop; see the pr()
                // beam_window_only argument.  Keys omitted when off =>
                // byte-identical.
                beam_window_only=beam_gate,
                beam_window_low=beam_window[0],
                beam_window_high=beam_window[1]),
            // Neutrino pattern recognition on the beam-coincident bundle.  The
            // beam_window gate (on cluster_t0 = matched flash time) replaces
            // uBooNE's single-main + beam_flash selection; companions are the
            // associated clusters sharing the main's matched_flash_gid.
            // dl_weights defaults to the uBooNE-trained SCN net = DL vertex ON
            // (docs/pr/4); '' falls back to the geometric vertex.
            tagger_check_neutrino: cm.tagger_check_neutrino(
                trackfitting_config_file=trackfitting_config_file,
                particle_dataset=wc.tn(particle_dataset),
                recombination_model=wc.tn(pdvd_recomb),
                perf=true,
                dl_weights=dl_weights,
                // The DL re-rank sub-knobs, pinned here rather than inherited: they
                // were inert while dl_weights was '' and went LIVE with the doc-pr/4
                // adoption, so SBND records the operating point it was validated at
                // (= the common/clus.jsonnet defaults as of 2026-07-30, hence the
                // compiled JSON is unchanged by this pinning).  min_accept and top_k
                // are threaded from the pr() args (defaults identical to the old
                // pinned literals => byte-identical when unset; doc pr/79).
                dl_vtx_rerank=dl_vtx_rerank,
                dl_vtx_top_k=dl_vtx_top_k,
                dl_vtx_min_accept_score=dl_vtx_min_accept_score,
                dl_vtx_score_scale=1000.0,
                beam_window_low=beam_window[0],
                beam_window_high=beam_window[1],
                nu_skip_cosmic_bundle_min_length=nu_skip_cosmic_bundle_min_length,
                // doc pr/94 Phase 2: the SAME switch that books the per-bundle
                // T_tagger/T_kine branches on tagger_output below also turns on the
                // per-bundle candidate loop here -- one row per in-beam-window
                // bundle needs both halves or the extra rows have nothing to carry.
                // nu_per_bundle_demoted_acts is wired straight from
                // evaluate_demoted_mains (not from its own TLA) so the
                // act_evaluated column can never drift from the admission gate the
                // cosmic taggers actually used.
                nu_per_bundle=nu_per_bundle,
                nu_per_bundle_demoted_acts=evaluate_demoted_mains,
                nu_per_bundle_min_length=nu_per_bundle_min_length,
                dir_weak_use_score=dir_weak_use_score,
                mip_dqdx_median=mip_dqdx_median,
                proton_dir_vote=proton_dir_vote,
                endpoint_trim_retry=endpoint_trim_retry,
                fit_vertex_min_seg_length=fit_vertex_min_seg_length,
                // doc pr/36 sec 10 (F1): same fiducial + margins as
                // tagger_check_{stm,tgm,fc} above -- one containment definition
                // across the stage.  Keys omitted when off.
                fiducial=(if neutrino_consistent_fv || cosmic_consistent_fv || nue_sp_consistent_fv then wc.tn(pdvd_pr_fv) else null),
                fv_tolerance=(if neutrino_consistent_fv || cosmic_consistent_fv || nue_sp_consistent_fv then pdvd_pr_fv_margins else []),
                // sbnd_xin/docs/74 G1/G2: cosmic_tagger() containment on the same
                // fiducial + margins.  Key omitted when off => byte-identical.
                // sbnd_xin/docs/75: nue/single-photon tagger containment on the
                // same fiducial (each site's own hardcoded tolerance -- see
                // NeutrinoTaggerNuE.cxx).  Key omitted when off => byte-identical.
                teb_min_len=teb_min_len,
                teb_min_arm=teb_min_arm,
                teb_min_arm_pts=teb_min_arm_pts,
                teb_stub_max=teb_stub_max,
                teb_accept_range=teb_accept_range,
                teb_rise_r1=teb_rise_r1,
                teb_rise_r2=teb_rise_r2,
                teb_abs_end_min=teb_abs_end_min,
                teb_dip_floor=teb_dip_floor,
                teb_score_cap_r1=teb_score_cap_r1,
                teb_score_cap_r2=teb_score_cap_r2,
                teb_turn_angle=teb_turn_angle,
                teb_turn_baseline=teb_turn_baseline,
                teb_turn_skirt=teb_turn_skirt,
                kink_dqdx_hot_ratio=kink_dqdx_hot_ratio,
              // doc 77 round 2: these five are read elsewhere in pr() too, so they stay
              // named parameters and join the knob bag here rather than at the job.
              knobs=tcn_knobs + {
                  [if cathode_x != null then 'cathode_x']: cathode_x,
                  [if cosmic_consistent_fv then 'cosmic_consistent_fv']: true,
                  [if mcs_enable then 'mcs_enable']: true,  // doc 80; sub-knobs arrive via tcn_knobs
                  [if mip_dqdx != null then 'mip_dqdx']: mip_dqdx,
                  [if neutrino_type_bitmask then 'neutrino_type_bitmask']: true,
                  [if nue_sp_consistent_fv then 'nue_sp_consistent_fv']: true,
              }),
            // NuMu / nue BDT scorers (UbooneNumuBDTScorer / UbooneNueBDTScorer,
            // geometry-free TaggerInfo consumers).  The weights are the
            // uBooNE-TRAINED XMLs from wire-cell-data uboone/weights/ -- the same
            // 35 files the qlport gate job books -- so on SBND the scores are
            // UNCALIBRATED (availability + relative ranking only; SBND retraining
            // is docs/pr/2 gap G1).  Must run after tagger_check_neutrino, nue
            // after numu.  Only compiled in when named in pipeline_names.
            local bdt_weights_dir = 'uboone/weights',
            numu_bdt_scorer: cm.numu_bdt_scorer(
                numu1_weights_xml=     bdt_weights_dir + '/numu_tagger1.weights.xml',
                numu2_weights_xml=     bdt_weights_dir + '/numu_tagger2.weights.xml',
                numu3_weights_xml=     bdt_weights_dir + '/numu_tagger3.weights.xml',
                cosmict10_weights_xml= bdt_weights_dir + '/cos_tagger_10.weights.xml',
                numu_xgboost_xml=      bdt_weights_dir + '/numu_scalars_scores_0923.xml',
                fast_xgb_forest=fast_xgb_forest),
            nue_bdt_scorer: cm.nue_bdt_scorer(
                mipid_weights_xml=       bdt_weights_dir + '/mipid_BDT.weights.xml',
                gap_weights_xml=         bdt_weights_dir + '/gap_BDT.weights.xml',
                hol_lol_weights_xml=     bdt_weights_dir + '/hol_lol_BDT.weights.xml',
                cme_anc_weights_xml=     bdt_weights_dir + '/cme_anc_BDT.weights.xml',
                mgo_mgt_weights_xml=     bdt_weights_dir + '/mgo_mgt_BDT.weights.xml',
                br1_weights_xml=         bdt_weights_dir + '/br1_BDT.weights.xml',
                br3_weights_xml=         bdt_weights_dir + '/br3_BDT.weights.xml',
                br3_3_weights_xml=       bdt_weights_dir + '/br3_3_BDT.weights.xml',
                br3_5_weights_xml=       bdt_weights_dir + '/br3_5_BDT.weights.xml',
                br3_6_weights_xml=       bdt_weights_dir + '/br3_6_BDT.weights.xml',
                stemdir_br2_weights_xml= bdt_weights_dir + '/stem_dir_br2_BDT.weights.xml',
                trimuon_weights_xml=     bdt_weights_dir + '/stl_lem_brm_BDT.weights.xml',
                br4_tro_weights_xml=     bdt_weights_dir + '/br4_tro_BDT.weights.xml',
                mipquality_weights_xml=  bdt_weights_dir + '/mipquality_BDT.weights.xml',
                pio_1_weights_xml=       bdt_weights_dir + '/pio_1_BDT.weights.xml',
                pio_2_weights_xml=       bdt_weights_dir + '/pio_2_BDT.weights.xml',
                stw_spt_weights_xml=     bdt_weights_dir + '/stw_spt_BDT.weights.xml',
                vis_1_weights_xml=       bdt_weights_dir + '/vis_1_BDT.weights.xml',
                vis_2_weights_xml=       bdt_weights_dir + '/vis_2_BDT.weights.xml',
                stw_2_weights_xml=       bdt_weights_dir + '/stw_2_BDT.weights.xml',
                stw_3_weights_xml=       bdt_weights_dir + '/stw_3_BDT.weights.xml',
                stw_4_weights_xml=       bdt_weights_dir + '/stw_4_BDT.weights.xml',
                sig_1_weights_xml=       bdt_weights_dir + '/sig_1_BDT.weights.xml',
                sig_2_weights_xml=       bdt_weights_dir + '/sig_2_BDT.weights.xml',
                lol_1_weights_xml=       bdt_weights_dir + '/lol_1_BDT.weights.xml',
                lol_2_weights_xml=       bdt_weights_dir + '/lol_2_BDT.weights.xml',
                tro_1_weights_xml=       bdt_weights_dir + '/tro_1_BDT.weights.xml',
                tro_2_weights_xml=       bdt_weights_dir + '/tro_2_BDT.weights.xml',
                tro_4_weights_xml=       bdt_weights_dir + '/tro_4_BDT.weights.xml',
                tro_5_weights_xml=       bdt_weights_dir + '/tro_5_BDT.weights.xml',
                nue_xgboost_xml=         bdt_weights_dir + '/XGB_nue_seed2_0923.xml',
                fast_xgb_forest=fast_xgb_forest),
            // PR-stage Magnify-tracking ROOT dump (docs/pr/3): fork of the uBooNE
            // writer reading the unnamed TrackFitting slot + PRGraph filled by
            // tagger_check_neutrino, with the two-TPC channel convention and
            // per-point APA from PR::Fit::paf.  Only active when named in
            // pipeline_names (the WireCellRoot plugin must be loaded by the job).
            // evt_out_prefix, not output_dir: a group running in one process
            // writes output_dir/pr_evt<ID>/tracking-pr.root, the same path the
            // per-event job writes.  Empty evt_subdir => unchanged.
            local tracking_pr_root = evt_out_prefix + 'tracking-pr.root',
            tracking_visitor: {
                type: 'PdvdPrMagnifyTrackingVisitor',
                name: 'pr',
                data: {
                    grouping: 'live',
                    output_filename: tracking_pr_root,
                    runNo: runNo,
                    subRunNo: subRunNo,
                    eventNo: eventNo,
                    anodes: [wc.tn(a) for a in anodes],
                    detector_volumes: wc.tn(dv),
                    dQdx_scale: 0.1,
                    dQdx_offset: -1000.0,
                    flag_skip_vertex: false,
                    // Readout length in ticks; only clamps T_bad_ch time ranges.
                    nticks: nticks,
                    // doc 87.  Key omitted when off => byte-identical config.
                    [if save_in_scope then 'save_in_scope']: true,
                },
            },
            // T_tagger/T_kine writer (UbooneTaggerOutputVisitor, reused as-is: it
            // is a pure TaggerInfo/KineInfo dump with no geometry).  Opens the
            // tracking file in UPDATE mode, so it must be named AFTER
            // tracking_visitor (and after both BDT scorers so the scores are set).
            // doc pr/36 sec 10.8 (F7): neutrino_type_bitmask books the
            // neutrino_type/I branch; key omitted when off => schema-identical.
            // doc pr/94 Phase 1: nu_per_bundle books the per-bundle identity +
            // per-activity cosmic-flag branches; key omitted when off =>
            // schema-identical.  Plumbing only -- nothing populates them yet.
            // doc 80: mcs_output derives from the SAME mcs_enable argument as
            // TaggerCheckNeutrino's computation key -- independent knobs would
            // allow branches full of -1 or computed-and-discarded results.
            tagger_output: cm.tagger_output(output_filename=tracking_pr_root,
                                            neutrino_type_bitmask=neutrino_type_bitmask,
                                            nu_per_bundle=nu_per_bundle,
                                            mcs_output=mcs_enable),
            // PR event-display calib dump (docs/pr/26): ONE self-contained JSON per
            // event carrying the PR-graph segments as polylines, the associated
            // track/shower points, the Steiner skeleton with its terminal flag, the
            // fitted 2-D charge per (apa,face,plane) and the dead regions -- i.e.
            // the union of what the Bee PR layers and the Magnify tracking file
            // carry, in one file the Bokeh viewer can read with the stdlib alone.
            // Read-only; mutates nothing.  Only active when named in
            // pipeline_names => compiled config byte-identical otherwise.
            // Must run AFTER tagger_check_neutrino (it reads the TrackFitting slot
            // and the PR graph that stage fills).
            pr_display: {
                type: 'PrDisplayDump',
                name: 'pr',
                data: {
                    grouping: 'live',
                    // The event id is in the NAME here, so a group run needs a
                    // conversion rather than the configured (constant) eventNo.
                    output_filename: evt_out_prefix + 'calib-pr-evt' + std.toString(eventNo) + '.json',
                    runNo: runNo,
                    subRunNo: subRunNo,
                    eventNo: eventNo,
                    anodes: [wc.tn(a) for a in anodes],
                    detector_volumes: wc.tn(dv),
                    // Same convention as the Bee track_fit layer and the Magnify
                    // writer; the dump stores RAW dQ and records these so the
                    // viewer can reproduce the Bee colouring if it wants to.
                    dQdx_scale: 0.1,
                    dQdx_offset: -1000.0,
                    nticks: nticks,
                    // doc pr/45.  C++ default false; mirrors the Bee shower_track
                    // layer knob.  Key omitted when off => byte-identical.
                    [if pseudo_shower_track_paint
                     then 'pseudo_shower_track_paint']: true,
                },
            },
        },
        local cm_pipeline = [cm_by_name[n] for n in pipeline_names],
        // The taggers' configs only name the recombination/particle-dataset
        // components; emit them (and the caller's LinterpFunctions etc. via
        // extra_uses) when a tagger is in the pipeline.
        local tagger_uses = (if std.member(pipeline_names, 'tagger_check_stm')
                             || std.member(pipeline_names, 'tagger_check_neutrino')
                             then [pdvd_recomb] + extra_uses else [])
                            + (if std.member(pipeline_names, 'tagger_check_tgm')
                               || std.member(pipeline_names, 'tagger_check_fc')
                               || (stm_consistent_fv
                                   && std.member(pipeline_names, 'tagger_check_stm'))
                               // doc pr/36 sec 10 (F1): the neutrino tagger also
                               // names pdvd_pr_fv when its consistent-FV knob is
                               // on.  Redundant in the default pipeline (TGM/FC
                               // already pull it in) => compiled config unchanged
                               // there; load-bearing only for a reduced pipeline.
                               || ((neutrino_consistent_fv || cosmic_consistent_fv || nue_sp_consistent_fv)
                                   && std.member(pipeline_names, 'tagger_check_neutrino'))
                               then [pdvd_pr_fv] else []),
        local bee_zip_path = evt_out_prefix + 'mabc-pr.zip',
        local mabc = g.pnode({
            type: 'MultiAlgBlobClustering',
            name: 'clus_pr',
            data: {
                inpath: 'pointtrees/%d',
                outpath: 'pointtrees/%d',
                perf: true,
                bee_dir: bee_dir,
                bee_zip: if pr_bee then bee_zip_path else '',  // doc 87
                bee_detector: 'protodunevd',
                initial_index: 0,
                use_config_rse: true,
                runNo: runNo,
                subRunNo: subRunNo,
                eventNo: eventNo,
                save_deadarea: true,
                dead_area_version: 2,
                save_opflash: false,
                anodes: [wc.tn(a) for a in anodes],
                detector_volumes: wc.tn(dv),
                cluster_id_order: 'tree',
                bee_points_sets: [
                    {
                        name: 'clustering',
                        detector: 'protodunevd',
                        algorithm: 'clustering',
                        pcname: '3d',
                        coords: clus_maker.t0cor_coords,
                        individual: false,
                        // Same filter semantics as clus_all_apa (default 1): the dump
                        // keys on the per-cluster runtime scope-filter flag, which is
                        // NOT persisted through the tarball -- it is re-established by
                        // running switch_scope at the head of the PR pipeline.  A
                        // pass-through job (pipeline_names=[]) therefore dumps
                        // nothing; the round-trip identity gate runs
                        // pipeline_names=['switch_scope'].
                    },
                    // Neutrino-PR layers, dumped after TaggerCheckNeutrino runs (the
                    // visitor key is the full type:name; prefix='pr').  Inert unless
                    // tagger_check_neutrino is in pipeline_names.  PRGraph fit points
                    // are already T0-corrected, hence plain x/y/z.
                    {
                        name: 'track_fit',
                        visitor: 'TaggerCheckNeutrino:pr',
                        grouping: 'live',
                        detector: 'protodunevd',
                        algorithm: 'track_fit',
                        pcname: '3d',            // not used for PRGraph dumps, but required
                        coords: ['x', 'y', 'z'],
                        individual: false,
                        dQdx_scale: 0.1,
                        dQdx_offset: -1000.0,
                        // C++ default false.  Append the PR-graph vertex fit points
                        // (real_cluster_id=-1) like the prototype's
                        // fill_skeleton_info_magnify rows (docs/pr/3).  Key present
                        // only when the PR visitor is in the pipeline => default
                        // production compiled config stays byte-identical.
                        [if std.member(pipeline_names, 'tagger_check_neutrino')
                         then 'include_vertex_points']: true,
                        // require_pr_graph: this layer is PR output.  Without it,
                        // an event where TaggerCheckNeutrino selects no candidate
                        // gets the WHOLE clustering dumped here in raw (un-T0-cor)
                        // coordinates -- see docs/pr/3 sec. 9.  C++ default false
                        // (the uBooNE steiner-bound sets rely on that fallback);
                        // key present only when the PR visitor is in the pipeline
                        // => default compiled config byte-identical.
                        [if std.member(pipeline_names, 'tagger_check_neutrino')
                         then 'require_pr_graph']: true,
                    },
                    {
                        name: 'shower_track',    // associated points: q=15000 shower, q=0 track
                        visitor: 'TaggerCheckNeutrino:pr',
                        grouping: 'live',
                        detector: 'protodunevd',
                        algorithm: 'shower_track',
                        pcname: '3d',
                        coords: ['x', 'y', 'z'],
                        individual: false,
                        use_associate_points: true,
                        // C++ default false.  Prototype per-particle sub-clustering:
                        // non-shower segments get real_cluster_id = cluster*1000+seg
                        // (NeutrinoID::fill_point_info) so Bee can color by particle
                        // (docs/pr/3).  Key present only when the PR visitor is in
                        // the pipeline => default compiled config byte-identical.
                        [if std.member(pipeline_names, 'tagger_check_neutrino')
                         then 'particle_ids']: true,
                        // doc pr/45.  C++ default false.  Muon-typed (+-13)
                        // pseudo-showers paint as track (q=0), matching the PF
                        // tree's "mu-" verdict from the same cached type.  Key
                        // omitted when off => byte-identical.
                        [if pseudo_shower_track_paint
                         then 'pseudo_shower_track_paint']: true,
                        // require_pr_graph: this layer is PR output.  Without it,
                        // an event where TaggerCheckNeutrino selects no candidate
                        // gets the WHOLE clustering dumped here in raw (un-T0-cor)
                        // coordinates -- see docs/pr/3 sec. 9.  C++ default false
                        // (the uBooNE steiner-bound sets rely on that fallback);
                        // key present only when the PR visitor is in the pipeline
                        // => default compiled config byte-identical.
                        [if std.member(pipeline_names, 'tagger_check_neutrino')
                         then 'require_pr_graph']: true,
                    },
                    {
                        name: 'vertices',        // PR graph vertices; main vertex q=15000
                        visitor: 'TaggerCheckNeutrino:pr',
                        grouping: 'live',
                        detector: 'protodunevd',
                        algorithm: 'vertices',
                        pcname: '3d',
                        coords: ['x', 'y', 'z'],
                        individual: false,
                        use_graph_vertices: true,
                        // require_pr_graph: this layer is PR output.  Without it,
                        // an event where TaggerCheckNeutrino selects no candidate
                        // gets the WHOLE clustering dumped here in raw (un-T0-cor)
                        // coordinates -- see docs/pr/3 sec. 9.  C++ default false
                        // (the uBooNE steiner-bound sets rely on that fallback);
                        // key present only when the PR visitor is in the pipeline
                        // => default compiled config byte-identical.
                        [if std.member(pipeline_names, 'tagger_check_neutrino')
                         then 'require_pr_graph']: true,
                    },
                ]
                // STM fit-trajectory layer, dumped after TaggerCheckSTM runs (doc
                // sbnd_xin/docs/40).  q = dQ*0.1 - 1000, same convention as the
                // PR track_fit layer.  Entry present only when save_stm_fit is on
                // => compiled config byte-identical otherwise.  Name must avoid
                // the substring '-track' (bee3 models.py filters such files).
                + (if save_stm_fit then [{
                    name: 'stm_fit',
                    visitor: 'TaggerCheckSTM:pr',
                    grouping: 'live',
                    detector: 'protodunevd',
                    algorithm: 'stm_fit',
                    pcname: 'stm_fit',
                    coords: ['x', 'y', 'z'],
                    individual: false,
                    dQdx_scale: 0.1,
                    dQdx_offset: -1000.0,
                }] else []),
                // Particle-flow Bee output ("mc" jsTree JSON), emitted once after
                // TaggerCheckNeutrino runs; inert when the visitor is not in the pipeline.
                bee_pf: [
                    {
                        name: 'mc',
                        visitor: 'TaggerCheckNeutrino:pr',
                        grouping: 'live',
                        // C++ defaults false/0 (legacy).  Prototype mc.json parity
                        // (docs/pr/3): TDatabasePDG-style names + integer MeV, and
                        // the WCReader::KeepMC display floors (5 MeV em, 10 MeV
                        // nucleon).  Keys present only when the PR visitor is in
                        // the pipeline => default compiled config byte-identical.
                        [if std.member(pipeline_names, 'tagger_check_neutrino')
                         then 'prototype_names']: true,
                        [if std.member(pipeline_names, 'tagger_check_neutrino')
                         then 'em_ke_min']: 5 * wc.MeV,
                        // doc pr/38: nucleon floor lowered from the prototype's
                        // 10 MeV (WCReader::KeepMC) to 3 MeV, owner decision
                        // 2026-08-05 -- sub-10-MeV protons attached at the
                        // neutrino vertex (18255-447477 3.2 MeV, 18255-52657
                        // 8 MeV) must show in the particle flow.
                        [if std.member(pipeline_names, 'tagger_check_neutrino')
                         then 'np_ke_min']: 3 * wc.MeV,
                        // doc pr/34 §10 port-fidelity knobs.  C++ defaults false;
                        // key omitted when off => byte-identical pre-knob config.
                        [if pf_track_main_cluster_only then 'pf_track_main_cluster_only']: true,
                        [if pf_track_bridged_clusters then 'pf_track_bridged_clusters']: true,   // doc pr/40 round 9 B2
                        [if pf_shower_vertex_barrier then 'pf_shower_vertex_barrier']: true,
                        [if pf_shower_parent_precedence then 'pf_shower_parent_precedence']: true,
                        [if pf_pi0_node_per_id then 'pf_pi0_node_per_id']: true,
                        [if pf_pdg_name_prototype_fallback then 'pf_pdg_name_prototype_fallback']: true,
                        // doc pr/38 Round 4.  C++ default false; key omitted when
                        // off => byte-identical pre-knob config.
                        [if pf_orphan_track_parentage then 'pf_orphan_track_parentage']: true,
                        // doc pr/65 round 3.  C++ default false; key omitted when
                        // off => byte-identical pre-knob config.
                        [if pf_orphan_audit_only then 'pf_orphan_audit_only']: true,
                        // doc pr/84 round 2 (F1/F2).  C++ defaults false (scalars
                        // 3 cm / 8 cm); keys omitted when off/null =>
                        // byte-identical pre-knob config.  Params in cm.
                        [if pf_direct_when_touching then 'pf_direct_when_touching']: true,
                        [if pf_touch_max != null then 'pf_touch_max']: pf_touch_max * wc.cm,
                        [if pf_pseudo_gap_from_main then 'pf_pseudo_gap_from_main']: true,
                        [if pf_unique_node_ids then 'pf_unique_node_ids']: true,
                        [if pf_drop_stray_satellites then 'pf_drop_stray_satellites']: true,
                        // doc pr/93 round 4; params in cm.
                        [if pf_orphan_confident_track then 'pf_orphan_confident_track']: true,
                        [if pf_orphan_track_min_cm != null then 'pf_orphan_track_min']: pf_orphan_track_min_cm * wc.cm,
                        [if pf_orphan_guard_freed then 'pf_orphan_guard_freed']: true,  // doc pr/123 r2
                        // doc pr/128; params in cm.
                        [if pf_orphan_near_cross_cluster then 'pf_orphan_near_cross_cluster']: true,
                        [if pf_orphan_near_gap_cm != null then 'pf_orphan_near_gap']: pf_orphan_near_gap_cm * wc.cm,
                        [if pf_orphan_near_min_len_cm != null then 'pf_orphan_near_min_len']: pf_orphan_near_min_len_cm * wc.cm,
                        [if pf_orphan_near_end_tol_cm != null then 'pf_orphan_near_end_tol']: pf_orphan_near_end_tol_cm * wc.cm,
                        [if pf_orphan_near_kink_deg != null then 'pf_orphan_near_kink_deg']: pf_orphan_near_kink_deg,
                        [if pf_conn4_near_candidate then 'pf_conn4_near_candidate']: true,
                        [if pf_conn4_near_gap_cm != null then 'pf_conn4_near_gap']: pf_conn4_near_gap_cm * wc.cm,
                        [if pf_track_owns_loose_vertex then 'pf_track_owns_loose_vertex']: true,
                    },
                ],
                pipeline: wc.tns(cm_pipeline),
            },
        }, nin=1, nout=1, uses=anodes + [dv, pcts] + cm_pipeline + tagger_uses),
        local sink = g.pnode({
            type: 'TensorFileSink',
            name: 'clus_pr',
            data: {
                outname: if tensor_outname == '' then 'trash-pr.tar.gz' else tensor_outname,
                prefix: 'clustering_',
                dump_mode: tensor_outname == '',
            },
        }, nin=1, nout=0),
        ret:: if dump then g.pipeline([mabc, sink]) else g.pipeline([mabc]),
    }.ret,
    detector_volumes(anodes):: clus_maker.detector_volumes(anodes),
}
