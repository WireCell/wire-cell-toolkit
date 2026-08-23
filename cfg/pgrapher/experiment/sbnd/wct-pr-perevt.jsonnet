// CANONICAL in-tree SBND pattern-recognition job — single source of truth.
// Promoted from sbnd_xin 2026-07-27 (sbnd_xin/docs/64_cfg-sync.md);
// sbnd_xin/wct-pr-perevt.jsonnet is now a one-line re-export of this file, so
// the working chain and the in-tree config cannot drift.  Runnable directly:
//   wire-cell -c pgrapher/experiment/sbnd/wct-pr-perevt.jsonnet --tla-...
// (WIRECELL_PATH must contain toolkit/cfg + wire-cell-data).
//
// The TLA defaults below are the SBND PRODUCTION operating point wherever a
// knob affects physics; pure diagnostic outputs (save_stm_fit, save_tensors,
// dl_weights) stay off.  Production additionally passes -stm-fit.
//
// trackfitting_config defaults to the canonical in-tree SBND parameter file
// pgrapher/experiment/sbnd/sbnd_track_fitting.json, which TaggerCheckSTM resolves
// through WIRECELL_PATH (Persist::resolve; an absolute path is passed through
// unchanged, so the runners' absolute paths keep working).  Passing '' selects the
// C++ preset defaults, which are uBooNE-hard-coded and never right for SBND.
//
// Per-event SBND pattern-recognition (PR) job — standalone, self-contained.
//
// Loads the persisted post-QL point-cloud tree written by the Q/L job
// (run_ql_evt.sh -save-pctree -> work/ql_evt<ID>/pctree-evt<ID>.tar.gz) and
// runs the PR-tail visitors on it inside a MultiAlgBlobClustering node.  With
// pipeline_names=[] it is instead the round-trip identity gate: the
// 'clustering' Bee layer of the output mabc-pr.zip must be content-identical
// to the Q/L job's mabc-all-apa.zip clustering layer.
// See sbnd/docs/sbnd-pattern-recognition.md.  Driven by run_pr_evt.sh /
// run_nusel_evt.sh.
//
// Usage (called from run_pr_evt.sh):
//   wire-cell \
//     --tla-str  input=work/ql_evt2/pctree-evt2.tar.gz \
//     --tla-code anode_indices='[0,1]' \
//     --tla-str  output_dir=work/pr_evt2 \
//     --tla-code run=0 --tla-code subrun=0 --tla-code event=2 \
//     --tla-str  reality=sim \
//     --tla-code 'pipeline_names=[]' \
//     -c wct-pr-perevt.jsonnet

local g = import 'pgraph.jsonnet';
local wc = import 'wirecell.jsonnet';
local tools_maker = import 'pgrapher/common/tools.jsonnet';

function(
    input          = '.',
    anode_indices  = [0, 1],
    output_dir     = '.',
    run            = 0,
    subrun         = 0,
    event          = 0,
    reality        = 'sim',
    // Same LAr TLAs as the Q/L job so the anode/params objects are identical.
    // DL/DT: sbndcode's production diffusion (wcsimsp_sbnd.fcl), restored
    // 2026-07-27, reverting the 6.5781/13.1349 retune of 2026-07-25 (docs/66).
    // Inert here in the same sense as lifetime below -- compiling this job at
    // either pair gives byte-identical JSON, zero 'DL' keys.  The diffusion the
    // PR/STM track fit actually uses lives in sbnd_track_fitting.json, which is
    // read at RUNTIME and is NOT inert.
    DL             = 4.0,
    DT             = 8.8,
    // lifetime: electron lifetime in ms.  35 = the SBND simparams.jsonnet value,
    // so this chain and the simulation now state the same number (owner,
    // 2026-07-27).  It was 6.0, an inherited placeholder from the first standalone
    // Q/L test (wcp-porting-img 655bd6a: "DL=6.2 DT=9.8 lifetime=6
    // driftSpeed=1.565") whose DL/DT half was corrected to the SBND physical
    // values in 9f498089 while lifetime was not.
    // The change is provably inert: this knob only feeds the sim Drifter's
    // attenuation and NO reco component reads it -- compiling this job with
    // lifetime=6 and lifetime=35 gives byte-identical JSON (zero 'lifetime' keys).
    // NO electron-lifetime / charge-attenuation correction is applied anywhere in
    // imaging, clustering, Q/L matching or PR; adding one would be a real change
    // needing a measured data lifetime and a dQ/dx revisit, not this knob.
    // See sbnd_xin/docs/64_cfg-sync.md sec 4.
    lifetime       = 35.0,
    driftSpeed     = 1.563,
    // PR visitors to run, by name (see clus_pr's cm_by_name in
    // cfg/pgrapher/experiment/sbnd/clus.jsonnet).  The default IS the SBND
    // production tagger chain (run_full1k_nusel.sh with -unmerge-assoc), so a
    // bare `wire-cell -c` run of this job exercises TGM/STM/FC.  Pass [] for the
    // pass-through round-trip identity gate (run_pr_evt.sh does), or a subset
    // like ['switch_scope','steiner','fiducialutils','tagger_check_stm'] for a
    // single-tagger demo.  Production additionally appends 'stm_magnify' when
    // save_stm_fit=true (the -stm-fit ROOT dump); it is left out of the default
    // because save_stm_fit defaults false.
    // Neutrino-PR chain (docs/pr/2-3): append
    // ['tagger_check_neutrino','numu_bdt_scorer','nue_bdt_scorer',
    //  'tracking_visitor','tagger_output'] -- ordering matters (BDTs after the
    // neutrino stage, nue after numu, tagger_output after tracking_visitor
    // because it opens tracking-pr.root in UPDATE mode).
    // protect_bundle + steiner_refresh (doc pr/23 sec 9): SBND PRODUCTION
    // DEFAULT since 2026-08-02 (owner, after the sec 8 fresh-tree gate: 0
    // event_label / nu_evaluated flips in 572 valfast events).  They sit
    // after tagger_check_fc (cosmic verdicts on UNSPLIT clusters,
    // prototype wire-cell-prod-stm.cxx:806) and before any appended
    // tagger_check_neutrino (prototype wire-cell-prod-nue.cxx:1322).
    // steiner_refresh (replace=false) MUST follow protect_bundle: it
    // rebuilds only the steiner products the split purged; a replace=true
    // second pass dangles the tagger-stage GraphAlgorithms (evt 54095).
    // NOT bit-identical to the pre-pr/23 chain by design; drop BOTH names
    // from pipeline_names for the pre-flip chain (runner:
    // SBND_PROTECT_BUNDLE=0).
    pipeline_names = ['switch_scope', 'unmerge_bundle', 'unmerge_assoc', 'steiner',
                      'fiducialutils', 'tagger_check_tgm', 'tagger_check_stm',
                      'tagger_check_fc', 'protect_bundle', 'steiner_refresh'],
    // TrackFitting parameter JSON, required whenever tagger_check_stm is in the
    // pipeline: the C++ preset defaults are uBooNE-hard-coded, never right for
    // SBND.  DEFAULT = the canonical in-tree file; TaggerCheckSTM resolves it
    // through WIRECELL_PATH (Persist::resolve), and an absolute path is passed
    // through unchanged, so the runners' absolute paths keep working.  '' selects
    // the uBooNE presets.
    trackfitting_config = 'pgrapher/experiment/sbnd/sbnd_track_fitting.json',
    // MIP dQ/dx scale in e/cm handed to TaggerCheckSTM.  56000 = the SBND value
    // (docs/48), matching the *DeDx tables now regenerated at 0.5 kV/cm.  Pass
    // 50000 (the MicroBooNE value, and the C++ default) to isolate the
    // reference-table change from the MIP-scale change in an A/B.
    mip_dqdx       = 56000,
    // Cathode kink veto (doc pr/20 Part II B0), both cm.  segment_search_kink's
    // four accept criteria are each guarded by para_angle > 7.5-15 deg, a guard
    // against isochronous imaging artifacts that is INVERTED at the cathode:
    // the crossing is drift-x dominated so para_angle is 61-78 deg, wide open,
    // while the ~2 cm transverse mismatch across the gap supplies the turn.
    // cathode_kink_xcut suppresses kink ACCEPT candidates within that distance
    // of x = cathode_x; the angle arithmetic itself is untouched.
    // null on both = C++ default 0 = OFF = the legacy kink search, with the keys
    // omitted (gate B0-2).  SBND DEFAULT ON (owner 2026-08-02 after the Part VI
    // Bee scan): 5 cm around x = 0, matching clus.jsonnet's pr() default so a
    // bare run is production.  NOT bit-identical -- doc pr/20 Part VI sec 1.
    // Pass -A cathode_kink_xcut=null -A cathode_x=null for the legacy search.
    cathode_x         = 0,
    cathode_kink_xcut = 5,
    // Wide-baseline cathode kink ACCEPT (doc pr/47 sec 8, O1), the converse
    // of the veto above: a genuine kink AT the crossing never reaches the
    // legacy accept thresholds (the gap/distortion suppresses every local
    // refl_angle window -- 52085's 33-38 deg junction reads ~23 deg), so a
    // fifth accept path fires when the skirt-excluded PCA turn angle across
    // the crossing is >= this many degrees.  SBND DEFAULT ON (owner
    // 2026-08-07): 25 deg, mid-gap in the bimodal census distribution (doc
    // pr/47 sec 6).  NOT bit-identical -- doc pr/47 sec 8.
    // Pass -A cathode_wide_kink_angle=null for the legacy search.  Skirt /
    // baseline null = C++ defaults 3 cm / 15 cm.
    cathode_wide_kink_angle    = 25,
    cathode_wide_kink_skirt    = null,
    cathode_wide_kink_baseline = null,
    // shower_topo_demote_len (cm, doc pr/25 sec 3): demote any segment the
    // topology test flags kShowerTopology whose geometric length exceeds this,
    // so it receives real track PID instead of the hard-coded pdg=11/score=100.
    // Rationale: the test's only measurement axis satisfies
    // dir_3.xhat = sin(angle-to-drift), so for a near-isochronous segment it
    // measures spread along the DRIFT axis, where points sit on a 0.313 cm
    // time-slice lattice -- against a 0.4 cm "large spread" cut.  86 of 91 long
    // firings across the 572-event valfast manifest carry no evidence above
    // that noise floor, and an owner hand-scan of all 10 such segments on a
    // selected nu-candidate main cluster found 10/10 tracks, 0 showers.
    // null = C++ default 0 = OFF = byte-identical.  SHIPS OFF: the flip is the
    // owner's call.  50 reproduces the scan-supported rule for 9 of the 10
    // scanned events; ~45 covers all 10 (evt 400504 measures 49.1 cm by
    // segment_track_length, the length the guard uses).
    // **SBND PRODUCTION DEFAULT ON at 50 cm, owner 2026-08-03**, after a
    // hand-scan of all 10 long shower-topology segments on a selected
    // nu-candidate main cluster returned 10/10 TRACKS, 0 showers.
    // Measured at 50 cm over 572 events: 23 long segments flip
    // shower->track across 22 events, 0 the other way; 17/572 events (3.0%)
    // move numu_score; evt 321107 goes pdg 11 -> 13, numu_score
    // -0.783 -> +0.317.  Score moves are MIXED -- evt 286353 drops
    // 2.023 -> -1.139 -- accepted by the owner with that caveat on record
    // (doc pr/25 sec 3.10).  50 covers 9 of the 10 scanned events; evt
    // 400504 measures 49.1 cm by segment_track_length (the length this
    // guard uses), so ~45 would cover all 10.
    // Set null for the legacy behavior (long segments stay shower-eligible).
    shower_topo_demote_len = 50,
    // ---- doc sbnd_xin/docs/pr/30 §11 port-fidelity knobs -----------------
    // SBND operating point for the four pr/30 findings.  This file is the
    // single source of the SBND operating point (doc 68), so a bare run
    // reproduces production exactly; diagnostic arms are produced by editing
    // a COPY of cfg/ (WIRECELL_PATH override) or the runner's env hooks.
    //   fit_exclusion (P1)            true  => the knobbed do_multi_tracking
    //                                          sites pass flag_exclusion, as 28
    //                                          of 30 live prototype sites do.
    //   graph_endpoint_strict (P8)    true  => PR::add_segment REFUSES a
    //                                          vertex/segment pair whose vertices
    //                                          are not at the segment's ends.
    //                                          The WARN + counter are always on.
    //   oov_prototype_parity (F2)     true  => a point outside every TPC votes
    //                                          the way the prototype's own helper
    //                                          answers, at all three sites.
    //   first_seg_local_pca (P2)      false => drop the toolkit-only local-PCA
    //                                          refinement of the first segment's
    //                                          endpoints (null/true = production).
    //   other_seg_relaxed_accept (P4) false => drop the toolkit-only
    //                                          0.72/15cm/1.05 acceptance clause
    //                                          (null/true = production).
    // **SBND PRODUCTION DEFAULT ON, owner 2026-08-20** (doc pr/98 §7: fits
    // equal-or-better in 11/12 top movers after the -1.0 sentinel fix; perf
    // rounds put the cost at ~1.15x median / 1.85x worst on nueCC48).  C++
    // default stays false; -A fit_exclusion=false (or SBND_FIT_EXCLUSION=false)
    // restores the pre-flip production path byte-exactly (doc pr/98 §10).
    fit_exclusion = true,
    graph_endpoint_strict = false,
    graph_endpoint_tol = null,
    // **SBND PRODUCTION DEFAULT ON, owner 2026-08-04** (doc pr/30 §12.10).
    // This is a port BUG, not tunable behaviour: at all three out-of-TPC guard
    // sites the toolkit's implicit vote is the OPPOSITE of what the
    // prototype's own helper returns for a point with no readout --
    //   modify_segment_isochronous  toolkit "connected" vs is_good_point FALSE
    //   examine_vertices_1p         toolkit "dead"      vs get_closest_dead_chs FALSE
    //   examine_vertices_3          toolkit "not unique" (the segment is then
    //                               REMOVED) vs a normally-evaluated distance
    // Flipped ON while it is provably free: measured byte-identical on 48/48
    // nueCC data events (work-pr30-baseHEAD vs work-pr30-f2on), because the
    // guards fire once in total on this manifest.  The population where they
    // WOULD fire -- the cathode region and the readout-window edges -- is
    // barely represented here, so landing it now costs nothing and avoids
    // having to A/B it later as a real behaviour change.
    // Set false for the pre-flip arm.
    oov_prototype_parity = true,
    first_seg_local_pca = null,
    other_seg_relaxed_accept = null,
    // shower_topo_proto_dir (doc pr/31 §11, F2 was P2)
    //                               true  => skip the stage-3
    //                                        segment_determine_shower_direction
    //                                        call, so a kShowerTopology segment
    //                                        keeps the direction
    //                                        segment_is_shower_topology set.
    // The prototype runs determine_dir_shower_topology here, which does not
    // touch flag_dir; its determine_shower_direction() is called from exactly
    // one place in the whole prototype tree and that place is stage 4
    // (NeutrinoID_track_shower.h:1532).  Deliberately NOT flipped on: unlike
    // oov_prototype_parity this is not a clear port bug -- the prototype's own
    // function carries a "// hack for now" comment and has both of its
    // direction blocks commented out, so the toolkit's 305-line PCA may be the
    // better physics.  OFF here = today's path = byte-identical; the arm that
    // measures it sets this true.
    shower_topo_proto_dir = false,
    // ---- doc pr/32 §11: the four kept findings of the stage-4 (neutrino
    // vertex identification) port audit.
    //
    // **ALL FOUR ARE SBND PRODUCTION DEFAULTS ON, owner 2026-08-04.**  Each is
    // a port BUG or a missing port, not a toolkit design choice, and all four
    // land byte-identical on the nueCC48 manifest: work-pr32r2-off48 vs
    // work-pr32r2-{f1on48,f2on48,f34on48,allon48} is 48/48 on pctree-pr
    // member-content hashes and identical on nusel-table.tsv /
    // nusel-events.tsv.  They are engaged, not inert -- see the PR32AUDIT
    // counters quoted per knob below.
    //
    // vertex_dir_use_fit_point (F1, was P1): calc_conflict_maps and
    //   compare_main_vertices_all_showers measure from the vertex fit SNAPPED
    //   to the nearest Steiner node; the prototype uses get_fit_pt(), the
    //   continuous fit (NeutrinoID_track_shower.h:1804/1808, :1468/:1480/
    //   :1504-1505/:1542/:1551/:1567/:1574).  Eleven expressions in two
    //   functions, against 23 expressions in the same file that already read
    //   `fit().valid() ? fit : wcpt`.  NOT byte-identical when on.
    // MEASURED: 5047 reads on 48 events, ALL with a valid fit;
    //   mean |fit - wcpt| 0.613 cm, max 11.41 cm.
    vertex_dir_use_fit_point = true,
    // shower_traj_recheck_parity (F2, was P3): the improve_vertex
    //   shower-trajectory recheck.  Prototype reads the STORED flag at both
    //   outer gates (NeutrinoID_improve_vertex.h:248, :287) and recomputes the
    //   inner test at the 10 cm default with the 50000-analog scale; the
    //   toolkit recomputes both gates and runs the inner test at 1.0 cm with
    //   the 43e3-analog.  Fixing the inner parameters ALONE makes the block
    //   dead code, so this knob moves all three together -- including letting
    //   segment_is_shower_trajectory clear kShowerTrajectory, which the
    //   prototype's is_shower_trajectory does on every call and the toolkit
    //   has never done.
    // MEASURED: the stored flag and a fresh 10 cm test
    //   disagree 22 times on 48 events; parity demotes 4 segments and runs the
    //   recheck body 4 times instead of 7.
    shower_traj_recheck_parity = true,
    // main_vertex_require_descriptor (F3, was P7): compare_main_vertices
    //   guards two of its six blocks on descriptor_valid(); an
    //   invalid-descriptor candidate reaches the argmax with score
    //   `0 + (0.5 if in FV) - conflicts/4` and can beat real candidates.
    //   Expected byte-identical -- the path looks unreachable -- but that is a
    //   control-flow argument, so the drop is counted (PR32AUDIT f3_dropped).
    // MEASURED: 0 of 2219 candidates dropped on 48
    //   events -- P7 is now MEASURED dead code, not argued dead code.
    main_vertex_require_descriptor = true,
    // main_vertex_candidate_flag (F4, was P12): set kMainCandidate on each
    //   per-cluster main-vertex candidate, the prototype's
    //   map_cluster_main_candidate_vertices (NeutrinoID_track_shower.h:1332).
    //   DIAGNOSTIC ONLY: only PrDisplayDump reads it, exactly as the prototype
    //   exposes it only to app-level tree fillers.
    // MEASURED: 3774 vertices flagged on 48 events.
    main_vertex_candidate_flag = true,
    // ---- doc pr/31 §12: the §10.12 topology/PID/direction port-fidelity
    // round.  All five fixes ship default-OFF in C++/cfg; the values here are
    // the SBND operating point.
    //
    // **ALL FIVE ARE SBND PRODUCTION DEFAULTS ON, owner 2026-08-04** (the
    // owner's standing instruction for this round: production keeps bug
    // fixes and improvements).  Each restores a prototype behaviour the port
    // dropped by accident, and every one is MEASURED NULL on the nueCC48
    // manifest: gate work-pr31r2-off48 vs work-pr32r2-allon48 48/48; per-knob
    // arms work-pr31r2-{f5on,f6on,f3on,f1onb,f4on}48 and the joint
    // work-pr31r2-allonb48 are 48/48 pctree-hash identical to off, with both
    // nusel TSVs, every T_tagger/T_kine leaf and every mabc-pr.zip member
    // unchanged.  Null does NOT mean inert: F1's preserve path fires on 47/48
    // events (its first implementation crashed there, doc pr/31 §12.5) -- the
    // rewritten 4-momenta are simply never consumed by any persisted output.
    // The flips therefore change no production number today and close five
    // latent divergence classes.
    //
    // cont_muon_dir3_30cm (F5, was P6): find_cont_muon_segment_nue's hoisted
    //   dir3 always at 30 cm as the prototype computes it, instead of falling
    //   back to the 15 cm dir1 for a short reference segment.  The reachable
    //   divergent case is "short reference segment, long neighbour".
    cont_muon_dir3_30cm = true,
    // track_comp_empty_abstain (F6, was P7): an empty dQ/dx comparison window
    //   ABSTAINS from the direction gate instead of confirming the
    //   orientation.  0.0 is the prototype's degenerate answer, verified by
    //   execution (zero-bin TH1F KolmogorovTest -> ks1 == ks2 == 0 ->
    //   eval_ks_ratio false).
    track_comp_empty_abstain = true,
    // shower_topo_reset (F3, was P13): segment_is_shower_topology clears
    //   kShowerTopology and dirsign at entry, before its early returns, as the
    //   prototype does (ProtoSegment.cxx:319-321) -- no stale flag survives a
    //   re-test.  Also closes shower_topo_proto_dir's stated residual.
    shower_topo_reset = true,
    // reclass_preserve_4mom (F1, was P1+P3a+P4): the 15 reclassification
    //   sites preserve the existing 4-momentum, recomputing only where the
    //   prototype's get_particle_4mom(3)>0 guard passes.  The only fix in this
    //   round that moves kine_reco_Enu directly.
    reclass_preserve_4mom = true,
    // dir_track_median_local (F4, was P8): determine_dir_track's median
    //   dQ/dx over the SAME local vector the PID receives (prototype
    //   nth_element), not the filtered helper rebuild.
    dir_track_median_local = true,
    // examine_showers_vertex_by_index (F7, was P5): DELIBERATELY DORMANT.
    //   Orders examine_all_showers' vertex pair by graph index before the
    //   asymmetric 165/150-degree branches -- A deterministic convention, not
    //   provably the prototype's.  Stays false pending pr/30 F4's
    //   find_vertices adjudication (that decision owns all three known
    //   order-sensitive callers).
    examine_showers_vertex_by_index = false,
    // Steiner TERMINAL filter fidelity -- doc pr/29 D1 and D12.
    //
    // **SBND PRODUCTION DEFAULT ON, owner 2026-08-04**: both are port BUGS, not
    // tunable behaviour.  The toolkit filter was strictly tighter than the WCP
    // prototype in two independent ways, and on SBND evt 388 that discarded
    // 47.7% of every cluster's Steiner terminals and left 24 clusters below the
    // two-terminal minimum -- i.e. with no Steiner tree at all.  With both on:
    // 20.0% and 3.  NOT bit-identical: PR segments 75 -> 88, vertices 118 ->
    // 131, numu_score -3.199 -> -2.166, kine_reco_Enu 2900.5 -> 2865.6 MeV,
    // and the event's pi0 candidate disappears (kine_pio_flag 1 -> 0) -- the
    // event label, TGM/STM/FC verdicts do not move.  Measured on ONE event
    // against a zero noise floor (same-binary repeat byte-identical); doc pr/29
    // sec.11.4/11.7.
    //
    // The C++ defaults stay OFF, so uBooNE/ICARUS/PDHD/PDVD are unaffected --
    // this is an SBND operating-point decision and lives here (doc 68: the
    // operating point is in cfg only, and a bare run IS production).
    // Set steiner_terminal_wire_tol=0 and steiner_terminal_adjacent_slice=false
    // for the pre-fix arm that every legacy comparison needs.
    //   steiner_terminal_wire_tol = 1     restores the prototype's one wire of
    //     slack on both sides of all three planes in the terminal containment
    //     test (PR3DCluster_steiner.h:285-290, :310-315, :336-341).  Under the
    //     toolkit's half-open [min,max) convention that is `>= min-1 && < max+1`
    //     -- NOT `<= max+1`, which would be asymmetric.  get_extreme_wcps shares
    //     the same C++ helper and is deliberately left at 0, matching its own
    //     prototype counterpart (PR3DCluster_path.h:111-119), which has no slack.
    //   steiner_terminal_adjacent_slice = true    makes the prototype's t+-1
    //     slice fallback resolve.  The map is keyed by blob->slice_index_min(),
    //     which is in TICKS, so the historical step of 1 names no real slice on
    //     SBND (nticks_live_slice = 4) and the branch is dead code.  ON, the
    //     step comes from the face's own nticks-per-slice, not a literal 4.
    // Both can only ADD terminals: each restores a way to PASS the filter.
    steiner_terminal_wire_tol = 1,
    steiner_terminal_adjacent_slice = true,
    // Steiner EDGE-WEIGHT charge fidelity (doc pr/29 D2).  ON at the owner's
    // instruction, same reasoning as the two above: the prototype weights
    // steiner edges with charges computed under the disable_dead_mix_cell this
    // chain actually passes (false, CreateSteinerGraph.cxx:234 ->
    // PR3DCluster_steiner.h:514,:521), and the toolkit dropped the argument at
    // the call so create_enhanced_steiner_graph's `= true` default won.  Not a
    // tunable: the toolkit was computing a charge the port did not intend.
    //   ON  -> only planes with a NONZERO charge value enter the Qs/Qt RMS.
    //   OFF -> all three planes enter and DEAD ones (charge_uncertainty > 1e10)
    //          are subtracted -- an independent predicate, so the two disagree
    //          on any live-but-zero or dead-but-nonzero plane.
    // Unlike the terminal-filter knobs this moves weights in EITHER direction,
    // so it can add or drop tree edges rather than only add terminals.
    // Set to false for the pre-fix arm.
    steiner_edge_charge_forward_dead_mix = true,
    // Isochronous first-segment endpoint finding (doc pr/24 round 2, SBND evt
    // 271851): for a long cluster whose quantile-trimmed drift-x extent is
    // small (a filled 2-D sheet), the first PR segment's endpoints come from
    // the sheet's principal axis instead of the wire-footprint boundary
    // metric (which degenerates to sheet-edge corners and produces the
    // two-edge-track fan), and the local-PCA endpoint refinement is skipped
    // on that branch.
    // **SBND PRODUCTION DEFAULT ON, owner 2026-08-03**, after the doc pr/24
    // sec 15 (round 3) scan.  Round 3 is what made it flippable: the endpoint
    // is the UNTRIMMED axial extreme (round 2's quantile-trimmed one left a tip
    // stub that find_other_segments claimed, planting a 0.9 deg "vertex" in the
    // straight track of evt 284794), and iso_endpoint_min_aspect requires real
    // 2-D sheet-ness so 1-D tracks never enter the branch -- evts 284794 and
    // 59899 are byte-identical to legacy again.  Evidence: evt 271851 vertex
    // 30.9 -> 2.4 cm from truth on the DL arm (the production chain), nue_score
    // -15 -> +4.30; evt 122660 nue -15 -> +4.30; 0/17 mcp1k events gain a
    // straight-through segment junction (round 2: 4/17); every non-firing event
    // stays archive-identical.  CAVEAT on record (doc pr/24 sec 15.6): on the
    // DIAGNOSTIC geometric-vertex arm evt 271851 goes 1.7 -> 6.6 cm and loses
    // the nue selection -- accepted by the owner, and the first thing to check
    // if the next valfast census surprises.
    // The four sizing knobs stay null = the C++ defaults (40 cm min length,
    // 25 cm max drift extent, 0.35 frac, 0.02 quantile, 4 cm diagnostic tube
    // radius, 0.12 min sheet aspect).
    // Set iso_endpoint=false for the legacy wire-footprint boundary endpoints
    // (SBND_ISO_ENDPOINT=0 on the run_pr_chain_batch.sh runner).
    // NOT bit-identical to the pre-flip chain -- a deliberate production change.
    iso_endpoint = true,
    iso_endpoint_min_length = null,
    iso_endpoint_max_xext = null,
    iso_endpoint_xext_frac = null,
    iso_endpoint_xext_quantile = null,
    iso_endpoint_tube_radius = null,
    iso_endpoint_min_aspect = null,
    // examine_vertices_3 extension-retraction guard (doc pr/24 sec 18, round
    // 5, SBND 18259-42280 / 18255-271851 / 18255-350186).  get_local_extension's
    // Hough-based direction estimate can return a point CLOSER to the far
    // endpoint than the vertex it started from -- worst at the axial extreme
    // of an isochronous sheet, which is exactly where iso_endpoint (above)
    // picks its seed, so an iso-picked trunk tip can lose 7.5-8.9 cm of its
    // delivered trajectory to a stage nominally meant to extend it further.
    // SBND production default since 2026-08-06 (owner: "you can put them on
    // for SBND running", doc pr/24 sec 19.1, gate work-pr24r6-off48 vs
    // work-pr24r5-off48 -- legacy restorable via SBND_V3_EXT_GUARD=0).
    // v3_extension_min_gain stays null = the C++ default (-1 cm); the key is
    // suppressed in the compiled config either way.
    v3_extension_guard = true,
    v3_extension_min_gain = null,
    // doc sbnd_xin/docs/pr/67 -- LOG-ONLY diagnostic probe for "the fitted
    // trajectory does not cover the image", worst in isochronous topologies
    // (owner cases 18264-137238 / 18259-42280 / 18345-21073 / 18255-58717).
    // Names the find_iso_first_segment_endpoints gate that rejected a cluster
    // (today only the aspect gate logs), reports get_local_extension's
    // perpendicular-to-drift no-op, the per-round find_other_segments census,
    // and what examine_end_ps_vec trims off each end.
    // C++ default false; key suppressed when off => byte-identical config.
    // Set SBND_TRAJ_COVER_PROBE=1 on run_pr_chain_batch.sh.
    traj_cover_probe = false,
    // doc pr/67 -- counterfactual ONLY, not a production setting.
    // find_proto_vertex's nrounds_find_other_tracks is hardcoded (2 for the
    // main cluster) with no config surface; > 0 overrides it so the owner's
    // "not sufficient rounds of branch searching" hypothesis can be measured.
    // C++ default 0 = keep the hardcoded 2.  A value > 0 CHANGES OUTPUT.
    // Set SBND_PR_FIND_OTHER_ROUNDS=<n>.
    pr_find_other_rounds = null,
    // protect_bundle stage knobs (doc pr/23): the PR-stage overclustering
    // protection (uboone's second graph examination).  The stage is in
    // pipeline_names by DEFAULT since the sec 9 production flip; the knobs
    // are inert only when it is dropped.  The cathode re-join
    // values are the SBND operating point in INTERNAL units (unlike
    // cathode_kink_xcut above, which is cm), matching clus.jsonnet's pr()
    // defaults so a bare run is production; nulls = prototype-faithful
    // (re-join pass disabled -- would break cathode crossers, doc pr/20).
    // **SBND PRODUCTION DEFAULT 'relaxed_strict_img_2d_rescue_long', owner
    // 2026-08-11 (doc pr/62; adds S7 on top of doc pr/57 round 6's S6+rescue,
    // preserved below).** S6's flood-fill kill test is capped at
    // `s6_dis_cap`=30cm (its cost tracks candidate distance, and its
    // `cell_budget` breaker fails CLOSED on exhaustion -- raising the cap
    // trades missed gaps for silent false-kills of large sparse real
    // objects).  Above 30cm the only surviving test was S1-S3's per-plane
    // *independent* radius query, which a track can pass with charge from
    // three unrelated nearby objects.  S7 (`Graphs::long_corridor_bad`,
    // `connect_graph_relaxed_strict.cxx`) closes exactly that gap: a
    // corridor-restricted BFS between the candidate's own two endpoint
    // lattice cells, O(D) not O(D^2) so the budget stays non-live, breaker
    // inverted to fail OPEN (abstain, not kill, on exhaustion) -- the
    // opposite of S6's posture, deliberately.  Operating point
    // min_gapped_planes=1 (mirrors S6's owner-instructed single-plane rule);
    // the more conservative min_gapped_planes=2 flavor
    // ('relaxed_strict_img_2d_rescue_long2') ships alongside, not selected
    // here.  117-event validation (48 nueCC + 19 NCpi0 + 50 PR-data): 17/117
    // movers, nusel (final neutrino selection) byte-identical in every
    // event at BOTH operating points (identical mover sets on this sample);
    // S7 kill rate on evaluated candidates 97.4% (184/189) -- higher than
    // doc pr/56 sec 8.4's own S6 ceiling (89% at 15-30cm), an aggressive
    // point, not a conservative one.  min_gapped_planes/gap_floor_cm are
    // UNFITTED -- no hand-scan labels exist yet in this distance band,
    // unlike S6's 899-label fit below.  Owner reviewed the 17-event Bee
    // before/after and flipped.  protect_bundle is the only consumer, so
    // clustering and TGM/STM/FC verdicts stay byte-identical for every
    // event S7 does not touch.  See doc pr/62.
    //
    // S6+rescue itself (doc pr/57 round 6, owner 2026-08-10, superseding
    // round 7's 'relaxed_strict_img'): adds S6, a per-plane 2D wire/tick
    // gap-kill OFF in every earlier default, plus a pure post-kill RESCUE
    // (`Graphs::two_d_rescue_ok`, killed -> kept only, never the reverse)
    // fitted to the owner's full 899-label hand scan (230 events / 847
    // component pairs).  Rescue recognizes dead-W-channel gaps, prolonged/
    // isochronous induction-plane signal on high-direction-consistency long
    // tracks, W allowed a tighter gap than U/V, and W-anchored two-view
    // footprint co-location.  Also repairs a dir-MST emission gap: a killed
    // closest candidate no longer silently drops a pair's dir-bridge
    // emissions when S6/rescue judged it should survive.  Final table on the
    // owner's hand scan: bad-separation (must reconnect) 124/127,
    // good-separation (must stay split) 154/156, the two misses being the
    // owner's own unsatisfiable label triangles.
    //
    // **SBND PRODUCTION DEFAULT '..._rescue_long_wtrack', owner 2026-08-11
    // (doc pr/64 round 4; adds the S6 W-plane long-track exception on top of
    // pr/62's S7, preserved above).** W (collection) has no excuse channel in
    // S6 -- `excused[3]={excuse_u,excuse_v,false}` -- so a genuinely straight
    // long track with a small W-only hole is broken by the one plane that
    // cannot be excused (18259-174224 plain 6-wire W charge hole,
    // 18255-276836 drift-parallel all-plane prolongation where only W
    // convicts, 18255-314507).  `Graphs::two_d_w_track_ok` revives an
    // S6-killed candidate iff W is the SOLE voting plane and both components
    // are long/thin/globally-collinear (Lmin>6cm, npmin>=50, Tmax<2cm,
    // angle<25deg, tightened by Tmax<1.7cm OR angle<6deg -- the owner-chosen
    // zero-regression point on the 899-label hand scan: protects evt122660's
    // good pair, forgoes the marginal 64959/172656 recoveries), or a dead-W
    // band explains the gap (dead_w>=3, npmin>=20, dis<3cm).  Revive-only,
    // composed after the shipped rescue, whose local-axis W branch is
    // byte-untouched; the exception uses the GLOBAL per-component PCA axis
    // (the local break-point axis is noisy exactly at a break -- 174224:
    // local 20.7deg vs global 4.2deg).  Validation (doc pr/64 round 4):
    // off-gates 0/66 movers on same-config fresh arms + 1000-event bare
    // rerun fully attributed to pr/65's flip; on-arm movers on the FULL
    // 1000-event PR-data sample are EXACTLY the three target events
    // (174224 fully rejoined, 276836/314507 partially -- their all-plane S7
    // breaks remain by design) plus evt122660's owner-OK 4/5/6 component
    // merge on nueCC48; the owner-labelled good pair 122660 13-16 stays
    // split on-arm; nusel byte-identical in all 1067 gated events.
    //
    // Legacy escapes: -A protect_graph_name=relaxed_strict_img_2d_rescue_long
    // (SBND_PROTECT_GRAPH=relaxed_strict_img_2d_rescue_long) restores
    // pre-pr/64 production (S7, no W exception);
    // =relaxed_strict_img_2d_rescue restores pre-pr/62
    // production (round 6, S6+rescue only, no S7);
    // =relaxed_strict_img_2d_rescue_long2 selects the conservative S7
    // min_gapped_planes=2 point; =relaxed_strict_img restores round 7;
    // =relaxed_strict restores round 6's predecessor; =relaxed restores
    // pre-round-6.
    protect_graph_name          = 'relaxed_strict_img_2d_rescue_long_wtrack',   // null => 'relaxed'
    // C++ default true (doc pr/23 ordering): a TGM/STM/lm-convicted in-window
    // main does not open its bundle for splitting.  null => key omitted.
    protect_skip_convicted      = null,
    // doc pr/94 round 3.  When true a convicted main still OPENS its bundle,
    // so the bundle's unconvicted members -- the secondary activity the
    // demoted-main fallback goes on to call the neutrino -- get the second
    // graph examination every member of an unconvicted bundle gets.  The
    // convicted cluster itself is still never split.  SBND 18255/395148: the
    // 198.9 cm secondary keeps a graph bridge that fits 17 cm of trajectory
    // through empty space because ClusteringProtectBundle logged
    // "OC53SKIP main ident=10 ... convicted STM=1 -- bundle not opened".
    // **SBND PRODUCTION DEFAULT ON since 2026-08-19 (owner flip, doc sec
    // 9.13).**  This is the knob that removes SBND 395148's 17 cm trajectory
    // excursion through empty space (fit points > 3 cm from any charge
    // 15 -> 0, worst 8.40 -> 0.83 cm).  null => C++ default false; pre-flip
    // arm: SBND_OPEN_CONVICTED_BUNDLES=0.
    protect_open_convicted_bundles = true,
    protect_cathode_x           = 0,
    protect_cathode_rejoin_xcut = 5 * wc.cm,
    protect_cathode_rejoin_dyz  = 4 * wc.cm,
    protect_cathode_rejoin_dis  = 8 * wc.cm,
    // Direction-agreement fallback for a dyz-only failure (doc pr/25 sec 1,
    // SBND evt 489327).  **SBND PRODUCTION DEFAULT ON, owner 2026-08-03.**
    // A cathode crosser whose two halves are collinear but whose transverse
    // travel exceeds cathode_rejoin_dyz -- which happens by construction as
    // a track tips toward the cathode plane -- is re-joined on direction
    // agreement instead of being left split.  Operating point from doc
    // pr/25 sec 1: measured 5/5 target re-joins, 1/1 junk pair still
    // rejected (89.2 deg direction mismatch, also fails dir_npts), 8/8
    // regression events byte- AND score-identical.
    // INTERNAL units except _angle (degrees) and _dir_npts (count).
    // Set null for the legacy dyz-only behavior.  Per doc 68 the SBND
    // operating point lives HERE only; clus.jsonnet's clus_pr()/pr()
    // function defaults stay null (the doc pr/23 both-files trap is about
    // an explicit null in the OTHER direction overriding a value set here).
    protect_cathode_rejoin_perp       = 3 * wc.cm,
    protect_cathode_rejoin_angle      = 20.0,
    protect_cathode_rejoin_dir_radius = 15 * wc.cm,
    protect_cathode_rejoin_dir_npts   = 20,
    // TaggerCheckSTM's containment gate (cluster_fc_check) uses the same fiducial
    // + margins as tagger_check_tgm / tagger_check_fc, so "contained" means one
    // thing across all three verdicts.  true = SBND production (docs/49); false
    // restores the pre-doc-49 FiducialUtils sensitive-volume fallback, which is
    // what the A/B compares against.  Runner flag: -no-stm-fv.
    stm_consistent_fv = true,
    // doc-63 round-1 STM acceptance guards (charge-desert one-objectness,
    // spike-not-ramp nu-vertex veto, eval ratio2 cap).  C++ default false;
    // key omitted when off => byte-identical compiled config.  DEFAULT TRUE
    // = SBND production as of doc 63 (owner 2026-07-26); pass false for a
    // pre-campaign A/B.  Runner flag: -stm-guards / SBND_STM_GUARDS=1.
    stm_accept_guards = true,
    // doc-63 round-2 muon-consistency guard on detect_proton (an end region
    // matching the muon hypothesis in shape and normalization is not called a
    // proton).  C++ default false; key omitted when off => byte-identical.
    // DEFAULT TRUE = SBND production as of doc 63 (owner 2026-07-26).
    // Runner flag: -stm-proton-guard / SBND_STM_PROTON_GUARD=1.
    stm_proton_muon_guard = true,
    // doc-63 round-3 cathode-truncation veto (a fitted stop within ~5 cm of
    // the CPA with no Bragg rise did not stop -- drift-boundary truncation).
    // C++ default false; key omitted when off => byte-identical.  DEFAULT
    // TRUE = SBND production as of doc 63 (owner 2026-07-26).  Runner flag:
    // -stm-cathode-guard / SBND_STM_CATHODE_GUARD=1.
    stm_cathode_guard = true,
    // doc-63 round-4a fix of the inverted face selection in the STM tagger's
    // dist_to_anode helper (the anode-clipped-TGM catch fired at the SBND
    // cathode).  C++ default false; key omitted when off => byte-identical.
    // DEFAULT TRUE = SBND production as of doc 63 (owner 2026-07-26, after
    // validation).  Runner flag: -stm-anode-fix / SBND_STM_ANODE_FIX=1.
    stm_anode_dist_fix = true,
    // doc-63 round-4b/4c second-track vetoes (long straight leftover past
    // the kink = V topology; long MIP-like other-track segment).  C++
    // default false; key omitted when off => byte-identical.  DEFAULT TRUE
    // = SBND production as of doc 63 (owner 2026-07-26, after validation).
    // Runner flag: -stm-track-guard / SBND_STM_TRACK_GUARD=1.
    stm_second_track_guard = true,
    // doc-63 round-5 stop-region vetoes.  5a deficit_guard: an end whose
    // last-5cm median dQ/dx is below 0.6 MIP with no Bragg rise in the last
    // 15 cm is a reconstruction truncation, not a stop.  5b
    // vertex_kink_guard: a sharp (>45 deg) turn within 12 cm of the stop
    // into a >2.2 MIP prong is a nu vertex plus proton, not a Bragg.  C++
    // defaults false; keys omitted when off => byte-identical.  DEFAULT
    // TRUE = SBND production as of doc 63 (owner 2026-07-26, after
    // validation).  Runner flags: -stm-deficit-guard / SBND_STM_DEFICIT_GUARD=1,
    // -stm-vertex-guard / SBND_STM_VERTEX_GUARD=1.
    stm_deficit_guard = true,
    stm_vertex_kink_guard = true,
    // doc-66 sec 12 diffusion-margin cut package: Michel-veto res_length
    // floor 6 -> 6.5 cm, detect_proton track_medium gate 1.0 -> 1.05,
    // block-B ks2 entry 0.05 -> 0.055, C1 peak clause 4.3 -> 4.1.  Restores
    // the four vetoes the 4.0/8.8 diffusion revert lost by hairline margins
    // (owner-adjudicated; full-1000-event sweep found zero collateral).  C++
    // defaults are the prototype constants; false omits the keys =>
    // byte-identical pre-package config.  The four values themselves are
    // module defaults in clus.jsonnet (pr(...) args stm_michel_res_cm etc.).
    // DEFAULT TRUE = SBND production as of doc 66 (owner 2026-07-27).
    // Runner flag: -no-stm-d66cuts / SBND_STM_D66CUTS=0.
    stm_d66_cuts = true,
    // Detector-extent literals of the PR chain (docs/pr/2 sec 2e(iv)).  The
    // uBooNE prototype's cosmic_tagger quotes its four "reaches the top" cuts as
    // 17 / 15 / 37 / 67 cm below its y = +117 cm top face (100/102/80/50 cm);
    // pr_y_top re-anchors the same offsets to a detector whose top is elsewhere.
    // DEFAULT 200 = SBND (sensitive |y| <= 199.965 cm) => 183/185/163/133 cm.
    // The uBooNE arm of an A/B is one flag: --tla-code 'pr_y_top=117'
    // (reproduces 100/102/80/50 exactly), paired with vertex_z_prior_scale=200.
    pr_y_top = 200.0,
    // Denominator (cm) of the upstream-z main-vertex penalty (z-min_z)/scale.
    // uBooNE 200 cm over a 1037 cm detector; DEFAULT 100 = SBND's 501 cm scaled
    // (200 x 501/1037 = 96.6).  See the clus.jsonnet arg comment for the
    // alternative reading that keeps 200.
    vertex_z_prior_scale = 100.0,
    // SSM beam-line reference directions [x,y,z] in the detector frame, feeding
    // the 8 ssm_*_angle_{target,absorber} BDT features (docs/pr/2 sec 2e(i)).
    // null = the C++ defaults = the prototype's uBooNE BNB-target
    // (0.46,0.05,0.885) / NuMI-absorber (0.33,0.75,-0.59) vectors, so the
    // compiled config is unchanged.  NO SBND VALUE EXISTS YET -- these are
    // reachable, not calibrated.  Diagnostic use: --tla-code
    // 'ssm_target_dir=[0,0,1]' makes ssm_*_angle_target identical to
    // ssm_*_angle_z, which is how the plumbing was verified.
    ssm_target_dir   = null,
    ssm_absorber_dir = null,
    // Charge -> kinetic-energy calibration constants of NeutrinoEnergyReco
    // (docs/pr/2 sec 2e(iii)), i.e. what turns the summed 2D charge of a shower
    // or track into kine_charge [MeV]:
    //   E = sum_p(w_p Q_p)/sum(w) / recom / fudge * kine_w_value * 1e-6 MeV
    // null = the C++ defaults = the uBooNE-tuned literals these keys replaced,
    // so the compiled config is unchanged.  NO SBND VALUE EXISTS FOR ANY OF
    // THEM -- reachable, not calibrated:
    //   kine_recom_factor        0.7   average recombination survival, track
    //   kine_fudge_factor        0.95  residual scale, track
    //   kine_shower_recom_factor 0.5   } same pair when the object is
    //   kine_shower_fudge_factor 0.8   } flagged shower-like
    //   kine_proton_recom_factor 0.35  |pdg|==2212 (fudge stays at the track one)
    //   kine_plane_weights  [0.25,0.25,1.0]  [U,V,W] charge-average weights
    //   kine_plane_asym_switch   0.04  drop the largest plane above this
    //                                  (median,max) relative asymmetry
    //   kine_w_value             23.6  argon W-value in eV
    // The recombination/fudge pairs are field- and calibration-dependent; the
    // plane weights and the asymmetry switch encode uBooNE induction-plane
    // quality.  Sensitivity check: --tla-code 'kine_recom_factor=0.6' scales
    // every kine_charge by 0.7/0.6.
    // 2026-07-30 (docs/pr/10 sec 6): the three RECOMBINATION factors now
    // default to the SBND transfer values -- the uBooNE empiricals scaled by
    // the table-integrated survival ratio (official uBooNE Box at 0.273
    // kV/cm -> doc-55 free-power SBND fit, C excluded).  The fudge factors
    // stay uBooNE (they absorb gain/lifetime normalization, which the C of
    // the fit carries -- keeping both would double-count).  Effect: Enu in
    // the T_kine tree drops 12-14% on nuecc48 172230/235435/444187.
    // null on any factor restores its uBooNE value.
    kine_fudge_factor        = null,
    kine_recom_factor        = 0.87,   // 0.70 x 1.249 (track)
    kine_shower_fudge_factor = null,
    kine_shower_recom_factor = 0.58,   // 0.50 x 1.169 (shower)
    kine_proton_recom_factor = 0.51,   // 0.35 x 1.453 (proton)
    kine_plane_weights       = null,
    kine_plane_asym_switch   = null,
    kine_w_value             = null,
    // doc pr/35 sec 10.2 (F1 = P1+P8): read the shower PDG live from the
    // start segment at the four fill_kine_tree sites (prototype kine.h:53
    // :67 :175 :187) instead of Shower's cached particle_type, whose refresh
    // path is incomplete.  C++ default false = the cached read.
    // SBND PRODUCTION DEFAULT ON (owner 2026-08-04, doc pr/35 sec 11): the
    // unconditional cached-vs-live counter measured 2/48 nueCC events with a
    // stale cache (13 vs live 11); knob-on moves ONLY the calib-pr kine block
    // on exactly those two (137238: type 13->11 + Enu -0.511 MeV; 469665:
    // type array only).  Gate work-pr35-off48 vs work-pr34-prod48 +
    // work-pr35-prod48 vs work-pr35-f1on48.  =false restores the cached read.
    kine_shower_pdg_live     = true,
    // Muon median-dQ/dx-vs-length envelope [c0, c1, pivot_cm, power]:
    //   dQ_dx_cut = c0 + c1*(pivot/L)^power   (a multiple of mip_dqdx_median)
    // used by nine tagger cuts (numu x2, vertex-finder, nue x4, ssm, cosmic).
    // DEFAULT = the docs/pr/10 sec 4 SBND fit (2026-07-30): SBND stopping-
    // muon table median (0.5 kV/cm, / mip_dqdx_median=48000) times the
    // uBooNE empirical/table margin g(L), same functional form.  Scales as
    // 1/mip_dqdx_median -- re-derive when 48000 becomes a measurement.
    // null = the C++ defaults = the prototype's empirical uBooNE refit
    // [0.8866, 0.9533, 18, 0.4234] (byte-identical pre-knob config; also
    // bit-identical when passed explicitly -- verified, docs/pr/10 sec 7).
    muon_dqdx_curve          = [0.8826, 1.0587, 18, 0.4745],
    // Recombination-model selection for BOTH taggers (STM + neutrino PR).
    // DEFAULT ON (owner 2026-07-30, docs/pr/10 sec 7): sbnd_power_recomb,
    // the free-power Modified Box fitted to SBND stopping tracks (docs/55
    // sec 7g canonical) -- every model-driven dQ/dx -> dE/dx conversion
    // (kinematics, PID energies) uses the measured SBND recombination.
    // false restores sbnd_box_recomb (A=1.0, B=0.255 at 0.5 kV/cm), the
    // byte-identical pre-pr/10 config.
    use_power_recomb         = true,
    // Single-photon stem dE/dx: DEFAULT ON (owner 2026-07-30) -- route
    // shw_sp_vec_{median,mean}_dedx through the configured recombination
    // model above, with sp_mean_dedx_cut = 2.23 MeV/cm, the physical-scale
    // transfer of the legacy 2.3 (docs/pr/2 sec 2e(i) third correctness
    // item; docs/pr/10 sec 5).  false/null restore the inline uBooNE-field
    // (0.273 kV/cm) inverse Box and the 2.3 literal.
    sp_dedx_use_recomb_model = true,
    sp_mean_dedx_cut         = 2.23,
    // dl_vtx_cut (mm) is a configurability thread only (docs/pr/2 sec
    // 7.4); null keeps the C++ 25.0 (= 2.5 cm) default.
    dl_vtx_cut               = null,
    // When set, re-save the (post-PR) tree to this TensorDM tarball.  Used by the
    // round-trip gate (re-save with pipeline_names=[] and compare member hashes
    // against the input tarball).
    save_tensors   = '',
    // SCN (DL) neutrino-vertex weights, WIRECELL_PATH-resolved.
    // DEFAULT = ON, the uBooNE-trained net (owner adopted 2026-07-30 on nueCC48
    // evt 18253/1/172230: the geometric vertex sat at the far end of a proton
    // track, the DL vertex moved it 9.7 cm onto the true interaction point --
    // docs/pr/4).  '' restores the geometric-only vertex, which is what every
    // identity gate must keep passing (CLAUDE.md M4).
    // The weights are still uBooNE-TRAINED -- SBND retraining is docs/pr/2 gap G3.
    // Inert unless 'tagger_check_neutrino' is in pipeline_names.
    // REQUIRES libpython preloaded in the job:
    //     LD_PRELOAD=$(python3 -c "import sysconfig;print(sysconfig.get_config_var('LIBDIR'))")/libpython3.11.so.1.0
    // Without it the SCN module import fails and TaggerCheckNeutrino silently
    // falls back to the geometric vertex after one WARN line -- always
    //     grep -c "DL vertex failed" <log>   # expect 0
    dl_weights     = 'uboone/scn_vtx/t48k-m16-l5-lr5d-res0.5-CP24.pth',
    // DL re-rank operating point (TaggerCheckNeutrino), threaded as TLAs for
    // A/B runs (doc pr/79).  min_accept 4.0 -> 10.0 adopted 2026-08-15
    // (owner, doc pr/79: +36/473 on the hand-scan live A/B); top_k=5 kept —
    // the pr/79 step-3 k=20 arm measured exactly zero realized gain (0
    // fixed / 0 regressed; admitted candidates all lose the rerank).
    // Pass 4.0 for the pre-flip arm.
    dl_vtx_min_accept_score = 10.0,
    dl_vtx_top_k   = 5,
    // doc pr/105: false = the legacy single-argmax DL branch (no composite
    // re-rank; top voxel snapped, dl_vtx_cut gate, else the traditional
    // vertex).  Default true = production; unset => byte-identical compiled
    // config.  Validation: --tla-code dl_vtx_rerank=false (or
    // SBND_DL_VTX_RERANK=false).
    dl_vtx_rerank = true,
    // Beam window [low, high) in us on cluster_t0 (= matched flash time) selecting
    // the bundle that gets neutrino PR.  [0,0] disables the gate (then
    // tagger_check_neutrino falls back to uBooNE single-main selection, which on
    // SBND picks an arbitrary main -- always set a window with tagger_check_neutrino).
    beam_window_us = [0.2, 2.2],
    // Run the PR tail (steiner + tagger_check_{tgm,stm,fc}) ONLY on the
    // beam-coincident bundle, i.e. the mains whose cluster_t0 is inside
    // beam_window_us (plus, for steiner, the companions sharing their
    // matched_flash_gid).  TRUE = SBND production default as of docs/56;
    // false restores the evaluate-every-bundle behavior with a compiled config
    // byte-identical to the pre-doc-56 one.  Inert when beam_window_us is empty.
    // Runner flag: -no-bwonly.
    beam_window_only = true,
    // nu_skip_cosmic / nu_skip_cosmic_bundle: TaggerCheckNeutrino refuses to run
    // neutrino PR on an in-window main already convicted as cosmic -- per-main
    // (TGM/STM/lm_flag) and, lifted to the whole flash bundle, per-bundle.
    // BOTH TRUE = the SBND production operating point (docs/pr/3 sec. 8,
    // docs/pr/16 sec. 7); these TLAs only surface what clus.jsonnet's pr()
    // already defaults to, so the compiled config is unchanged when they are
    // not passed.  Setting them false is a DIAGNOSTIC, not an operating point:
    // it makes PR run on a convicted cosmic, which is how you get track_fit /
    // shower_track / vertices layers for an event that legitimately has no
    // neutrino candidate (e.g. SBND evt 116962, sole in-window main TGM).
    // Runner env: SBND_NU_SKIP_COSMIC=0|1.
    nu_skip_cosmic = true,
    nu_skip_cosmic_bundle = true,
    // Enable the ported check_neutrino_candidate veto in tagger_check_tgm so
    // in-beam-window bundles may be tagged TGM (C++ default false; key
    // omitted when off => byte-identical pre-port config).
    tgm_neutrino_candidate = true,
    // Require charge along the chord between the two extreme points before a
    // pair may tag TGM (C++ default false; key omitted when off =>
    // byte-identical pre-fix config).  Guards against the QL flash-time merge
    // grafting a detached fragment into the tagged cluster -- see
    // docs/29_tgm-chord-charge.md.
    tgm_chord_charge = true,
    // How the chord-charge guard measures support (C++ default "chord"; key
    // omitted then => byte-identical).  "chord" samples the STRAIGHT segment
    // between the extremes -- it falsely rejects curved tracks (evt285185
    // cluster 16: a continuous 480 cm top->anode crosser bows 10 cm off its
    // chord and lost a real TGM, doc 31).  "path" requires a piecewise charge
    // path through the cluster's own points with no jump > chord_max_gap.
    // The -chord runner flag passes "path" (SBND_TGM_CHORD_MODE overrides).
    tgm_chord_mode = 'path',
    // Find the extreme points PER connected component and union them, instead
    // of 8 global extremes over the whole flash-merged cluster (C++ default
    // false; key omitted when off).  Must be used WITH tgm_chord_charge -- see
    // docs/29_tgm-chord-charge.md.  The -chord runner flag sets both.
    tgm_component_extremes = true,
    // Let a component shorter than component_min_length still donate its
    // extremes when it is path-connected (30 cm-step charge path) to a
    // component that passed the length cut, so a fragmented genuine track END
    // keeps its wall exit (evt286681 cluster 7) while detached specks stay
    // dropped (C++ default false; key omitted when off => byte-identical).
    // Only meaningful WITH tgm_component_extremes.  Runner flag: -rescue.
    tgm_component_rescue = true,
    // Rescued-end pairs must ALSO pass the straight-chord support test even
    // in path mode (C++ default false; key omitted when off => byte-identical
    // doc-32 rescue behavior).  Path mode alone lets a rescued speck pair
    // across TWO merged cosmics through an L-shaped charge detour (evt288727
    // cluster 6, doc 33).  Only meaningful WITH tgm_component_rescue.
    // Runner flag: -rescue-chord.
    tgm_rescue_chord = true,
    // A pair may tag TGM only when at least one end lies in the cluster's
    // MAIN charge component -- the largest 30 cm path component; a cathode
    // crosser is ONE such component (C++ default false; key omitted when off
    // => byte-identical).  A merged-in fragment that is itself through-going
    // otherwise tags the whole bundle on its own within-component pair,
    // which the chord guard deliberately allows and the nu-candidate veto
    // cannot protect (it walks the pair's own path): evt289343 cluster 9,
    // in-beam bundle tagged TGM by a 26 cm corner-clipping cosmic fragment
    // 450 cm from the main track (doc 36).  Runner flag: -main-pair.
    tgm_main_pair = true,
    // How tgm_main_pair identifies the main (C++ default "path"; key omitted
    // then => byte-identical doc-36 behavior).  "path" = largest 30 cm path
    // component (proxy; wrong if a merged-in cosmic outweighs the main).
    // "real" = per-blob real_cluster_main flash-merge provenance -- exact,
    // needs a pctree saved with run_ql_evt.sh -save-rcid; falls back to the
    // proxy on old tarballs (doc 38).  Runner flag: -main-pair-real.
    tgm_main_pair_mode = 'real',
    // Exempt a flag_demoted_main cluster from main_pair_rejects (C++ default
    // false; key omitted when off => byte-identical).  With tgm_main_pair on,
    // a demoted main's own real_cluster_main/path-component provenance is
    // all-zero by construction post-carve, so the guard vetoed EVERY
    // demoted-main pair unconditionally, before any CASE-A/CASE-B geometry
    // ran (doc pr/25, SBND evt 320029/18255-1: cluster 30, a 37 cm
    // corner-clipping shape with both ends on a boundary).
    //
    // SBND DEFAULT ON (owner 2026-08-02, doc pr/25 sec 2): small-group
    // measurement (15 events: the target + the pr/20 Part I 14-event
    // demoted-main census) showed the fix works on the target (cluster 30
    // correctly convicts and its charge is excluded via the already-shipped
    // skip_cosmic_companions guard) but ALSO surfaces a pre-existing,
    // unrelated check_tgm CASE-A gap (the "2 extreme groups" branch has no
    // length floor), letting a handful of wall-hugging demoted-main debris
    // specks (8-14 points) pick up a spurious tag.  Owner judged the
    // measured impact small enough to accept (one event moved
    // kine_reco_Enu_MeV by +3.6 out of ~800, i.e. ~0.4%) rather than block
    // on a separate co-requisite fix to check_tgm.  NOT bit-identical.
    // Set false to restore the pre-flip behavior (runner:
    // SBND_TGM_EXEMPT_DEMOTED_MAIN=0).
    tgm_exempt_demoted_main = true,
    // Downstream-z (z ~ 500 cm face) inset of the TGM/FC fiducial box, in cm
    // (default 3 = byte-identical legacy margin).  Shared by tagger_check_tgm
    // and tagger_check_fc so containment keeps one meaning.  Runner flag:
    // -fvz <cm>.
    tgm_fv_zmax_margin = 5,
    // Downstream-z inset used by check_tgm's CASE-A INTERIOR-support tests
    // (chord midpoints + waypoint re-check) when > 0, in cm (default 0 =
    // OFF, key omitted => byte-identical; the interior tests then share
    // tgm_fv_zmax_margin).  Makes the doc-32 widening endpoint-only so a
    // wall-hugging corner clipper keeps its midpoint support (doc 35,
    // evt287517 cluster 16 / evt289805 cluster 9).  TGM only; FC untouched.
    // Runner flag: -fvzi <cm>.
    tgm_fv_zmax_margin_interior = 3,
    // Drift-x and vertical-y insets of the TGM/FC fiducial box, in cm, both
    // faces symmetric (defaults 2 / 2.5 = byte-identical legacy margins).
    // Shared by tagger_check_tgm and tagger_check_fc, same as
    // tgm_fv_zmax_margin.  Runner flags: -fvx <cm> / -fvy <cm>.
    tgm_fv_x_margin = 2.5,
    tgm_fv_y_margin = 3,
    // Persist the per-pass STM track fits (C++ default false; key omitted
    // when off => byte-identical): cluster PCs stm_fit/stm_pass/stm_eval, a
    // Bee 'stm_fit' layer in mabc-pr.zip, and (when 'stm_magnify' is added
    // to pipeline_names) tracking-stm.root for Magnify-tracking-SBND.  Also
    // gates loading the WireCellRoot plugin.  Runner flag: -stm-fit.
    save_stm_fit = false,
    // ---- doc pr/34 §10: particle-flow (Bee mc tree) port-fidelity knobs ----
    // **ALL FIVE ARE SBND PRODUCTION DEFAULTS ON, owner 2026-08-04** (doc
    // pr/34 §11; C++ defaults stay false).  Display-only, MEASURED: the ON
    // arms move ONLY mabc-pr.zip::data/0/0-mc.json (F2 5/48 evts, F3 6/48,
    // F1/F4/F5 null) while pctree-pr, the other six mabc members and both
    // nusel TSVs are 48/48 byte-identical (gate work-pr34-off48 vs
    // work-pr31r2-prod48; joint work-pr34-allon48; bare==prod
    // work-pr34-prod48).  F4's merged pi0 home = HIGHEST-ENERGY daughter's
    // parent (owner decision 2026-08-04, deliberately not prototype
    // first-writer-wins).  Runner env overrides: SBND_PF_* (tri-state).
    pf_track_main_cluster_only = true,
    // doc pr/40 round 9 B2: let the PF track BFS traverse nv-bridged
    // clusters.  C++ default false; key omitted when off => byte-identical.
    // SBND PRODUCTION DEFAULT ON since 2026-08-18 (with the bridge above).
    pf_track_bridged_clusters = true,
    pf_shower_vertex_barrier = true,
    pf_shower_parent_precedence = true,
    pf_pi0_node_per_id = true,
    pf_pdg_name_prototype_fallback = true,
    // doc pr/38 Round 4 pf_orphan_track_parentage: graph-faithful parentage
    // for barrier-orphaned PF track segments (a guard-excluded muon that
    // continues an EM arm attaches under that shower's leaf; a proton behind
    // an orphaned pi+ chains under the pi+; SBND 18255-142421 segs
    // 7011/7012/7018).  Display-only: moves ONLY mc.json.  C++ default
    // false.  SBND PRODUCTION DEFAULT ON (gate work-pr44-off48 vs
    // pr43_cleanhead_ref48b 48/48 + nusel 0-diff; on-census 1/48 nueCC48 +
    // 4/19 ncpi0 events, every move an attributed re-parent -- doc pr/38 R4
    // + doc pr/44).  Runner env: SBND_PF_ORPHAN_TRACK_PARENTAGE.
    pf_orphan_track_parentage = true,
    // doc pr/65 round 3 (rung 3) -- audit-only orphan net: a still-unclaimed
    // segment is named in the log instead of fabricated as a root-level PF
    // particle (owner requirement: never add an arbitrary particle to the PF
    // root).  Display-only: moves ONLY mc.json.  SBND PRODUCTION ON -- owner
    // flip 2026-08-11 with rung 1 above (which claims the fragments first;
    // audit fired 0 unclaimed on all 48 nueCC gate events).  Legacy escape:
    // -A pf_orphan_audit_only=false (or SBND_PF_ORPHAN_AUDIT_ONLY=false).
    // C++ default stays false.
    pf_orphan_audit_only = true,
    // doc pr/84 round 2 (F1/F2) -- vertex-touching pseudo-parent suppression
    // and the remote-gap anchor.  Display-only: moves ONLY mc.json (gate
    // pr84r2_disp_gate PASS=512 FAIL=0; nusel scores 0 movers; no label
    // flips).  SBND PRODUCTION ON 2026-08-17: F1 suppressed 32 spurious
    // gamma/neutron carriers in 512 events (owner evts 283713 x2, 316025,
    // 407280 + round-1's 65289/347129/169626/174752); F2 re-anchored 106
    // remote carriers to draw their real gaps.  pf_touch_max stays null =
    // C++ 3 cm.  Rung 2 (pf_touch_cross_main) stays OFF: the F1.0 probe
    // showed Flags::main_cluster is NOT set on the event-body cluster at PF
    // writer time (evt 64921 deferred -- root cause is vertex determination,
    // doc pr/52 territory).  Runner env: SBND_PF_DIRECT_WHEN_TOUCHING,
    // SBND_PF_TOUCH_MAX, SBND_PF_TOUCH_CROSS_MAIN, SBND_PF_TOUCH_CROSS_MAX,
    // SBND_PF_PSEUDO_GAP_FROM_MAIN.
    pf_direct_when_touching = true,
    pf_touch_max = null,
    pf_touch_cross_main = false,
    pf_touch_cross_max = null,
    pf_pseudo_gap_from_main = true,
    // doc pr/84 round 3 (G1) -- guarantee unique jsTree node ids in mc.json.
    // Bee keys its PF tree model by node id, so a repeated id is invalid
    // input: on SBND 394532 a node and its own descendant both carried id
    // 8033 and selecting it blanked the whole PF panel (owner report
    // 2026-08-18).  Display-only: moves ONLY mc.json, and only when a
    // collision exists.  Kept as the standing invariant even with
    // shower_dedup_start_seg on (which removes the known source).  Runner
    // env: SBND_PF_UNIQUE_NODE_IDS.  C++ default stays false.
    // SBND PRODUCTION ON 2026-08-18 (owner flip): with S1 also on it fired 0
    // times on the 24-event gate manifest -- it is carried as the standing
    // invariant, not as a fix.  Legacy escape: -A pf_unique_node_ids=false.
    pf_unique_node_ids = true,
    // doc sbnd_xin/docs/pr/92: mirror the kine-side stray-satellite drop in
    // the Bee PF tree (inert while kine_drop_stray_satellites is off -- the
    // dropped-id set is then empty).  C++ default false.  Runner env:
    // SBND_PF_DROP_STRAY_SATELLITES.
    // SBND PRODUCTION DEFAULT ON since 2026-08-18 (doc pr/92, with
    // kine_drop_stray_satellites above).
    pf_drop_stray_satellites = true,
    // doc pr/93 round 4 -- see the round-4 block below.  SBND PRODUCTION ON
    // (owner round 2026-08-18; validation in that block).
    pf_orphan_confident_track = true,
    pf_orphan_track_min_cm = null,  // null => C++ default 50cm
    pf_track_owns_loose_vertex = true,
    // doc pr/38: F2's ON-behavior was CORRECTED in place (owner decision
    // 2026-08-05, no new knobs): the barrier now excludes each shower's
    // start vertex (prototype map_vtx_segs parity, WCShower.cxx:547) so
    // main-track junctions stay traversable, and the same knob emits
    // BFS-unreached non-shower main-cluster segments as root-level leaves
    // (prototype flat-loop mc_mother=0 parity).  Display-only: moves ONLY
    // mc.json (SBND 18255-219295/489330/56243 protons + vertex proton stubs
    // with the 3 MeV np_ke_min floor).
    // How 'unmerge_bundle' (when named in pipeline_names) identifies the main
    // sub-cluster of a flash-merged bundle (C++ default "real").  "real" =
    // per-blob real_cluster_main/real_cluster_id provenance -- exact, needs a
    // pctree saved with run_ql_evt.sh -save-rcid; a cluster without it is
    // left UNSPLIT (warning), never proxied.  "component" = longest connected
    // component -- a clustering decision (breaks cathode crossers), opt-in
    // only.  Runner flags: -unmerge / -unmerge-comp (doc 45).
    unmerge_bundle_mode = 'real',
    // restore_demoted_mains (doc pr/20 Part I P2; C++ default false = OFF):
    // when the OUTER (flash) un-merge splits a bundle, additionally tag each
    // split-off part that was ITSELF a matched bundle main before the merge
    // with flag_demoted_main -- read from the per-blob "real_cluster_was_main"
    // array, so the Q/L stage must have run with save_bundle_main_provenance
    // (else the visitor warns and flags nothing).  The part keeps
    // flag_associated_cluster and does NOT get flag_main_cluster back.
    //
    // SBND DEFAULT ON since doc pr/20 Part X.  Fails CLOSED: on a legacy Q/L
    // tree without the array the visitor warns once per cluster and flags
    // nothing, so the chain still runs -- but P3/P4 below are then inert and
    // the run is NOT the production operating point.  Regenerate the Q/L tree.
    // Set false to restore the pre-flip behaviour (runner:
    // SBND_RESTORE_DEMOTED_MAINS=0).
    restore_demoted_mains = true,
    // require_provenance (doc pr/23 sec 4.2; C++ default false = the historical
    // warn-and-skip).  SBND DEFAULT TRUE: with restore_demoted_mains on, a
    // pctree with NO 'real_cluster_was_main' array is a LEGACY pre-pr/20 tree
    // and the job ABORTS loudly instead of silently reproducing pre-pr/20
    // behaviour (the mistake documented in doc pr/23 sec 4.2).  Intentional
    // legacy-tree runs (e.g. valfast's pinned PR-tail hubs) declare it:
    // runner SBND_REQUIRE_WASMAIN=0.
    require_provenance = true,
    // evaluate_demoted_mains (doc pr/20 Part I P3; C++ default false = OFF):
    // let TaggerCheckTGM / STM / FC evaluate a cluster carrying
    // flag_demoted_main.  Inert unless restore_demoted_mains above is on, since
    // nothing else sets that flag.
    //
    // SBND DEFAULT ON since doc pr/20 Part X.  On its own this produces a
    // verdict and acts on nothing -- gated at 1000 events (doc pr/20 Part IX
    // sec 4): 4 differing cells vs the knob-off arm, all `kine_reco_Enu` last
    // digit, and 3 of the 4 events also differ in the A/A control with the
    // identical value change.  Set false to restore the pre-flip behaviour
    // (runner: SBND_EVAL_DEMOTED_MAINS=0).
    evaluate_demoted_mains = true,
    // skip_cosmic_companions / cosmic_companion_min_length (doc pr/20 Part I
    // P4; C++ defaults false / 0 = OFF): drop a TGM- or STM-tagged COMPANION
    // from the neutrino's other_clusters, so its charge cannot reach the
    // particle-flow tree or kine_reco_Enu.  A tagged companion SHORTER than the
    // floor (cm) stays in regardless, so a mis-tagged neutrino daughter is
    // never silently lost.  Inert without evaluate_demoted_mains above.
    //
    // SBND DEFAULT ON since doc pr/20 Part X.  **This changes reconstructed
    // physics** -- 14 of 1000 mcp1k events move (doc pr/20 Part IX sec 5-6),
    // zero beam labels change, and the motivating event 18255/59003 goes
    // kine_reco_Enu 1202.5 -> 841.0 MeV.  nueCC48 is untouched: 0 of its 69
    // admitted demoted mains is convicted (gate PI-8, run at floor 0, which by
    // the monotonicity of `L >= floor` discharges it at every floor).
    //
    // THE FLOOR IS NOT DERIVED FROM THE IMPACT DISTRIBUTION.  Part IX sec 6b:
    // length does not predict impact (the four inert drops span 3.3-158.7 cm)
    // and a floor does not stop the vertex relocations (the largest is at
    // 39.5 cm).  15 cm is a tagger-plausibility guard -- every conviction in
    // 1000 events came from STM, never TGM, and evt 73004's companion is an
    // 11.9 cm "stopping muon" that stops 177 cm from the cathode.  Do not
    // re-derive this number from a length histogram; it is not in there.
    // Set false / null to restore the pre-flip behaviour (runner:
    // SBND_SKIP_COSMIC_COMPANIONS=0).
    skip_cosmic_companions = true,
    cosmic_companion_min_length = 15,
    // nu_fallback_demoted_mains (sbnd_xin/docs/73 sec 12, round 3; C++
    // default false = OFF): when the primary loop selects NO candidate
    // (every in-window main convicted or vetoed), consider DEMOTED mains
    // (flag_demoted_main) under the same window / cosmic / bundle-veto
    // gates.  Never runs when a main-cluster candidate exists.  Motivated
    // by SBND data evt 65289: the cathode-rescue join legitimately became
    // an STM and the examined, untagged 88.9 cm former main was never a
    // candidate.  Inert unless restore_demoted_mains is on; pairs with
    // evaluate_demoted_mains so candidates carry tagger verdicts.
    // SBND PRODUCTION ON since 2026-08-17 (owner flip on the docs/73
    // sec 12 validation: recovers 65289 and four census events, one
    // accepted re-admission 398690 -- sec 12.8).  Escape:
    // SBND_NU_FALLBACK_DEMOTED=0 omits the key (pre-round-3 baseline).
    nu_fallback_demoted_mains = true,
    // sp_photon_flag (doc pr/26 sec. 8.2; C++ default false = OFF): store the
    // single-photon tagger's verdict in TaggerInfo::photon_flag, the way
    // prototype NeutrinoID.cxx:271 does.  The port already runs
    // singlephoton_tagger() and fills its ~90 shw_sp_* BDT features -- only the
    // return value was discarded, leaving the uBooNE tagger ntuple's
    // photon_flag branch a constant 0.
    //
    // SBND DEFAULT ON (owner 2026-08-03).  This changes ONE written output
    // branch and nothing else: photon_flag has no readers anywhere in clus/,
    // root/ or cfg/ -- its only other occurrence is
    // UbooneTaggerOutputVisitor.cxx:1078's SCALAR_BR -- so no reconstruction,
    // no selection and no other tagger field moves.  Measured on the 7-event
    // gate (doc pr/26 sec. 9.3): 1215 of 1216 T_tagger branch-values identical,
    // the one that moves is photon_flag 0 -> 1 on evt 172230.
    // Set false to restore the pre-fix gap (runner: SBND_SP_PHOTON_FLAG=0).
    sp_photon_flag = true,
    // ---- doc sbnd_xin/docs/pr/36 sec 11 tagger-stage knobs -----------------
    // SBND PRODUCTION DEFAULTS (owner 2026-08-04, doc pr/36 sec 11): the
    // gate for this stage is TWO artifacts (sec 10.9): T_tagger in
    // tracking-pr.root (leaf compare) PLUS the tagger block of
    // calib-pr-evt<ID>.json (PR_EXTRA_STAGES=pr_display) -- match_isFC is NOT
    // booked in T_tagger.  Gate labels: work-pr35-prod48 vs work-pr36-off48
    // 48/48 everything; work-pr36-allon48 (env-forced) == work-pr36-prod48
    // (bare).  Escape hatches: SBND_NEUTRINO_CONSISTENT_FV=0 etc.
    //
    // F1 ON: match_isFC recomputed on the same fiducial + margins the
    // STM/TGM/FC taggers use (mirrors stm_consistent_fv above).  Measured
    // 6/48 nueCC events flip contained->exiting (54095 74544 137238 168596
    // 268784 360535), numu_score moves on all six, nue_score only on 137238;
    // no nusel verdict flips.  =false restores the FiducialUtils fallback.
    neutrino_consistent_fv = true,
    // sbnd_xin/docs/74 G1/G2 ON (owner 2026-08-20 "please implement these
    // fixes"): cosmic_tagger() containment on the same sbnd_pr_fv + margins
    // as TGM/STM/FC (its inside_fv lambda + the flag-1 vertex test; the
    // FiducialUtils fallback has no wall inset and excludes the CPA slab,
    // leaving 5 of the 10 cosmict checks structurally near-dead on SBND).
    // Measured (doc 74 sec 9): OFF gate 234/234 byte-identical; ON census
    // 117 nu-MC events -> 2 false cosmic tags cleared (52085 flag_1 = the
    // CPA-band artifact at x=-1.09 cm, 48895 flag_2), numu_score moves on
    // 12 events, nue_score never.  =false restores the FiducialUtils
    // fallback.
    cosmic_consistent_fv = true,
    // sbnd_xin/docs/75 ON (owner 2026-08-20, "things are good, turn them on
    // for SBND production"): routes the SAME sbnd_pr_fv + margins into the
    // nue/single-photon taggers' containment tests (angular_cut,
    // shower_to_wall, bad_reconstruction_2/_2_sp) -- the identical
    // zero-margin FiducialUtils inconsistency doc 74 fixed for cosmic_tagger,
    // found in the nue/SP tagger family by the doc 75 FV audit.  OFF gate
    // 286/286 archives + 143/143 events byte-identical; ON census 2/48 flag
    // flips (one inert, one fix-direction) + 2/48 small nue_score moves on
    // nueCC48, 1/19 on NCpi0 (a known F1-margin event), 0/50 on numu-50, no
    // ADVERSE (doc 75 sec 6).  =false restores the FiducialUtils fallback.
    nue_sp_consistent_fv = true,
    // F3 OFF (owner 2026-08-04): single-photon SCE gate.  Vacuous today (no
    // SBND SCE helper; clus_geom_helper is ''), proven zero-movement with
    // the knob forced; OFF keeps kine and single-photon independently
    // gateable when a helper lands.
    sp_sce_correction = false,
    // F4 ON: graph-index-ordered tagger accumulation sets (M4 house rule,
    // prototype n/a).  48/48 byte-identical on nueCC48 -- value is run-to-run
    // stability under a different address layout.
    tagger_ordered_segment_sets = true,
    // F5 ON: prototype wcpt-identity stem-endpoint rule (18 sites).  16/18
    // sites exact-match every call; the two neither-match firings (268067,
    // 350186, site 17) measured at 4.8-95 cm = the indirect-shower
    // population where the prototype also picks back.  Movement: those two
    // events, 4 shw_sp_lol_* ntuple fields each, scores untouched.
    stem_endpoint_wcpt_parity = true,
    // F6 ON: broken_muon_id counts distinct cluster IDs (prototype
    // semantics).  48/48 byte-identical; id==pointer count everywhere on
    // this sample (counter f6_id_ptr_disagree is the doc-53 tripwire).
    broken_muon_cluster_id_count = true,
    // F7 ON: neutrino_type verdict bitmask + its T_tagger branch (/I).
    // Knob-on diff = exactly the one new branch on 47/47 PR events.
    neutrino_type_bitmask = true,
    // doc pr/94 Phase 6 -- per-bundle identity + per-activity cosmic-flag
    // T_tagger/T_kine rows: ONE ROW PER IN-BEAM-WINDOW FLASH BUNDLE instead of
    // one per event, so a cosmic-convicted activity can no longer take a
    // co-bundled neutrino candidate down with it (SBND 18255/395148).
    // **SBND PRODUCTION DEFAULT ON since 2026-08-19 (owner flip after the
    // round-3 Bee scan, bee/pr94r3 -- doc sec 9.13).**  Consumers that read
    // one number per event must go through scripts/pr94_rows.py
    // primary_index(); row 0 is the longest selected activity, i.e. the
    // candidate the pre-pr/94 chain would itself have chosen.  Known open
    // issue carried into production, owner-accepted: an activity can be
    // REPORTED under a bundle it was not matched to, because
    // ClusteringExamineBundles' 80 ns flash-t0 merge overwrites its
    // matched_flash_gid (doc sec 9.8/9.9, mcp2k evt 73038).  Pre-flip arm:
    // SBND_NU_PER_BUNDLE=0.
    nu_per_bundle = true,
    // doc pr/94 Phase 5b round 2 -- the dot guard.  Length floor (cm) for a
    // per-bundle candidate, exempting the legacy event-wide winner (so the
    // row the legacy chain reports can never be floored away).  Without it,
    // per-bundle mode promoted sub-cm blobs: on 1000 mcp1k events round 1
    // added 143 candidates, ALL with a seed under 5 cm and 87 reconstructing
    // to 100-149 MeV -- a dot fitted as a muon at rest; round 2 (this
    // exemption) took that to 0.  KNOWN OPEN ISSUE at this value (doc pr/94
    // sec 9.8, unrelated to the floor): a per-bundle candidate can inherit
    // a bright cosmic-mate's flash gid via ClusteringExamineBundles' flash-t0
    // merge (80 ns window, no spatial check) despite being independently
    // matched to a DIFFERENT, weaker flash of its own -- confirmed on mcp2k
    // evt 73038.  STILL NOT FIXED, and knowingly carried into production by
    // the 2026-08-19 flip (owner decision after the round-3 Bee scan: the
    // display half was fixed, the bookkeeping half was declined -- doc sec
    // 9.9/9.13).  null => C++ default 0 = no floor.
    nu_per_bundle_min_length = 15,
    // doc pr/94 round 3.  Give the SELECTED neutrino candidate the
    // main-cluster PR treatment for the duration of its own pass, even when it
    // is a demoted main.  The PR chain reads main-ness from Flags::main_cluster
    // (NeutrinoPatternBase.cxx:2797, NeutrinoVertexFinder.cxx:3450,
    // NeutrinoTrackShowerSep.cxx:2013), which ClusteringUnmergeBundle clears on
    // a demoted main -- so today a fallback-selected or per-bundle candidate
    // silently loses examine_vertices_3, improve_vertex +
    // fix_maps_shower_in_track_out, main_cluster_initial_pair_vertices,
    // break_two_end_dqdx and the main-branch endpoint ordering.  NOT gated on
    // nu_per_bundle -- the legacy fallback path has the same defect.
    // **SBND PRODUCTION DEFAULT ON since 2026-08-19 (owner flip, doc sec
    // 9.13).**  C++ default false.  Pre-flip arm: SBND_NU_SELECTED_AS_MAIN=0.
    nu_selected_as_main = true,
    // sbnd_xin/docs/75 ON (owner 2026-08-20, "things are good, turn them on
    // for SBND production"): closes a gap in nu_selected_as_main's own
    // guard -- the DL/SCN vertex path can move Flags::main_cluster onto a
    // DIFFERENT cluster in the candidate's bundle mid-pass, which the
    // narrow guard's restore does not undo.  Snapshots/restores the whole
    // {main_cluster} u other_clusters set instead.  OFF gate byte-identical.
    // A flip-equivalence check (doc sec 9) found this fires much more often
    // than the original per-round census sampled: 16/143 standard-sample
    // events + 2/26 enriched-manifest events, NOT confined to promoted
    // candidates as first assumed.  Of those 16, 15 correct only the
    // persisted main-cluster flag (zero vertex/score change); one (evt
    // 37112) moves nue_score by 0.03 on an already-nonselecting score, a
    // known chronically boundary-sensitive event.  No ADVERSE (doc 75 sec
    // 9.3).  =false restores the narrow-guard-only path.
    nu_selected_as_main_snapshot_all = true,
    // ---- doc sbnd_xin/docs/pr/33 sec 11 EM-shower-clustering knobs, ALL ON
    // (owner 2026-08-05; see the sbnd clus.jsonnet clus_pr arg comments).
    // Gate labels: work-pr33-base48 (clean-HEAD binary) vs work-pr33-off48
    // 48/48 everything; per-knob arms vs off48; work-pr33-allon48
    // (env-forced 8) == work-pr33-prod48 (bare).  NOTE: baselines predating
    // 2026-08-05 (work-pr35/36-*) are NOT comparable -- they were produced
    // by a stale-object binary (doc pr/33 sec 11.2).  Escape hatches:
    // SBND_DAUGHTER_COUNT_PROTO_MAIN_VERTEX=0 etc (run_pr_chain_batch.sh).
    //
    // F1a/F1b ON: prototype calculate_num_daughter_tracks callee at both
    // sites.  Measured NULL 48/48 (main-vertex site: 3 calls, 2 events
    // value-differ, proton-skip verdict flips 0; examine site: 40 calls,
    // 2 value-differ, no cut crossed).
    daughter_count_proto_main_vertex = true,
    daughter_count_proto_examine_showers = true,
    // F2a/F2b/F2c ON: read the PDG off the prototype's object (4 sites
    // start-segment, 1 inverted site shower-type, 2 sites exact ==13).
    // Each measured NULL 48/48; all counter disagreements live in evt
    // 137238 (sites :170 x1, :525 x5, :1247 x1) and never cross a decision.
    shower_pdg_from_start_segment = true,
    shower_pdg_from_shower_type = true,
    shower_pdg_exact_muon_test = true,
    // F3 ON: shared pi0-id allocation stream across the two finders.
    // Measured NULL 48/48 (finders fire 10-with/1-without, overlap 0;
    // ssmsp_* branches stable = the sec 10.10 scoping holds).
    pi0_id_shared_allocator = true,
    // F4 ON: is_shower gains the prototype's abs(pdg)==11 disjunct at the
    // center-point site (147 firings / 35 events).  Moves 17/48 events in
    // shower-derived channels only (T_kine pio block, shw_sp_*/mgo_*/sig_*,
    // nue_score on 4 events, Bee mc.json + shower_track); NO nusel verdict
    // flips; pctree + steiner + proj stable 48/48.
    shower_flag_pdg_electron = true,
    // F5 ON: shower_less same-index tie-break by stable shower id (house
    // rule, prototype n/a).  The fallback IS reachable (2 hits, evt 235435,
    // the one id_pi0_without_vertex event) but knob-on is byte-identical
    // 48/48 -- the counter stays as the doc-53-style tripwire.
    shower_less_id_tiebreak = true,
    // doc pr/39: exclude a shower's own start vertex from the end_point
    // farthest-vertex search (prototype map_vtx_segs parity).  SBND
    // production default since 2026-08-06 (owner: "turn it on for SBND"
    // after reviewing the recovery numbers, 32/83 -> 0/83 reversed showers
    // across the 9-event NCpi0 gate) -- exactly the v3_extension_guard
    // idiom (doc pr/24 sec 19.1).  No C++ change --
    // m_shower_endpoint_exclude_start_vertex{false} is still the library
    // default; only the SBND operating point flips it on.
    shower_endpoint_exclude_start_vertex = true,
    // doc pr/91 round 1 F1: the same end_point search must also skip a
    // node NO member segment of the shower touches.  set_start_vertex()
    // calls add_vertex(), so a conn-2/3 shower's view carries a foreign
    // cluster's vertex; the exclusion above hides it only while that
    // shower owns it, and Shower::add_shower imports it into an absorber
    // where it then wins the farthest-vertex search.  Measured on
    // 169626/174752/347129/394532 with shower_dedup_start_seg ON: 6 orphan
    // imports, 5 wrong end points (394532's 30 MeV and 66 MeV showers end
    // on each other).  C++ default false.
    //
    // SBND PRODUCTION ON since 2026-08-18 -- owner flip after scanning the
    // pr/91 sec 2b Bee pair ("This is good, you can flip is on this for SBND
    // production for now").  before e81dfbf9-3801-46a4-ad61-31b5511127f1 /
    // after 04107fda-6306-44e0-aa8c-12a352e32235.  Gates on the 24-event
    // pr/84 r3 manifest: knob-off vs shipped production PASS=48/48; knob-on
    // moves 5 events and ONLY their mabc-pr.zip member data/0/0-mc.json --
    // every pctree tarball byte-identical, nusel tsv byte-identical (zero
    // score/label movers), kine_reco_Enu and kine_energy_particle unchanged.
    // Six showers move and `end` is the only field that changes on any of
    // them.  No C++ change -- m_shower_endpoint_skip_orphan_vtx{false} is
    // still the library default; only the SBND operating point flips it on.
    //
    // Footprint caveat kept on the record: the 24 events were chosen for this
    // defect, and nueCC 168596 -- the ONE mover we did not hand-pick, with
    // zero shower_dedup_start_seg firings, its orphan left behind by
    // id_pi0_with_vertex's pi0 re-seat -- shows the population is not the
    // dedup's.  1/20 non-hand-picked events moved; that is a footprint signal,
    // not a rate.  A standard-manifest gate is still owed.
    shower_endpoint_skip_orphan_vtx = true,
    // doc sbnd_xin/docs/pr/91 round 3: Shower::complete_structure_with_
    // start_segment's flood-fill frontier test switches from view MEMBERSHIP
    // (!has_node) to VISITED (Shower::m_walked_nodes, the prototype's
    // map_vtx_segs equivalent).  Fixes the mechanism behind SBND nueCC
    // 168596's spurious pi0: the 2039 MeV electron's former start vertex
    // 14027 walls off a 4.74 cm proton stub and, past it, a 7.7 cm electron
    // stub 96% inside this shower's own point cloud, left as a separate
    // shower that pairs into a fake pi0.  C++ default false.
    //
    // SBND PRODUCTION ON since 2026-08-18 -- owner-authorized round-3 flip
    // after the validation gate: 67-event manifest (nueCC48 + NCpi0 19),
    // fresh binary, knob-off vs knob-on.  V1 (knob off vs a genuine
    // pre-change build, git-stashed and rebuilt) PASS on evt 168596 +
    // 2 spot events, 6/6 archives byte-identical.  V2 (knob-on vs knob-off,
    // all 67 events): nuecc48 PASS=95/96 -- the ONE mover is 168596's
    // mabc-pr.zip, nothing else; ncpi0 (the pi0-veto sample) PASS=38/38,
    // zero movers -- no genuine two-gamma pi0 was disturbed. Both samples'
    // merged nusel-table.tsv are byte-identical (zero score/label movers,
    // zero nu-candidate status changes) and every pctree tarball matches.
    // 168596 mechanism, confirmed by probe (SHOWER_WALK_DEBUG): the ONLY
    // re-expansion across all 67 events is shower_id=2 re-walking vertex
    // 14027 via segment 14093 -- exactly the round-2 diagnosis.  Shower
    // count 17->16 (the 7.7cm electron stub is absorbed, not orphaned);
    // pio_mass 114.06 -> no pi0 pairs (all -1.0 sentinel); kine_reco_Enu
    // 2331.29 -> 2324.26 MeV (-0.3%, the borrowed-charge correction).
    // See sbnd_xin/docs/pr/91_em-shower-clustering-round1.md sec 11.
    shower_walk_visited_parity = true,
    // doc sbnd_xin/docs/pr/40 sec 17 (2026-08-06): track (proton/pion/muon)
    // mis-identified as electron.  F1 restores prototype-faithful PID
    // persistence (segment_determine_dir_track); F2/F3 spare a segment
    // whose own median dQ/dx is decisively proton- or muon-like from an
    // unconditional track-to-electron reclassification, at the wholesale-
    // conversion sites (F2) and the shower-topology test (F3).  C++
    // defaults all false.  SBND PRODUCTION DEFAULT since 2026-08-06 (owner:
    // "Flip to SBND ON if gates pass" -- gates G1-G7 all passed): gate
    // work-pr40-off48 vs work-pr40-base48 (clean-source reference) 48/48
    // byte-identical; work-pr40-on48 fixes 8 of the 9 owner-reported cases
    // to a hadron pdg (388->13, 74544/267597/269774/423981/489330->2212/
    // 211, 433451->13), the 9th (evt 256587) correctly left untouched
    // because its own median dQ/dx (1.26x MIP) sits in the deliberately
    // conservative gap between the muon and proton thresholds, not a bug;
    // population moves 3/48 events beyond the 9, all confined to the
    // nusel-table stmfit field, zero label (verdict) flips; wcdoctest-clus
    // 97/97.
    track_pid_persist_dqdx = true,
    shower_reclass_dqdx_guard = true,
    shower_topo_dqdx_guard = true,
    // doc sbnd_xin/docs/pr/40 round 2: two follow-on defects measured on the Bee
    // display of the pr/40 fix above.  F4 track_pid_persist_4mom: the
    // non-free-end persistence branch stored a rest-mass-only 4-momentum
    // (E=mass, zero momentum), so Aux::ParticleInfo::kinetic_energy()
    // (=E-mass) read exactly ZERO for every segment it rescued (SBND evt
    // 174637 seg 9050, Bee PF node "mu- 0 MeV"); segment_cal_4mom has no
    // free-end dependence, so it is now called unconditionally instead.
    // SBND PRODUCTION DEFAULT ON: G1 48/48 byte-identical, G2a fixes the
    // reported case (0 -> 86 MeV), G4 census 1 -> 0 zero-energy PF nodes.
    // F6 reclass_never_computed_ke_floor: a related reclass_pinfo
    // negative-KE construction (KE = -mass on the never-computed path),
    // fixed alongside.  SBND PRODUCTION DEFAULT ON: verified by direct unit
    // test (the negative-KE precondition was not observed on the 48-event
    // census either arm, so this flip rests on the doctest, not population
    // evidence -- see doc pr/40 round 2).
    // F5 shower_proton_daughter_pion: an electron cannot father a proton --
    // set_default_shower_particle_info relabels pion (211) instead of
    // electron when the candidate segment emanates from the neutrino
    // vertex and its far end is a PID'd, charge-confirmed proton (SBND evt
    // 256587 seg 11079, daughter segment 11080 at 3.72x MIP).  Round 2 found
    // the override was reverted end-to-end by a downstream writer
    // (Shower::update_particle_type, PRShower.cxx); round 3 closed that with
    // a matching guard threaded into update_particle_type itself.  SBND
    // PRODUCTION DEFAULT ON (doc pr/40 round 3): G1 48/48 byte-identical,
    // evt 256587 seg 11079 now 211 end-to-end, population census 2/2209
    // segments move, zero nusel-table.tsv verdict diff.  See
    // porting_dictionary.md and doc pr/40 round 3.
    track_pid_persist_4mom = true,
    shower_proton_daughter_pion = true,
    reclass_never_computed_ke_floor = true,
    // doc sbnd_xin/docs/pr/40 round 4 -- two follow-on defects measured on the
    // Bee display of round 2/3's F5 fix.  F7 shower_proton_daughter_pion_dissolve:
    // F5 relabelled the pdg but left the shower flags set, so the pion was
    // still wrapped as a Shower (its proton daughter was pre-claimed into
    // the shower's segment set and never got its own particle-flow node;
    // the pi+ Bee node's displayed endpoint was the shower's, not the
    // segment's own).  F8 muon_multi_proton_pion: a muon segment whose far
    // (non-neutrino-vertex) end is a multi-proton (>=2, charge-confirmed)
    // hadronic vertex is relabelled pion; no propagation across a degree-2
    // kink to the sibling muon segment further from the vertex.
    // SBND PRODUCTION DEFAULT ON (doc pr/40 round 4): G1 48/48 events byte-
    // identical (96/96 archives) vs a git-stash clean-HEAD reference; G2a/G2b
    // both owner-reported cases fixed exactly as specified; G3 nusel-table.tsv
    // diff 0 lines across the 48-event manifest; G4 census exactly 5 segments
    // move, all attributed; G5 wcdoctest-clus 99/99 test cases.
    shower_proton_daughter_pion_dissolve = true,
    muon_multi_proton_pion = true,
    // doc sbnd_xin/docs/pr/40 round 5 -- muon mis-identified as electron,
    // three owner-reported Bee cases (84229, 54341, 55715), three
    // independent mechanisms.  F9 track_pid_persist_dqdx_electron_guard:
    // narrows F1 so it no longer rescues an undirected (no free end)
    // electron guess that would otherwise poison a NEIGHBORING segment's
    // flag_shower_in test.  F10 shower_connect_main_vertex_straight_guard:
    // the main-vertex EM-shower selection heuristic excludes a candidate
    // whose own start segment is long and straight.  F11
    // shower_traj_straight_guard: segment_is_shower_trajectory gets the
    // same straightness veto F3 already gave segment_is_shower_topology's
    // dQ/dx.
    // SBND PRODUCTION DEFAULT ON (doc pr/40 round 6, together with F12/F14
    // below): round 5 measured these three alone do NOT reach the displayed
    // outcome (G2 failed, F11 regressed evt 55715 seg 15005); round 6's
    // boundary-level F12/F14 close that gap, and the round-6 gates cover the
    // five as a set -- G1 48/48 events (96/96 archives) byte-identical off,
    // G2 all three owner cases fixed, G3 nusel-table.tsv 0-line diff, G4
    // census attributed (26/48 events restructure, 42 pdg moves dominated by
    // the intended 11->13 recoveries, long-muon reassembly intact 2==2).
    track_pid_persist_dqdx_electron_guard = true,
    shower_connect_main_vertex_straight_guard = true,
    shower_traj_straight_guard = true,
    // doc sbnd_xin/docs/pr/40 round 6 -- the boundary-level fixes round 5's
    // G2 measurement demanded.  F12 shower_absorb_track_guard: the shower
    // flood-fill no longer absorbs a confidently PID'd straight non-electron
    // track (long-muon pseudo-showers exempt).  F14
    // michel_stem_muon_rescue: the Michel stopping-muon rescue reaches a
    // confident-direction proton-called stem at a multi-prong stopping
    // vertex.  SBND PRODUCTION DEFAULT ON (same gate set as above).
    shower_absorb_track_guard = true,
    // F13 shower_connect_protected_pion_guard: MEASURED DEAD as shaped
    // (its motivating segment is already pdg 2212 by candidate-selection
    // time -- see NeutrinoPatternBase.h), kept as a documented negative
    // result, NEVER flip.
    shower_connect_protected_pion_guard = false,
    michel_stem_muon_rescue = true,
    // doc sbnd_xin/docs/pr/74 round 2 P1: examine_direction's
    // flag_shower_in cascade relabels a downstream |pdg|==13/pdg==0 segment
    // electron with no length ceiling and no charge test (prototype has the
    // identical unconditional branch, NeutrinoID_track_shower.h:2004).
    // When on, refuse the relabel for a segment BOTH long
    // (> shower_in_max_len, C++ 40 cm) AND MIP-like (median dQ/dx <
    // shower_in_mip_hi x MIP median, C++ 1.3).  C++ default false.  Key
    // omitted when off => byte-identical pre-pr/74 config.
    // SBND PRODUCTION DEFAULT ON since 2026-08-13 (doc pr/74 round 2:
    // off-gate 0/117, on-census 7/117 movers all attributed, nusel 0/117).
    shower_in_cascade_guard = true,
    shower_in_max_len = null,
    shower_in_mip_hi = null,
    // doc sbnd_xin/docs/pr/40 round 9 -- the rounds-7+8 straight-track PID
    // guard family + the B2 cross-cluster bridge (SBND 286906/409546 +
    // 54629/320865 classes; must-not-touch 521075).  Five guards decline an
    // unconditional pdg-11 write on a straight long track
    // (segment_is_straight_long_track, continuation-aware at the two
    // shower-connect sites); the bridge replaces the cross-cluster conn-2
    // electron fabrication with a real 2-point graph bridge when the
    // steiner-cloud gap is < shower_nv_bridge_max_gap (C++ 1.8 cm).
    // C++ defaults ALL false (scalars 25 deg / 1.8 cm live in C++).  Keys
    // omitted when off/null => byte-identical pre-round-9 config.
    // SBND PRODUCTION DEFAULT ON since 2026-08-18 (doc pr/40 round 9:
    // V1 pre-change-binary gate PASS 6/6 archives; V2 OFF-vs-ON over
    // 45-Bee + nueCC48 + NCpi0-19 (100 events): 30 mabc-only movers, 0
    // pctree movers, nusel tables byte-identical 3/3 samples; census
    // straight-long e- 63 -> 42; bridges fired on 286906/409546/407280
    // only, 521075 untouched; V3 bare-run composition gate PASS).
    shower_connect_from_vertices_straight_guard = true,
    shower_connect_start_seg_straight_guard = true,
    examine_direction_dirsign_shower_in_guard = true,
    daughter_shower_angle_reclass_straight_guard = true,
    shower_topo_reexam_straight_guard = true,
    sfv_kink_max = null,
    shower_nv_bridge_track = true,
    shower_nv_bridge_max_gap = null,
    // doc sbnd_xin/docs/pr/97 D1 -- deterministic main_pi sentinel.  OFF
    // pending an owner flip: turning it on changes which vertex an
    // other-cluster shower attaches to on every event where the overall main
    // vertex lives outside main_cluster.  SBND_MAIN_PI_INIT=true flips it.
    shower_nv_main_pi_init = false,
    // doc sbnd_xin/docs/pr/92 -- drop stray satellite showers (overclustered
    // cosmics / second neutrinos) from kine_reco_Enu and (via
    // pf_drop_stray_satellites below) the Bee PF tree.  Candidates: conn-2/3
    // showers in NON-main clusters above 20 MeV, not pi0-paired, not within
    // 8 cm of a main-cluster attachment; dropped when the fresh shower axis
    // is > 60 deg off the attachment vertex (Arm A), the attachment is
    // > 90 cm away or outside the main cluster AND the axis is >= 45 deg
    // off the main vertex (Arm B), or the start segment is the collinear
    // continuation of an out-of-shower straight long track (Arm C).
    // C++ defaults false (scalars 20 MeV / 8 cm / 60 deg / 45 deg / 90 cm /
    // 30 cm / 25 deg live in C++).  Keys omitted when off/null =>
    // byte-identical pre-pr/92 config.  Runner env:
    // SBND_KINE_DROP_STRAY_SATELLITES + SBND_KINE_SAT_* scalars.
    // SBND PRODUCTION DEFAULT ON since 2026-08-18 (doc pr/92: V1 pre/post
    // binary gate PASS 8/8 archives; probe round over mcp1k-33 + nueCC48 +
    // NCpi0-19 (100 events, 1202 satellite candidates): 51 drops in 30
    // events, all 8 target offenders dropped, every must-keep sentinel
    // retained; V2 OFF-vs-ON: movers = exactly those 30 events, mabc-only,
    // nusel selections byte-identical 3/3 samples; 350935 -449 MeV cosmic
    // gamma gone, 321371 -98 MeV cosmic neutron/mu gone, 389538 -1362 MeV
    // second neutrino incl. the 955 MeV neutron->proton gone).
    kine_drop_stray_satellites = true,
    kine_sat_min_energy = null,
    kine_sat_prox_max = null,
    kine_sat_angle_bad = null,
    kine_sat_angle_main = null,
    kine_sat_far_dis = null,
    kine_sat_axis_dis_cut = null,
    kine_sat_cont_kink = null,
    // pr/92 round 2 (owner retune): direction arms restricted to TRACK-like
    // satellites; EM-shower-like satellites drop only when > em_far_dis
    // from the main vertex with a folded main-vertex angle failure.
    // C++ defaults 3 / 150 cm.
    kine_sat_track_max_nseg = null,
    kine_sat_em_far_dis = null,
    // doc sbnd_xin/docs/pr/74 round 2 P2: the F14 Michel rescue accepts ANY
    // shower-like sibling at the stem's far vertex; on a nueCC event that
    // sibling is the EM shower trunk and the rescue paints a muon at the
    // neutrino vertex (18255-90055).  When on, additionally require the
    // graph beyond the far vertex to be Michel-sized (total track length <
    // michel_stem_max_far_len, C++ 40 cm).  C++ default false.  Key omitted
    // when off => byte-identical pre-pr/74 config.
    // SBND PRODUCTION DEFAULT ON since 2026-08-13 (same gate set).
    michel_stem_michel_check = true,
    michel_stem_max_far_len = null,
    // doc sbnd_xin/docs/pr/74 round 2 K4: shower formation walks outward
    // from the main vertex and starts a shower at the first shower-like
    // segment; the track stem it walked PAST is structurally excluded
    // (18255-90055 trunk 11045; 18255-469665 chain 15003/15001; both
    // prototype-shared gaps).  When on, a post-pass absorbs the chain from
    // each substantial EM shower's attach vertex back toward the main
    // vertex while segments are short (< stem_backfill_max_len, C++ 30 cm)
    // and not charge-hot (median dQ/dx < stem_backfill_mip_hi x MIP, C++
    // 3.5 -- a Bragg proton stops the walk).  C++ defaults false/30/3.5/40
    // (+ stem_backfill_mip_lo 0.75 and a both-endpoint junction guard, the
    // round-2 iteration that closed the stranded-PF-orphan class).
    // Keys omitted when off/null => byte-identical pre-pr/74 config.
    // SBND PRODUCTION DEFAULT ON since 2026-08-13 (same gate set).
    shower_stem_backfill = true,
    stem_backfill_max_len = null,
    stem_backfill_mip_lo = null,
    stem_backfill_mip_hi = null,
    stem_backfill_min_shower_len = null,
    // doc sbnd_xin/docs/pr/74 round 2 K5 = pr/65's deferred rung 2: promote
    // a graph-unreachable, unclaimed main-cluster segment (18306-142421 seg
    // 7013, 41.9 cm / 266 MeV, PF-invisible today) through the prototype's
    // own connection_type=3 pseudo-gamma path.  C++ defaults false/10 cm.
    // Keys omitted when off/null => byte-identical pre-pr/74 config.
    // SBND PRODUCTION DEFAULT ON since 2026-08-13 (same gate set; closes
    // pr/65's "0 unclaimed" gap on NC-pi0).
    shower_conn3_unreachable = true,
    conn3_unreachable_min_len = null,
    // doc pr/84 round 2 (F3) -- bridge disconnected main-cluster components
    // whose closest approach to the reachable side is within this radius (cm)
    // BEFORE clustering_points, so they classify conn-1 naturally;
    // shower_conn3_unreachable above stays the backstop for wider gaps.
    // SBND PRODUCTION ON at 1 (cm) 2026-08-17: 10 fires / 512 events, all
    // sub-cm-family, 9 score movers all target-family (283713 enu 1513->2034
    // = the stranded 567 MeV muon rejoins; 66272 rescues an invisible
    // pi+ -> proton branch; nueCC 168596 consolidates its 1929 MeV primary,
    // 30->17 nodes), ZERO nu-candidate label flips, ncpi0 untouched.
    // Deliberately NOT the swept 3 cm: at 1-3 cm the bridges turn
    // speculative -- nueCC 38856's two 2.5-2.9 cm bridges fragmented its
    // 1244 MeV electron and flipped nue 3.25 -> -3.45 (recorded in doc
    // pr/84 sec 13 as the tuning evidence).  Runner env:
    // SBND_CONN3_STITCH_MAX.
    conn3_stitch_max = 1,
    // doc pr/84 round 3 (S1) -- one shower per start segment.  Two
    // PR::Showers can be built on the SAME start segment (attributed with
    // WCT_SHOWER_CREATE_DEBUG: shower_clustering_in_other_clusters picks a
    // start segment with no claim check at all, and the K5 conn3_unreachable
    // branch reads a map_segment_in_shower that update_shower_maps only
    // refreshes at the end of the function).  The twin renders a duplicate PF
    // node AND is counted twice in kine_energy_particle -- SBND 394532
    // kine_reco_Enu 352.2 MeV vs 255.5 de-duplicated.  Changes physics
    // output, so it is a knob: the group collapses onto its most directly
    // connected member, which absorbs the others' segments, and only that
    // shower's kinematics are recomputed.  Runner env:
    // SBND_SHOWER_DEDUP_START_SEG.  C++ default stays false.
    // SBND PRODUCTION ON 2026-08-18 (owner flip: "this is a clear bug ...
    // aim to fix the underlying issue").  Gate manifest = the 24 round-2 Bee
    // events: knob-off byte-identical to production (pr83r3_hash_gate
    // PASS=48/48), knob-on moves EXACTLY the 4 twin events and inside them
    // only mabc-pr.zip::data/0/0-mc.json plus the kine/tagger blocks -- every
    // pctree tarball, every other event and every nusel label byte-identical.
    // 6 absorptions: 169626 (Enu 825.4->737.3), 174752 (188.7->176.0),
    // 347129 (700.8->571.9), 394532 (352.2->248.2); no nu-candidate flips.
    // The 492-event round-2 census is what bounds the population: exactly
    // these 4 events carry twins.  Legacy escape:
    // -A shower_dedup_start_seg=false.
    shower_dedup_start_seg = true,
    // doc sbnd_xin/docs/pr/74 round 4 K6 shower_traj_michel_stem: a stopping
    // muon that emits a Michel electron at the neutrino vertex is
    // reconstructed as ONE EM shower, because track/shower separation flags
    // the muon kShowerTrajectory on its wiggliness and the ordinary track PID
    // is never consulted (and, measured, ABSTAINS when it is).  18255-506746
    // seg 21048: 20.4 cm at 1.57x MIP, far vertex degree 2, 17.5 cm terminal
    // downstream, 93 deg kink -> reported by the owner as "107 MeV electron
    // is muon + Michel".  When on, that one segment is demoted to a stopping
    // muon so the Michel becomes its own shower.  C++ defaults
    // false/15cm/45cm/1.3x/40cm/40deg.  Keys omitted when off/null =>
    // byte-identical pre-round-4 config.
    // SBND PRODUCTION DEFAULT ON since 2026-08-13 (doc pr/74 round 4 gates:
    // off-gate 0/117 archives + 0/117 pctree member hashes; knob footprint
    // 1/117 events, 506746 only, 0-mc.json + the paint layer; nusel 0/117;
    // dangling PF roots gained 0/117; nu-vertex 0/117; flip and escape gates
    // byte-exact both ways).
    shower_traj_michel_stem = true,
    michel_stem_traj_min_len = null,
    michel_stem_traj_max_len = null,
    michel_stem_traj_mip_lo = null,
    michel_stem_traj_max_far_len = null,
    michel_stem_traj_min_kink_deg = null,
    // doc pr/44 shower_long_muon_keep_type: a MULTI-segment long-muon
    // pseudo-shower (cached type 13 at the in_main_cluster seed) keeps its
    // muon start segment -- the update_particle_type majority vote there is
    // a toolkit-only addition (18f09178) absent from the prototype, and it
    // relabels a pure muon chain's start segment 13 -> 11 (SBND
    // 18255-142421: fake "e- 163 MeV" paired into the pi0).  C++ default
    // false.  SBND PRODUCTION DEFAULT ON (same gate set as
    // pf_orphan_track_parentage above; fires on NO nueCC48 event, on ncpi0
    // only 142421 where it restores the owner-truth three-prong vertex --
    // doc pr/44).  Runner env: SBND_SHOWER_LONG_MUON_KEEP_TYPE.
    shower_long_muon_keep_type = true,

    // doc pr/40 round 10 shower_bragg_protect_start_segment: spares a
    // segment from examine_all_showers' cluster-wide "every non-shower
    // segment here becomes electron" reclassification
    // (NeutrinoTrackShowerSep.cxx) when it is longer than 20 cm, carries a
    // confident (<1.0) Bragg/dE-dx-template PID score -- the population
    // segment_dqdx_spares_electron_reclass's flat median-dQ/dx ratio test
    // (F2, doc pr/40) cannot reach (ratio in [1.2,1.75], the gap between
    // its "clean MIP" and "proton-like" spares) -- AND sits in the MAIN
    // interaction cluster.  C++ default false.  SBND PRODUCTION DEFAULT ON
    // (2026-08-18): SBND 18255-314507 seg 17002 (32.3 cm, xMIP 1.57x)
    // restored to muon (was mislabelled e- 151 MeV).  The is_main_cluster
    // restriction was added same-day after owner Bee review of the
    // initial fix's second mover (18255-259542) identified it as a
    // genuine photon by topology -- its spared segment sat in a separate
    // SATELLITE cluster (124, disjoint from the main interaction), where a
    // good Bragg score is not reliable evidence (a photon's early
    // conversion stem can score well against the muon template before the
    // cascade visibly multiplies); satellite clusters already have their
    // own dedicated classifier (kine_drop_stray_satellites, doc pr/92)
    // that correctly kept 259542 as EM.  Population census (nueCC48 48/48
    // + ncpi0 19/19 + 31-event mcp1k electron-misID subset), post-
    // restriction: ZERO nueCC48 movers, ZERO ncpi0 movers, ONE mcp1k mover
    // (314507, the motivating case).  Runner env:
    // SBND_SHOWER_BRAGG_PROTECT_START_SEGMENT.
    shower_bragg_protect_start_segment = true,

    // doc pr/93 round 3 -- four knobs for the "electron is really tracks /
    // hadronic-pi0 shower" family (SBND 18255-55595/348471/69314/292643/
    // 315167).  C++ defaults false; key-suppressed when off => compiled
    // config byte-identical.  SBND PRODUCTION ON (owner round 2026-08-18):
    // rescues 4/5 owner events (55595 fake 458 MeV e- gone; 348471 ->
    // p 719 MeV; 292643 -> pi 162 MeV; 315167 1046.7 -> 164.3 MeV EM stub
    // with the 150.7cm score-0.10 proton kept separate and kept proton;
    // 69314 named residual -- its pion stamp carries the unscored-100
    // sentinel at 38.4cm, below any floor that does not regress real
    // electrons).  100-event validation (mcp1k-33 + nueCC48 + NCpi0):
    // nueCC48 movers 2/48, BOTH documented-defect/adjudicated (137238
    // pr/74 shape-B pencil excluded; 46363 13.5cm conn-2 satellite);
    // NCpi0 movers 1/19 (285567, pr/74 shape-A fake e- -> proton, fix
    // direction); mcp1k movers 8/33, ALL in the pr/40 mis-ID census, all
    // rescue-direction.  Manifest straight-long pdg-11 census rows
    // 41 -> 35.  The shared 50cm floor (shower_pid_guard_min_len, C++
    // default) is REQUIRED: un-floored attribution arms regressed 9+17/48
    // nueCC48 events (real 22-47cm electron stems carry confident
    // 0.11-0.64 proton/muon template scores).  Runner envs
    // SBND_SHOWER_RECLASS_CASE_B_DQDX_GUARD / SBND_SHOWER_ACCEPT_PID_GUARD
    // / SBND_SHOWER_VOTE_TRACK_PID_COUNTS / SBND_SHOWER_CONE_ABSORB_GUARD
    // / SBND_SHOWER_PID_GUARD_MIN_LEN.
    shower_reclass_case_b_dqdx_guard = true,
    shower_accept_pid_guard = true,
    shower_pid_guard_min_len = null,  // null => C++ default 50cm
    shower_vote_track_pid_counts = true,
    shower_cone_absorb_guard = true,

    // doc pr/93 round 4 -- PF-hierarchy fine-tunes + 137238 cross-cluster
    // muon.  C++ defaults false/null; key-suppressed when off => compiled
    // config byte-identical.  SBND PRODUCTION ON (owner round 2026-08-18):
    //   shower_detach_track_stem: peel the main-cluster track prefix off a
    //     track-headed shower and re-root the EM remainder at the prefix's
    //     far vertex, conn 2 (348471: "proton 308 -> pi0 113 (g355+g20) +
    //     g74 + g11", Enu 2075.7 -> 1090.8 = the round-3 proton-mass-on-
    //     charge-aggregate regression repaired; 292643: "pi+ 88 -> mu- 58 ->
    //     4 gammas", Enu 950.4 -> 1073.6).
    //   kine_count_orphan_tracks + pf_orphan_confident_track (below):
    //     count/emit confident straight-long main-cluster orphan tracks
    //     freed by shower_cone_absorb_guard (315167: PF gains root
    //     "proton 595 MeV", Enu 722.1 -> 1326.0 = +595.3 KE +8.6 binding).
    //   straight_cont_cross_cluster + sccc_bridge_body: demote a main-vertex
    //     shower-trajectory stem that is the cross-cluster continuation of a
    //     straight long track across a pr/57 W-gap split, and bridge the
    //     body into the PF/kine chain (137238: "e- 152 MeV" -> "mu- 60 ->
    //     bridge -> mu- 211 + mu- 65" with delta rays as small EM leaves,
    //     Enu 1087.1 -> 1101.4).  Owner's angle-conditioned tiers; base
    //     retuned 5->6cm / 15->18deg (137238's body measures g=5.68cm
    //     K=17.0deg in the fitted tangents); aligned tier 12cm/7.5deg at
    //     C++ defaults.  The demoted stem + bridged cluster are shielded
    //     from pass-2 seeding, from_vertices Step-3 analysis, and the
    //     examine_showers retarget (all sccc-scoped sets).
    //   pf_track_owns_loose_vertex: a vertex the track BFS walked a real
    //     segment to is not claimable by a root shower whose only tie is
    //     the loose fill_sets view (69314: the muon's 67 MeV e- + 18 MeV
    //     gamma chain re-parent from the 595 MeV shower to the muon;
    //     render-only, Enu byte-unchanged).
    // Validation 2026-08-18 (work-pr93r4-{off2,on2}-{mcp1k,nuecc48,ncpi0} +
    // ctrl pair): OFF gate PASS 200/200 vs round-3 production; movers
    // exactly the 4 targets + 447477 (carrier renumber, nil) + 137238;
    // NCpi0 0/19 touched; pr/57 negative controls (61579/55715 byte-
    // identical; 21073/84229/122660 untouched) clean; runtime/RSS
    // unchanged.  Runner envs SBND_SHOWER_DETACH_TRACK_STEM /
    // SBND_KINE_COUNT_ORPHAN_TRACKS / SBND_STRAIGHT_CONT_CROSS_CLUSTER /
    // SBND_SCCC_BRIDGE_BODY / SBND_PF_ORPHAN_CONFIDENT_TRACK /
    // SBND_PF_TRACK_OWNS_LOOSE_VERTEX (+ SBND_SCCC_* numerics).
    shower_detach_track_stem = true,
    // doc pr/99 round 2 -- projective-ghost member drop inside shower
    // membership (395148).  C++ defaults false/0.7/0.25/10cm = legacy.
    shower_ghost_member_drop = true,  // SBND PRODUCTION ON 2026-08-20 (doc pr/99 round 2)
    shower_ghost_overlap_frac = null,
    shower_ghost_dqdx_ratio = null,
    shower_ghost_min_len = null,
    // doc pr/99 round 3 -- C1 kine-charge cell-ownership dedup + C1b
    // prototype cloud-rebuild parity (168596 Enu double count) + A5
    // hadronic-shower re-type (315167/395148 labels).  C++ defaults
    // false/null = legacy.  Key omitted when off => byte-identical
    // pre-round-3 config.
    kine_charge_dedup = true,    // SBND PRODUCTION ON 2026-08-20 (doc pr/99 round 3; owner-adjudicated full winner-take-all)
    kine_charge_rebuild = true,  // SBND PRODUCTION ON 2026-08-20 (doc pr/99 round 3; prototype cloud-rebuild parity)
    // doc pr/101 -- Enu accounting round: K1 shower<->track cell ownership
    // (37112), K2 paper mass/binding rules (proton shower +938 MeV), K3
    // hadronic-shower KE = sum dE/dx, K4 long-muon range (0 dQdx / 1 range /
    // 2 range with dead-end + ratio fallback), K5 main-vertex member guard.
    // C++ defaults false/0/null = legacy.  Key omitted when off =>
    // byte-identical pre-pr/101 config.
    kine_charge_track_ctx = true,    // SBND PRODUCTION ON 2026-08-20 (doc pr/101 K1; gate 234 PASS, 0 nue/numu flips from this knob)
    kine_mass_rules = true,          // SBND PRODUCTION ON 2026-08-20 (doc pr/101 K2; latent on the 117-evt manifest)
    kine_hadronic_dqdx = true,       // SBND PRODUCTION ON 2026-08-20 (doc pr/101 K3; owner: object-level)
    kine_long_muon_mode = 2,         // SBND PRODUCTION ON 2026-08-20 (doc pr/101 K4; range w/ dead-end + ratio fallback)
    kine_long_muon_ratio_lo = null,  // C++ default 0.3
    kine_long_muon_ratio_hi = null,  // C++ default 0.5
    kine_mainvtx_used_guard = true,  // SBND PRODUCTION ON 2026-08-20 (doc pr/101 K5; latent on the manifest)
    shower_hadronic_tag = true,  // SBND PRODUCTION ON 2026-08-20 (doc pr/99 round 3, A5)
    shower_hadronic_min_len = null,   // cm; C++ default 10
    shower_hadronic_scan_len = null,  // cm; C++ default 30
    shower_hadronic_bin = null,       // cm; C++ default 3
    shower_hadronic_r_cyl = null,     // cm; C++ default 8
    shower_hadronic_r_core = null,    // cm; C++ default 1.2
    shower_hadronic_growth_max = 0.7,    // calibrated (109-shower roster; protected primaries min 2.32)
    shower_hadronic_growth_bragg = null, // C++ default 1.2
    shower_hadronic_bragg_ratio = null,  // C++ default 3.0
    shower_hadronic_stem_ratio = 2.8,    // calibrated (reaches 395148 at C++ stem 2.98; gammas ~2)
    kine_count_orphan_tracks = true,
    kine_orphan_track_min = null,  // null => C++ default 50cm
    straight_cont_cross_cluster = true,
    sccc_bridge_body = true,
    sccc_max_gap = 6,          // cm; base tier (C++ default 5; 137238 g=5.68)
    sccc_kink_max = 18,        // deg; base tier (C++ default 15; 137238 K=17.0)
    sccc_gap_aligned = null,   // null => C++ default 12cm (aligned tier)
    sccc_kink_tight = null,    // null => C++ default 7.5deg (aligned tier)

    // doc pr/43 round 2 -- three PID-consistency knobs for the remaining
    // owner cases (18255: 54351 / 56463 / 57661).  K1
    // single_muon_proton_chain_veto: the vertex muon selection's proton veto
    // walks the bounded degree-2 continuation chain (a muon cannot end in a
    // proton), demoting the chain to pion and re-picking.  K2
    // single_muon_long_muon_claim: a long-muon accumulation chain claims the
    // vertex muon slot so a second pdg-13 arm demotes to pion.  K3
    // pid_flag_reconcile: late reconciliation pass (post shower clustering,
    // pre taggers) -- forced-electron terminal rescue behind a main-vertex
    // proton chain + stale shower-flag/wrapper cleanup on confirmed tracks.
    // C++ defaults false.  SBND PRODUCTION ON (2026-08-07): flipped per the
    // owner's "flip if clean" policy after the round-2 gates -- G1 48/48
    // byte-identical off, G2 3/3 owner shapes, G3 zero movement on all four
    // nueCC48 on-arms (K1/K2/K3/all 0/48, nusel 0-diff) and a single fully
    // attributed K2 move on ncpi0-19 (142421 stub 7023 mu- 4 -> pi+ 4 MeV).
    // Gate arms: work-pr43r2-{off48,onK1-48,onK2-48,onK3-48,onall-48,
    // cleanref48,off19n,onall19n}.
    single_muon_proton_chain_veto = true,
    single_muon_long_muon_claim = true,
    pid_flag_reconcile = true,
    // doc pr/45 (18255-56463 follow-ups).  other_seg_empty_2d_guard: in
    // find_other_segments a segment with zero fit points on the queried face
    // returns the -1.0 empty-2D-tree sentinel, which counts as "distance
    // zero" and lets ONE far-TPC segment tag the entire near face as
    // explained -- the 30 cm isochronous tail beyond segment 14006 could
    // never seed a residual component.  pseudo_shower_track_paint: the Bee
    // shower_track layer painted the 411 cm cathode-crossing muon red
    // because paint reads shower MEMBERSHIP while the PF tree reads the
    // shower's cached type (mu-); muon-typed (+-13) pseudo-showers now paint
    // as track.  C++ defaults false.  SBND PRODUCTION ON (2026-08-07):
    // flipped per the owner's standing "flip if clean" policy -- G1 96/96 +
    // 38/38 byte-identical off; guard: nueCC48 0/48 (archive-level),
    // mcp1k-200 3/200 all attributed (cathode-crossing clusters, new/
    // re-formed segments: 276836, 404684, 407280), nusel 0-diff; paint:
    // display-only, 0 PF movers, 2/48 nueCC48 shower_track q flips
    // (137238, 400474), nusel 0-diff.  Gate arms work-pr45-{off48,onA48,
    // onB48,off19n,m200on,m200onA,m200off}.  Doc pr/45.
    other_seg_empty_2d_guard = true,
    pseudo_shower_track_paint = true,
    // doc pr/46 (18255-55595).  long_muon_stub_bridge: the long-muon chain
    // walk's 15 deg junction-angle test fails on a short (< 6 cm) vertex
    // stub's noisy fitted direction, so a broken long muon behind such a
    // stub never enters the long-muon category (55595: 2.4 cm stub at
    // 30.5-35.4 deg vs its 192.9 cm MIP continuation; the muon-slot logic
    // then crowns a 2.7 cm sibling "mu- 48 MeV" and demotes the gateway
    // stub to pi+).  When on, a < 6 cm incoming segment bridges to a
    // > 35 cm MIP (< 1.3) continuation up to 45 deg unless the junction
    // has another track-like arm > 10 cm (multiple substantial outgoing
    // tracks = genuine hadronic vertex; delta rays/fragments do not veto).
    // C++ default false.  SBND PRODUCTION ON (2026-08-07): flipped per the
    // owner's standing "flip if clean" policy -- G1 off byte-identical
    // 96/96 + 38/38 vs clean-source base arms; knob-on ARCHIVE-identical on
    // nueCC48 (96/96) and ncpi0 (38/38); full-1k census 3/1000 movers, all
    // the designed stub-bridge assembly (55595 2.4cm+181.5cm -> mu- 455;
    // 73004 4.4cm+88.8cm -> mu- 226; 349945 3.1cm+172.8cm -> mu- 464,
    // chain traces nsegs=2 all), nusel 0-diff everywhere; must-not-merge
    // controls hold (66118 angle 70-82 > 45 cap; 61461 multiplicity veto;
    // 284145 hadronic vertex).  Gate arms work-pr46-{base48,base19n,off48,
    // off19n,on48,on19n,oncase,m1koffb,m1konb}.  Doc pr/46.
    long_muon_stub_bridge = true,
    // doc pr/48 (18255-51513/56211/57903/59335/57485) -- back-to-back track
    // fixes.  two_end_break: the two-end residual-range break pass (nu
    // vertex mid-segment on a single unbroken track, dQ/dx rising at BOTH
    // ends -- two stopping ends imply a junction at the interior dip; a
    // simple-topology + both-ends-in-FV gated break at the dip / turn-max
    // index, template-score accepted).  kink_walk_dqdx_stop +
    // kink_break_protect: the 59335 fixes (a correct C4 kink accept walks
    // past its own evidence, then examine_vertices_4's unconditional < 2 cm
    // floor erases the break).  C++ defaults false.  SBND DEFAULT ON (owner
    // round, doc pr/48 sec 9): off-gates work-pr48-base48 vs -off48c 48/48 +
    // -base19n vs -off19nc 19/19 member-identical, nusel byte-identical;
    // 1k footprint 69/1000 movers, 0/1000 nusel diffs, every mover classified
    // and every new TEB break examined (sec 9.6).  Pass -A two_end_break=false
    // etc. (or the SBND_TWO_END_BREAK / SBND_KINK_WALK_DQDX_STOP /
    // SBND_KINK_BREAK_PROTECT runner envs) to restore the legacy path.
    // teb_* operating point rides the C++ defaults.
    two_end_break = true,
    // doc pr/90 round 2 -- two refinements of the two_end_break pass, both
    // C++ default 0 = legacy.  teb_turn_min_arm_frac: route R2's turn argmax
    // PREFERS an index whose PCA arms can each span this fraction of
    // teb_turn_baseline when one clears teb_turn_angle on its own; the
    // legacy argmax (starved near-end candidates included) is the fallback
    // (320865: a starved 4-pt/1.9 cm near-end arm outscored the true 33 deg
    // corner by 5 deg of PCA jitter; genuine corners at 4-5 cm from an end
    // keep their legacy break -- pr/90 sec 8.6).  teb_second_max (cm): the
    // entry gate tolerates extra >stub prongs when exactly one segment
    // exceeds this cap (172832/61681: a second 11-13 cm prong made n_long=2
    // and the strict gate declined; measured NEGATIVE on its own motivating
    // events, pr/90 sec 8.5 -- stays OFF).  null = keys suppressed =>
    // byte-identical pre-fix config.  Escapes: SBND_TEB_TURN_MIN_ARM_FRAC /
    // SBND_TEB_SECOND_MAX runner envs (or -A).
    //
    // teb_turn_min_arm_frac SBND PRODUCTION ON 2026-08-17 (owner request,
    // pr/90 sec 8.8-8.9 gates: knobs-off 1067/1067 byte-identical vs the
    // harv3 production arms on the shipping binary (work-*-pr90off2);
    // knob-on live A/B 1066/1067 identical, the single mover IS the target
    // event 320865, toward, 37.29 -> 1.22 cm from the hand-scan click, zero
    // ADVERSE (pr90_movers.py vs TAGS_HARV3); nueCC48 48/48 + NCpi0 19/19
    // untouched.  0.4 (not the sec-6 illustrative 0.7) per the sec-8.3
    // census: excludes every starved arm (<= 6.1 cm) with >2x margin while
    // keeping the genuine 18.4 cm-arm break of evt 172942.
    teb_turn_min_arm_frac = 0.4,
    teb_second_max = null,
    // doc pr/90 round 4 (sec 9.5 D1/D3/D4) -- three knobs for the round-3
    // residual classes, C++ defaults false/0 = legacy.
    // teb_chain_topology: when n_long > 1, admit iff the cluster's segment
    // graph is a simple path ("still a line, no 3-track vertex") and the
    // candidate is the unique longest segment; chain-admitted candidates go
    // to route R3 only.  teb_r3_turn (deg) / teb_r3_hot (x mip median):
    // R3 breaks at the largest 10 cm-baseline local turn that carries a
    // vertex-activity spot within +-2 cm, refined to the activity maximum
    // (172832: t10 plateau 19-23.5 deg + 2.50x MIP at the click vs t35
    // 18.3 < 25; 61681: t10 54 deg + 3.1x MIP).  teb_bragg_veto_turn (deg):
    // veto an accepted R2 break below this turn when its short-arm end is
    // not Bragg-consistent (peak >= 2x MIP AND hot extent <= peak cm/MIP;
    // sec 9.4b owner calibration: kills 26.5-27.4 deg vs keeps >= 32.5).
    // false/null = keys suppressed => byte-identical pre-fix config.
    // Escapes: SBND_TEB_CHAIN_TOPOLOGY / SBND_TEB_R3_TURN / SBND_TEB_R3_HOT
    // / SBND_TEB_BRAGG_VETO_TURN runner envs (or -A).
    //
    // teb_chain_topology / teb_r3_* STAY OFF: the D1+D3 live A/B was net
    // NEGATIVE (pr/90 sec 10.6: 19 ADVERSE vs 6 toward on harv3 labels,
    // two cosmict flips) -- no local scalar separates a true interior
    // junction from an energetic delta ray.
    //
    // teb_bragg_veto_turn SBND PRODUCTION ON 2026-08-17 (owner request,
    // pr/90 sec 10.8-10.9 gates: knobs-off 1000/1000 byte-identical vs
    // production work-mcp1k-pr90on2, nueCC48/NCpi0 48/48+19/19; D4-only
    // A/B movers EXACTLY the five vetoed near-end breaks; labels: 291064
    // 159.36 -> 0.00 toward, 64503 sanctioned sec 9.0; 349461/278420
    // mid-track breaks byte-identical under the 15 cm near-end scope).
    teb_chain_topology = false,
    teb_r3_turn = null,
    teb_r3_hot = null,
    teb_bragg_veto_turn = 30.0,
    kink_walk_dqdx_stop = true,
    kink_break_protect = true,
    // doc pr/49 (18255-57441) -- cross-cluster projection-ghost deweighting
    // in the trajectory fit's 2D charge association (V-plane ghost: an
    // unrelated cluster's real charge aliases with the fitted track in one
    // view only and detours the fit; a live candidate cell outside the
    // fitted cluster's OWN blob coverage that sits INSIDE an OUT-OF-SCOPE
    // cluster's -- round 3: one with no segment in the current fit -- keeps
    // its measurement at x0.1 weight; cells covered by nobody keep full
    // weight, so dead-channel single-view charge stays usable by design).
    // **SBND PRODUCTION DEFAULT ON since the owner flip 2026-08-08** (doc
    // pr/49 round 3: off-gates off48d/off50d 48/48+50/50, census 42/48+28/50
    // all sentinel-gated, nusel-events 0/98; known open item: 172230-class
    // near-vertex robustness, doc pr/49).  C++ default -1 = off; value =
    // wire/slice tolerance in cells (0 = strict, the validated operating
    // point).  Companion numerics (ghost_dis 0 = scope-only, weight 0.1)
    // ride the C++ defaults.  Pass -A fit_blob_coverage=-1 (or the
    // SBND_FIT_BLOB_COVERAGE=-1 runner env) to restore the legacy path for
    // an A/B.
    fit_blob_coverage = 0,
    // doc pr/50 -- suspend the deweighting above during find_proto_vertex
    // (partition stage) so the recursive break partition forms on legacy
    // fits; all later fitting stages (main-vertex determination,
    // improve_vertex, final trajectory + dQ/dx) keep the ghost protection.
    // C++ default false = pr/49 behavior; false omits the key =>
    // byte-identical.  Validation: -A fit_blob_coverage_defer=true (or the
    // SBND_FIT_BLOB_COVERAGE_DEFER=true runner env).
    fit_blob_coverage_defer = false,
    // doc pr/50 -- main-vertex kink-consistency snap (172230-class
    // near-vertex robustness).  SBND PRODUCTION ON (owner flip after the
    // pr/50 Bee hand-scan: 172230 vertex recovered to 2.6 mm, 4/98
    // sample-wide firings all scanned; the defer knob stays OFF).  Numerics
    // (vks_* null) ride the C++ defaults.  Pass -A vertex_kink_snap=false
    // (or the SBND_VERTEX_KINK_SNAP=false runner env) to restore the
    // legacy path for an A/B.
    vertex_kink_snap = true,
    vks_radius = null,
    vks_min_dis = null,
    vks_angle = null,
    vks_margin = null,
    vks_collinear = null,
    vks_skirt = null,
    vks_baseline = null,
    vks_min_arm = null,
    vks_fit_miss = null,
    vks_hot_ratio = null,
    // doc pr/85 -- carry the old vertex's arms through the snap residual
    // when the snap arc is below this (cm): the residual edge old-vertex ->
    // corner IS the pr/85 "interposed stub" (mode 1a-VIA), left behind by
    // the snap above.  C++ default 0 = off.  DEFAULT OFF pending the pr/85
    // implementation-round gates.  Validation: -A vks_carry_prong=<cm> (or
    // the SBND_VKS_CARRY_PRONG runner env).  null/0 omits the key =>
    // byte-identical.
    vks_carry_prong = null,
    // doc pr/104 -- main-vertex junction snap: re-point the main vertex to
    // a nearby multi-prong track junction when the selected vertex has no
    // track prong of its own, or the joint line fit of both vertices' prongs
    // lands on the junction (18255-405707/65289/66712/282072/345633).
    // C++ default false.  SBND PRODUCTION ON 2026-08-21 (doc pr/104 sec 6:
    // OFF gates PASS 30/38/96/2000, zero ADVERSE movers on mcp1k/NCpi0,
    // nueCC48 nue-selected 39->40, owner pre-authorized "same as before").
    // Numerics (vjs_* null) ride the C++ defaults (5 cm / 3 cm / 2 / 150 deg /
    // 0.5 cm / 1.0 cm / min_move 1.0 cm).  Rollback: -A vertex_junction_snap=false
    // (or the SBND_VERTEX_JUNCTION_SNAP=false runner env).
    vertex_junction_snap = true,
    vjs_radius = null,
    vjs_min_arm = null,
    vjs_min_prongs = null,
    vjs_collinear = null,
    vjs_fit_margin = null,
    vjs_fit_rms = null,
    // doc pr/104: also arbitrate a main vertex that vertex_kink_snap itself
    // created (405707/65289/345633 are kink-snap products); C++ default false.
    // SBND PRODUCTION ON 2026-08-21 (doc pr/104 sec 6, validated together with
    // vertex_junction_snap; 5 of the 13 final fires are kink-snap products).
    vjs_override_kink_snap = true,
    vjs_min_move = null,          // C++ default 1.0 cm
    // esva_ignore_empty_2d (sbnd_xin/docs/73 sec 12, round 3; C++ default
    // false = legacy): eliminate_short_vertex_activities case 5 treats the
    // empty-2D-index sentinel (-1: the pre-existing segment has no points in
    // the query point's APA -- possible only on a cathode-crossing cluster)
    // as "no information" instead of "covered within 0.45 cm".  Legacy
    // vacuously deletes the cross-cathode junction segment the rescue
    // creates (SBND data evt 78242: 132.5 cm fitted muon segment removed,
    // 71 cm track_fit hole, far half then absorbed into an EM shower via
    // absorb_unreachable_main).  Single-APA clusters are unreachable by this
    // knob.  SBND PRODUCTION ON since 2026-08-17 (owner flip on the
    // docs/73 sec 12 validation).  Escape: SBND_ESVA_IGNORE_EMPTY_2D=0
    // omits the key (pre-round-3 baseline).
    esva_ignore_empty_2d = true,
    // doc pr/51 -- main-vertex graph audit (near-vertex graph-shape repair:
    // duplicate-corridor merge / charge-less-bridge removal / micro-stub
    // absorb + re-seat / one refit; 131357 / 268067 / 360535 / 142421 /
    // 285567).  SBND PRODUCTION ON 2026-08-15 (doc pr/85 sec 10 gates:
    // knob-off 1022/1022 byte-identical; 57-label census stubs 95->42,
    // interposed prongs 26->17, zero new; every >1 cm mover adjudicated
    // against the owner's hand labels lands toward-or-on the click; nueCC48
    // recovers 2 nue events, loses none; runtime/RSS neutral).  Owner
    // pre-authorized the flip conditional on those gates.  Escape for A/B:
    // -A main_vertex_graph_audit=false (or SBND_MAIN_VERTEX_GRAPH_AUDIT).
    // Operating point (doc pr/85 sec 10.4-10.6):
    //   mvga_reseat_angle = 0  -- re-seat DISABLED: both label-set re-seat
    //     firings moved the vertex OFF the owner's click (280017 1.64 cm,
    //     314838 0.60 cm);
    //   mvga_stub = 2.5        -- fitted by the sec 10.5 sweep;
    //   mvga_dup_frac = 0.8, mvga_stub_pts = 3 -- close the marginal-overlap
    //     absorb class (nfit=4 overlap=0.75), the whole of the adverse-mover
    //     population in sec 10.6.
    main_vertex_graph_audit = true,
    mvga_radius = null,
    mvga_dup_tol = null,
    mvga_dup_frac = 0.8,
    mvga_dup_angle = null,
    mvga_bridge_mip = null,
    mvga_reconnect = null,
    mvga_stub = 2.5,
    mvga_stub_pts = 3,
    mvga_reseat_angle = 0,
    // doc pr/51 round 3 -- op3 satellite-anchor radius (cm): reaches
    // terminal micro-stubs at satellite vertices near, not just at, the
    // main vertex (142421's 7082/7023, 285567's residual).  SBND PRODUCTION
    // ON 2026-08-15 with the audit flip above (doc pr/85 sec 10).
    mvga_satellite = 3.0,
    // doc pr/85 -- op3 interposed-stub absorb at the main-vertex anchor:
    // opens op3's terminal-only degree(far)==1 line to INTERPOSED stubs
    // (far vertex carries the real prongs, degree >= 2 -- the whole of the
    // pr/85 mode-1a-VIA class, 21/462 events), splicing the far prongs
    // through the stub onto the main vertex.  Inert unless
    // main_vertex_graph_audit.  SBND PRODUCTION ON 2026-08-15 (doc pr/85
    // sec 10 gates; 15 firings across the three samples, every adjudicated
    // one <= 0.46 cm vertex motion).  Escape: -A mvga_interposed=false (or
    // SBND_MVGA_INTERPOSED); angle rides the C++ default 150 deg.
    mvga_interposed = true,
    // doc pr/86 -- the straight track that misses the vertex (the owner's
    // named priority): an interposed segment of 2.5-10 cm carries two
    // thirds of the Class-B defect, above the mvga_stub=2.5 splice ceiling.
    // SBND PRODUCTION ON 2026-08-15 (doc pr/86 implementation round gates:
    // knob-off 2134/2134 archives byte-identical vs the pr/85 flip arms;
    // full-sample Class-B cases 90->48 and orphans 118->82 with ZERO events
    // worse; every mover <= 0.28 cm, zero adverse vs the hand labels; the
    // four pr/85 sec 10.6 adverse-mover events at exactly 0.00; nueCC48 +4
    // nue recoveries (163543, 268784, 400474, 423981) vs 1 named loss
    // (122660 -- a Class-B repair at a vertex 4.95 cm off the owner's
    // click, i.e. a pre-existing mis-pick the splice entrenched; owner
    // hand-check); runtime/RSS neutral).  Owner pre-authorized the flip
    // conditional on validation.  Escapes: SBND_MVGA_INTERPOSED_LEN,
    // SBND_MVGA_INTERPOSED_ANGLE, SBND_MVGA_INTERPOSED_DEG1,
    // SBND_MVGA_SAT_DUP_FRAC (or the -A TLAs).
    // P1 -- interposed-splice candidate ceiling, cm.  Widens ONLY the
    // splice; the terminal absorb keeps mvga_stub=2.5.  10 = the measured
    // edge of the interposed-segment length distribution (sec 10.3).
    mvga_interposed_len = 10,
    // P2 -- far-end collinearity gate, deg (C++ default 150).  The Class-A
    // pile-up sits at 120-150 deg and both live declines were 139.8/130.1;
    // 130 fitted by the corner sweep (l10a130: orphans 102->68 vs 92 at
    // 150, zero adverse movers).
    mvga_interposed_angle = 130,
    // P4 -- satellite-anchor op3 overlap threshold.  All four pr/85
    // sec 10.6 adverse absorbs were main-anchor d=0.00; 0.7 restores the
    // pre-pr/85 threshold at satellites only (evt30504's wanted absorb).
    mvga_sat_dup_frac = 0.7,
    // P1b -- interposed splice at degree-1 main anchors (26 of the 86
    // Class-B cases sit there; fixes 349945 and 168432).  The terminal
    // absorb still requires degree >= 2.
    mvga_interposed_deg1 = true,
    // doc pr/86 sec 15 (round 2) -- the round-1 splice repaired the graph
    // but not the picture: carry_prong_execute concatenates wcpt chains and
    // the op4 refit seeds from wcpts, so the rendered trajectory never
    // moved (owner: 349945's direct track still missing, 38856's turn
    // still there).  Round 2 makes the repair visible: R1 re-derives the
    // spliced near-anchor stretch straight (es2 recipe: 0.6 cm samples,
    // is_good_point charge veto, Steiner snap), R2 late-collapses degree-2
    // zigzag junctions near the main vertex, with op3<->op3.5 interleave +
    // created-stub charge-gated splice (the 349945 chain: collapse, then
    // splice the op1-created elbow, then straighten).  SBND PRODUCTION ON
    // 2026-08-16 (sec 15 gates: knobs-off 2134/2134 byte-identical on the
    // shipping binary; kink>=60 near-vertex band 60->43 full-sample
    // [pre-round-1 was 39]; zigzag ratio>=1.5 12->6; 349945 approach 3
    // hops/ratio 2.42 -> 1 hop/1.16; ZERO adverse movers 3 samples; zero
    // nue-score sign changes; pr/85 four adverse events untouched;
    // runtime/RSS neutral).  Escapes: SBND_MVGA_SPLICE_STRAIGHTEN,
    // SBND_MVGA_APPROACH_COLLAPSE, SBND_MVGA_STRAIGHTEN_RADIUS (or -A).
    // R1 -- straighten reach (cm) past the dissolved junction.
    mvga_splice_straighten = 5,
    // R2 -- op3.5 junction-collapse radius (cm) around the main vertex.
    mvga_approach_collapse = 15,
    // R1/R2 charge-veto radius (cm; C++ default 0 = prototype 0.2, which
    // vetoes every long-stub straighten -- the charge ridge deviates
    // ~1.6 cm from the straight chord; 1.0 is grid-validated, uniquely
    // required by 349945 + 122660, zero adverse anywhere).
    mvga_straighten_radius = 1.0,
    // doc pr/83 r3 -- the duplicate-corridor round, SBND PRODUCTION ON
    // (owner flip 2026-08-17): op1 unscoped (-1), op1 threshold 0.7 for
    // >=10cm pairs (C++ length gate; 390842 guard), post-op3 dup pass
    // (class A), abandoned-cluster dup audit (Mechanism C + losing-candidate
    // orphans, 350935/359980).  Gates: knob-off 1022/1022 byte-identical;
    // census 17->0 (511-evt) and 14->0 (mcp2k), zero new findings.
    // mvga_carry_max stays null (not needed; class A cleared without it).
    // Escapes: SBND_MVGA_OP1_RADIUS, SBND_MVGA_OP1_DUP_FRAC,
    // SBND_MVGA_OP1_POST, SBND_MVGA_CARRY_MAX, SBND_SWAP_ORPHAN_DUP_AUDIT
    // (or -A).  null/false omit the keys => byte-identical legacy path.
    mvga_op1_radius = -1,
    mvga_op1_dup_frac = 0.7,
    mvga_op1_post = true,
    mvga_carry_max = null,
    swap_orphan_dup_audit = true,
    // doc pr/83 r4 -- projective duplicate collapse at the main vertex,
    // SBND PRODUCTION ON (owner flip 2026-08-18, 4-event Bee scan
    // approved): a 1-track-1-shower stem split into two 3D tracks that
    // overlap in >= 2 of 3 wire views; the charge-starved member reads
    // stem dQ/dx ratio 0.08-0.28 final, 0.33-0.48 at mvga time -- hence
    // ratio 0.55, margin over the measured 0.47/0.48 and far from
    // MIP-parity two-prongs (geometry gates alone: zero false pairs in
    // 559 events).  Gates: knob-off 1024/1024 byte-identical; projective
    // census 4->0 (511) + 2->0 (mcp2k), zero new; r3 census stays 0;
    // pr/86 census identical; movers = exactly the census events.
    // Escapes: SBND_MVGA_PROJ_DUP_FRAC, SBND_MVGA_PROJ_DQDX_RATIO
    // (or -A).  null omits the keys => byte-identical legacy.
    mvga_proj_dup_frac = 0.7,
    mvga_proj_dqdx_ratio = 0.55,
    // doc pr/83 r4b (284206): op1-proj's own angle ceiling, SBND
    // PRODUCTION ON (owner request 2026-08-18): the residual stem pair
    // reads 22 deg -- just over op1's shared 20 -- and its ratio 0.52
    // already passes 0.55, so the angle was the only blocker.  Small-set
    // check: 284206 stem 3 -> 1 track (duplicate 161 MeV electron gone);
    // all 7 other knob-affected events byte-identical at 25.
    // Escape: SBND_MVGA_PROJ_ANGLE (or -A).  null omits => byte-identical.
    mvga_proj_angle = 25,
    // doc pr/99 round 2 -- op3.5 approach-collapse guards + op1-post charge
    // second-opinion.  C++ defaults 0/false = legacy.  null/false => key
    // omitted => byte-identical pre-fix config.
    // SBND PRODUCTION ON 2026-08-20 (doc pr/99 round 2; owner pre-authorized
    // flip on validation PASS).  ac_veto_radius stays OFF: 0.2 cm measured
    // ADVERSE (kills the 349945 design case -- re-confirms pr/86 Stage A's
    // deliberate 1.0 cm relax).  Ghost thresholds ride the C++ defaults
    // (overlap 0.7 / dqdx 0.25 / min_len 10 cm).
    mvga_ac_veto_radius = null,
    mvga_ac_chord_max = 30,
    mvga_ac_no_cascade = true,
    // doc pr/103 (SBND 18255-405707) -- mvga op0 pass-through split: a prong of
    // a junction J near the main vertex whose fitted interior passes through the
    // main vertex is split there (the J-M piece becomes the connecting stub,
    // so op1 no longer deletes the connector as a "duplicate" and op3 can
    // re-anchor J's other prong onto the vertex).  C++ default 0 = off;
    // passthru_tol C++ default 1.0 cm.  null => key omitted => byte-identical.
    // SBND PRODUCTION ON 2026-08-21 (doc pr/103 sec 7; owner pre-authorized the
    // flip on validation PASS: OFF gate 2000/2000+96/96+38/38+30/30, zero ADVERSE
    // movers, no nueCC48 nue change).  passthru_tol rides the C++ default 1.0 cm.
    mvga_passthru = 4,
    mvga_passthru_tol = null,
    // doc pr/103 -- op3 interposed splice: when the far-angle gate declines at
    // the main anchor (the 3-track vertex split in two by a 0.6-4 cm stub:
    // 65289/345633/400856/287517), try the per-prong charge-verified straighten
    // splice instead.  C++ default false; false => key omitted => byte-identical.
    mvga_interposed_fallback = true,   // SBND PRODUCTION ON 2026-08-21 (doc pr/103 sec 7)
    // doc pr/103 sec 6: measured-far_angle floor for the fallback (deg).  C++
    // default 0 = only "measured"; validated production value 45 (hairpins /
    // back-folds below it only shorten tracks).  null => key omitted.
    mvga_interposed_fallback_min_angle = 45,   // SBND PRODUCTION ON 2026-08-21 (doc pr/103 sec 6 ledger)
    mvga_dup_starved_asym = 0.55,
    mvga_dup_starved_mip = 0.8,
    mvga_dup_starved_span = 0.5,
    // doc pr/51 (18255-506746) -- DL rerank cross-cluster swap guard: with
    // the guard on, an accepted DL vertex can never swap the main cluster
    // (506746: one confident uBooNE-net voxel, s_dl = +576, moved the
    // vertex 28 cm onto a non-flash-matched cluster).  DEFAULT OFF pending
    // owner review.  Validation: -A dl_vtx_swap_guard=true (or the
    // SBND_DL_VTX_SWAP_GUARD runner env).  false omits the key =>
    // byte-identical.
    dl_vtx_swap_guard = false,
    // doc pr/106 sec 10 -- exclusion-free charge cloud for the DL vertex net:
    // with fit_exclusion on, the cells equidistant from two prongs are dropped
    // from both and the net's vertex voxel is starved (nueCC DL-alone 34 -> 42
    // /47 with exclusion off).  This refits each cluster once without
    // exclusion, feeds that cloud to the net, and restores the exclusion fit.
    // C++ default false => key omitted => byte-identical.  Validation:
    // --tla-code dl_vtx_cloud_no_exclusion=true (or SBND_DL_VTX_CLOUD_NO_EXCLUSION).
    dl_vtx_cloud_no_exclusion = false,
    // doc pr/112 sec 11 -- the DUAL CHAIN: a second, exclusion-free PR pass
    // suggests the neutrino vertex and the production (exclusion-ON) pass
    // decides.  dl_vtx_dual_chain C++ default false => key omitted =>
    // byte-identical (also the retrain-era off switch).  dual_chain_mode
    // "snap" (OFF final vertex snapped to the nearest production candidate,
    // accepted iff d <= dual_chain_transfer_max cm) | "voxels" (OFF top-K
    // replaces the ON inference in the rerank) | "union" (pooled top-K, plus
    // dual_chain_vtx_weight * max(0, 1-d/D) proximity term).
    // dual_chain_transfer=false is the PROBE (pass runs, agreement flag
    // recorded, nothing moves); NOT transfer_max=0, which transfers.
    // Runner env: SBND_DL_VTX_DUAL_CHAIN / SBND_DUAL_CHAIN_{MODE,TRANSFER,
    // TRANSFER_MAX,ALLOW_CLUSTER_SWAP,VTX_WEIGHT}.
    dl_vtx_dual_chain = false,
    dual_chain_mode = null,
    dual_chain_transfer = false,
    dual_chain_transfer_max = null,
    dual_chain_allow_cluster_swap = null,
    dual_chain_vtx_weight = null,
    // doc pr/107 -- dQ/dx fit keeps every trajectory point (prototype
    // parity).  do_multi_tracking's toolkit-only third form_map_graph pass
    // (before dQ_dx_multi_fit) re-applies the zero-quantity drop to the final
    // trajectory; with fit_exclusion on that deletes the junction-adjacent
    // points whose cells update_association stripped (442 points over 47
    // nueCC48 events vs 86 with exclusion off), which starves the DL vertex
    // net's input cloud.  The prototype fits dQ/dx on every trajectory point.
    // C++ default false => key omitted => byte-identical.  Validation:
    // --tla-code dqdx_fit_keep_all_points=true (or SBND_DQDX_FIT_KEEP_ALL_POINTS).
    dqdx_fit_keep_all_points = false,
    // doc pr/89 Arm C (C2) -- rule-1 outgoing-prong topology term in the DL
    // rerank composite: s_topo = w * (frac - center), vote-gated.  C++
    // defaults 0/0 = term never computed; null omits the keys =>
    // byte-identical.  The offline C1 replay selected weight 3.0, center 0
    // (+12/924, pr/89 sec 11.5).  Validation: --tla-code dl_vtx_topo_weight=3.0
    // (or the SBND_DL_VTX_TOPO_WEIGHT runner env).  DEFAULT OFF pending the
    // live A/B.
    dl_vtx_topo_weight = null,
    dl_vtx_topo_center = null,
    // doc pr/51 round 3 -- apply the traditional main-vertex path's cluster
    // swap decision instead of silently discarding it (a latent bug: the
    // decision fires today but never reaches the caller).  DEFAULT OFF
    // pending owner review.  Validation: -A main_vertex_swap_apply=true (or
    // the SBND_MAIN_VERTEX_SWAP_APPLY runner env).  false omits the key =>
    // byte-identical.
    main_vertex_swap_apply = false,
    // doc pr/51 round 4 -- diagnostic-only rough-path probe for the
    // near-vertex short-cut investigation (path-COST, not graph-shape).
    // No graph/fit/segment content is ever changed; every line is TRACE.
    // Validation: -A rough_path_probe=true (or the SBND_ROUGH_PATH_PROBE
    // runner env).  false omits the key => byte-identical.
    rough_path_probe = false,
    // doc pr/51 round 5 -- steiner gap penalty, the H1 short-cut fix the
    // round-4 probe validated: when scale > 0, do_rough_path routes on a
    // per-cluster "steiner_graph_gap" flavor whose edges are re-weighted
    // w' = w * (1 + scale * unsupported_fraction), so a gap-spanning chord
    // no longer beats following the charge around a corner.  C++ defaults:
    // scale 0 (OFF), dead_alpha 0.25, min_edge 0.5 cm, sample_step 0.3 cm,
    // point_radius 0.2 cm.  SBND PRODUCTION ON at scale 2.0 (2026-08-12,
    // owner pre-authorized on clean validation): off-gates 0/117
    // byte-identical, on-census 101/117 display-level movers with nusel
    // 0/117, wall/RSS deltas negligible; targets 131357/234638/268067/
    // 506746 fixed (work-pr51r5-* arms, doc pr/51 round 5).  Legacy escape
    // for A/B: -A steiner_gap_penalty=0 (or SBND_STEINER_GAP_PENALTY=0)
    // restores the pre-flip production bare, byte-exact (C++ gate is
    // scale <= 0; flip/escape gates 0/117 both directions).
    steiner_gap_penalty = 2.0,
    sgp_dead_alpha = null,
    sgp_min_edge = null,
    sgp_sample_step = null,
    sgp_point_radius = null,
    // doc pr/51 round 6 -- weak-charge deficit term on the same gap flavor
    // (residuals 18259-131357 3-track V, 18255-506746 branch turn: chords
    // that are image-supported but charge-poor, invisible to the
    // unsupported-fraction penalty at any scale).  C++ defaults:
    // weak_scale 0 (OFF), weak_qref 2000 (calc_charge_wcp charge units).
    // SBND PRODUCTION ON at (5.0, 6000) (2026-08-12): off-gates 0/117
    // byte-identical; nusel/nusel-table 0/117 at both grid scales; both
    // round-6 targets fixed (131357 clean 2-track V, MV 0.94 cm from the
    // true corner; 506746 long track reaches the vertex).  DISCLOSED COST
    // (doc pr/51 round 6 census): nu-vertex moves >10 cm on 25/115 and
    // kine_reco_Enu by >100 MeV on ~44/115 -- same footprint class as the
    // accepted round-5 flip measured identically (17 and 30), majority the
    // same bistable events; every >10cm mover is in the pr51r6 Bee sets.
    // Escape for A/B: SBND_SGP_WEAK_SCALE=0 (or -A sgp_weak_scale=0)
    // restores round-5 production byte-exact; note steiner_gap_penalty=0
    // kills BOTH terms (pre-round-5 legacy).
    sgp_weak_scale = 5.0,
    sgp_weak_qref = 6000.0,

    // doc pr/73 -- per-edge DEBUG sentinel for the steiner_graph_gap scan
    // (endpoints, midpoint, w, bad, both recovered vertex charges,
    // deficit).  Log-only diagnostic; answers "where are the penalized
    // edges" rather than "how many".  C++ default false = legacy; false
    // omits the key => byte-identical compiled config.
    // Validation: -A sgp_edge_probe=true (or SBND_SGP_EDGE_PROBE=true).
    sgp_edge_probe = false,
    // doc sbnd_xin/docs/pr/73 round 2, fix F3a -- bound what the round-6
    // penalty may do to the ROUTE.  When >= 0 and the gap flavor is in use,
    // do_rough_path also routes on the untouched base flavor and keeps the
    // BASE route whenever the penalized one strays further than this many cm
    // from it.  C++ default -1 = off (unbounded, today's behaviour); null
    // omits the key => byte-identical compiled config.
    // NB 0 is a meaningful cap, so the off-test is `< 0`, not `<= 0`.
    // SBND PRODUCTION ON at 3 cm -- owner flip 2026-08-13, doc pr/73 sec 9.
    // Fixes the ISO case 18255-57903: the corridor returns from two segments
    // (path/chord 1.326 & 1.252, bow 2.52 & 2.94 cm, a 69 deg hairpin) to one
    // 47.4 cm segment at 1.026 / 1.86 cm / jitter 0.081 -- the pre-round-5
    // answer to three decimals, charge coverage included.
    // The owner's stated criterion for this round was "keep the round 5/6
    // improvements, fix the ISO case", and both halves are PROVEN: 18259-131357
    // and 18255-506746, the two events round 6 was built to fix, are
    // byte-identical with the guard on.
    // Footprint 46/117 events, nusel 0/117 flips; the nu-vertex movers are
    // under hand scan (docs/pr/pr73f3a-movers.index.txt).  Open item: +2
    // dangling PF roots on 18255-285567.
    // Escape for A/B: -A sgp_max_sep=-1 (or SBND_SGP_MAX_SEP=-1).
    sgp_max_sep = 3.0,
    // doc sbnd_xin/docs/pr/83 -- orient break_segment() splits to the wcpt
    // path (find_vertices) instead of boost source/target.  On a reversed
    // graph edge (routinely produced by the examine_vertices re-routes) the
    // legacy slice hands each child the WRONG half of the parent's
    // wcpts/fits: vertex fit points land 67-118 cm from their wcpts, both
    // children's fitted trajectories stack onto one arm, and the fit-vs-wcpt
    // divergence seeds runaway vertex-activity bridges (mcp1k
    // 283040/59899/72586 -- "many overlapping tracks on one long track").
    // Prototype parity: WCP resolves start/end by wcpt-index equality tested
    // in both orientations (NeutrinoID_proto_vertex.h:595-601) and cannot
    // cross.  C++ default false.
    // Escape for A/B: -A break_seg_orient=false (or SBND_BREAK_SEG_ORIENT=false).
    // SBND PRODUCTION ON 2026-08-15 (doc pr/83 gates: knob-off 1022/1022
    // byte-identical; on-arm footprint 25->12 findings, 0 new; both >1cm
    // vertex movers land on the owner's hand labels).  Owner pre-authorized
    // the flip conditional on those gates.
    break_seg_orient = true,
    // doc sbnd_xin/docs/pr/75 -- record, per event, HOW the neutrino vertex
    // was chosen: compare_main_vertices' additive score per candidate, the DL
    // top-K voxels, the seven rerank composite terms, and the accept route
    // (doc pr/52 sec 1).  PrDisplayDump emits it as "vertex_scoreboard" so a
    // hand scan can rank candidates and become acceptance-tuning input.
    // Pure recording -- no decision reads it.  C++ default false = legacy;
    // false omits the key => byte-identical compiled config.
    // Validation: -A vertex_scoreboard=true (or SBND_VERTEX_SCOREBOARD=true;
    // the driver turns it on automatically for PR_EXTRA_STAGES=pr_display).
    vertex_scoreboard = false,
    // doc sbnd_xin/docs/pr/79 sec 10 -- live-feature harvest for DL-vertex
    // training: the exact live SCN input cloud (pre-voxelization) + the
    // traditional-path per-candidate features (proton topology, z prior,
    // degree, FV, conflicts, fit chi2, global-scorer rows).  Recording only;
    // REQUIRES vertex_scoreboard (the runner auto-enables it).  C++ default
    // false = legacy; false omits the key => byte-identical compiled config.
    // Validation: -A dl_vtx_harvest=true (or SBND_DL_VTX_HARVEST=true).
    dl_vtx_harvest = false,
    // doc pr/51 round 7 -- robust vertex fit (dynamic per-leg re-seat-free
    // direction windows for MyFCN, disagreement-gated, relaxed prior on
    // substituted 2-leg vertices).  C++ defaults: robust false, main_only
    // true, min_len 10, rin_margin 2, rout_frac 0.5, rout_min 9,
    // rout_max 18, angle 20, min_pts 5, min_aniso 3, prior_range 1
    // (lengths cm, angle deg); satellites below ride the C++ defaults.
    // SBND PRODUCTION ON (2026-08-12), owner instruction "if the validation
    // pass, turn on the knob for SBND" after the 117-evt validation: off-gate
    // 0/117 byte-identical, nusel 0/117 at ON, TOTAL knob footprint 3/117
    // events -- 18255-57903 (self-confirmed hairpin vertex, census legs
    // 74/44 deg, moves 4.62 cm onto the muon line), 18255-56982 (35 deg leg,
    // 2.23 cm), 423981 (0.08 cm micro-adjust, >=3-leg tight prior) -- zero
    // nu-vtx movers >10 cm, zero |dEnu| > 100 MeV, wall/RSS at noise.
    // Escape for A/B: SBND_MVFIT_ROBUST=false (or -A mvfit_robust=false)
    // restores round-6 production byte-exact (escape gate 0/117).
    mvfit_robust = true,
    mvfit_main_only = null,
    mvfit_min_len = null,
    mvfit_rin_margin = null,
    mvfit_rout_frac = null,
    mvfit_rout_min = null,
    mvfit_rout_max = null,
    mvfit_angle = null,
    mvfit_min_pts = null,
    mvfit_min_aniso = null,
    mvfit_prior_range = null,
    // doc pr/54 -- keep well-supported isolated residual segments in
    // find_other_segments (18255-142421 "missing gammas": a separated EM
    // shower of the main cluster is fit and then silently discarded because
    // neither endpoint touches the existing graph).  SBND PRODUCTION ON --
    // owner flip 2026-08-09 after the Bee before/after hand-scan of 142421
    // ("this was a bug"): bare production now IS the validated
    // work-pr54-on142421 arm config.  Legacy escape for A/B: -A
    // other_seg_keep_isolated=false (or SBND_OTHER_SEG_KEEP_ISOLATED=false)
    // restores the pre-flip production bare, byte-exact (doc pr/54 SS12
    // flip proofs).  C++ floors 25 points / 3 cm ride as defaults (null
    // omits their keys); C++ knob default itself stays false.
    other_seg_keep_isolated = true,
    other_seg_keep_isolated_min_points = null,
    other_seg_keep_isolated_min_length = null,
    // doc pr/102 P1 -- OR-disjuncts on the pr/54 keep above: min_nnf admits
    // a candidate below the 25-terminal floor when its number_not_faked is
    // at least this (pr/67 sec 10.3 S1; 18255-70084 read nnf 10/12 vs pr/54
    // noise at 0); len_admit (cm) admits any candidate at least this long
    // (pr/102 sec 4 B1: 67-145 cm candidates dropped at 16-23 terminals).
    // C++ defaults 0 = off; null omits the keys => byte-identical.
    //
    // len_admit: SBND PRODUCTION ON at 30 -- owner authorization 2026-08-20
    // ("implement P1 ... if validation pass, turn them on by default"),
    // validation doc pr/102 sec 8: nuecc48 sweep 2/48 changed ZERO nue
    // flips, mcp1k 6/1000 changed zero ADVERSE movers, census fixes the
    // 145.5/67.1 cm dropped tracks (18255-292384/387850).  Named
    // hand-check item: 18255-284794 numu 2.12 -> -1.28 (its real 71.4 cm
    // second track is recovered; vertex unmoved).  Legacy escape:
    // -A other_seg_keep_isolated_len_admit=0 (or SBND_OSEG_LEN_ADMIT=0...
    // any empty env leaves this default) -- C++ knob default stays 0.
    //
    // min_nnf: STAYS OFF -- validation FAILED at 4 (nueCC48 nue ledger
    // -4/+1) and carries one named nue loss at 8 (389538 4.3 -> -15 vs
    // 30504 gain); doc pr/102 sec 8.3.  Owner hand-scan before any flip.
    other_seg_keep_isolated_min_nnf = null,
    other_seg_keep_isolated_len_admit = 30.0,
    // doc pr/102 P2 -- 3-D uncovered-charge radius (cm): imaged charge
    // farther than this from EVERY existing fitted trajectory cannot be
    // 2-D-tagged in find_other_segments step 1 and counts toward
    // number_not_faked at step 8 / re-eval (the B2 nnf=0 fragmentation
    // family).  C++ default 0 = off; null omits the key => byte-identical.
    other_seg_uncover_3d = null,
    // doc pr/67 round 3 (S2) -- isochronous-snap size gate in cm, the first
    // clause of the guard at NeutrinoOtherSegments.cxx:721.  The machinery
    // behind it (modify_vertex/segment_isochronous) is the only thing that
    // ATTACHES an isochronously displaced branch to its parent, and at the
    // legacy 10 cm it never saw the doc pr/67 branches (dir_mag 4.70 / 4.38 /
    // 4.34 cm on 18264-137238, 18259-42280, 18345-21073).
    //
    // SBND PRODUCTION ON at 4.0 -- owner flip 2026-08-12 after the Bee
    // before/after hand-scan of the three targets plus all eight collateral
    // neutrino-vertex movers ("based on my scan, the overall is positive").
    // Bare production now IS the validated work-pr67f-on{48,19,50} arm config.
    // Legacy escape for A/B: -A iso_snap_min_dir_mag=10.0 (or
    // SBND_ISO_SNAP_MIN_DIR_MAG=10.0) restores the pre-flip production bare,
    // byte-exact (doc pr/67 sec 11.12 flip proof).  C++ knob default itself
    // stays 10.0.
    //
    // Owner-accepted cost, doc pr/67 sec 11.7: this moves the reconstructed
    // neutrino vertex on 30 of 117 census events, 9 by more than 10 cm (max
    // 82.4 cm).  Do NOT read a vertex mover here as a regression without
    // re-reading that section first.
    iso_snap_min_dir_mag = 4.0,
    // doc pr/65 round 3 -- offer graph-unreachable main-cluster segments
    // (the disconnected components other_seg_keep_isolated above creates) to
    // the shower absorbers by relaxing the cluster()==main_cluster guards to
    // a main_vertex-reachability test.  SBND PRODUCTION ON -- owner flip
    // 2026-08-11 (doc pr/65 round 3): bare production IS the validated
    // work-pr65-on48 arm config.  Gate work-pr65-base48 vs work-pr65-off48
    // 0/48 archives + 0/48 pctree + nusel 0/48; on-census 21/48 movers, all
    // sentinel-attributed, 0 unclaimed left for the orphan net (evt 54095:
    // 12/12 fragments absorbed, shower 2334 -> 2633 MeV, 9 -> 3 root nodes).
    // Legacy escape for A/B: -A shower_absorb_unreachable_main=false (or
    // SBND_SHOWER_ABSORB_UNREACHABLE_MAIN=false).  C++ default stays false.
    shower_absorb_unreachable_main = true,
    // doc pr/59 round 2 -- per-cluster orphaned-associate_points rescue
    // (18255-142421 seg 20; 116944-71372 segs 19052/19053/136199).  SBND
    // PRODUCTION ON -- owner flip 2026-08-10 after the work-pr59r2-on19 gate
    // (12/12 known orphans rescued, zero manufactured-orphan regressions,
    // nusel-table.tsv byte-identical across the 19-event manifest) and the
    // Bee before/after hand-scan of 71372/142421: bare production now IS the
    // validated work-pr59r2-on19 arm config.  Legacy escape for A/B: -A
    // assoc_full_recluster=false (or SBND_ASSOC_FULL_RECLUSTER=false)
    // restores the pre-flip production bare, byte-exact (doc pr/59 sec 9
    // gate proofs).  C++ knob default itself stays false.
    assoc_full_recluster = true,
    // doc pr/64 round 7 -- reassign, instead of drop, an association point
    // that loses (or never enters) clustering_points_segments' Stage-C ghost
    // removal to a segment in the SAME cluster that actually wins the global
    // 2D projection contest.  SBND PRODUCTION ON 2026-08-11, per owner
    // pre-authorization ("if validation passed, turn this on") conditioned on
    // the 48-event nueCC off-gate/mover census (doc pr/64 round 7):
    // nusel-table.tsv byte-identical 48/48 (selection untouched); mabc-pr.zip
    // moves on 47/48 (the 48th has no PF reconstruction at all); every moved
    // event's kine_reco_Enu goes UP (mean +10.7 MeV, max +51.9 MeV on this
    // sample) -- monotonic and one-directional by construction (the rescue
    // only ever adds a previously-silently-dropped point, never removes one),
    // consistent with correcting a systematic under-count, not a random
    // perturbation.  NOTE: this does NOT resolve the owner's originally
    // reported 18259-18625 blob at (142.1,78.3,176.5) -- round 7 traced that
    // specific loss to a DIFFERENT, not-yet-fixed mechanism (a segment that
    // legitimately won those points is later removed from the graph with no
    // re-association pass); see doc pr/64 round 7 "still open" section.
    // Legacy escape for A/B: -A assoc_reassign_orphans=false (or
    // SBND_ASSOC_REASSIGN_ORPHANS=false) restores the pre-flip production
    // bare, byte-exact.  C++ knob default itself stays false.
    assoc_reassign_orphans = true,
    // doc pr/64 round 8 -- clear a merge survivor's associate_points when
    // examine_structure_final_1/_1p/_3 (called unconditionally inside
    // determine_main_vertex) deletes a segment that had non-empty
    // associate_points, so pr/59's reassociate_cluster_orphans any_orphan
    // trigger correctly re-fires instead of leaving the survivor's stale
    // cloud in place.  Root-caused the round 7 "still open" item: the
    // 18259-18625 blob at (142.1,78.3,176.5) is exactly this -- a segment
    // legitimately wins those points, is then absorbed into its main-vertex
    // neighbor by examine_structure_final_1p, and its points are discarded
    // with no re-derivation.  SBND PRODUCTION ON 2026-08-12, per owner
    // confirmation of the round 8 validation: nusel-table.tsv byte-identical
    // 48/48 on the nueCC48 sample (selection untouched); mabc-pr.zip moves
    // on only 3/48 (54095/174637/389538, a much narrower footprint than
    // round 7's Stage-C mechanism); every mover's kine_reco_Enu goes UP only
    // (+4.377/+2.246/+0.602 MeV) -- monotonic, matching the mechanism
    // (clearing only ever lets pr/59 add points back, never removes any);
    // the reported evt 18625 blob itself is now resolved (0 -> 8 bbox
    // points recovered).  wcdoctest-clus 176/176.  Legacy escape for A/B:
    // -A assoc_clear_on_merge=false (or SBND_ASSOC_CLEAR_ON_MERGE=false)
    // restores the pre-flip production bare, byte-exact.  C++ knob default
    // itself stays false.
    assoc_clear_on_merge = true,
    // doc sbnd_xin/docs/pr/72 round 2 -- guard examine_structure_3 against
    // merging a genuine near-vertex track stub into an unrelated
    // shower/track trunk (18255-196649: a real 6.28cm stub, terminal at the
    // true neutrino vertex, silently absorbed into a 33cm shower trunk).
    // SBND PRODUCTION ON 2026-08-12, owner flip: a 480-point grid scan over
    // the 48 nueCC + 19 NC pi0 + 50 PR data blast radius picked stub_max=7cm/
    // len_ratio=2.0 (ang3_min=15deg/ang_ratio=1.0/require_terminal=true
    // already correct) as the operating point that suppresses exactly ONE
    // junction in the whole 117-event population -- 196649's own target --
    // with an empty residual (misses) list; on-arm full-reconstruction
    // validation confirmed 1/117 archive-level movers, 0/117 selection
    // flips, and bundle-level TGM/STM/FC/lm byte-identical for 196649
    // itself (evidence the main vertex determination is unperturbed).
    // wcdoctest-clus 180/180.  Legacy escape for A/B: -A
    // es3_stub_guard=false (or SBND_ES3_STUB_GUARD=0) restores the pre-flip
    // production bare, byte-exact.  C++ knob default itself stays false;
    // key suppressed in the compiled config when off => byte-identical.
    // Numeric sub-parameters stay null (component's own C++ default, fitted
    // from the same 117-event census -- see doc pr/72 round 2).
    es3_stub_guard = true,
    es3sg_stub_max = null,
    es3sg_len_ratio = null,
    es3sg_ang3_min = null,
    es3sg_ang_ratio = null,
    es3sg_require_terminal = null,
)
    local base = import 'pgrapher/experiment/sbnd/simparams.jsonnet';
    local params = base {
        lar: super.lar {
            DL:          DL * wc.cm2 / wc.s,
            DT:          DT * wc.cm2 / wc.s,
            lifetime:    lifetime * wc.ms,
            drift_speed: driftSpeed * wc.mm / wc.us,
        },
    };

    local tools_all = tools_maker(params);
    local anodes  = [tools_all.anodes[i] for i in anode_indices];

    local source = g.pnode({
        type: 'TensorFileSource',
        name: 'pr_pctree',
        data: { inname: input, prefix: 'clustering_' },
    }, nin=0, nout=1);

    local clus_mod = import 'clus.jsonnet';
    local clus_maker = clus_mod(
        output_dir=output_dir,
        runNo=run,
        subRunNo=subrun,
        eventNo=event,
        reality=reality,
    );
    // dQ/dx (e/cm, SBND 0.5 kV/cm) + range (detector-agnostic) LinterpFunctions.
    // Resolved through WIRECELL_PATH (was a relative '../particle_dataset.jsonnet'
    // when this job lived under sbnd_xin, which only resolved from the real,
    // non-symlinked working dir).
    local pds = (import 'pgrapher/experiment/sbnd/particle_dataset.jsonnet')();
    local pr = clus_maker.pr(anodes, dump=true,
                             pipeline_names=pipeline_names,
                             tensor_outname=save_tensors,
                             trackfitting_config_file=trackfitting_config,
                             particle_dataset=pds.particle_dataset,
                             extra_uses=pds.all,
                             dl_weights=dl_weights,
                             dl_vtx_min_accept_score=dl_vtx_min_accept_score,
                             dl_vtx_top_k=dl_vtx_top_k,
                             dl_vtx_rerank=dl_vtx_rerank,
                             beam_window=[t * wc.us for t in beam_window_us],
                             beam_window_only=beam_window_only,
                             nu_skip_cosmic=nu_skip_cosmic,
                             nu_skip_cosmic_bundle=nu_skip_cosmic_bundle,
                             tgm_neutrino_candidate=tgm_neutrino_candidate,
                             tgm_chord_charge=tgm_chord_charge,
                             tgm_chord_mode=tgm_chord_mode,
                             tgm_component_extremes=tgm_component_extremes,
                             tgm_component_rescue=tgm_component_rescue,
                             tgm_rescue_chord=tgm_rescue_chord,
                             tgm_main_pair=tgm_main_pair,
                             tgm_main_pair_mode=tgm_main_pair_mode,
                             tgm_fv_zmax_margin=tgm_fv_zmax_margin,
                             tgm_fv_zmax_margin_interior=tgm_fv_zmax_margin_interior,
                             tgm_fv_x_margin=tgm_fv_x_margin,
                             tgm_fv_y_margin=tgm_fv_y_margin,
                             save_stm_fit=save_stm_fit,
                             pf_track_main_cluster_only=pf_track_main_cluster_only,
                             pf_track_bridged_clusters=pf_track_bridged_clusters,
                             pf_shower_vertex_barrier=pf_shower_vertex_barrier,
                             pf_shower_parent_precedence=pf_shower_parent_precedence,
                             pf_pi0_node_per_id=pf_pi0_node_per_id,
                             pf_pdg_name_prototype_fallback=pf_pdg_name_prototype_fallback,
                             pf_orphan_track_parentage=pf_orphan_track_parentage,
                             pf_orphan_audit_only=pf_orphan_audit_only,
                             pf_direct_when_touching=pf_direct_when_touching,
                             pf_touch_max=pf_touch_max,
                             pf_touch_cross_main=pf_touch_cross_main,
                             pf_touch_cross_max=pf_touch_cross_max,
                             pf_pseudo_gap_from_main=pf_pseudo_gap_from_main,
                             pf_unique_node_ids=pf_unique_node_ids,
                             pf_drop_stray_satellites=pf_drop_stray_satellites,
                             pf_orphan_confident_track=pf_orphan_confident_track,
                             pf_orphan_track_min_cm=pf_orphan_track_min_cm,
                             pf_track_owns_loose_vertex=pf_track_owns_loose_vertex,
                             unmerge_bundle_mode=unmerge_bundle_mode,
                             restore_demoted_mains=restore_demoted_mains,
                             require_provenance=require_provenance,
                             evaluate_demoted_mains=evaluate_demoted_mains,
                             tgm_exempt_demoted_main=tgm_exempt_demoted_main,
                             skip_cosmic_companions=skip_cosmic_companions,
                             cosmic_companion_min_length=cosmic_companion_min_length,
                             nu_fallback_demoted_mains=nu_fallback_demoted_mains,
                             sp_photon_flag=sp_photon_flag,
                             mip_dqdx=mip_dqdx,
                             stm_consistent_fv=stm_consistent_fv,
                             stm_accept_guards=stm_accept_guards,
                             stm_proton_muon_guard=stm_proton_muon_guard,
                             stm_cathode_guard=stm_cathode_guard,
                             stm_anode_dist_fix=stm_anode_dist_fix,
                             stm_second_track_guard=stm_second_track_guard,
                             stm_deficit_guard=stm_deficit_guard,
                             stm_vertex_kink_guard=stm_vertex_kink_guard,
                             stm_d66_cuts=stm_d66_cuts,
                             // Same offsets below the top face as clus.jsonnet's
                             // pr() defaults, re-anchored to pr_y_top.
                             cosmic_y_top_main=pr_y_top - 17,
                             cosmic_y_top_strict=pr_y_top - 15,
                             cosmic_y_top_loose=pr_y_top - 37,
                             cosmic_y_small_piece=pr_y_top - 67,
                             cathode_x=cathode_x,
                             cathode_kink_xcut=cathode_kink_xcut,
                             cathode_wide_kink_angle=cathode_wide_kink_angle,
                             cathode_wide_kink_skirt=cathode_wide_kink_skirt,
                             cathode_wide_kink_baseline=cathode_wide_kink_baseline,
                             shower_topo_demote_len=shower_topo_demote_len,
                             fit_exclusion=fit_exclusion,
                             graph_endpoint_strict=graph_endpoint_strict,
                             graph_endpoint_tol=graph_endpoint_tol,
                             oov_prototype_parity=oov_prototype_parity,
                             first_seg_local_pca=first_seg_local_pca,
                             other_seg_relaxed_accept=other_seg_relaxed_accept,
                             shower_topo_proto_dir=shower_topo_proto_dir,
                             vertex_dir_use_fit_point=vertex_dir_use_fit_point,
                             shower_traj_recheck_parity=shower_traj_recheck_parity,
                             main_vertex_require_descriptor=main_vertex_require_descriptor,
                             main_vertex_candidate_flag=main_vertex_candidate_flag,
                             cont_muon_dir3_30cm=cont_muon_dir3_30cm,
                             track_comp_empty_abstain=track_comp_empty_abstain,
                             shower_topo_reset=shower_topo_reset,
                             reclass_preserve_4mom=reclass_preserve_4mom,
                             dir_track_median_local=dir_track_median_local,
                             examine_showers_vertex_by_index=examine_showers_vertex_by_index,
                             steiner_terminal_wire_tol=steiner_terminal_wire_tol,
                             steiner_terminal_adjacent_slice=steiner_terminal_adjacent_slice,
                             steiner_edge_charge_forward_dead_mix=steiner_edge_charge_forward_dead_mix,
                             iso_endpoint=iso_endpoint,
                             iso_endpoint_min_length=iso_endpoint_min_length,
                             iso_endpoint_max_xext=iso_endpoint_max_xext,
                             iso_endpoint_xext_frac=iso_endpoint_xext_frac,
                             iso_endpoint_xext_quantile=iso_endpoint_xext_quantile,
                             iso_endpoint_tube_radius=iso_endpoint_tube_radius,
                             iso_endpoint_min_aspect=iso_endpoint_min_aspect,
                             v3_extension_guard=v3_extension_guard,
                             traj_cover_probe=traj_cover_probe,
                             pr_find_other_rounds=pr_find_other_rounds,
                             v3_extension_min_gain=v3_extension_min_gain,
                             protect_graph_name=protect_graph_name,
                             protect_skip_convicted=protect_skip_convicted,
                             protect_open_convicted_bundles=protect_open_convicted_bundles,
                             protect_cathode_x=protect_cathode_x,
                             protect_cathode_rejoin_xcut=protect_cathode_rejoin_xcut,
                             protect_cathode_rejoin_dyz=protect_cathode_rejoin_dyz,
                             protect_cathode_rejoin_dis=protect_cathode_rejoin_dis,
                             protect_cathode_rejoin_perp=protect_cathode_rejoin_perp,
                             protect_cathode_rejoin_angle=protect_cathode_rejoin_angle,
                             protect_cathode_rejoin_dir_radius=protect_cathode_rejoin_dir_radius,
                             protect_cathode_rejoin_dir_npts=protect_cathode_rejoin_dir_npts,
                             vertex_z_prior_scale=vertex_z_prior_scale,
                             ssm_target_dir=ssm_target_dir,
                             ssm_absorber_dir=ssm_absorber_dir,
                             kine_fudge_factor=kine_fudge_factor,
                             kine_recom_factor=kine_recom_factor,
                             kine_shower_fudge_factor=kine_shower_fudge_factor,
                             kine_shower_recom_factor=kine_shower_recom_factor,
                             kine_proton_recom_factor=kine_proton_recom_factor,
                             kine_plane_weights=kine_plane_weights,
                             kine_plane_asym_switch=kine_plane_asym_switch,
                             kine_w_value=kine_w_value,
                             kine_shower_pdg_live=kine_shower_pdg_live,
                             neutrino_consistent_fv=neutrino_consistent_fv,
                             cosmic_consistent_fv=cosmic_consistent_fv,
                             nue_sp_consistent_fv=nue_sp_consistent_fv,
                             sp_sce_correction=sp_sce_correction,
                             tagger_ordered_segment_sets=tagger_ordered_segment_sets,
                             stem_endpoint_wcpt_parity=stem_endpoint_wcpt_parity,
                             broken_muon_cluster_id_count=broken_muon_cluster_id_count,
                             neutrino_type_bitmask=neutrino_type_bitmask,
                             nu_per_bundle=nu_per_bundle,
                             nu_per_bundle_min_length=nu_per_bundle_min_length,
                             nu_selected_as_main=nu_selected_as_main,
                             nu_selected_as_main_snapshot_all=nu_selected_as_main_snapshot_all,
                             daughter_count_proto_main_vertex=daughter_count_proto_main_vertex,
                             daughter_count_proto_examine_showers=daughter_count_proto_examine_showers,
                             shower_pdg_from_start_segment=shower_pdg_from_start_segment,
                             shower_pdg_from_shower_type=shower_pdg_from_shower_type,
                             shower_pdg_exact_muon_test=shower_pdg_exact_muon_test,
                             pi0_id_shared_allocator=pi0_id_shared_allocator,
                             shower_flag_pdg_electron=shower_flag_pdg_electron,
                             shower_less_id_tiebreak=shower_less_id_tiebreak,
                             shower_endpoint_exclude_start_vertex=shower_endpoint_exclude_start_vertex,
                             shower_endpoint_skip_orphan_vtx=shower_endpoint_skip_orphan_vtx,
                             shower_walk_visited_parity=shower_walk_visited_parity,
                             track_pid_persist_dqdx=track_pid_persist_dqdx,
                             shower_reclass_dqdx_guard=shower_reclass_dqdx_guard,
                             shower_topo_dqdx_guard=shower_topo_dqdx_guard,
                             track_pid_persist_4mom=track_pid_persist_4mom,
                             shower_proton_daughter_pion=shower_proton_daughter_pion,
                             shower_proton_daughter_pion_dissolve=shower_proton_daughter_pion_dissolve,
                             muon_multi_proton_pion=muon_multi_proton_pion,
                             track_pid_persist_dqdx_electron_guard=track_pid_persist_dqdx_electron_guard,
                             shower_connect_main_vertex_straight_guard=shower_connect_main_vertex_straight_guard,
                             shower_traj_straight_guard=shower_traj_straight_guard,
                             shower_absorb_track_guard=shower_absorb_track_guard,
                             shower_connect_protected_pion_guard=shower_connect_protected_pion_guard,
                             michel_stem_muon_rescue=michel_stem_muon_rescue,
                             shower_in_cascade_guard=shower_in_cascade_guard,
                             shower_in_max_len=shower_in_max_len,
                             shower_in_mip_hi=shower_in_mip_hi,
                             shower_connect_from_vertices_straight_guard=shower_connect_from_vertices_straight_guard,
                             shower_connect_start_seg_straight_guard=shower_connect_start_seg_straight_guard,
                             examine_direction_dirsign_shower_in_guard=examine_direction_dirsign_shower_in_guard,
                             daughter_shower_angle_reclass_straight_guard=daughter_shower_angle_reclass_straight_guard,
                             shower_topo_reexam_straight_guard=shower_topo_reexam_straight_guard,
                             sfv_kink_max=sfv_kink_max,
                             shower_nv_bridge_track=shower_nv_bridge_track,
                             shower_nv_main_pi_init=shower_nv_main_pi_init,
                             shower_nv_bridge_max_gap=shower_nv_bridge_max_gap,
                             kine_drop_stray_satellites=kine_drop_stray_satellites,
                             kine_sat_min_energy=kine_sat_min_energy,
                             kine_sat_prox_max=kine_sat_prox_max,
                             kine_sat_angle_bad=kine_sat_angle_bad,
                             kine_sat_angle_main=kine_sat_angle_main,
                             kine_sat_far_dis=kine_sat_far_dis,
                             kine_sat_axis_dis_cut=kine_sat_axis_dis_cut,
                             kine_sat_cont_kink=kine_sat_cont_kink,
                             kine_sat_track_max_nseg=kine_sat_track_max_nseg,
                             kine_sat_em_far_dis=kine_sat_em_far_dis,
                             michel_stem_michel_check=michel_stem_michel_check,
                             michel_stem_max_far_len=michel_stem_max_far_len,
                             shower_stem_backfill=shower_stem_backfill,
                             stem_backfill_max_len=stem_backfill_max_len,
                             stem_backfill_mip_lo=stem_backfill_mip_lo,
                             stem_backfill_mip_hi=stem_backfill_mip_hi,
                             stem_backfill_min_shower_len=stem_backfill_min_shower_len,
                             shower_conn3_unreachable=shower_conn3_unreachable,
                             conn3_unreachable_min_len=conn3_unreachable_min_len,
                             conn3_stitch_max=conn3_stitch_max,
                             shower_dedup_start_seg=shower_dedup_start_seg,
                             shower_traj_michel_stem=shower_traj_michel_stem,
                             michel_stem_traj_min_len=michel_stem_traj_min_len,
                             michel_stem_traj_max_len=michel_stem_traj_max_len,
                             michel_stem_traj_mip_lo=michel_stem_traj_mip_lo,
                             michel_stem_traj_max_far_len=michel_stem_traj_max_far_len,
                             michel_stem_traj_min_kink_deg=michel_stem_traj_min_kink_deg,
                             shower_long_muon_keep_type=shower_long_muon_keep_type,
                             shower_bragg_protect_start_segment=shower_bragg_protect_start_segment,
                             shower_reclass_case_b_dqdx_guard=shower_reclass_case_b_dqdx_guard,
                             shower_accept_pid_guard=shower_accept_pid_guard,
                             shower_pid_guard_min_len=shower_pid_guard_min_len,
                             shower_vote_track_pid_counts=shower_vote_track_pid_counts,
                             shower_cone_absorb_guard=shower_cone_absorb_guard,
                             shower_detach_track_stem=shower_detach_track_stem,
                             shower_ghost_member_drop=shower_ghost_member_drop,
                             shower_ghost_overlap_frac=shower_ghost_overlap_frac,
                             shower_ghost_dqdx_ratio=shower_ghost_dqdx_ratio,
                             shower_ghost_min_len=shower_ghost_min_len,
                             kine_charge_dedup=if kine_charge_dedup == null then false else kine_charge_dedup,
                             kine_charge_rebuild=if kine_charge_rebuild == null then false else kine_charge_rebuild,
                             kine_charge_track_ctx=if kine_charge_track_ctx == null then false else kine_charge_track_ctx,
                             kine_mass_rules=if kine_mass_rules == null then false else kine_mass_rules,
                             kine_hadronic_dqdx=if kine_hadronic_dqdx == null then false else kine_hadronic_dqdx,
                             kine_long_muon_mode=kine_long_muon_mode,
                             kine_long_muon_ratio_lo=kine_long_muon_ratio_lo,
                             kine_long_muon_ratio_hi=kine_long_muon_ratio_hi,
                             kine_mainvtx_used_guard=if kine_mainvtx_used_guard == null then false else kine_mainvtx_used_guard,
                             shower_hadronic_tag=if shower_hadronic_tag == null then false else shower_hadronic_tag,
                             shower_hadronic_min_len=shower_hadronic_min_len,
                             shower_hadronic_scan_len=shower_hadronic_scan_len,
                             shower_hadronic_bin=shower_hadronic_bin,
                             shower_hadronic_r_cyl=shower_hadronic_r_cyl,
                             shower_hadronic_r_core=shower_hadronic_r_core,
                             shower_hadronic_growth_max=shower_hadronic_growth_max,
                             shower_hadronic_growth_bragg=shower_hadronic_growth_bragg,
                             shower_hadronic_bragg_ratio=shower_hadronic_bragg_ratio,
                             shower_hadronic_stem_ratio=shower_hadronic_stem_ratio,
                             kine_count_orphan_tracks=kine_count_orphan_tracks,
                             kine_orphan_track_min=kine_orphan_track_min,
                             straight_cont_cross_cluster=straight_cont_cross_cluster,
                             sccc_bridge_body=sccc_bridge_body,
                             sccc_max_gap=sccc_max_gap,
                             sccc_kink_max=sccc_kink_max,
                             sccc_gap_aligned=sccc_gap_aligned,
                             sccc_kink_tight=sccc_kink_tight,
                             single_muon_proton_chain_veto=single_muon_proton_chain_veto,
                             single_muon_long_muon_claim=single_muon_long_muon_claim,
                             pid_flag_reconcile=pid_flag_reconcile,
                             other_seg_empty_2d_guard=other_seg_empty_2d_guard,
                             long_muon_stub_bridge=long_muon_stub_bridge,
                             two_end_break=two_end_break,
                             teb_turn_min_arm_frac=teb_turn_min_arm_frac,
                             teb_second_max=teb_second_max,
                             teb_chain_topology=teb_chain_topology,
                             teb_r3_turn=teb_r3_turn,
                             teb_r3_hot=teb_r3_hot,
                             teb_bragg_veto_turn=teb_bragg_veto_turn,
                             kink_walk_dqdx_stop=kink_walk_dqdx_stop,
                             kink_break_protect=kink_break_protect,
                             fit_blob_coverage=fit_blob_coverage,
                             fit_blob_coverage_defer=fit_blob_coverage_defer,
                             vertex_kink_snap=vertex_kink_snap,
                             vks_radius=vks_radius,
                             vks_min_dis=vks_min_dis,
                             vks_angle=vks_angle,
                             vks_margin=vks_margin,
                             vks_collinear=vks_collinear,
                             vks_skirt=vks_skirt,
                             vks_baseline=vks_baseline,
                             vks_min_arm=vks_min_arm,
                             vks_fit_miss=vks_fit_miss,
                             vks_hot_ratio=vks_hot_ratio,
                             vks_carry_prong=vks_carry_prong,
                             vertex_junction_snap=vertex_junction_snap,
                             vjs_radius=vjs_radius,
                             vjs_min_arm=vjs_min_arm,
                             vjs_min_prongs=vjs_min_prongs,
                             vjs_collinear=vjs_collinear,
                             vjs_fit_margin=vjs_fit_margin,
                             vjs_fit_rms=vjs_fit_rms,
                             vjs_override_kink_snap=vjs_override_kink_snap,
                             vjs_min_move=vjs_min_move,
                             esva_ignore_empty_2d=esva_ignore_empty_2d,
                             main_vertex_graph_audit=main_vertex_graph_audit,
                             mvga_radius=mvga_radius,
                             mvga_dup_tol=mvga_dup_tol,
                             mvga_dup_frac=mvga_dup_frac,
                             mvga_dup_angle=mvga_dup_angle,
                             mvga_bridge_mip=mvga_bridge_mip,
                             mvga_reconnect=mvga_reconnect,
                             mvga_stub=mvga_stub,
                             mvga_stub_pts=mvga_stub_pts,
                             mvga_reseat_angle=mvga_reseat_angle,
                             mvga_satellite=mvga_satellite,
                             mvga_interposed=mvga_interposed,
                             mvga_interposed_angle=mvga_interposed_angle,
                             mvga_interposed_len=mvga_interposed_len,
                             mvga_sat_dup_frac=mvga_sat_dup_frac,
                             mvga_interposed_deg1=mvga_interposed_deg1,
                             mvga_splice_straighten=mvga_splice_straighten,
                             mvga_approach_collapse=mvga_approach_collapse,
                             mvga_straighten_radius=mvga_straighten_radius,
                             mvga_op1_radius=mvga_op1_radius,
                             mvga_op1_dup_frac=mvga_op1_dup_frac,
                             mvga_op1_post=mvga_op1_post,
                             mvga_carry_max=mvga_carry_max,
                             swap_orphan_dup_audit=swap_orphan_dup_audit,
                             mvga_proj_dup_frac=mvga_proj_dup_frac,
                             mvga_proj_dqdx_ratio=mvga_proj_dqdx_ratio,
                             mvga_proj_angle=mvga_proj_angle,
                             mvga_ac_veto_radius=mvga_ac_veto_radius,
                             mvga_ac_chord_max=mvga_ac_chord_max,
                             mvga_ac_no_cascade=mvga_ac_no_cascade,
                             mvga_passthru=mvga_passthru,
                             mvga_passthru_tol=mvga_passthru_tol,
                             mvga_interposed_fallback=mvga_interposed_fallback,
                             mvga_interposed_fallback_min_angle=mvga_interposed_fallback_min_angle,
                             mvga_dup_starved_asym=mvga_dup_starved_asym,
                             mvga_dup_starved_mip=mvga_dup_starved_mip,
                             mvga_dup_starved_span=mvga_dup_starved_span,
                             dl_vtx_swap_guard=dl_vtx_swap_guard,
                             dl_vtx_cloud_no_exclusion=dl_vtx_cloud_no_exclusion,
                             dl_vtx_dual_chain=dl_vtx_dual_chain,
                             dual_chain_mode=dual_chain_mode,
                             dual_chain_transfer=dual_chain_transfer,
                             dual_chain_transfer_max=dual_chain_transfer_max,
                             dual_chain_allow_cluster_swap=dual_chain_allow_cluster_swap,
                             dual_chain_vtx_weight=dual_chain_vtx_weight,
                             dqdx_fit_keep_all_points=dqdx_fit_keep_all_points,
                             dl_vtx_topo_weight=dl_vtx_topo_weight,
                             dl_vtx_topo_center=dl_vtx_topo_center,
                             main_vertex_swap_apply=main_vertex_swap_apply,
                             rough_path_probe=rough_path_probe,
                             steiner_gap_penalty=steiner_gap_penalty,
                             sgp_dead_alpha=sgp_dead_alpha,
                             sgp_min_edge=sgp_min_edge,
                             sgp_sample_step=sgp_sample_step,
                             sgp_point_radius=sgp_point_radius,
                             sgp_weak_scale=sgp_weak_scale,
                             sgp_weak_qref=sgp_weak_qref,
                             sgp_edge_probe=sgp_edge_probe,
                             sgp_max_sep=sgp_max_sep,
                             break_seg_orient=break_seg_orient,
                             vertex_scoreboard=vertex_scoreboard,
                             dl_vtx_harvest=dl_vtx_harvest,
                             mvfit_robust=mvfit_robust,
                             mvfit_main_only=mvfit_main_only,
                             mvfit_min_len=mvfit_min_len,
                             mvfit_rin_margin=mvfit_rin_margin,
                             mvfit_rout_frac=mvfit_rout_frac,
                             mvfit_rout_min=mvfit_rout_min,
                             mvfit_rout_max=mvfit_rout_max,
                             mvfit_angle=mvfit_angle,
                             mvfit_min_pts=mvfit_min_pts,
                             mvfit_min_aniso=mvfit_min_aniso,
                             mvfit_prior_range=mvfit_prior_range,
                             other_seg_keep_isolated=other_seg_keep_isolated,
                             other_seg_keep_isolated_min_points=other_seg_keep_isolated_min_points,
                             other_seg_keep_isolated_min_length=other_seg_keep_isolated_min_length,
                             other_seg_keep_isolated_min_nnf=other_seg_keep_isolated_min_nnf,
                             other_seg_keep_isolated_len_admit=other_seg_keep_isolated_len_admit,
                             other_seg_uncover_3d=other_seg_uncover_3d,
                             iso_snap_min_dir_mag=iso_snap_min_dir_mag,
                             shower_absorb_unreachable_main=shower_absorb_unreachable_main,
                             assoc_full_recluster=assoc_full_recluster,
                             assoc_reassign_orphans=assoc_reassign_orphans,
                             assoc_clear_on_merge=assoc_clear_on_merge,
                             es3_stub_guard=es3_stub_guard,
                             es3sg_stub_max=es3sg_stub_max,
                             es3sg_len_ratio=es3sg_len_ratio,
                             es3sg_ang3_min=es3sg_ang3_min,
                             es3sg_ang_ratio=es3sg_ang_ratio,
                             es3sg_require_terminal=es3sg_require_terminal,
                             pseudo_shower_track_paint=pseudo_shower_track_paint,
                             reclass_never_computed_ke_floor=reclass_never_computed_ke_floor,
                             muon_dqdx_curve=muon_dqdx_curve,
                             use_power_recomb=use_power_recomb,
                             sp_dedx_use_recomb_model=sp_dedx_use_recomb_model,
                             sp_mean_dedx_cut=sp_mean_dedx_cut,
                             dl_vtx_cut=dl_vtx_cut);

    local graph = g.intern(
        innodes=[source],
        outnodes=[pr],
        edges=[g.edge(source, pr, 0, 0)],
    );

    local app = {
        type: 'Pgrapher',
        data: { edges: g.edges(graph) },
    };

    // WireCellRoot hosts SbndMagnifyTrackingVisitor and the PR output
    // components (SbndPrMagnifyTrackingVisitor, UbooneTaggerOutputVisitor,
    // the Uboone BDT scorers).  It used to be loaded only when one of those
    // was in the job, so the default compiled config stayed byte-identical.
    // That is no longer possible: master commit 8d069234 ("add sce and
    // dumper", merged here in 87ada3d5) put SCEFieldTH3 -- which lives in
    // WireCellRoot -- unconditionally into the DetectorVolumes `uses` of
    // pgrapher/experiment/sbnd/clus.jsonnet, so EVERY job importing that
    // module instantiates it.  Without the plugin the job aborts at
    // configure time with "No factory for class SCEFieldTH3" (the Q/L job
    // never saw this because it always loads WireCellRoot).  The component
    // is inert here -- both SBND realities set use_sce=false -- so loading
    // the plugin only pays the SCE-map open, and no reconstruction changes.
    // The narrower alternative is to gate that `uses` entry on use_sce in
    // clus.jsonnet, which also moves the Q/L job's compiled config.
    local needs_root = true;

    local cmdline = {
        type: 'wire-cell',
        data: {
            plugins: ['WireCellGen', 'WireCellPgraph', 'WireCellAux', 'WireCellSio',
                      'WireCellSigProc', 'WireCellImg', 'WireCellClus']
                     + (if needs_root then ['WireCellRoot'] else []),
            apps: ['Pgrapher'],
        },
    };

    [cmdline] + g.uses(graph) + [app]
