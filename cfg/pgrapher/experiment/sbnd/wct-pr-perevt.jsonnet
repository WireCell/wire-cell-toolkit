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
    // SBND operating point for the four pr/30 findings.  ALL FIVE ARE AT THE
    // LEGACY VALUE: this file is the single source of the SBND operating point
    // (doc 68), so a bare run reproduces production exactly, and the knob-on
    // arms are produced by editing a COPY of cfg/ (WIRECELL_PATH override),
    // never by flipping a default here.
    //   fit_exclusion (P1)            true  => the 27 knobbed do_multi_tracking
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
    fit_exclusion = false,
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
    // protect_bundle stage knobs (doc pr/23): the PR-stage overclustering
    // protection (uboone's second graph examination).  The stage is in
    // pipeline_names by DEFAULT since the sec 9 production flip; the knobs
    // are inert only when it is dropped.  The cathode re-join
    // values are the SBND operating point in INTERNAL units (unlike
    // cathode_kink_xcut above, which is cm), matching clus.jsonnet's pr()
    // defaults so a bare run is production; nulls = prototype-faithful
    // (re-join pass disabled -- would break cathode crossers, doc pr/20).
    protect_graph_name          = null,   // null => 'relaxed'
    // C++ default true (doc pr/23 ordering): a TGM/STM/lm-convicted in-window
    // main does not open its bundle for splitting.  null => key omitted.
    protect_skip_convicted      = null,
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
    pf_shower_vertex_barrier = true,
    pf_shower_parent_precedence = true,
    pf_pi0_node_per_id = true,
    pf_pdg_name_prototype_fallback = true,
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
                             pf_shower_vertex_barrier=pf_shower_vertex_barrier,
                             pf_shower_parent_precedence=pf_shower_parent_precedence,
                             pf_pi0_node_per_id=pf_pi0_node_per_id,
                             pf_pdg_name_prototype_fallback=pf_pdg_name_prototype_fallback,
                             unmerge_bundle_mode=unmerge_bundle_mode,
                             restore_demoted_mains=restore_demoted_mains,
                             require_provenance=require_provenance,
                             evaluate_demoted_mains=evaluate_demoted_mains,
                             tgm_exempt_demoted_main=tgm_exempt_demoted_main,
                             skip_cosmic_companions=skip_cosmic_companions,
                             cosmic_companion_min_length=cosmic_companion_min_length,
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
                             protect_graph_name=protect_graph_name,
                             protect_skip_convicted=protect_skip_convicted,
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
                             sp_sce_correction=sp_sce_correction,
                             tagger_ordered_segment_sets=tagger_ordered_segment_sets,
                             stem_endpoint_wcpt_parity=stem_endpoint_wcpt_parity,
                             broken_muon_cluster_id_count=broken_muon_cluster_id_count,
                             neutrino_type_bitmask=neutrino_type_bitmask,
                             daughter_count_proto_main_vertex=daughter_count_proto_main_vertex,
                             daughter_count_proto_examine_showers=daughter_count_proto_examine_showers,
                             shower_pdg_from_start_segment=shower_pdg_from_start_segment,
                             shower_pdg_from_shower_type=shower_pdg_from_shower_type,
                             shower_pdg_exact_muon_test=shower_pdg_exact_muon_test,
                             pi0_id_shared_allocator=pi0_id_shared_allocator,
                             shower_flag_pdg_electron=shower_flag_pdg_electron,
                             shower_less_id_tiebreak=shower_less_id_tiebreak,
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
