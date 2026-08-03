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
                             unmerge_bundle_mode=unmerge_bundle_mode,
                             restore_demoted_mains=restore_demoted_mains,
                             require_provenance=require_provenance,
                             evaluate_demoted_mains=evaluate_demoted_mains,
                             skip_cosmic_companions=skip_cosmic_companions,
                             cosmic_companion_min_length=cosmic_companion_min_length,
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
                             protect_graph_name=protect_graph_name,
                             protect_skip_convicted=protect_skip_convicted,
                             protect_cathode_x=protect_cathode_x,
                             protect_cathode_rejoin_xcut=protect_cathode_rejoin_xcut,
                             protect_cathode_rejoin_dyz=protect_cathode_rejoin_dyz,
                             protect_cathode_rejoin_dis=protect_cathode_rejoin_dis,
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
    // the Uboone BDT scorers); loaded only when something in the job needs it
    // so the default compiled config stays byte-identical.
    local needs_root = save_stm_fit
        || std.member(pipeline_names, 'numu_bdt_scorer')
        || std.member(pipeline_names, 'nue_bdt_scorer')
        || std.member(pipeline_names, 'tracking_visitor')
        || std.member(pipeline_names, 'tagger_output');

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
