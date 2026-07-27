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
    DL             = 6.5781,
    DT             = 13.1349,
    // lifetime: electron lifetime in ms.  INERT in this chain -- it only feeds
    // the sim Drifter's attenuation, and no reco component reads it (grep the
    // compiled config: zero 'lifetime' keys).  NO lifetime / charge-attenuation
    // correction is applied anywhere in imaging, clustering, Q/L matching or PR.
    // The 6.0 is an inherited placeholder from the first standalone Q/L test
    // (wcp-porting-img 655bd6a: "DL=6.2 DT=9.8 lifetime=6 driftSpeed=1.565"); the
    // DL/DT of that same triple were later corrected to the SBND physical values
    // (9f498089), lifetime never was because nothing consumes it.  SBND
    // simparams.jsonnet says 35 ms.  See sbnd_xin/docs/64_cfg-sync.md sec 4.
    lifetime       = 6.0,
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
    pipeline_names = ['switch_scope', 'unmerge_bundle', 'unmerge_assoc', 'steiner',
                      'fiducialutils', 'tagger_check_tgm', 'tagger_check_stm',
                      'tagger_check_fc'],
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
    // When set, re-save the (post-PR) tree to this TensorDM tarball.  Used by the
    // round-trip gate (re-save with pipeline_names=[] and compare member hashes
    // against the input tarball).
    save_tensors   = '',
    // SCN vertex weights (WIRECELL_PATH-resolved, e.g.
    // 'uboone/scn_vtx/t48k-m16-l5-lr5d-res0.5-CP24.pth').  '' = geometric vertex.
    // NOTE: only uBooNE-trained weights exist; SBND use is an untuned demo.
    dl_weights     = '',
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
                             mip_dqdx=mip_dqdx,
                             stm_consistent_fv=stm_consistent_fv,
                             stm_accept_guards=stm_accept_guards,
                             stm_proton_muon_guard=stm_proton_muon_guard,
                             stm_cathode_guard=stm_cathode_guard,
                             stm_anode_dist_fix=stm_anode_dist_fix,
                             stm_second_track_guard=stm_second_track_guard,
                             stm_deficit_guard=stm_deficit_guard,
                             stm_vertex_kink_guard=stm_vertex_kink_guard);

    local graph = g.intern(
        innodes=[source],
        outnodes=[pr],
        edges=[g.edge(source, pr, 0, 0)],
    );

    local app = {
        type: 'Pgrapher',
        data: { edges: g.edges(graph) },
    };

    local cmdline = {
        type: 'wire-cell',
        data: {
            plugins: ['WireCellGen', 'WireCellPgraph', 'WireCellAux', 'WireCellSio',
                      'WireCellSigProc', 'WireCellImg', 'WireCellClus']
                     // WireCellRoot hosts SbndMagnifyTrackingVisitor; loaded
                     // only with the knob so the compiled config stays
                     // byte-identical when off.
                     + (if save_stm_fit then ['WireCellRoot'] else []),
            apps: ['Pgrapher'],
        },
    };

    [cmdline] + g.uses(graph) + [app]
