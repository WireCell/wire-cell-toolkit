// SBND imaging + clustering + charge/light (Q/L) matching, faithfully following
// Xin's standalone chain (sbnd_xin/wct-clus-matching-standalone.jsonnet +
// wct-img-all.jsonnet) but reading the artROOT DIRECTLY instead of dumped files:
//   - charge: wclsCookedFrameSource (recob::Wire)        [vs ClusterFileSource npz]
//   - light : wclsOpFlashSource per TPC (ITensorSet)     [vs TensorFileSource tar.gz]
//
// Imaging uses the TOOLKIT img.jsonnet 'multi-3view' with full_deghost=true: one
// per-anode pipe emits the live (active, multi_active 2-view-recovering) clusters
// on port 0 and the dead (masked, multi_masked_2view) clusters on port 1 -- the
// same live+dead views Xin produces with wct-img-all.jsonnet.  These feed the
// toolkit clus.jsonnet per_apa (PointTreeBuilding live/dead -> MABC).
//
// Matching uses the canonical FlashTensorToOpticalPCs + QLMatching (WireCellMatch),
// JOINT over both APAs (one node, premerged all_apa), exactly as Xin's default.
//
// All MABC nodes (per-APA + all-APA) write into ONE shared Bee zip: mabc.zip.
//
// Requires the toolkit cfg + photodet on WIRECELL_PATH -- source sbnd/setup-ap.sh.
// Nothing under sbnd_xin is imported or modified (we import the in-tree canonical
// pgrapher/experiment/sbnd/{img,clus,qlmatching,cathode_fiducial}.jsonnet directly).

local g = import "pgraph.jsonnet";
local wc = import "wirecell.jsonnet";

local tools_maker = import 'pgrapher/common/tools.jsonnet';

local reality = std.extVar('reality');

// The reco chain's reality config (use_sce, pos_offset_on) is grouped inside
// the toolkit clus maker keyed by `reality` -- see clus.jsonnet `reco` local.

// The follow-up tail -- labeler_truth -> STM/TGM/FC PR taggers -> labeler_tagger -> dump --
// is assembled below in this entry config (see the "follow-up tail" section);
// clus.jsonnet now produces only the clustering+matching all-APA MABC.  Requires
// the "WireCellAIML" plugin + BOTH "wclsTensorSetLabeler:labeler_truth" and
// ":labeler_tagger" inputers in the fcl (both fcls have them); the labelers run
// in BOTH realities (sim attaches nu truth + truth_per_track + truth Bee sets +
// nugraph-with-truth, data attaches RSE + input-only HDF5).  The
// tagger_stm/tgm/fc Bee sets are emitted in both realities.

// Canonical SBND simparams (toolkit).  drift_speed (1.563 mm/us) flows into
// QLMatching from here; no DL/DT/lifetime/driftSpeed extVars are needed.
local params = import 'pgrapher/experiment/sbnd/simparams.jsonnet';

local tools_all = tools_maker(params);
local tools = tools_all {anodes: [tools_all.anodes[n] for n in [0, 1]]};
local nanode = std.length(tools.anodes);

// ---- artROOT inputs (must match the names used in the fcl inputers) ----
local wcls_input = {
    sigs: g.pnode({
        type: 'wclsCookedFrameSource',
        name: 'sigs',
        data: {
            nticks: params.daq.nticks,
            frame_scale: 50,                              // scale up input recob::Wire
            summary_scale: 50,                            // scale up input summary
            frame_tags: ["orig"],
            recobwire_tags: std.extVar('recobwire_tags'),
            trace_tags: std.extVar('trace_tags'),
            summary_tags: std.extVar('summary_tags'),
            input_mask_tags: std.extVar('input_mask_tags'),
            output_mask_tags: std.extVar('output_mask_tags'),
        },
    }, nin=0, nout=1),
    // wclsOpFlashSource emits an ITensorSet [nflashes x (npmts+1)] (npmts=312 ->
    // 313 cols) -- exactly the opflash matrix FlashTensorToOpticalPCs reads on
    // port 1, so it feeds flash_attach directly (no opflash file dump/read).
    opflashes0: g.pnode({
        type: 'wclsOpFlashSource',
        name: 'tpc0',
        data: { art_tag: std.extVar('opflash0_input_label') },
    }, nin=0, nout=1),
    opflashes1: g.pnode({
        type: 'wclsOpFlashSource',
        name: 'tpc1',
        data: { art_tag: std.extVar('opflash1_input_label') },
    }, nin=0, nout=1),
    opflashes: [wcls_input.opflashes0, wcls_input.opflashes1],
};

// ---- frame fanout: one sigs frame -> per-anode frame with the per-anode
// 'gauss<id>'/'wiener<id>' trace tags the imaging pre_proc expects ----
local fanout_rules = [
    {
        frame: { '.*': 'orig%d' % tools.anodes[n].data.ident },
        trace: {
            gauss: 'gauss%d' % tools.anodes[n].data.ident,
            wiener: 'wiener%d' % tools.anodes[n].data.ident,
        },
    }
    for n in std.range(0, nanode - 1)
];
local frame_fan = g.pnode({
    type: 'FrameFanout',
    name: 'sig_fanout',
    data: { multiplicity: nanode, tag_rules: fanout_rules },
}, nin=1, nout=nanode);

// ---- imaging: toolkit multi-3view + full_deghost (live=port0, dead=port1) ----
local img = import 'pgrapher/experiment/sbnd/img.jsonnet';
local img_maker = img();
local img_pipes = [
    img_maker.per_anode(tools.anodes[n], 'multi-3view', add_dump=false, full_deghost=true)
    for n in std.range(0, nanode - 1)
];

// ---- clustering: toolkit per_apa, shared Bee ----
local clus = import 'pgrapher/experiment/sbnd/clus.jsonnet';
// rse_from_ident: each frame's tensor ident carries the real event id, so the
// Bee display is labelled with the true event number (matches Xin's chain).
// rse_from_metadata=true: every MABC prefers the run/subrun/event carried in
// its input tensor-set metadata, stamped by the wclsTensorSetMetadataAttacher
// nodes below.  rse_from_ident stays on as the fallback for any node whose
// input carries no metadata (C++ precedence: metadata > ident > config), so
// nothing regresses if an attacher is removed.
local clus_maker = clus(rse_from_ident=true, rse_from_metadata=true, reality=reality);
// Single shared Bee sink: every MABC node (per-APA + all-APA) writes into this
// one zip instead of one zip per node.
local bee_shared = {
    type: 'BeeSink',
    name: 'mabc_shared',
    data: { outname: 'mabc.zip', initial_index: 0 },
};
// ---- RSE metadata attachers (larwirecell wclsTensorSetMetadataAttacher) ----
// The WCT graph cannot otherwise know the art run/subrun: a tensor ident is the
// EVENT number alone (CookedFrameSource -> SimpleFrame(event.event(), ...)),
// which is why MABC's rse_from_ident writes run=subrun=0.  These nodes stamp
// runNo/subRunNo/eventNo into the tensor-set metadata; each downstream MABC
// picks them up via rse_from_metadata and forwards them to its own output.
//
// They are O(1): the output shares the input's ITensor vector by pointer, with
// no point-cloud round trip, so instantiating three costs microseconds.
//
// Placement -- one immediately upstream of every MABC whose RSE matters:
//   rse_apa0/1  after PointTreeBuilding, before the per-APA MABCs
//   rse_all_apa after the joint matcher, before clus_all_apa
// clus_pr needs none: labeler_truth sits directly upstream and stamps the same
// three keys from its own art::Event.
// Each MUST also be listed in the fcl 'inputers' or it never sees an
// art::Event; it warns and passes through rather than stamping 0/0/0.
local rse_attach = [
    g.pnode({
        type: 'wclsTensorSetMetadataAttacher',
        name: 'rse_apa%d' % n,
        data: {},
    }, nin=1, nout=1)
    for n in std.range(0, nanode - 1)
];
local rse_all_apa = g.pnode({
    type: 'wclsTensorSetMetadataAttacher',
    name: 'rse_all_apa',
    data: {},
}, nin=1, nout=1);

// save_assoc_id=true: the per-APA clustering_isolated pass writes the isolated-
// grouping provenance (assoc_cluster_id/assoc_cluster_main); this flag
// homogenizes it so it SURVIVES serialization out of the per-APA MABC (the
// serializer silently drops heterogeneous perblob keys).  Without it the
// provenance is born here but dropped before matching_joint, so the all-APA
// save_assoc_cluster_id has nothing to preserve and the pr_node's unmerge_assoc
// finds nothing to split ("unmerged 0 main clusters") -- leaving the tagger
// main un-refined and mis-verdicted (MC 32-10-10 cluster 22: STM=0 un-split vs
// STM=1 after the isolated split, matching sbnd_xin's 2-step).  Pairs with the
// all_apa save_real/assoc_cluster_id below (real_cluster_id is instead created
// by all_apa's own examine_bundles, so it needs no per-APA flag).
local clus_pipes = [
    clus_maker.per_apa(tools.anodes[n], dump=false, bee_sink=bee_shared, save_assoc_id=true,
                       // Clustering-stage performance variants, part of the SBND
                       // production operating point but configured on the
                       // clustering entry points rather than through pr(), so the
                       // generated PR operating point cannot carry them -- the
                       // issue-17 gate is what finds them missing.
                       eb_fast=true, po_fast=true, dg_fast=true,
                       // spliced between PointTreeBuilding and the per-APA MABC
                       pre_mabc=rse_attach[n])
    for n in std.range(0, nanode - 1)
];

// ---- Q/L matching: canonical FlashTensorToOpticalPCs + joint QLMatching ----
local semimodel_file = 'semi-analytical-sbnd.json';
local pmt_nl = true;

local qlm = (import 'pgrapher/experiment/sbnd/qlmatching.jsonnet')(params);
local cathode_fv = (import 'pgrapher/experiment/sbnd/cathode_fiducial.jsonnet')(
    cx=0.5*wc.cm, cy=0.5*wc.cm, cz=0.5*wc.cm);
// cluster tree (port 0) + opflash matrix (port 1) -> tree with flash PCs attached
local flash_attach = [qlm.flash_attach(n) for n in std.range(0, nanode - 1)];
// ONE joint QLMatching node fed by both APAs' flash_attach outputs; it matches
// each APA and merges the per-APA trees, feeding the all-APA MABC directly.
local matching_joint = qlm.matching_joint(
    tools.anodes, clus_maker.detector_volumes(tools.anodes),
    reality, semimodel_file, cathode_fiducial=cathode_fv.tn, pmt_nl=pmt_nl);

// all-APA clustering + matching ONLY (joint matcher already merged ->
// premerged=true; dump=false -> no terminal sink here, the tail below consumes
// the MABC tensor output).  Bee (clustering-global etc.) still goes to the
// shared mabc.zip.
// save_real_cluster_id / save_assoc_cluster_id: attach the flash-merge (and
// isolated-grouping) provenance to the grouping so the pr_node's unmerge_bundle
// / unmerge_assoc can actually split the flash-merged bundles.  WITHOUT these,
// unmerge logs "no flash-merge provenance ... not split" and every cluster stays
// merged across the cathode -> the taggers judge the wrong (merged) geometry
// (e.g. nue evt 46363 cluster 19 spans both TPCs -> FC=false, vs the compact
// post-split main -> FC=true in sbnd_xin's 2-step chain).
// PR operating point -- RETIRED AS A SWITCH by doc sbnd_xin/120 sec 4.  All
// three former modes now produce the same graph, because the SBND production
// operating point lives in pr()'s OWN DEFAULTS
// (cfg/pgrapher/experiment/sbnd/clus.jsonnet) and this chain gets it simply by
// calling pr().  sbnd/pr-operating-point.jsonnet is deleted.
//
// Why each mode stopped meaning anything:
//
//   'sync' applied the GENERATED mirror.  Doc 120 sec 2 measured what that
//   mirror actually delivered: exact for the 242 knobs it carried, and silently
//   missing 25 that were added to the PR job after its last hand-run -- 17 pf_*
//   Bee knobs, six tracking-pr.root record keys, and nu_bundle_flash_group,
//   which is physics.  Its generator exists in neither repository.  A single
//   source cannot omit; a mirror can, and no check on the mirror can see it.
//
//   'preflip' claimed to reproduce the pre-2026-08-29 behaviour, in which this
//   chain called clus_maker.pr() with 7 arguments and inherited clus.jsonnet's
//   defaults.  That claim had ALREADY expired before doc 120 touched it: the arm
//   reaches its values by INHERITING defaults, so every pr() default that moved
//   since 2026-08-29 moved the arm with it -- the doc-118 trajectory flip on
//   2026-09-20 and doc 119's proj_pad flip on 2026-09-21 among them.
//   Reproducing the issue-16 MC and issue-18 data campaigns means checking out
//   the tree at that date, not selecting a mode whose meaning drifts.
//
//   'bare' was the regeneration baseline for the generator that no longer runs.
//
// Pinned to 'sync' rather than read, so that the two all-APA operating-point
// members below -- which still test it, and which no pr() default can reach --
// take the production value for EVERY caller.  That closes the mode out
// completely: an fcl that still sets pr_operating_point (the -preflip variant
// does) keeps compiling and now gets production, which is the whole point of
// the retirement.  The extVar is no longer consulted at all.
local pr_operating_point = 'sync';

local clus_all_apa = clus_maker.all_apa(tools.anodes, dump=false, eb_fast=true,
                                        bee_sink=bee_shared, premerged=true,
                                        save_real_cluster_id=true,
                                        save_assoc_cluster_id=true,
                                        // The two all-APA members of the SBND PR
                                        // operating point.  They live on
                                        // clus_all_apa(), not pr(), so no pr()
                                        // default can reach them -- and
                                        // two_chain_gate.py cannot see them
                                        // either, because the local chain builds
                                        // its all-APA stage in a DIFFERENT job
                                        // (wct-clus-matching-perevt.jsonnet) that
                                        // the gate does not compile.  Doc
                                        // sbnd_xin/120 sec 6 item 3 carries that
                                        // as an open item.  The test below is now
                                        // always true (the local is pinned to
                                        // 'sync'), so every caller gets these.
                                        save_bundle_main_provenance=
                                            pr_operating_point == 'sync',
                                        bee_flash_pred_min=
                                            if pr_operating_point == 'sync' then 0 else null);

// ==== follow-up tail (moved out of clus.jsonnet into this entry config) ====
//   clus_all_apa (MABC) -> pr_node (STM/TGM/FC taggers) -> labeler -> tail_dump

// Beam gate shared by the tagger PR pass AND the labeler's tagger-Bee cluster_id
// encoding (which uses it to tell "beam-window candidate" mains from out-of-
// window mains) -- ONE source so the two never drift.
local beam_window = [0.2 * wc.us, 2.2 * wc.us];

// PR tagger pass: SBND production operating point (clus_maker.pr defaults),
// running the nusel tagger visitors so each cluster carries flag_STM/TGM/FC for
// the labeler to read.  dump=false -> a pass-through tensor node.
local pds = (import 'pgrapher/experiment/sbnd/particle_dataset.jsonnet')();
// enable_tracking_root: emit tracking-pr.root (tracking_visitor + tagger_output).
//
// MUST be false when running more than one event per lar process.  Both
// writers use ONE fixed filename -- SbndPrMagnifyTrackingVisitor opens it
// RECREATE and UbooneTaggerOutputVisitor opens it UPDATE, every event -- so
// `lar -n 5` leaves an 8 KB file with T_tagger and T_kine ABSENT and Trun
// holding only the last event.  It is silent: rc=0 and every other output is
// fine.  One event per process (the issue-11 harness pattern) is safe.
//
// CONSEQUENCE OF TURNING IT OFF: tagger_output is what writes T_tagger and
// T_kine, so you also lose numu_score / nue_score / kine_reco_Enu.  The Bee
// side is unaffected -- the merged "mc" node still carries Enu and both
// scores, and the tagger_* layers still carry the cosmic verdicts.
//
// Set via the fcl:
//   physics.producers.sig2img.wcls_main.params.enable_tracking_root: "false"
local enable_tracking_root =
    if std.extVar('enable_tracking_root') == 'false' then false else true;



local pr_pipeline_names =
    ['switch_scope', 'unmerge_bundle', 'unmerge_assoc', 'steiner',
     'fiducialutils', 'tagger_check_tgm', 'tagger_check_stm', 'tagger_check_fc',
     'protect_bundle', 'steiner_refresh', 'tagger_check_neutrino',
     'numu_bdt_scorer', 'nue_bdt_scorer']
    + (if enable_tracking_root then ['tracking_visitor', 'tagger_output'] else []);

// doc sbnd_xin/120 sec 4: ONE call, no mode switch and no mirror.  pr()'s own
// defaults ARE the SBND production operating point now, so naming no knobs here
// is what gets production -- and a knob added to the PR job can no longer fail
// to reach this chain, because there is no second list for it to be missing
// from.  Gate: two_chain_gate.py, which compiles this file and the local PR job
// and requires their shared components to agree key for key.
local pr_node =
        clus_maker.pr(
        tools.anodes, dump=false,
        // Share the ONE Bee zip.  Without this the PR display layers
        // (track_fit / shower_track / vertices -- the only place the fitted
        // trajectory, per-particle shower/track colouring and reconstructed
        // vertices appear) go to mabc-pr.zip, which the bulk harness deletes.
        // Sharing renames this node's two colliding sets to clustering-pr /
        // mc-pr; see clus.jsonnet.
        bee_sink=bee_shared,
        // FULL 15-stage SBND production PR chain, matching sbnd_xin's
        // run_pr_chain_batch.sh PIPELINE string exactly (docs/5-pr-chain-in-1step).
        // Ordering is load-bearing, not stylistic:
        //   * protect_bundle + steiner_refresh sit AFTER the cosmic taggers (uboone
        //     takes cosmic verdicts on UNSPLIT clusters, wire-cell-prod-stm.cxx:806)
        //     and BEFORE tagger_check_neutrino (Protect_Over_Clustering exists only
        //     in the nue executable, wire-cell-prod-nue.cxx:1322).  SBND production
        //     default since 2026-08-02.
        //   * steiner_refresh must immediately follow protect_bundle: it rebuilds
        //     only the steiner products the split purged.
        //   * nue_bdt_scorer after numu_bdt_scorer.
        //   * tagger_output after tracking_visitor -- it reopens tracking-pr.root
        //     in UPDATE mode.
        // The last four need the WireCellRoot plugin (both fcls load it).
        pipeline_names=['switch_scope', 'unmerge_bundle', 'unmerge_assoc', 'steiner',
                        'fiducialutils', 'tagger_check_tgm', 'tagger_check_stm', 'tagger_check_fc',
                        'protect_bundle', 'steiner_refresh', 'tagger_check_neutrino',
                        'numu_bdt_scorer', 'nue_bdt_scorer']
                       + (if enable_tracking_root then ['tracking_visitor', 'tagger_output'] else []),
        particle_dataset=pds.particle_dataset, extra_uses=pds.all, beam_window=beam_window);
        // doc sbnd_xin/120 sec 4.5: `iso_endpoint=true` used to be passed here, to
        // match the 2-step job (apc doc pr/24 sec 16).  It is NOT a pr() parameter
        // -- it is a tcn_knobs bag key -- so this call could never have compiled.
        // It survived because it sat in the 'preflip' branch, which production
        // never took and jsonnet therefore never evaluated.  The value is now
        // delivered by pr()'s own tcn_knobs default (clus.jsonnet:2107), which is
        // where it always belonged.

// wclsTensorSetLabeler (larwirecell "WireCellAIML" plugin).  Deps come from
// clus_maker primitives; reads raw (non-t0-corrected) blob coords, so
// drift_speed/time_offset/tick MUST match the BlobSampler.
local lbl_dv = clus_maker.detector_volumes(tools.anodes, face='');
local lbl_pcts = clus_maker.pc_transforms(lbl_dv);
local lbl_fv = clus_maker.fiducial_box();

// TWO labeler instances, one either side of pr_node (docs/5-pr-chain-in-1step
// sec 6.7).  The labeler does four separable jobs and only the tagger Bee sets
// need the PR verdicts, so truth is attached UPSTREAM of the PR stage:
//
//   clus_all_apa -> labeler_truth -> pr_node -> labeler_tagger -> tail_dump
//
// Why it matters: the PR stage SPLITS clusters (switch_scope, unmerge_bundle,
// unmerge_assoc, protect_bundle) and the MABC renumbers idents 1..N after every
// step, so a labeler placed after it attaches truth to a PR-configuration-
// dependent clustering.  Upstream, the truth labels are stable against every PR
// knob.  They still survive the splits: trackid is written into each BLOB
// node's "scalar" local PC, while switch_scope erases only the CLUSTER-level
// "perblob" dataset, and separate() moves blob nodes wholesale.
//
// The labeler passes every non-live tensor through untouched (it re-serializes
// only pointtrees/N/live), so the dead grouping still reaches pr_node's
// fiducialutils.
local labeler_common = {
    inpath: 'pointtrees/%d',
    grouping: 'live',
    reality: reality,
    anodes: [wc.tn(anode) for anode in tools.anodes],
    detector_volumes: wc.tn(lbl_dv),
    pc_transforms: wc.tn(lbl_pcts),
    drift_speed: clus_maker.drift_speed,
    time_offset: clus_maker.time_offset,
    tick: 0.5 * wc.us,
    // shift the priorSCE (true) depos true->reco before association.
    sce_field: wc.tn(clus_maker.sce_field_fwd),
    sce_correction: true,
    n_sample_truth_depo_sce: 1,
    pf_ke_min: 10 * wc.MeV,
    pf_fiducial: wc.tn(lbl_fv),
    pf_nu_only: true,
    truth_tracks_nu_only: true,
    bee_sink: wc.tn(bee_shared),
    // tagger_stm/tgm/fc Bee sets use the SAME corrected coords as
    // clustering_global (data x_t0cor/y_cor/z_cor, sim x_t0cor/y/z).
    tagger_coords: clus_maker.bee_coords,
    // beam gate for the 4-case tagger cluster_id (0 non-main, 1 out-of-window
    // main, 2 in-window untagged, 3 in-window tagged); matches pr's window.
    beam_window: beam_window,
};
local labeler_uses =
    tools.anodes + [lbl_dv, lbl_pcts, clus_maker.sce_field_fwd, lbl_fv, bee_shared];

// (A) BEFORE pr_node: owns the truth association (writes the per-blob trackid),
// the truth_per_track tensor, the truth/sed Bee sets and the nugraph HDF5.
// Emits NO tagger sets -- the flags do not exist yet, so they would be all-zero
// AND would collide with (B)'s entries in the shared zip (the tagger set names
// are hardcoded, not derived from bee_algorithm).
local labeler_truth = g.pnode({
    type: 'wclsTensorSetLabeler',
    name: 'labeler_truth',
    data: labeler_common {
        label_blobs: true,
        // NO 'pf': the truth particle tree is not written here.  It is handed to
        // pr_node via pf_metadata_key and comes back out as ONE merged 'mc',
        // with the reco flow hung under a "reco nu" node -- Bee shows only one
        // particle tree per event, so this is the only way to see both.
        bee_sets: ['truth', 'sed'],
        pf_metadata_key: 'bee_pf_truth',
        hdf5_output: true,
    },
}, nin=1, nout=1, uses=labeler_uses);

// (B) AFTER pr_node: emits ONLY the tagger verdict Bee sets, which need
// flag_STM/flag_TGM/flag_FC/lm_flag from the PR stage.
// label_blobs=false is REQUIRED, not an optimization: the trackid write-back is
// not reality-gated, so a second labeling pass would stamp trackid=-1 over (A)'s
// labels in data mode, and would pay the (expensive) depo projection twice in
// sim.  With it false the projection is skipped and the existing per-blob
// trackid is read instead.
local labeler_tagger = g.pnode({
    type: 'wclsTensorSetLabeler',
    name: 'labeler_tagger',
    data: labeler_common {
        label_blobs: false,
        bee_sets: ['tagger'],
        hdf5_output: false,
    },
}, nin=1, nout=1, uses=labeler_uses);

// Terminal sink for the labeled pctree tarball (dump_mode=false -> real tensors,
// the historical labeler deliverable).  Bee -> mabc.zip and nugraph HDF5 are
// written upstream by the MABC/labeler, so this tarball is just a byproduct.
local tail_dump = g.pnode({
    type: 'TensorFileSink',
    name: 'clus_all_apa',
    data: {
        outname: 'trash-all-apa.tar.gz',
        prefix: 'clustering_',
        dump_mode: false,
    },
}, nin=1, nout=0);

// ---- assemble the graph ----
//  sigs ─FrameFanout─┬─ img[0] ─(live,dead)─ clus[0] ─ flash_attach[0] ─┐
//                    └─ img[1] ─(live,dead)─ clus[1] ─ flash_attach[1] ─┤
//  opflash[n] ───────────────────────────────────────(port1) flash_attach[n]
//        flash_attach[n] ─(port n)─ matching_joint ─ rse_all_apa ─ clus_all_apa (all-APA MABC)
//        clus_all_apa ─ labeler_truth ─ pr_node (taggers) ─ labeler_tagger ─ tail_dump
local graph = g.intern(
    innodes=[wcls_input.sigs],
    centernodes=[frame_fan] + img_pipes + clus_pipes + flash_attach
                + [matching_joint] + wcls_input.opflashes
                + [rse_all_apa, clus_all_apa, labeler_truth, pr_node, labeler_tagger],
    outnodes=[tail_dump],
    edges=
        [g.edge(wcls_input.sigs, frame_fan, 0, 0)]
        + [g.edge(frame_fan, img_pipes[n], n, 0) for n in std.range(0, nanode - 1)]
        + [g.edge(img_pipes[n], clus_pipes[n], 0, 0) for n in std.range(0, nanode - 1)]   // live / active
        + [g.edge(img_pipes[n], clus_pipes[n], 1, 1) for n in std.range(0, nanode - 1)]   // dead / masked
        + [g.edge(clus_pipes[n], flash_attach[n], 0, 0) for n in std.range(0, nanode - 1)]
        + [g.edge(wcls_input.opflashes[n], flash_attach[n], 0, 1) for n in std.range(0, nanode - 1)]
        + [g.edge(flash_attach[n], matching_joint, 0, n) for n in std.range(0, nanode - 1)]
        + [g.edge(matching_joint, rse_all_apa, 0, 0)]
        + [g.edge(rse_all_apa, clus_all_apa, 0, 0)]
        + [g.edge(clus_all_apa, labeler_truth, 0, 0)]
        + [g.edge(labeler_truth, pr_node, 0, 0)]
        + [g.edge(pr_node, labeler_tagger, 0, 0)]
        + [g.edge(labeler_tagger, tail_dump, 0, 0)]
);

local app = {
  type: 'Pgrapher',
  data: { edges: g.edges(graph) },
};

local cmdline = {
    type: "wire-cell",
    data: {
        // wcls reads plugins from the FCL; this list is for standalone use.
        plugins: ["WireCellGen", "WireCellPgraph", "WireCellAux", "WireCellSio",
                  "WireCellSigProc", "WireCellImg", "WireCellRoot", "WireCellTbb",
                  "WireCellClus", "WireCellMatch"],
        apps: ["Pgrapher"]
    }
};

[cmdline] + g.uses(graph) + cathode_fv.configs + [app]
