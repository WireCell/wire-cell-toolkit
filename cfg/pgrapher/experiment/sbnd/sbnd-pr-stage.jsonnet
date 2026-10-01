// The SBND pattern-recognition (PR) stage, as ONE definition shared by
//   wcls-img-clus-matching-pr-lib.jsonnet  (the LArSoft 1-step chain: PR runs in-process
//                                            after clus_all_apa + labeler_truth)
//   wct-pr.jsonnet                          (step 2 of the 2-step split: PR re-run
//                                            standalone from the step-1 ITensorSet tar)
// so the two cannot drift: both build their clus_pr MABC by calling node() below with
// their own clus maker and Bee sink.  The operating point itself is pr()'s own defaults
// (cfg/pgrapher/experiment/sbnd/clus.jsonnet, doc sbnd_xin/120 sec 4) -- no knob is set here.
local wc = import 'wirecell.jsonnet';

// dQ/dx (e/cm, SBND 0.5 kV/cm) + range LinterpFunctions.
local pds = (import 'pgrapher/experiment/sbnd/particle_dataset.jsonnet')();

{
    // Beam gate shared by the tagger PR pass AND the labeler's tagger-Bee cluster_id
    // encoding (which uses it to tell "beam-window candidate" mains from out-of-
    // window mains) -- ONE source so the two never drift.
    beam_window: [0.2 * wc.us, 2.2 * wc.us],

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
    // The last four need the WireCellRoot plugin.
    // enable_tracking_root=false drops tracking_visitor + tagger_output, i.e. no
    // tracking-pr.root, hence no T_tagger / T_kine (numu_score / nue_score /
    // kine_reco_Enu); the Bee "mc" node keeps Enu and both scores regardless.
    // tagger_bee=true appends the tagger verdict Bee sets (TaggerBeeVisitor) --
    // the standalone step-2 job, which has no art-side labeler_tagger to write them.
    pipeline_names(enable_tracking_root=true, tagger_bee=false)::
        ['switch_scope', 'unmerge_bundle', 'unmerge_assoc', 'steiner',
         'fiducialutils', 'tagger_check_tgm', 'tagger_check_stm', 'tagger_check_fc',
         'protect_bundle', 'steiner_refresh', 'tagger_check_neutrino',
         'numu_bdt_scorer', 'nue_bdt_scorer']
        + (if enable_tracking_root then ['tracking_visitor', 'tagger_output'] else [])
        + (if tagger_bee then ['tagger_bee'] else []),

    // The clus_pr MABC as a pass-through tensor node (dump=false).  bee_sink: the
    // shared IBeeSink the PR display layers (clustering-pr / track_fit /
    // shower_track / vertices / mc) go to; it must be non-null for the production
    // Bee content (the "mc" node's merge of the upstream truth tree, the BDT scores
    // in its text) -- see clus.jsonnet pr().
    // dl_vtx_dump (ai-helper issue 35): record the DL-vertex network calls into
    // tracking-pr.root (T_dlvtx_call / T_dlvtx_cloud); recording only, default off.
    node(clus_maker, anodes, bee_sink, enable_tracking_root=true, tagger_bee=false, dl_vtx_dump=false)::
        clus_maker.pr(anodes, dump=false, bee_sink=bee_sink, dl_vtx_dump=dl_vtx_dump,
                      pipeline_names=$.pipeline_names(enable_tracking_root, tagger_bee),
                      particle_dataset=pds.particle_dataset, extra_uses=pds.all,
                      beam_window=$.beam_window),
}
