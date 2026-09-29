// STEP 1 of the SBND 2-step chain (ai-helper issue 33), LArSoft (wcls):
//   recob::Wire + recob::OpHit -> imaging -> per-APA clustering -> Q/L matching
//   -> all-APA clustering (clus_all_apa) -> labeler_truth
//   -> wclsTruthInformationAttacher:truth -> qlpctree.tar.gz
// The tar holds, per event, the matched bundles (live + dead groupings with the
// CTPC / dead-wind / flash / light / provenance point clouds), the set metadata
// (RSE, the labeler's nu_* keys) and, on MC, the truth_nu / truth_pf tables.
// Step 2 re-runs pattern recognition on it standalone:
//   wire-cell -c pgrapher/experiment/sbnd/wct-pr.jsonnet --tla-str input=qlpctree.tar.gz ...
// Light: the SBND standalone production light (hit-rebuilt flashes + the scenario-1
// light gate), as wcls-img-clus-matching-xin-hits.jsonnet.  Everything up to and
// including the truth attacher is the 1-step chain's own graph (the shared lib).
(import 'pgrapher/experiment/sbnd/wcls-img-clus-matching-xin-lib.jsonnet')(
    flash_source='hits', xtpc_sc1_light_gate=true, xtpc_sc1_overpred_max=2.9, stage='ql')
