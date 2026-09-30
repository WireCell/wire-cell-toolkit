// OBSOLETE since 2026-09-30 (ai-helper issue 33): production is the 2-step chain,
// wcls-img-clus-matching.jsonnet (step 1, LArSoft) + wct-pr.jsonnet (step 2, standalone PR),
// validated identical to this 1-step job.  Kept runnable -- the whole chain in one lar
// process -- as pgrapher/experiment/sbnd/obsolete/wcls-img-clus-matching-pr-hits.jsonnet (renamed from wcls-img-clus-matching-xin-hits.jsonnet).
//
// SBND LArSoft 1-step job -- imaging + clustering + Q/L + the PR tail -- with the
// SBND standalone PRODUCTION light of 2026-09-25 (doc sbnd_xin/123 sec 17 + 19):
// flashes rebuilt from the reco1 PMT OpHits (wclsOpHitSource:tpc<N> ->
// SBNDOpFlashFinder:tpc<N>, the art-side twin of SBNDReco1OpHitSource ->
// SBNDOpFlashFinder in sbnd_xin/wct-reco1-dump.jsonnet) and the QLMatching
// scenario-1 light gate + over-prediction ceiling (wct-clus-matching-perevt.jsonnet
// xtpc_sc1_light_gate=true, xtpc_sc1_overpred_max=2.9).  Everything else is the
// reco1-flash job, wcls-img-clus-matching-pr-flash.jsonnet.
//
// Requires in the fcl: inputers wclsOpHitSource:tpc0/tpc1 (instead of
// wclsOpFlashSource:tpc0/tpc1), params ophit0_input_label / ophit1_input_label
// (the recob::OpHit product, "ophitpmt" in reco1), and "WireCellFlash" in plugins.
// See wcp-porting-validation sbnd/wcls-img-clus-matching-pr-hits.fcl.

(import 'pgrapher/experiment/sbnd/wcls-img-clus-matching-pr-lib.jsonnet')(
    flash_source='hits',
    xtpc_sc1_light_gate=true,
    xtpc_sc1_overpred_max=2.9)
