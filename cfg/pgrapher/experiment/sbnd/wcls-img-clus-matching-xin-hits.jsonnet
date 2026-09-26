// SBND LArSoft 1-step job -- imaging + clustering + Q/L + the PR tail -- with the
// SBND standalone PRODUCTION light of 2026-09-25 (doc sbnd_xin/123 sec 17 + 19):
// flashes rebuilt from the reco1 PMT OpHits (wclsOpHitSource:tpc<N> ->
// SBNDOpFlashFinder:tpc<N>, the art-side twin of SBNDReco1OpHitSource ->
// SBNDOpFlashFinder in sbnd_xin/wct-reco1-dump.jsonnet) and the QLMatching
// scenario-1 light gate + over-prediction ceiling (wct-clus-matching-perevt.jsonnet
// xtpc_sc1_light_gate=true, xtpc_sc1_overpred_max=2.9).  Everything else is the
// reco1-flash job, wcls-img-clus-matching-xin.jsonnet.
//
// Requires in the fcl: inputers wclsOpHitSource:tpc0/tpc1 (instead of
// wclsOpFlashSource:tpc0/tpc1), params ophit0_input_label / ophit1_input_label
// (the recob::OpHit product, "ophitpmt" in reco1), and "WireCellFlash" in plugins.
// See wcp-porting-validation sbnd/wcls-img-clus-matching-xin-hits.fcl.

(import 'pgrapher/experiment/sbnd/wcls-img-clus-matching-xin-lib.jsonnet')(
    flash_source='hits',
    xtpc_sc1_light_gate=true,
    xtpc_sc1_overpred_max=2.9)
