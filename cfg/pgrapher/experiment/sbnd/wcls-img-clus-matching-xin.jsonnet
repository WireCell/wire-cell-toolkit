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

//
// The SBND LArSoft 1-step job with the reco1 recob::OpFlash light (wclsOpFlashSource).
// The implementation is wcls-img-clus-matching-xin-lib.jsonnet (a function); this
// file compiles byte-identically to the pre-2026-09-26 monolithic job.  The
// hit-flash variant is wcls-img-clus-matching-xin-hits.jsonnet.

(import 'pgrapher/experiment/sbnd/wcls-img-clus-matching-xin-lib.jsonnet')()
