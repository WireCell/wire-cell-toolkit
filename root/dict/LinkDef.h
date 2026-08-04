// ROOT dictionary for WireCellRoot.
//
// WHY THIS FILE EXISTS -- please read before deleting it.
//
// These three nested-STL pragmas are what let TTree::Branch() create the
// T_proj_data branches written by root/'s three Magnify-tracking visitors:
//
//   SbndPrMagnifyTrackingVisitor.cxx    (SBND neutrino-PR chain, tracking-pr.root)
//   SbndMagnifyTrackingVisitor.cxx      (SBND STM chain,          tracking-stm.root)
//   UbooneMagnifyTrackingVisitor.cxx    (uBooNE,                  tracking-*.root)
//
// Each declares std::vector<std::vector<int>> members (channel, time_slice,
// charge, charge_err, charge_pred) and branches them.  ROOT refuses to create a
// branch for an STL collection that has no *compiled* CollectionProxy: it
// reports on the error stream, returns a null branch WITHOUT throwing, and the
// tree silently loses five of its six branches.  A TClass that is merely
// interpreted (GetState() == kInterpreted) still returns a non-null
// GetCollectionProxy(), so that is not a usable check -- only a dictionary
// generated here makes the class kHasTClassInit.
//
// HISTORY / MERGE HAZARD.  This file previously also carried the SBND reco1
// art-file mirror classes, and was deleted wholesale by upstream 5f684887
// ("moved to standalone wire-cell-sbnd-reco1", 2026-07-28) when that code left
// the toolkit.  The three std:: pragmas below were collateral damage -- they
// belong to root/, not to reco1 -- and every tracking-*.root written between
// that merge and the commit restoring this file has an empty T_proj_data.  The
// failure is silent in both directions, so if a future merge deletes this file
// again, nothing will fail loudly.  See sbnd_xin/docs/pr/26 sec 5.1.
//
// Keep this dictionary limited to std:: collections that root/ itself branches.

#ifdef __ROOTCLING__
#pragma link C++ class vector < vector<int> > +;
#pragma link C++ class vector < vector<float> > +;
#pragma link C++ class vector < vector<double> > +;
#endif
