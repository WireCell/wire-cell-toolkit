// SBND reco1 art-file mirror classes (bare-ROOT read of recob::Wire /
// recob::OpFlash / timing products; see WireCellRoot/SBNDReco1Products.h).
// Included OUTSIDE the __ROOTCLING__ guard so rootcling sees the class
// declarations when generating the dictionary.
#include "../inc/WireCellRoot/SBNDReco1Products.h"

#ifdef __ROOTCLING__
#pragma link C++ class vector < vector<int> > +;
#pragma link C++ class vector < vector<float> > +;
#pragma link C++ class vector < vector<double> > +;

// SBND reco1 mirror classes.  Pragmas mirror the TFile::MakeProject
// LinkDef for the same file (only the subset the reco1 sources read).
#pragma link C++ class lar::range_t<unsigned long>+;
#pragma link C++ class lar::sparse_vector<float>+;
#pragma link C++ class lar::sparse_vector<float>::datarange_t+;
#pragma link C++ class vector<lar::sparse_vector<float>::datarange_t>+;
#pragma link C++ class recob::Wire+;
#pragma link C++ class vector<recob::Wire>+;
#pragma link C++ class recob::OpFlash+;
#pragma link C++ class vector<recob::OpFlash>+;
#pragma link C++ class recob::OpHit+;
#pragma link C++ class vector<recob::OpHit>+;
#pragma link C++ class raw::ptb::Trigger+;
#pragma link C++ class raw::ptb::ChStatus+;
#pragma link C++ class raw::ptb::Feedback+;
#pragma link C++ class raw::ptb::Misc+;
#pragma link C++ class raw::ptb::WordIndex+;
#pragma link C++ class raw::ptb::sbndptb+;
#pragma link C++ class vector<raw::ptb::sbndptb>+;
#pragma link C++ class sbnd::timing::DAQTimestamp+;
#pragma link C++ class vector<sbnd::timing::DAQTimestamp>+;
#pragma link C++ class artdaq::detail::RawEventHeader+;
#pragma link C++ class art::Hash<2>+;
#pragma link C++ class art::RunID+;
#pragma link C++ class art::SubRunID+;
#pragma link C++ class art::EventID+;
#pragma link C++ class art::Timestamp+;
#pragma link C++ class art::EventAuxiliary+;
#pragma link C++ class art::RunAuxiliary+;
#pragma link C++ class art::EDProduct+;
#pragma link C++ class art::Wrapper<vector<recob::Wire> >+;
#pragma link C++ class art::Wrapper<vector<recob::OpFlash> >+;
#pragma link C++ class art::Wrapper<vector<recob::OpHit> >+;
#pragma link C++ class art::Wrapper<vector<int> >+;
#pragma link C++ class art::Wrapper<vector<double> >+;
#pragma link C++ class art::Wrapper<vector<raw::ptb::sbndptb> >+;
#pragma link C++ class art::Wrapper<vector<sbnd::timing::DAQTimestamp> >+;
#pragma link C++ class art::Wrapper<artdaq::detail::RawEventHeader>+;
#endif
