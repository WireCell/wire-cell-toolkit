/** Shared bare-ROOT reading helpers for SBND reco1 art files.
 *
 * Gotchas learned probing the files (do not "simplify" these away):
 *  - NEVER call tree->GetEntry(i): the Events tree holds ~2500 branches
 *    including DAQ fragment products with no dictionary; deserializing
 *    them crashes.  Always read per-branch: branch->GetEntry(i).
 *  - Read art products through the TOP-LEVEL wrapper branch
 *    ("<product>.") with an art::Wrapper<...> object address.  Setting a
 *    typed vector address on the ".obj" member sub-branch silently reads
 *    nothing.
 *
 * Only the SBNDReco1*Source components use this header.
 */
#ifndef WIRECELLROOT_SBNDRECO1READER
#define WIRECELLROOT_SBNDRECO1READER

#include "WireCellRoot/SBNDReco1Products.h"
#include "WireCellUtil/Exceptions.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"

#include <memory>
#include <string>

namespace WireCell {
    namespace Root {
        namespace SBNDReco1 {

            // Read one art product (wrapper) from the Events tree at the
            // given entry.  Returns false if the branch does not exist.
            template <typename WrapperType>
            bool read_wrapper(TTree* events, const std::string& branch_name, Long64_t entry,
                              WrapperType& wrapper)
            {
                TBranch* br = events->GetBranch(branch_name.c_str());
                if (!br) return false;
                WrapperType* ptr = &wrapper;
                br->SetAddress(&ptr);
                br->GetEntry(entry);
                br->ResetAddress();
                return true;
            }

            inline bool read_event_aux(TTree* events, Long64_t entry, art::EventAuxiliary& aux)
            {
                TBranch* br = events->GetBranch("EventAuxiliary");
                if (!br) return false;
                art::EventAuxiliary* ptr = &aux;
                br->SetAddress(&ptr);
                br->GetEntry(entry);
                br->ResetAddress();
                return true;
            }

            // Locate the tree entry with the given run/subrun/event
            // (subrun < 0 matches any).  Returns -1 when not found.
            inline Long64_t find_entry(TTree* events, int run, int subrun, int event)
            {
                art::EventAuxiliary aux;
                const Long64_t nent = events->GetEntries();
                for (Long64_t ent = 0; ent < nent; ++ent) {
                    if (!read_event_aux(events, ent, aux)) return -1;
                    if ((int) aux.id_.subRun_.run_.run_ != run) continue;
                    if (subrun >= 0 && (int) aux.id_.subRun_.subRun_ != subrun) continue;
                    if ((int) aux.id_.event_ != event) continue;
                    return ent;
                }
                return -1;
            }

            // Begin time of the given run from the Runs tree, or an
            // invalid (0/0) timestamp when unavailable.
            inline art::Timestamp run_begin_time(TFile* tfile, unsigned int run)
            {
                art::Timestamp invalid;
                TTree* runs = dynamic_cast<TTree*>(tfile->Get("Runs"));
                if (!runs) return invalid;
                TBranch* br = runs->GetBranch("RunAuxiliary");
                if (!br) return invalid;
                art::RunAuxiliary raux;
                art::RunAuxiliary* rptr = &raux;
                br->SetAddress(&rptr);
                art::Timestamp found;
                for (Long64_t ent = 0; ent < runs->GetEntries(); ++ent) {
                    br->GetEntry(ent);
                    if (raux.id_.run_ == run) {
                        found = raux.beginTime_;
                        break;
                    }
                }
                br->ResetAddress();
                return found;
            }

            // The raw DAQ header timestamp.  The file carries one
            // RawEventHeader branch per event builder (DAQEVB01..07,
            // *P2); exactly one is present per event.  Scan branches
            // matching the label prefix and return the present one.
            // Returns 0 when none found.
            inline uint64_t raw_header_timestamp(TTree* events, Long64_t entry,
                                                 const std::string& prefix =
                                                     "artdaq::detail::RawEventHeader_daq_RawEventHeader_")
            {
                TObjArray* branches = events->GetListOfBranches();
                for (int i = 0; i < branches->GetEntries(); ++i) {
                    TBranch* br = (TBranch*) branches->At(i);
                    const std::string name = br->GetName();
                    if (name.compare(0, prefix.size(), prefix) != 0) continue;
                    art::Wrapper<artdaq::detail::RawEventHeader> wrapper;
                    art::Wrapper<artdaq::detail::RawEventHeader>* wptr = &wrapper;
                    br->SetAddress(&wptr);
                    br->GetEntry(entry);
                    br->ResetAddress();
                    if (wrapper.present) {
                        return wrapper.obj.timestamp;
                    }
                }
                return 0;
            }

            struct FilePtr {
                TFile* file{nullptr};
                TTree* events{nullptr};
                ~FilePtr()
                {
                    if (file) {
                        file->Close();
                        delete file;
                    }
                }
            };

            inline void open_events(FilePtr& fp, const std::string& filename, const std::string& who)
            {
                fp.file = TFile::Open(filename.c_str(), "READ");
                if (!fp.file || fp.file->IsZombie()) {
                    THROW(IOError() << errmsg{who + ": failed to open " + filename});
                }
                fp.events = dynamic_cast<TTree*>(fp.file->Get("Events"));
                if (!fp.events) {
                    THROW(IOError() << errmsg{who + ": no Events tree in " + filename});
                }
            }

        }  // namespace SBNDReco1
    }  // namespace Root
}  // namespace WireCell

#endif
