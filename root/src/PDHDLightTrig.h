/** Private helper shared by the PDHD light-data ROOT sources.
 *
 * Reads the trigoff/trigger_offset tree of the temporary PDHD light
 * ROOT files and selects the trigger candidate for one event.  See
 * flash/docs/design.md for the time conventions.
 */
#ifndef WIRECELLROOT_PDHDLIGHTTRIG
#define WIRECELLROOT_PDHDLIGHTTRIG

#include "WireCellUtil/Exceptions.h"

#include "TFile.h"
#include "TTree.h"

#include <cmath>
#include <cstdint>

namespace WireCell::Root::PDHDLight {

    struct Trigger {
        int run{-1}, subrun{-1}, event{-1}, tc_type{-1};
        uint64_t rd_timestamp{0}, tc_time{0};
        double offset_us{0};
        int ncandidates{0};
    };

    // Select the trigger candidate of (run, event) whose offset_us is
    // closest to nominal_offset_us.  subrun < 0 matches any subrun.
    inline Trigger read_trigger(TFile* file, int run, int subrun, int event,
                                double nominal_offset_us)
    {
        TTree* tree = dynamic_cast<TTree*>(file->Get("trigoff/trigger_offset"));
        if (!tree) {
            THROW(IOError() << errmsg{"PDHDLight: no trigoff/trigger_offset tree in input file"});
        }
        Int_t b_run, b_subrun, b_event, b_tctype;
        ULong64_t b_rd, b_tc;
        Double_t b_offset;
        tree->SetBranchAddress("run", &b_run);
        tree->SetBranchAddress("subrun", &b_subrun);
        tree->SetBranchAddress("event", &b_event);
        tree->SetBranchAddress("tc_type", &b_tctype);
        tree->SetBranchAddress("rd_timestamp", &b_rd);
        tree->SetBranchAddress("tc_time_candidate", &b_tc);
        tree->SetBranchAddress("offset_us", &b_offset);

        Trigger best;
        double best_dist = -1;
        const Long64_t nent = tree->GetEntries();
        for (Long64_t ent = 0; ent < nent; ++ent) {
            tree->GetEntry(ent);
            if (b_run != run or b_event != event) continue;
            if (subrun >= 0 and b_subrun != subrun) continue;
            ++best.ncandidates;
            const double dist = std::abs(b_offset - nominal_offset_us);
            if (best_dist < 0 or dist < best_dist) {
                best_dist = dist;
                best.run = b_run;
                best.subrun = b_subrun;
                best.event = b_event;
                best.tc_type = b_tctype;
                best.rd_timestamp = b_rd;
                best.tc_time = b_tc;
                best.offset_us = b_offset;
            }
        }
        tree->ResetBranchAddresses();
        if (best.ncandidates == 0) {
            THROW(ValueError() << errmsg{"PDHDLight: no trigger candidate for requested run/event"});
        }
        return best;
    }

}  // namespace WireCell::Root::PDHDLight

#endif
