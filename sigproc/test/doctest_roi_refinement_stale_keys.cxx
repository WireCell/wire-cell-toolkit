// ROI_refinement::CleanUpInductionROIs must not leave a deleted loose ROI
// behind as a key of the loose->tight containment map (contained_rois).
//
// Such a stale key outlives its ROI.  Refinement of the U plane finishes
// before V on the same ROI_refinement object, so a fresh V SignalROI that
// malloc places at the freed U address inherits U's tight ROIs -- a
// heap-dependent decision that showed up as run-to-run V-plane-only
// differences in SBND DetSim gauss/wiener (wcp-porting-img sbnd_xin/docs/127).
//
// ROI_refinement itself defaults to the legacy behaviour and
// set_erase_stale_contained(true) enables the fix; both are pinned here.
// OmnibusSigProc (its only user) turns the fix on by default via its
// r_erase_stale_contained knob, pinned by the last test case.  The stale key is only ever
// compared, never dereferenced, so the check below is safe on freed memory.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Array.h"
#include "WireCellUtil/Waveform.h"

#include "WireCellSigProc/OmnibusSigProc.h"

#include "../src/ROI_formation.h"
#include "../src/ROI_refinement.h"

#include <algorithm>
#include <cmath>

using namespace WireCell;
using namespace WireCell::SigProc;

namespace {

    const int nwire = 4;   // wires per plane
    const int nticks = 200;

    // Keys of contained_rois that are not a live loose ROI of plane 0 (U).
    size_t stale_u_keys(ROI_refinement& rr)
    {
        size_t nstale = 0;
        auto& loose = rr.get_u_rois();
        for (const auto& kv : rr.get_contained_rois()) {
            bool live = false;
            for (const auto& lst : loose) {
                if (std::find(lst.begin(), lst.end(), kv.first) != lst.end()) {
                    live = true;
                    break;
                }
            }
            if (!live) ++nstale;
        }
        return nstale;
    }

    // One isolated, sub-threshold U loose ROI [100,130] on wire 1 that
    // contains one tight ROI [110,120], on a triangular pulse peaking at 10
    // at tick 115 (SignalROI subtracts the straight line between its end
    // bins, so the pulse must fall toward the window edges).  The tight ROI
    // passes load_data's rms*th_factor cut (~8.3 > 1*3) but the loose ROI is
    // far below the induction cleanup thresholds (500/1000), so
    // CleanUpInductionROIs(0) deletes it.
    void load_one_bad_u_roi(ROI_refinement& rr, ROI_formation& rf)
    {
        rf.get_uplane_rms().assign(nwire, 1.0);
        rf.get_self_rois(1).push_back({110, 120});
        rf.get_loose_rois(1).push_back({100, 130});

        Array::array_xxf r_data = Array::array_xxf::Zero(nwire, nticks);
        for (int t = 109; t <= 121; ++t) r_data(1, t) = 10.0 * (1.0 - std::abs(t - 115) / 6.0);
        rr.load_data(0, r_data, rf);
    }
}  // namespace

TEST_CASE("ROI_refinement CleanUpInductionROIs contained_rois stale keys")
{
    Waveform::ChannelMaskMap cmm;

    SUBCASE("setup: the loose ROI is keyed before cleanup")
    {
        ROI_formation rf(cmm, nwire, nwire, nwire, nticks);
        ROI_refinement rr(cmm, nwire, nwire, nwire);
        load_one_bad_u_roi(rr, rf);
        REQUIRE(rr.get_u_rois().at(1).size() == 1);
        REQUIRE(rr.get_contained_rois().size() == 1);
        CHECK(stale_u_keys(rr) == 0);
    }

    SUBCASE("legacy (default): the deleted ROI stays keyed")
    {
        ROI_formation rf(cmm, nwire, nwire, nwire, nticks);
        ROI_refinement rr(cmm, nwire, nwire, nwire);
        load_one_bad_u_roi(rr, rf);
        rr.CleanUpInductionROIs(0);
        REQUIRE(rr.get_u_rois().at(1).empty());  // the ROI was deleted
        CHECK(stale_u_keys(rr) == 1);             // ...but its key survives
    }

    SUBCASE("erase_stale_contained: no key outlives its ROI")
    {
        ROI_formation rf(cmm, nwire, nwire, nwire, nticks);
        ROI_refinement rr(cmm, nwire, nwire, nwire);
        rr.set_erase_stale_contained(true);
        load_one_bad_u_roi(rr, rf);
        rr.CleanUpInductionROIs(0);
        REQUIRE(rr.get_u_rois().at(1).empty());
        CHECK(stale_u_keys(rr) == 0);
        CHECK(rr.get_contained_rois().empty());
    }
}

TEST_CASE("OmnibusSigProc r_erase_stale_contained defaults to true")
{
    WireCell::SigProc::OmnibusSigProc osp;
    auto cfg = osp.default_configuration();
    REQUIRE(cfg.isMember("r_erase_stale_contained"));
    CHECK(cfg["r_erase_stale_contained"].asBool() == true);
}
