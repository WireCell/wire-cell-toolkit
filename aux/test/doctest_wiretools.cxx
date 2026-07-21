#include "WireCellUtil/doctest.h"
#include "WireCellAux/Testing.h"

#include "WireCellAux/WireTools.h"


using namespace WireCell;
using namespace WireCell::Aux;

// Want to avoid making it more than once as it takes some time relative to what else is being tested.
static IAnodePlane::pointer get_anode()
{
    static IAnodePlane::pointer anode0 = WireCell::Testing::anodes("pdsp")[0];
    return anode0;
}


DOCTEST_TEST_SUITE("WireTools") {

    DOCTEST_TEST_CASE("wan_ordered_channels and nwires_wpid") {
        spdlog::info("Testing WireTools::wan_ordered_channels and nwires_wpid");

        auto anode = get_anode();
        
        const int nplanes = 3;
        for (int iplane_index = 0; iplane_index < nplanes; ++iplane_index) {
        
            std::vector<int> wpid_nums;
            size_t nwires_expected = 0;
            size_t nchans_expected = 0;
            std::vector<size_t> nwires_per_face_expected;
            for (auto iface : anode->faces()) {
                auto iplane = iface->planes()[iplane_index];
                DOCTEST_REQUIRE(iplane);
                const auto& iwires = iplane->wires();

                auto wpid = iwires[0]->planeid();
                int wpid_num = wpid.ident();
                DOCTEST_REQUIRE(wpid_num > 0);
                wpid_nums.push_back(wpid_num);
                nwires_expected += iwires.size();
                nwires_per_face_expected.push_back(iwires.size());

                nchans_expected += iplane->channels().size();
            }                

            DOCTEST_REQUIRE( wpid_nums.size() > 0 );
            DOCTEST_REQUIRE( nwires_expected > 0 );
            DOCTEST_REQUIRE( nchans_expected > 0 );

            auto chans_by_wan = WireTools::wan_ordered_channels(anode, wpid_nums);

            DOCTEST_REQUIRE( nchans_expected == chans_by_wan.size() );
            DOCTEST_REQUIRE( nwires_per_face_expected == WireTools::nwires_wpid(anode, wpid_nums));

        }
    }

    DOCTEST_TEST_CASE("wire_channels") {
        spdlog::info("Testing WireTools::wire_channels");
        auto anode = get_anode();

        const int nplanes = 3;
        std::vector<int> max_segments = {3,3,1};
        for (int iplane_index = 0; iplane_index < nplanes; ++iplane_index) {
        
            IChannel::vector all_chans;
            IWire::vector all_wires;

            for (auto iface : anode->faces()) {
                auto iplane = iface->planes()[iplane_index];
                DOCTEST_REQUIRE(iplane);

                const auto& iwires = iplane->wires();
                all_wires.insert(all_wires.end(), iwires.begin(), iwires.end());

                const auto& ichans = iplane->channels();
                all_chans.insert(all_chans.end(), ichans.begin(), ichans.end());
            }

            
            IChannel::vector wc_all = WireTools::wire_channels(all_wires, all_chans);
        
            DOCTEST_CHECK(wc_all.size() == all_wires.size());
        
            for (size_t i=0; i<all_wires.size(); ++i) {
                // Check wire_channels result
                DOCTEST_CHECK(wc_all[i] != nullptr);
                DOCTEST_CHECK(wc_all[i]->ident() == all_wires[i]->channel());
            
            }

            // 2. Test with a subset of channels (e.g., skip the first one)
            IChannel::vector subset_chans;
            if (all_chans.size() > 1) {
                subset_chans.insert(subset_chans.end(), all_chans.begin() + 1, all_chans.end());
            }
        
            IChannel::vector wc_subset = WireTools::wire_channels(all_wires, subset_chans);

            DOCTEST_CHECK(wc_subset.size() == all_wires.size());

            // Multiple wires may be on the omitted channel
            int nwires_no_channel=0;
            for (size_t i=0; i<all_wires.size(); ++i) {
                auto ich = wc_subset[i];

                if (ich == nullptr) {
                    ++nwires_no_channel;
                    continue;
                }


                DOCTEST_REQUIRE(ich != nullptr);
                DOCTEST_CHECK(ich->ident() == all_wires[i]->channel());
            
            }

            // We got rid of one channel so we must have at least one wire that
            // lacks a channel but no more than the max number of wire segments
            // a channel may have in the plane.
            DOCTEST_CHECK(nwires_no_channel > 0);
            DOCTEST_CHECK(nwires_no_channel <= max_segments[iplane_index]);
        }
    }
}
