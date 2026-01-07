#include "WireCellSpng/Testing.h"
#include "WireCellSpng/CellBasis.h"
#include "WireCellSpng/RayGridOG.h"
#include "WireCellSpng/Util.h"

#include "WireCellAux/Testing.h"


#include "WireCellIface/IWirePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IWire.h"
#include "WireCellIface/IChannel.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/RayHelpers.h"

#include <cmath>
#include <memory>
#include <algorithm> // For std::reverse check

using namespace WireCell;
using namespace WireCell::SPNG;
using namespace torch::indexing;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Return the maximum distance between columns given tensors of shape (N, 2).
torch::Tensor rgc_dist_stat(torch::Tensor a, torch::Tensor b)
{
    auto diff = a - b;
    auto dist = torch::norm(diff, 2, 1);
    return torch::stack({dist.min(),dist.mean(),dist.max()});
}

// Want to avoid making it more than once as it takes some time relative to what else is being tested.
static IAnodePlane::pointer get_anode()
{
    static IAnodePlane::pointer anode0 = WireCell::Testing::anodes("pdsp")[0];
    return anode0;
}

DOCTEST_TEST_SUITE("CellBasis") {


    DOCTEST_TEST_CASE("wire_endpoints") {
        spdlog::info("Testing CellBasis::wire_endpoints");

        auto anode = get_anode();

        spdlog::info("got anode");

        IWire::vector wires = anode->face(0)->plane(0)->wires();

        torch::Tensor endpoints = CellBasis::wire_endpoints(wires);

        spdlog::info("found endpoints {}", to_string(endpoints));

        int64_t nwires = wires.size();
        DOCTEST_REQUIRE(nwires == endpoints.size(0));

        for (int64_t index = 0; index < nwires; ++index) {
            auto ray = wires[index]->ray();

            DOCTEST_REQUIRE(ray.first.z() == endpoints[index][0][0].item<double>());
            DOCTEST_REQUIRE(ray.first.y() == endpoints[index][0][1].item<double>());
            DOCTEST_REQUIRE(ray.second.z() == endpoints[index][1][0].item<double>());
            DOCTEST_REQUIRE(ray.second.y() == endpoints[index][1][1].item<double>());
        }           

    }
    DOCTEST_TEST_CASE("cell_basis") {
        spdlog::info("Testing CellBasis::cell_basis with standard geometry");

        auto anode = get_anode();

        spdlog::info("got anode");

        auto iface = anode->face(0);
        spdlog::info("Making cell basis:");
        torch::Tensor basis = CellBasis::cell_basis(iface);
        spdlog::info("cell basis: {}", to_string(basis));

        const auto& iwires_u = iface->plane(0)->wires();
        const auto& iwires_v = iface->plane(1)->wires();
        const auto& iwires_w = iface->plane(2)->wires();

        size_t ncells = basis.size(0);
        DOCTEST_REQUIRE(ncells > 0);

        spdlog::info("testing lower bounds");

        // All must be non-negative
        DOCTEST_REQUIRE( ( basis >= 0 ).all().item<bool>() );

        spdlog::info("testing per-plane upper bounds");

        auto byplane = basis.t();

        DOCTEST_REQUIRE( ( byplane[0] < iwires_u.size() ).all().item<bool>() );
        DOCTEST_REQUIRE( ( byplane[1] < iwires_v.size() ).all().item<bool>() );
        DOCTEST_REQUIRE( ( byplane[2] < iwires_w.size() ).all().item<bool>() );

        spdlog::info("making spng raygrid from og:");
        auto rg_views = to_spng_views(iface->raygrid());
        SPNG::RayGrid::Coordinates coords(rg_views);
        
        const auto u_layer = torch::tensor({2+0}, torch::kLong);
        const auto v_layer = torch::tensor({2+1}, torch::kLong);
        const auto w_layer = torch::tensor({2+2}, torch::kLong);


        spdlog::info("re finding ray crossings");

        // Assure pair-wise crossings are mutually close
        auto uv_pts = coords.ray_crossing(u_layer, byplane[0],
                                          v_layer, byplane[1]);
        auto vw_pts = coords.ray_crossing(v_layer, byplane[1],
                                          w_layer, byplane[2]);
        auto wu_pts = coords.ray_crossing(w_layer, byplane[2],
                                          u_layer, byplane[0]);

        spdlog::info("uv_pts: {}", to_string(uv_pts));

        const auto uvw_stats = rgc_dist_stat(uv_pts, vw_pts);
        const auto vwu_stats = rgc_dist_stat(vw_pts, wu_pts);
        const auto wuv_stats = rgc_dist_stat(wu_pts, uv_pts);

        spdlog::info("uvw_stats: {}", to_string(uvw_stats));

        std::vector<std::string> stats = {"min", "mean", "max"};
        /// These bounds are merely TOFU.  The max bounds seem a bit large to me for PDSP.
        std::vector<double> lo_stat = {0.4*units::mm, 2*units::mm, 10*units::mm};
        std::vector<double> hi_stat = {0.9*units::mm, 4*units::mm, 20*units::mm};
        for (int ind = 0; ind<3; ++ind) {
            const double uvw = uvw_stats[ind].item<double>();
            const double vwu = vwu_stats[ind].item<double>(); 
            const double wuv = wuv_stats[ind].item<double>();

            spdlog::info("stat {}: |uv-vw|={} |vw-wu|={} |wu-uv|={}",
                         stats[ind], uvw, vwu, wuv);
            
            const double lo = lo_stat[ind];
            const double hi = hi_stat[ind];
            DOCTEST_CHECK(uvw < hi);
            DOCTEST_CHECK(vwu < hi);
            DOCTEST_CHECK(wuv < hi);
            DOCTEST_CHECK(uvw > lo);
            DOCTEST_CHECK(vwu > lo);
            DOCTEST_CHECK(wuv > lo);
        }
    }

    DOCTEST_TEST_CASE("wan_ordered_channels and nwires_wpid") {
        spdlog::info("Testing CellBasis::wan_ordered_channels and nwires_wpid");

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

            auto chans_by_wan = CellBasis::wan_ordered_channels(anode, wpid_nums);

            DOCTEST_REQUIRE( nchans_expected == chans_by_wan.size() );
            DOCTEST_REQUIRE( nwires_per_face_expected == CellBasis::nwires_wpid(anode, wpid_nums));

        }
    }

    DOCTEST_TEST_CASE("channel_idents") {
        spdlog::info("Testing CellBasis::channel_idents");
        auto anode = get_anode();

        auto face0 = anode->face(0);
        auto plane0 = face0->plane(0);
        
        IChannel::vector chans;
        for (const auto& ich : plane0->channels()) {
            chans.push_back(ich);
        }
        
        torch::Tensor idents_ten = CellBasis::channel_idents(chans);
        
        DOCTEST_CHECK(idents_ten.sizes().vec() == std::vector<int64_t>{(int64_t)chans.size()});
        DOCTEST_CHECK(idents_ten.dtype() == torch::kInt);
        
        for (size_t i=0; i<chans.size(); ++i) {
            DOCTEST_CHECK(idents_ten[i].item<int>() == chans[i]->ident());
        }
    }

    DOCTEST_TEST_CASE("wire_channels and wire_channel_index") {
        spdlog::info("Testing CellBasis::wire_channels and wire_channel_index");
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

            
            IChannel::vector wc_all = CellBasis::wire_channels(all_wires, all_chans);
            torch::Tensor wci_all = CellBasis::wire_channel_index(all_wires, all_chans);
        
            DOCTEST_CHECK(wc_all.size() == all_wires.size());
            DOCTEST_CHECK(wci_all.size(0) == (int64_t)all_wires.size());
            DOCTEST_CHECK(wci_all.dtype() == torch::kLong);
        
            for (size_t i=0; i<all_wires.size(); ++i) {
                // Check wire_channels result
                DOCTEST_CHECK(wc_all[i] != nullptr);
                DOCTEST_CHECK(wc_all[i]->ident() == all_wires[i]->channel());
            
                // Check wire_channel_index result
                int64_t index = wci_all[i].item<int64_t>();
                DOCTEST_CHECK(index >= 0);
                DOCTEST_CHECK(all_chans[index]->ident() == all_wires[i]->channel());
            }

            // 2. Test with a subset of channels (e.g., skip the first one)
            IChannel::vector subset_chans;
            if (all_chans.size() > 1) {
                subset_chans.insert(subset_chans.end(), all_chans.begin() + 1, all_chans.end());
            }
        
            IChannel::vector wc_subset = CellBasis::wire_channels(all_wires, subset_chans);
            torch::Tensor wci_subset = CellBasis::wire_channel_index(all_wires, subset_chans);

            DOCTEST_CHECK(wc_subset.size() == all_wires.size());
            DOCTEST_CHECK(wci_subset.size(0) == (int64_t)all_wires.size());

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
            
                int64_t index = wci_subset[i].item<int64_t>();
                DOCTEST_CHECK(index >= 0);
                DOCTEST_CHECK(subset_chans[index]->ident() == all_wires[i]->channel());
            
            }

            // We got rid of one channel so we must have at least one wire that
            // lacks a channel but no more than the max number of wire segments
            // a channel may have in the plane.
            DOCTEST_CHECK(nwires_no_channel > 0);
            DOCTEST_CHECK(nwires_no_channel <= max_segments[iplane_index]);
        }
    }
}

