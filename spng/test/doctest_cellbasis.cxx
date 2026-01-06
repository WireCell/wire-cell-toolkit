#include "WireCellSpng/Testing.h"
#include "WireCellSpng/CellBasis.h"
#include "WireCellSpng/RayGridOG.h"
#include "WireCellSpng/Util.h"

#include "WireCellAux/Testing.h"


#include "WireCellIface/IWirePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IAnodePlane.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/RayHelpers.h"

#include <cmath>
#include <memory>

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


    // // Expected calculations (15 cells total)
    
    // // U indices
    // auto expected_u = torch::tensor({0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2}, torch::kLong);
    
    // // V indices (U0: [-2, 2], U1: [-2, 2], U2: [-3, 1])
    // auto expected_v = torch::tensor({-2, -1, 0, 1, 2, -2, -1, 0, 1, 2, -3, -2, -1, 0, 1}, torch::kLong);
    
    // // W indices (i_w = (i_u - i_v) / 1.732, rounded, clamped [0, 4])
    // // Raw rounded W indices: {1, 1, 0, -1, -1, 2, 1, 1, 0, -1, 3, 2, 2, 1, 1}
    // // Clamped [0, 4]
    // auto expected_w = torch::tensor({1, 1, 0, 0, 0, 2, 1, 1, 0, 0, 3, 2, 2, 1, 1}, torch::kLong);
    
    // torch::Tensor expected_basis = torch::stack({expected_u, expected_v, expected_w}, 1);

    // DOCTEST_CHECK(basis.sizes().vec() == std::vector<int64_t>{15, 3});
    // DOCTEST_CHECK(basis.dtype() == torch::kLong);
    
    // DOCTEST_CHECK(basis.equal(expected_basis));
}

