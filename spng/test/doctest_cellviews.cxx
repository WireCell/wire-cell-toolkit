/**
 * Unit test for CellViews MP2/MP3 computation.
 *
 * Strategy: build cell-basis channel mappings with the same functions
 * CellViews uses internally, then synthesize U/V/W input images where
 * individual cells are selectively activated, and verify the MP2 and MP3
 * output images satisfy the expected boolean properties.
 *
 * Three scenarios per test cell, each on a distinct tick to avoid cross-talk:
 *   tick+0: all three channels active  → MP3 fires, MP2 silent (active & ~active impossible)
 *   tick+1: only V and W active        → MP2 fires for U view, MP3 silent for U
 *   tick+2: only U and W active        → MP2 fires for V view, MP3 silent for V
 *
 * A global mutual-exclusivity check verifies that no output channel has both
 * MP2 and MP3 set at the same tick in any output view.
 */

#include "WireCellSpng/Testing.h"
#include "WireCellSpng/CellViews.h"
#include "WireCellSpng/CellBasis.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

#include "WireCellAux/Testing.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IWire.h"
#include "WireCellIface/IChannel.h"

using namespace WireCell;
using namespace WireCell::SPNG;

// Cache the anode to avoid redundant plugin loads across test cases.
static IAnodePlane::pointer get_anode()
{
    static IAnodePlane::pointer anode = WireCell::Testing::anodes("pdsp")[0];
    return anode;
}

DOCTEST_TEST_SUITE("CellViews") {

    DOCTEST_TEST_CASE("mp2 and mp3 basic consistency") {
        spdlog::info("Testing CellViews MP2/MP3 consistency");

        auto anode = get_anode();

        // Use both faces as CellViews does by default.
        const std::vector<int> face_idents = {0, 1};

        // Build channels_per_view in the same face order that CellViews::configure
        // uses.  The resulting channel index space must match the input tensors.
        std::vector<IChannel::vector> channels_per_view(3);
        for (int view = 0; view < 3; ++view) {
            for (int fid : face_idents) {
                auto iface = anode->face(fid);
                const auto& ichans = iface->planes()[view]->channels();
                channels_per_view[view].insert(channels_per_view[view].end(),
                                               ichans.begin(), ichans.end());
            }
        }

        std::vector<int64_t> nchan(3);
        for (int view = 0; view < 3; ++view) {
            nchan[view] = (int64_t)channels_per_view[view].size();
        }
        spdlog::info("channels per view: U={} V={} W={}",
                     nchan[0], nchan[1], nchan[2]);

        // Compute the channel-basis for face 0, replicating what CellViews does
        // internally.  We use this to know which channel index triples correspond
        // to real physical cells so we can set up meaningful input pixels.
        auto iface0 = anode->face(face_idents[0]);
        torch::Tensor wire_basis = CellBasis::cell_basis(iface0);

        std::vector<torch::Tensor> w2c;
        for (int view = 0; view < 3; ++view) {
            IWire::vector wires = iface0->planes()[view]->wires();
            w2c.push_back(CellBasis::wire_channel_index(wires, channels_per_view[view]));
        }
        torch::Tensor chan_basis = CellBasis::index(wire_basis, w2c);
        // chan_basis shape: (ncells, 3) – columns are [u_idx, v_idx, w_idx]

        int64_t ncells = chan_basis.size(0);
        DOCTEST_REQUIRE(ncells > 0);
        spdlog::info("face 0 cell basis has {} cells", ncells);

        // Sample at most 50 cells spread across the basis to keep the test fast.
        const int64_t max_test_cells = 50;
        const int64_t cell_step = std::max(int64_t(1), ncells / max_test_cells);

        // Three ticks per test cell (see module doc above).
        const int64_t ticks_per_cell = 3;
        const int64_t ntick = ((ncells + cell_step - 1) / cell_step) * ticks_per_cell;

        // Create boolean input tensors, all False.
        std::vector<torch::Tensor> uvw(3);
        for (int view = 0; view < 3; ++view) {
            uvw[view] = torch::zeros({nchan[view], ntick}, torch::kBool);
        }

        // Helper: set a pixel in a 2D bool tensor.
        auto set_pix = [](torch::Tensor& t, int64_t ch, int64_t tick) {
            t.accessor<bool, 2>()[ch][tick] = true;
        };

        // Record each test cell's channel indices and base tick offset.
        struct TestCell {
            int64_t u_ch, v_ch, w_ch, base_tick;
        };
        std::vector<TestCell> test_cells;

        int64_t base_tick = 0;
        for (int64_t ci = 0; ci < ncells; ci += cell_step) {
            int64_t u_ch = chan_basis[ci][0].item<int64_t>();
            int64_t v_ch = chan_basis[ci][1].item<int64_t>();
            int64_t w_ch = chan_basis[ci][2].item<int64_t>();

            // Skip any cell with an out-of-range channel index (sentinel -1).
            if (u_ch < 0 || u_ch >= nchan[0]) continue;
            if (v_ch < 0 || v_ch >= nchan[1]) continue;
            if (w_ch < 0 || w_ch >= nchan[2]) continue;

            // tick+0: all three channels active (exercise MP3).
            set_pix(uvw[0], u_ch, base_tick + 0);
            set_pix(uvw[1], v_ch, base_tick + 0);
            set_pix(uvw[2], w_ch, base_tick + 0);

            // tick+1: V and W active only (exercise MP2 for U view).
            set_pix(uvw[1], v_ch, base_tick + 1);
            set_pix(uvw[2], w_ch, base_tick + 1);

            // tick+2: U and W active only (exercise MP2 for V view).
            set_pix(uvw[0], u_ch, base_tick + 2);
            set_pix(uvw[2], w_ch, base_tick + 2);

            test_cells.push_back({u_ch, v_ch, w_ch, base_tick});
            base_tick += ticks_per_cell;
        }

        DOCTEST_REQUIRE(test_cells.size() > 0);
        spdlog::info("built {} test cells, {} total ticks", test_cells.size(), ntick);

        // ----------------------------------------------------------------
        // Instantiate and configure CellViews.
        // ----------------------------------------------------------------
        CellViews cv;
        auto cfg = cv.default_configuration();

        // Anode registered by Testing::anodes("pdsp").
        cfg["anode"] = "AnodePlane:0";

        // face_idents default is {0,1} which matches our channel layout.

        // Request output for all three views.
        {
            Json::Value jv(Json::arrayValue);
            jv.append(0); jv.append(1); jv.append(2);
            cfg["out_views"] = jv;
        }

        // cell_views default is ["mp2","mp3"] → mp2 at feature dim 0, mp3 at 1.
        cv.configure(cfg);

        // ----------------------------------------------------------------
        // Build input tensor set and run the filter.
        // ----------------------------------------------------------------
        ITorchTensor::vector tv;
        for (int view = 0; view < 3; ++view) {
            tv.push_back(std::make_shared<SimpleTorchTensor>(uvw[view]));
        }
        auto input_ts = std::make_shared<SimpleTorchTensorSet>(
            0, Json::Value(), tv);

        auto output_ts = cv.filter_tensor(input_ts);
        DOCTEST_REQUIRE(output_ts != nullptr);

        auto otens = output_ts->tensors();
        DOCTEST_REQUIRE((int64_t)otens->size() == 3);  // one per out_view

        // Retrieve output tensors.  Shape is (2, nchan[view], ntick):
        //   dim 0: feature index – 0=mp2, 1=mp3 (per default cell_views order)
        //   dim 1: channel index
        //   dim 2: tick index
        torch::Tensor out[3];
        for (int view = 0; view < 3; ++view) {
            out[view] = (*otens)[view]->tensor();
            DOCTEST_REQUIRE(out[view].dim() == 3);
            DOCTEST_REQUIRE(out[view].size(0) == 2);
            DOCTEST_REQUIRE(out[view].size(1) == nchan[view]);
            DOCTEST_REQUIRE(out[view].size(2) == ntick);
        }

        const int mp2 = 0;
        const int mp3 = 1;

        // Helper: read a single pixel from a (2, nchan, ntick) bool tensor.
        auto get = [](const torch::Tensor& t, int feat, int64_t ch, int64_t tick) {
            return t.accessor<bool, 3>()[feat][ch][tick];
        };

        // ----------------------------------------------------------------
        // Per-cell checks.
        // ----------------------------------------------------------------
        int nfail = 0;
        for (const auto& tc : test_cells) {
            // --- tick+0: all three channels active ---
            // MP3 must fire for U, V and W channels.
            if (!get(out[0], mp3, tc.u_ch, tc.base_tick + 0)) {
                spdlog::error("MP3 not set for U ch={} tick={}", tc.u_ch, tc.base_tick);
                ++nfail;
            }
            DOCTEST_CHECK(get(out[0], mp3, tc.u_ch, tc.base_tick + 0));

            if (!get(out[1], mp3, tc.v_ch, tc.base_tick + 0)) {
                spdlog::error("MP3 not set for V ch={} tick={}", tc.v_ch, tc.base_tick);
                ++nfail;
            }
            DOCTEST_CHECK(get(out[1], mp3, tc.v_ch, tc.base_tick + 0));

            if (!get(out[2], mp3, tc.w_ch, tc.base_tick + 0)) {
                spdlog::error("MP3 not set for W ch={} tick={}", tc.w_ch, tc.base_tick);
                ++nfail;
            }
            DOCTEST_CHECK(get(out[2], mp3, tc.w_ch, tc.base_tick + 0));

            // MP2 must NOT fire when the channel itself is active.
            // (~active) & ... = false & ... = false, so MP2 cannot be set.
            DOCTEST_CHECK(!get(out[0], mp2, tc.u_ch, tc.base_tick + 0));
            DOCTEST_CHECK(!get(out[1], mp2, tc.v_ch, tc.base_tick + 0));
            DOCTEST_CHECK(!get(out[2], mp2, tc.w_ch, tc.base_tick + 0));

            // --- tick+1: V and W active, U not active ---
            // MP2 for U must fire (the other two planes have activity).
            if (!get(out[0], mp2, tc.u_ch, tc.base_tick + 1)) {
                spdlog::error("MP2 not set for U ch={} tick={}", tc.u_ch, tc.base_tick + 1);
                ++nfail;
            }
            DOCTEST_CHECK(get(out[0], mp2, tc.u_ch, tc.base_tick + 1));

            // MP3 must NOT fire for U (U channel is not active, so Uc=false).
            DOCTEST_CHECK(!get(out[0], mp3, tc.u_ch, tc.base_tick + 1));

            // --- tick+2: U and W active, V not active ---
            // MP2 for V must fire.
            if (!get(out[1], mp2, tc.v_ch, tc.base_tick + 2)) {
                spdlog::error("MP2 not set for V ch={} tick={}", tc.v_ch, tc.base_tick + 2);
                ++nfail;
            }
            DOCTEST_CHECK(get(out[1], mp2, tc.v_ch, tc.base_tick + 2));

            // MP3 must NOT fire for V (V channel is not active, so Vc=false).
            DOCTEST_CHECK(!get(out[1], mp3, tc.v_ch, tc.base_tick + 2));
        }

        spdlog::info("per-cell checks complete, {} failures", nfail);

        // ----------------------------------------------------------------
        // Global mutual-exclusivity check.
        //
        // For any channel c and tick t, uvw[view][c,t] is a single Boolean.
        // If it is True, Uc=True for every cell mapping to c, making ~Uc=False
        // and therefore MP2 = (~Uc)&... = False.  Conversely, if it is False,
        // Uc=False making MP3 = Uc&... = False.  Hence MP2 and MP3 can never
        // both be True for the same (channel, tick) in any view.
        // ----------------------------------------------------------------
        for (int view = 0; view < 3; ++view) {
            auto mp3_img = out[view][mp3];  // (nchan[view], ntick)
            auto mp2_img = out[view][mp2];  // (nchan[view], ntick)
            auto both = mp3_img & mp2_img;
            bool any_both = both.any().item<bool>();
            if (any_both) {
                spdlog::error("view={}: MP3 and MP2 both True at the same pixel", view);
            }
            DOCTEST_CHECK(!any_both);
        }

        spdlog::info("mutual exclusivity check complete");
    }

}  // DOCTEST_TEST_SUITE("CellViews")
