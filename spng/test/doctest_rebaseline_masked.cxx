/** Tests for rebaseline_masked() / rebaseline_masked_rows() and the
 * SPNGMaskRebaseliner node.  See Rebaseline.h for the algorithm.
 *
 * The first test case reproduces the bug that motivated this code: with a
 * tight ROI mask, rebaseline_zero() anchors its linear baseline on the first
 * and last samples *inside* the ROI, which sit on the pulse tails, so it
 * removes real charge.  Anchoring outside the ROI recovers it.
 */
#include "WireCellSpng/Rebaseline.h"
#include "WireCellSpng/MaskRebaseliner.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/Testing.h"
#include "WireCellSpng/Util.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"

using namespace WireCell;
using namespace WireCell::SPNG;

namespace {

    bool tensors_are_close(const torch::Tensor& a, const torch::Tensor& b, double tolerance = 1e-4) {
        if (a.sizes().vec() != b.sizes().vec()) {
            spdlog::error("Tensor size mismatch: {} vs {}", to_string(a), to_string(b));
            return false;
        }
        const bool ok = a.cpu().allclose(b.cpu(), tolerance, tolerance);
        if (!ok) {
            auto diff = (a.cpu() - b.cpu()).abs();
            spdlog::error("tensor mismatch: max|a-b|={} at flat index {} (a={} b={})",
                          diff.max().item<double>(), diff.argmax().item<int64_t>(),
                          a.cpu().flatten()[diff.argmax()].item<double>(),
                          b.cpu().flatten()[diff.argmax()].item<double>());
            if (a.numel() <= 64) {
                spdlog::error("a: {}\nb: {}", to_string(a, true), to_string(b, true));
            }
        }
        return ok;
    }

    MaskRebaselineOptions opts(int64_t pad, int64_t width = 1, int64_t dim = -1) {
        MaskRebaselineOptions o;
        o.anchor_pad = pad;
        o.anchor_width = width;
        o.dim = dim;
        return o;
    }

    torch::Tensor t1d(std::vector<float> v) {
        return torch::tensor(v, torch::kFloat32);
    }

    /// A Gaussian pulse on a flat zero baseline, nt ticks long, centered at t0.
    torch::Tensor gauss_pulse(int64_t nt, double t0, double sigma, double amp) {
        auto t = torch::arange(nt, torch::kFloat32);
        return amp * torch::exp(-0.5 * torch::pow((t - t0) / sigma, 2));
    }

    /// Random mask rows with ROIs of length >= 2 separated by gaps >= 2, and
    /// gaps at the tensor edges of 0 or >= 2 samples, as rebaseline_zero()
    /// requires to see the same ROIs (a single zero is a "crossing" to it,
    /// also at the tensor edge).
    torch::Tensor random_mask(int64_t nrows, int64_t nt, int64_t seed) {
        torch::manual_seed(seed);
        auto mask = torch::zeros({nrows, nt}, torch::kBool);
        auto acc = mask.accessor<bool, 2>();
        for (int64_t r = 0; r < nrows; ++r) {
            int64_t t = (r % 3 == 0) ? 0 : 2 + (r * 7) % 11;
            while (t < nt) {
                const int64_t len = 2 + (r * 13 + t * 3) % 9;
                int64_t end = std::min(nt, t + len);
                if (nt - end == 1) end = nt;       // no single trailing zero
                for (int64_t k = t; k < end; ++k) acc[r][k] = true;
                t = end + 2 + (r + t) % 7;         // gap >= 2
            }
        }
        return mask;
    }
}

TEST_SUITE("spng rebaseline_masked") {

    TEST_CASE("tight ROI mask: outside anchors keep the charge") {
        // Pulse sigma=3 ticks, amplitude 1000, tight mask |t-t0| <= 2 sigma so the
        // ROI edge samples sit at exp(-2) ~ 13.5% of peak, like SPNG's DNN ROIs.
        const int64_t nt = 64;
        auto wave = gauss_pulse(nt, 32.0, 3.0, 1000.0).unsqueeze(0);
        auto t = torch::arange(nt, torch::kFloat32);
        auto mask = (torch::abs(t - 32.0) <= 6.0).unsqueeze(0);
        const double q_masked = (wave * mask).sum().item<double>();

        // The legacy algorithm on the masked waveform: anchors on the ROI edge samples.
        auto legacy = rebaseline_zero(wave * mask);
        const double q_legacy = legacy.sum().item<double>() / q_masked;
        CHECK(q_legacy < 0.92);   // documents the bug: ~10% of the charge removed

        // Anchors 4 ticks outside the ROI (on the pulse tails, ~1% of peak).
        auto fixed = rebaseline_masked(wave, mask, opts(4, 1));
        const double q_fixed = fixed.sum().item<double>() / q_masked;
        CHECK(q_fixed > 0.99);
        // Output region is still the tight mask.
        CHECK(torch::all((fixed != 0) == mask).item<bool>());
    }

    TEST_CASE("pad=0 width=1 reproduces rebaseline_zero on masked input") {
        // Fixtures from doctest_rebaseline.cxx "spng rebaseline_zero" whose ROIs are
        // bounded by >= 2 zeros (or the tensor edge).
        std::vector<torch::Tensor> inputs = {
            t1d({0, 0, 5, 10, 5, 0, 0}),
            t1d({5, 10, 5, 0, 0}),
            t1d({0, 0, 5, 10, 5, 0, 0, 3, 6, 3, 0, 0}),
            t1d({1, 2, 3, 4, 5}),                       // whole-tensor ROI
            t1d({0, 0, 0, 0}),                          // all zeros
            t1d({0, 0, 2, -8, 2, 0, 0}),                // negative excursion
        };
        for (auto& input : inputs) {
            auto mask = input != 0;
            auto expected = rebaseline_zero(input);
            auto actual = rebaseline_masked(input, mask, opts(0, 1));
            CHECK(tensors_are_close(actual, expected));
        }
        // 2D, dim=1 (the production layout: rows are channels).
        auto in2d = torch::stack({t1d({0, 0, 5, 10, 5, 0, 0}), t1d({3, 6, 3, 0, 0, 7, 7})});
        CHECK(tensors_are_close(rebaseline_masked(in2d, in2d != 0, opts(0, 1, 1)),
                                rebaseline_zero(in2d, 1)));
        // dim=0 exercises the permutation path.
        auto in2dT = in2d.transpose(0, 1).contiguous();
        CHECK(tensors_are_close(rebaseline_masked(in2dT, in2dT != 0, opts(0, 1, 0)),
                                rebaseline_zero(in2dT, 0)));

        // Randomized: many rows, positive waveform so no accidental zeros.
        auto mask = random_mask(64, 512, 42);
        auto wave = torch::rand({64, 512}) * 100.0 + 1.0;
        auto masked = wave * mask;
        CHECK(tensors_are_close(rebaseline_masked(wave, mask, opts(0, 1)),
                                rebaseline_zero(masked), 1e-5));
    }

    TEST_CASE("min_roi_size") {
        // Single-tick ROI: legacy leaves it untouched (end-start=0 <= min_roi_size-1).
        auto input = t1d({0, 0, 7, 0, 0});
        CHECK(tensors_are_close(rebaseline_masked(input, input != 0, opts(0)), input));
        // min_roi_size=0: the one-sample "line" is the sample itself -> 0.
        auto o = opts(0);
        o.min_roi_size = 0;
        CHECK(tensors_are_close(rebaseline_masked(input, input != 0, o), torch::zeros_like(input)));
    }

    TEST_CASE("anchor window placement") {
        // wave index: 0  1  2  3  4  5  6  7  8  9 10 11
        auto wave = t1d({1, 2, 3, 4, 9, 9, 9, 4, 3, 2, 1, 0});
        auto mask = torch::zeros({12}, torch::kBool);
        mask.slice(0, 4, 7).fill_(true);                // ROI = ticks 4..6

        SUBCASE("pad=2 width=1: anchors at ticks 2 and 8") {
            // a = wave[2] = 3, b = wave[8] = 3 -> flat baseline 3.
            auto out = rebaseline_masked(wave, mask, opts(2, 1));
            CHECK(tensors_are_close(out, t1d({0, 0, 0, 0, 6, 6, 6, 0, 0, 0, 0, 0})));
        }
        SUBCASE("pad=1 width=3: means of ticks 1..3 and 7..9") {
            // a = mean(2,3,4) = 3, b = mean(4,3,2) = 3.
            auto out = rebaseline_masked(wave, mask, opts(1, 3));
            CHECK(tensors_are_close(out, t1d({0, 0, 0, 0, 6, 6, 6, 0, 0, 0, 0, 0})));
        }
        SUBCASE("sloped baseline") {
            auto w2 = t1d({0, 0, 0, 2, 9, 9, 9, 6, 0, 0, 0, 0});
            // pad=1: anchors at ticks 3 and 7 (a=2, b=6).  The baseline is the line
            // through the two anchor positions, i.e. over the dilated ROI 3..7:
            // baseline(t) = 2 + 4*(t-3)/4 -> t=4:3, 5:4, 6:5.
            auto out = rebaseline_masked(w2, mask, opts(1, 1));
            CHECK(tensors_are_close(out, t1d({0, 0, 0, 0, 6, 5, 4, 0, 0, 0, 0, 0})));
        }
        SUBCASE("keep_pad=true outputs the dilated region") {
            auto o = opts(2, 1);
            o.keep_pad = true;
            auto out = rebaseline_masked(wave, mask, o);
            // dilated ROI = ticks 2..8, anchors at its edges (ticks 2 and 8, both 3).
            CHECK(tensors_are_close(out, t1d({0, 0, 0, 1, 6, 6, 6, 1, 0, 0, 0, 0})));
        }
    }

    TEST_CASE("edges and neighbours") {
        SUBCASE("ROI touching the tensor edge falls back to the in-ROI edge sample") {
            auto wave = t1d({5, 10, 5, 0, 0, 0, 0, 0, 0});
            auto mask = wave != 0;
            auto out = rebaseline_masked(wave, mask, opts(4, 1));
            // Dilated ROI is clipped to 0..6: a = wave[0] = 5 (the ROI's own edge
            // sample), b = wave[6] = 0.  baseline(t) = 5 - 5*t/6 -> t=1: 25/6, t=2: 10/3.
            CHECK(tensors_are_close(out, t1d({0, 10 - 25.0f / 6, 5 - 10.0f / 3, 0, 0, 0, 0, 0, 0})));
            // Mirror image at the far edge.
            auto rwave = wave.flip(0);
            auto rout = rebaseline_masked(rwave, rwave != 0, opts(4, 1));
            CHECK(tensors_are_close(rout, out.flip(0)));
        }
        SUBCASE("ROIs closer than 2*pad merge") {
            // Two ROIs with a 3-tick gap, pad=4: the dilated regions overlap, so a
            // single baseline spans both; anchors are 4 ticks outside the outer edges.
            auto wave = t1d({0, 0, 0, 0, 0, 0, 4, 8, 4, 0, 0, 0, 4, 8, 4, 0, 0, 0, 0, 0, 0});
            auto mask = wave != 0;
            auto out = rebaseline_masked(wave, mask, opts(4, 1));
            // anchors: tick 2 (0) and tick 18 (0) -> baseline 0 -> output == wave.
            CHECK(tensors_are_close(out, wave));
        }
        SUBCASE("width larger than the gap is clipped to the gap") {
            // gap of 2 zeros between ROIs; width 5, pad 0: lo window of the second
            // ROI is clipped to [prev_end+1, s] = ticks 5..7 -> mean(0,0,4).
            auto wave = t1d({4, 8, 4, 0, 0, 0, 0, 4, 8, 4, 0, 0});
            auto mask = wave != 0;
            auto out = rebaseline_masked(wave, mask, opts(0, 5));
            // first ROI (0..2): lo window clipped to tick 0 -> a=4; hi window ticks
            // 2..6 -> mean(4,0,0,0,0)=0.8; baseline = 4 + (0.8-4)*(t-0)/2.
            // second ROI (7..9): lo window ticks 3..7 -> mean(0,0,0,0,4)=0.8;
            // hi window ticks 9..11 -> mean(4,0,0)=4/3; baseline = 0.8 + (4/3-0.8)*(t-7)/2.
            auto expect = t1d({0, 8 - 2.4f, 4 - 0.8f, 0, 0, 0, 0,
                               4 - 0.8f, 8 - (0.8f + (4.0f / 3 - 0.8f) * 0.5f), 4 - 4.0f / 3, 0, 0});
            CHECK(tensors_are_close(out, expect));
        }
        SUBCASE("three ROIs per row, multiple rows") {
            auto mask = random_mask(8, 128, 7);
            auto wave = torch::rand({8, 128}) + 1.0;
            auto out = rebaseline_masked(wave, mask, opts(2, 1));
            // Only masked samples are nonzero and the result is finite.
            CHECK(torch::all((out != 0) <= mask).item<bool>());
            CHECK(torch::all(torch::isfinite(out)).item<bool>());
        }
    }

    TEST_CASE("row selection") {
        auto wave = torch::stack({t1d({0, 0, 5, 10, 5, 0, 0}),
                                  t1d({0, 0, 5, 10, 5, 0, 0}),
                                  t1d({0, 0, 5, 10, 5, 0, 0})});
        auto mask = wave != 0;
        auto out = rebaseline_masked_rows(wave, mask, {{1, 2}}, opts(0, 1));
        CHECK(tensors_are_close(out[0], wave[0]));
        CHECK(tensors_are_close(out[1], t1d({0, 0, 0, 5, 0, 0, 0})));
        CHECK(tensors_are_close(out[2], wave[2]));
        // empty = all rows
        auto all = rebaseline_masked_rows(wave, mask, {}, opts(0, 1));
        CHECK(tensors_are_close(all, rebaseline_zero(wave, -1)));
        // bad ranges
        CHECK_THROWS_AS(rebaseline_masked_rows(wave, mask, {{2, 1}}, opts(0, 1)), ValueError);
        CHECK_THROWS_AS(rebaseline_masked_rows(wave, mask, {{0, 4}}, opts(0, 1)), ValueError);
    }

    TEST_CASE("remove_negative and mask dtypes") {
        auto wave = t1d({0, 0, 5, 1, 5, 0, 0});
        auto o = opts(0, 1);
        o.remove_negative = true;
        auto out = rebaseline_masked(wave, wave != 0, o);
        CHECK(tensors_are_close(out, t1d({0, 0, 0, 0, 0, 0, 0})));   // 1-5 = -4 clamped
        // float and int masks behave like bool
        auto fmask = (wave != 0).to(torch::kFloat32);
        auto imask = (wave != 0).to(torch::kInt32);
        auto ref = rebaseline_masked(wave, wave != 0, opts(1, 1));
        CHECK(tensors_are_close(rebaseline_masked(wave, fmask, opts(1, 1)), ref));
        CHECK(tensors_are_close(rebaseline_masked(wave, imask, opts(1, 1)), ref));
        // shape mismatch raises
        CHECK_THROWS_AS(rebaseline_masked(wave, torch::ones({3}, torch::kBool), opts(1, 1)), ValueError);
    }

    TEST_CASE("cpu direct path matches the vectorized path") {
        auto mask = random_mask(48, 300, 5);
        // Add some ROIs touching the tensor edges and some 1-sample gaps.
        mask.index_put_({0, torch::indexing::Slice(0, 5)}, true);
        mask.index_put_({1, torch::indexing::Slice(295, 300)}, true);
        mask.index_put_({2, torch::indexing::Slice(10, 20)}, true);
        mask.index_put_({2, torch::indexing::Slice(21, 30)}, true);
        auto wave = torch::rand({48, 300}) * 100.0 - 20.0;
        for (int64_t pad : {0, 1, 3, 4, 7}) {
            for (int64_t width : {1, 2, 5}) {
                for (bool keep : {false, true}) {
                    for (int64_t minsz : {1, 0, 3}) {
                        auto o = opts(pad, width);
                        o.keep_pad = keep;
                        o.min_roi_size = minsz;
                        auto direct = rebaseline_masked(wave, mask, o);
                        o.force_vectorized = true;
                        auto vec = rebaseline_masked(wave, mask, o);
                        INFO("pad=" << pad << " width=" << width << " keep_pad=" << keep << " min_roi_size=" << minsz);
                        CHECK(tensors_are_close(direct, vec, 1e-4));
                    }
                }
            }
        }
        // rows variant too
        auto o = opts(4, 2);
        auto direct = rebaseline_masked_rows(wave, mask, {{3, 20}, {30, 48}}, o);
        o.force_vectorized = true;
        CHECK(tensors_are_close(direct, rebaseline_masked_rows(wave, mask, {{3, 20}, {30, 48}}, o), 1e-4));
    }

    TEST_CASE("cuda parity") {
        if (!torch::cuda::is_available()) {
            spdlog::warn("CUDA not available, skipping cuda parity test");
            return;
        }
        auto mask = random_mask(64, 512, 99);
        auto wave = torch::rand({64, 512}) * 100.0 + 1.0;
        for (auto o : {opts(0, 1), opts(4, 1), opts(2, 3)}) {
            auto cpu = rebaseline_masked(wave, mask, o);
            auto gpu = rebaseline_masked(wave.cuda(), mask.cuda(), o);
            CHECK(gpu.is_cuda());
            CHECK(tensors_are_close(gpu.cpu(), cpu, 1e-4));
        }
    }

    TEST_CASE("spng maskrebaseliner node") {
        auto wave = torch::stack({t1d({0, 0, 5, 10, 5, 0, 0}),
                                  t1d({1, 2, 3, 4, 9, 9, 9})});    // row 1 has no ROI
        auto mask = torch::stack({t1d({0, 0, 1, 1, 1, 0, 0}),
                                  t1d({0, 0, 0, 0, 0, 0, 0})}) != 0;

        Configuration wmd, mmd;
        wmd["datapath"] = "/traces/gauss/v0";
        wmd["tag"] = "gauss";
        mmd["datapath"] = "/traces/Threshold/v0";
        mmd["tag"] = "roi";

        MaskRebaseliner::input_vector inv = {
            std::make_shared<SimpleTorchTensor>(wave, wmd),
            std::make_shared<SimpleTorchTensor>(mask, mmd)};

        MaskRebaseliner node;
        auto cfg = node.default_configuration();
        CHECK(cfg["anchor_pad"].asInt() == 4);
        cfg["anchor_pad"] = 0;
        cfg["datapath_format"] = "/traces/MaskRebaseliner/{tag}";
        node.configure(cfg);

        ITorchTensor::pointer out;
        CHECK(node(inv, out));
        REQUIRE(out);
        CHECK(tensors_are_close(out->tensor(), torch::stack({t1d({0, 0, 0, 5, 0, 0, 0}),
                                                             torch::zeros({7})})));
        CHECK(out->metadata()["tag"].asString() == "signal");
        CHECK(out->metadata()["datapath"].asString() == "/traces/MaskRebaseliner/signal");

        // enabled=false -> wave*mask
        cfg["enabled"] = false;
        node.configure(cfg);
        CHECK(node(inv, out));
        REQUIRE(out);
        CHECK(tensors_are_close(out->tensor(), wave * mask));

        // rows restrict the work
        cfg["enabled"] = true;
        cfg["rows"] = Json::arrayValue;
        Json::Value r(Json::arrayValue); r.append(1); r.append(2);
        cfg["rows"].append(r);
        node.configure(cfg);
        CHECK(node(inv, out));
        CHECK(tensors_are_close(out->tensor(), wave * mask));   // row 0 untouched, row 1 empty

        // EOS passes through
        MaskRebaseliner::input_vector eos = {nullptr, nullptr};
        CHECK(node(eos, out));
        CHECK(out == nullptr);

        // factory registration
        PluginManager::instance().add("WireCellSpng");
        auto icfg = Factory::lookup<IConfigurable>("SPNGMaskRebaseliner", "test");
        CHECK(icfg);
    }
}
