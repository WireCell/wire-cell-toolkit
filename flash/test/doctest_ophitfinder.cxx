#include "WireCellUtil/doctest.h"

#include "WireCellFlash/OpHitFinder.h"

#include <cmath>

using namespace WireCell;

static Configuration default_algo()
{
    Configuration algo;
    algo["adc_threshold"] = 3.0;
    algo["nsigma_threshold"] = 1.0;
    algo["tail_adc_threshold"] = 1.0;
    algo["tail_nsigma_threshold"] = 1.0;
    algo["end_adc_threshold"] = 1.0;
    algo["end_nsigma_threshold"] = 1.0;
    algo["min_pulse_width"] = 1;
    algo["num_presample"] = 2;
    algo["num_postsample"] = 2;
    return algo;
}

static Configuration split_algo(bool enable)
{
    Configuration algo = default_algo();
    algo["split_enable"] = enable;
    algo["split_min_prominence"] = 0.4;
    algo["split_min_peak"] = 3.0;
    algo["split_min_separation"] = 2;
    return algo;
}

// Two pulses riding a common decaying scintillation tail: the second
// peaks before the first's tail has fallen back to baseline, so the
// valley between them never reaches the tail threshold and
// sliding_window merges them into ONE pulse.
static std::vector<short> two_on_tail()
{
    std::vector<short> wf(200, 0);
    auto add = [&](int c, double amp, double tau) {
        for (int i = c; i < 200; ++i) {
            wf[i] = std::max<short>(wf[i], short(amp * std::exp(-(i - c) / tau)));
        }
    };
    add(52, 80.0, 20.0);  // pulse 1: peak 80 at t=52
    add(85, 50.0, 20.0);  // pulse 2: peak 50 at t=85, on pulse 1's ~16-ADC tail
    return wf;
}

TEST_CASE("ophitfinder sliding window single pulse")
{
    // Flat baseline at 10 with a triangular pulse.
    std::vector<short> wf(100, 10);
    const std::vector<short> pulse = {15, 30, 60, 100, 70, 40, 20, 12};
    for (size_t i = 0; i < pulse.size(); ++i) wf[40 + i] = pulse[i];

    auto pulses = Flash::OpHitFinder::sliding_window(wf, 10.0, 0.0, default_algo());
    REQUIRE(pulses.size() == 1);
    const auto& p = pulses[0];
    CHECK(p.t_max == 43);          // peak position
    CHECK(p.peak == doctest::Approx(90.0));  // 100 - baseline
    CHECK(p.t_start <= 40);
    CHECK(p.t_end >= 46);
    // Area is the baseline-subtracted sum over the pulse extent (the
    // pre-samples only contribute their positive part).
    double expected = 0;
    for (short v : pulse) expected += v - 10;
    CHECK(p.area == doctest::Approx(expected).epsilon(0.1));
}

TEST_CASE("ophitfinder sliding window two pulses")
{
    std::vector<short> wf(200, 0);
    for (int i = 0; i < 5; ++i) {
        wf[50 + i] = 50 - 10 * std::abs(i - 2);  // peak 50 at 52
        wf[120 + i] = 25 - 5 * std::abs(i - 2);  // peak 25 at 122
    }
    auto pulses = Flash::OpHitFinder::sliding_window(wf, 0.0, 0.0, default_algo());
    REQUIRE(pulses.size() == 2);
    CHECK(pulses[0].t_max == 52);
    CHECK(pulses[0].peak == doctest::Approx(50.0));
    CHECK(pulses[1].t_max == 122);
    CHECK(pulses[1].peak == doctest::Approx(25.0));
}

TEST_CASE("ophitfinder sliding window quiet waveform")
{
    std::vector<short> wf(100, 5);
    auto pulses = Flash::OpHitFinder::sliding_window(wf, 5.0, 0.0, default_algo());
    CHECK(pulses.empty());
}

TEST_CASE("ophitfinder split disabled keeps the merged pulse")
{
    auto wf = two_on_tail();
    auto pulses = Flash::OpHitFinder::sliding_window(wf, 0.0, 0.0, default_algo());
    REQUIRE(pulses.size() == 1);  // the two pulses merge

    auto subs = Flash::OpHitFinder::split_pulse(wf, 0.0, pulses[0], split_algo(false));
    REQUIRE(subs.size() == 1);    // disabled -> verbatim pass-through
    CHECK(subs[0].t_start == pulses[0].t_start);
    CHECK(subs[0].t_end == pulses[0].t_end);
    CHECK(subs[0].t_max == pulses[0].t_max);
    CHECK(subs[0].peak == doctest::Approx(pulses[0].peak));
    CHECK(subs[0].area == doctest::Approx(pulses[0].area));
}

TEST_CASE("ophitfinder split enabled separates the two pulses")
{
    auto wf = two_on_tail();
    auto pulses = Flash::OpHitFinder::sliding_window(wf, 0.0, 0.0, default_algo());
    REQUIRE(pulses.size() == 1);

    auto subs = Flash::OpHitFinder::split_pulse(wf, 0.0, pulses[0], split_algo(true));
    REQUIRE(subs.size() == 2);
    CHECK(subs[0].t_max == 52);
    CHECK(subs[0].peak == doctest::Approx(80.0));
    CHECK(subs[1].t_max == 85);
    CHECK(subs[1].peak == doctest::Approx(50.0));
    // The sub-windows partition the parent window: no gap, no overlap.
    CHECK(subs[0].t_start == pulses[0].t_start);
    CHECK(subs[1].t_end == pulses[0].t_end);
    CHECK(subs[1].t_start == subs[0].t_end + 1);
    // Each area is the unclamped baseline-subtracted sum of its sub-window.
    double a1 = 0;
    for (int i = subs[1].t_start; i <= subs[1].t_end; ++i) a1 += wf[i];
    CHECK(subs[1].area == doctest::Approx(a1));
}

TEST_CASE("ophitfinder split single peak is a no-op when enabled")
{
    std::vector<short> wf(100, 10);
    const std::vector<short> pulse = {15, 30, 60, 100, 70, 40, 20, 12};
    for (size_t i = 0; i < pulse.size(); ++i) wf[40 + i] = pulse[i];

    auto pulses = Flash::OpHitFinder::sliding_window(wf, 10.0, 0.0, default_algo());
    REQUIRE(pulses.size() == 1);
    auto subs = Flash::OpHitFinder::split_pulse(wf, 10.0, pulses[0], split_algo(true));
    CHECK(subs.size() == 1);  // one peak -> nothing to split
}
