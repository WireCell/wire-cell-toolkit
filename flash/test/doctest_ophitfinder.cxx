#include "WireCellUtil/doctest.h"

#include "WireCellFlash/OpHitFinder.h"

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
