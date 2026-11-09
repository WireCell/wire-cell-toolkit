// Tests of the TPG integer arithmetic against values and logic transcribed
// from dunetrigger TriggerSim/TPAlgTools (TPAlgTPCTool.hh and
// TPAlgTPCSimpleThreshold.hh at dunetrigger develop b65e250).

#include "WireCellTrig/TPGMath.h"

#include "WireCellUtil/doctest.h"

#include <algorithm>
#include <limits>
#include <random>
#include <vector>

using namespace WireCell::Trig;

TEST_CASE("trig tpgmath avx2_divide")
{
    CHECK(avx2_divide(100, 10) == 10);
    CHECK(avx2_divide(0, 10) == 0);
    CHECK(avx2_divide(-100, 10) == -10);
    // 1/10 * 2^15 truncates to 3276 so the result rounds low at times.
    CHECK(avx2_divide(15, 10) == 1);
    CHECK(avx2_divide(16, 10) == 2);
    CHECK(avx2_divide(32767, 10) == 3276);
}

TEST_CASE("trig tpgmath mode_pedestal")
{
    CHECK(mode_pedestal(std::vector<int16_t>{}) == 0);
    CHECK(mode_pedestal(std::vector<int16_t>{900, 901, 901, 902}) == 901);
    // A tie goes to the smallest value.
    CHECK(mode_pedestal(std::vector<int16_t>{902, 902, 900, 900}) == 900);
}

TEST_CASE("trig tpgmath frugal pedestal")
{
    FrugalPedestal fp;
    fp.accum_limit = 2;
    fp.pedestal = 100;
    // Three samples above moves accum past the limit and the pedestal up one.
    CHECK(fp(105) == 100);
    CHECK(fp(105) == 100);
    CHECK(fp(105) == 101);
    CHECK(fp.accum == 0);
    // Equal samples leave it alone.
    CHECK(fp(101) == 101);
    CHECK(fp.accum == 0);
    // Converges to a constant input.
    for (int i = 0; i < 1000; ++i) fp(50);
    CHECK(fp.pedestal == 50);
}

TEST_CASE("trig tpgmath running sums saturate")
{
    RunningSum rs;
    for (int i = 0; i < 100; ++i) rs(1000);
    CHECK(rs.sum > 0);
    RunningSum big;
    big.r_value = 10;
    big.scale_factor = 10;
    for (int i = 0; i < 100; ++i) big(30000);
    CHECK(big.sum == std::numeric_limits<int16_t>::max());

    // Below INT16_MIN tpglibs wraps and dunetrigger clamps.
    RunningSum wrap;
    wrap.r_value = 10;
    wrap.scale_factor = 10;
    RunningSum clamp = wrap;
    clamp.saturate_low = true;
    for (int i = 0; i < 100; ++i) {
        wrap(-30000);
        clamp(-30000);
    }
    CHECK(clamp.sum == std::numeric_limits<int16_t>::min());
    CHECK(wrap.sum != std::numeric_limits<int16_t>::min());

    AbsRunningSum ars;
    CHECK(ars.r_value == 8);
    ars.r_value = 10;
    ars.scale_factor = 10;
    for (int i = 0; i < 100; ++i) ars(-30000);
    CHECK(ars.sum == std::numeric_limits<int16_t>::max());
}

namespace {
    // tpglibs NaiveRunSumProcessor::process() for one channel.
    int16_t ref_tpglibs_rs(int16_t& rs, int16_t x, int16_t r, int16_t s)
    {
        int32_t scaled_rs = avx2_divide(rs, 10);
        scaled_rs *= r;
        int32_t scaled_signal = avx2_divide(x, 10);
        scaled_signal *= s;
        int32_t intermediate = scaled_signal + scaled_rs;
        rs = std::min(intermediate, (int32_t) INT16_MAX);
        return rs;
    }
    // dunetrigger TPAlgTPCRunningSum::run_sum().
    int16_t ref_dunetrigger_rs(int16_t& rs, int16_t sample, int r, int s)
    {
        const int16_t tmp = r * avx2_divide(rs, 10) + avx2_divide(sample, 10) * s;
        if (sample > 0 && rs > 0 && tmp < 0)
            rs = std::numeric_limits<int16_t>::max();
        else if (sample < 0 && rs < 0 && tmp > 0)
            rs = std::numeric_limits<int16_t>::min();
        else
            rs = tmp;
        return rs;
    }
}  // namespace

TEST_CASE("trig tpgmath running sum matches tpglibs and dunetrigger")
{
    std::mt19937 rng(54321);
    std::uniform_int_distribution<int> dist(-8000, 8000);
    RunningSum daq, dt;
    dt.saturate_low = true;
    int16_t ref_daq = 0, ref_dt = 0;
    for (int i = 0; i < 100000; ++i) {
        const int16_t x = dist(rng);
        REQUIRE(daq(x) == ref_tpglibs_rs(ref_daq, x, 9, 5));
        REQUIRE(dt(x) == ref_dunetrigger_rs(ref_dt, x, 9, 5));
    }
}

TEST_CASE("trig tpgmath threshold hitter minimum length")
{
    ThresholdHitter th;
    th.threshold = 10;
    th.min_samples_over_threshold = 2;
    ThresholdHit hit;
    const std::vector<int16_t> wf{0, 50, 0, 50, 50, 0};
    std::vector<ThresholdHit> hits;
    for (size_t i = 0; i < wf.size(); ++i) {
        if (th(i, wf[i], hit)) hits.push_back(hit);
    }
    REQUIRE(hits.size() == 1);
    CHECK(hits[0].sample_start == 3);
    CHECK(hits[0].samples_over_threshold == 2);
}

TEST_CASE("trig tpgmath threshold hitter")
{
    ThresholdHitter th;
    th.threshold = 10;
    ThresholdHit hit;
    const std::vector<int16_t> wf{0, 11, 30, 20, 10, 0, 50};
    std::vector<ThresholdHit> hits;
    for (size_t i = 0; i < wf.size(); ++i) {
        if (th(1000 + i, wf[i], hit)) hits.push_back(hit);
    }
    REQUIRE(hits.size() == 1);
    CHECK(hits[0].sample_start == 1001);
    CHECK(hits[0].samples_over_threshold == 3);
    CHECK(hits[0].samples_to_peak == 1);
    CHECK(hits[0].adc_integral == 61);
    CHECK(hits[0].adc_peak == 30);
    // The last sample opened a hit that only flush() closes.
    CHECK(th.open());
    REQUIRE(th.flush(1007, hit));
    CHECK(hit.sample_start == 1006);
    CHECK(hit.samples_over_threshold == 1);
    CHECK(!th.open());
    CHECK(!th.flush(1008, hit));
}

namespace {
    // Transcription of TPAlgTPCSimpleThreshold::process_waveform() with the
    // geometry lookup replaced by a given threshold and times kept in ticks.
    struct RefTP {
        int64_t time_start, time_over_threshold, time_peak;
        uint32_t adc_integral;
        uint16_t adc_peak;
    };
    std::vector<RefTP> reference_simple_threshold(const std::vector<short>& adcs, int16_t threshold,
                                                  int accum_limit, int64_t start)
    {
        std::vector<RefTP> out;
        int16_t pedestal = mode_pedestal(adcs);
        int16_t accum = 0;
        bool prev_was_over = false;
        uint32_t hit_charge = 0;
        uint16_t hit_tover = 0, hit_peak_adc = 0, hit_peak_time = 0;
        for (size_t i_t = 0; i_t < adcs.size(); ++i_t) {
            int16_t sample = adcs[i_t];
            if (sample > pedestal) ++accum;
            if (sample < pedestal) --accum;
            if (accum > accum_limit) { ++pedestal; accum = 0; }
            if (accum < -1 * accum_limit) { --pedestal; accum = 0; }
            sample -= pedestal;
            bool is_over = sample > threshold;
            if (is_over) {
                hit_charge += sample;
                if (sample > hit_peak_adc) {
                    hit_peak_adc = (uint16_t) sample;
                    hit_peak_time = hit_tover;
                }
                ++hit_tover;
            }
            if (prev_was_over && !is_over) {
                out.push_back({start + (int64_t) (i_t - hit_tover), hit_tover,
                               start + (int64_t) (i_t - hit_tover + hit_peak_time), hit_charge, hit_peak_adc});
                hit_charge = 0;
                hit_tover = 0;
                hit_peak_adc = 0;
                hit_peak_time = 0;
            }
            prev_was_over = is_over;
        }
        return out;
    }
}  // namespace

TEST_CASE("trig tpgmath simple threshold matches dunetrigger")
{
    std::mt19937 rng(12345);
    std::normal_distribution<double> noise(0.0, 4.0);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    const int16_t threshold = 20;
    const int accum_limit = 10;
    const int64_t start = 5000;

    std::vector<short> adcs(6000);
    for (size_t i = 0; i < adcs.size(); ++i) {
        adcs[i] = (short) (900 + noise(rng));
    }
    // sprinkle unipolar pulses
    for (size_t i = 100; i + 20 < adcs.size(); i += 137) {
        const double amp = 100 * uni(rng);
        for (int k = 0; k < 10; ++k) adcs[i + k] += (short) (amp * (k < 5 ? k : 10 - k) / 5.0);
    }

    auto ref = reference_simple_threshold(adcs, threshold, accum_limit, start);
    REQUIRE(ref.size() > 10);

    FrugalPedestal fp;
    fp.accum_limit = accum_limit;
    fp.pedestal = mode_pedestal(adcs);
    ThresholdHitter th;
    th.threshold = threshold;
    std::vector<ThresholdHit> got;
    ThresholdHit hit;
    for (size_t i = 0; i < adcs.size(); ++i) {
        int16_t sample = adcs[i];
        sample -= fp(sample);
        if (th(start + i, sample, hit)) got.push_back(hit);
    }
    REQUIRE(got.size() == ref.size());
    for (size_t i = 0; i < got.size(); ++i) {
        CHECK(got[i].sample_start == ref[i].time_start);
        CHECK(got[i].samples_over_threshold == ref[i].time_over_threshold);
        CHECK(got[i].sample_start + got[i].samples_to_peak == ref[i].time_peak);
        CHECK(got[i].adc_integral == ref[i].adc_integral);
        CHECK(got[i].adc_peak == ref[i].adc_peak);
    }
}
