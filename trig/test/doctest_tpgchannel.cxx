// Tests of the per-channel TPG chain.  The reference is a transcription of
// dunetrigger TPAlgTPCAbsRunningSum::process_waveform() (develop b65e250)
// with the geometry lookup replaced by a given threshold.

#include "WireCellTrig/TPGChannel.h"

#include "WireCellUtil/doctest.h"

#include <cstdlib>
#include <limits>
#include <random>
#include <vector>

using namespace WireCell::Trig;

namespace {

    std::vector<int16_t> make_waveform(size_t nsamples, unsigned seed)
    {
        std::mt19937 rng(seed);
        std::normal_distribution<double> noise(0.0, 4.0);
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        std::vector<int16_t> adcs(nsamples);
        for (auto& adc : adcs) adc = (int16_t) (900 + noise(rng));
        // bipolar pulses, as on an induction plane
        for (size_t i = 100; i + 20 < adcs.size(); i += 137) {
            const double amp = 150 * uni(rng);
            for (int k = 0; k < 6; ++k) adcs[i + k] += (int16_t) (amp * (k < 3 ? k : 6 - k) / 3.0);
            for (int k = 6; k < 12; ++k) adcs[i + k] -= (int16_t) (amp * (k < 9 ? k - 6 : 12 - k) / 3.0);
        }
        return adcs;
    }

    struct RefHit {
        int64_t start, tover, peak;
        uint32_t charge;
        uint16_t peak_adc;
    };

    void frugal(int16_t sample, int16_t& pedestal, int16_t& accum, int limit)
    {
        if (sample > pedestal) ++accum;
        if (sample < pedestal) --accum;
        if (accum > limit) { ++pedestal; accum = 0; }
        if (accum < -1 * limit) { --pedestal; accum = 0; }
    }

    std::vector<RefHit> reference_abs_running_sum(const std::vector<int16_t>& adcs, int16_t threshold)
    {
        const int accum_limit = 10, scale_factor = 5, r_value = 8, stabilisation_period = 1000;
        std::vector<RefHit> out;
        int16_t pedestal = mode_pedestal(adcs), accum = 0, running_sum = 0;
        int16_t abs_pedestal = 0, abs_accum = 0;
        bool prev_was_over = false;
        uint32_t hit_charge = 0;
        uint16_t hit_tover = 0, hit_peak_adc = 0, hit_peak_time = 0;
        auto abs_run_sum = [&](int16_t sample) {
            const int16_t tmp = r_value * avx2_divide(running_sum, 10) + avx2_divide(std::abs(sample), 10) * scale_factor;
            running_sum = tmp < 0 ? std::numeric_limits<int16_t>::max() : tmp;
        };
        for (int jj = 0; jj < stabilisation_period; ++jj) {
            int16_t sample = adcs[jj];
            frugal(sample, pedestal, accum, accum_limit);
            sample -= pedestal;
            abs_run_sum(sample);
            sample = running_sum;
            frugal(sample, abs_pedestal, abs_accum, accum_limit);
        }
        for (size_t i_t = 0; i_t < adcs.size(); ++i_t) {
            int16_t sample = adcs[i_t];
            frugal(sample, pedestal, accum, accum_limit);
            sample -= pedestal;
            abs_run_sum(sample);
            sample = running_sum;
            frugal(sample, abs_pedestal, abs_accum, accum_limit);
            sample -= abs_pedestal;
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
                out.push_back({(int64_t) (i_t - hit_tover), hit_tover, hit_peak_time, hit_charge, hit_peak_adc});
                hit_charge = 0;
                hit_tover = 0;
                hit_peak_adc = 0;
                hit_peak_time = 0;
            }
            prev_was_over = is_over;
        }
        return out;
    }

    TPGConfig abs_running_sum_config(int16_t threshold)
    {
        TPGConfig cfg;
        cfg.pedestal_init = TPGConfig::PedestalInit::mode;
        cfg.filter = TPGConfig::Filter::abs_running_sum;
        cfg.r_value = 8;
        cfg.scale_factor = 5;
        cfg.filter_pedestal = true;
        cfg.threshold = threshold;
        return cfg;
    }
}  // namespace

TEST_CASE("trig tpgchannel abs running sum matches dunetrigger")
{
    const int16_t threshold = 80;
    const auto adcs = make_waveform(6000, 777);
    const auto ref = reference_abs_running_sum(adcs, threshold);
    REQUIRE(ref.size() > 10);

    TPGChannel tpg(abs_running_sum_config(threshold));
    tpg.settle(adcs, 1000);
    std::vector<ThresholdHit> got;
    tpg(adcs, 0, got);
    REQUIRE(got.size() == ref.size());
    for (size_t i = 0; i < got.size(); ++i) {
        CHECK(got[i].sample_start == ref[i].start);
        CHECK(got[i].samples_over_threshold == ref[i].tover);
        CHECK(got[i].samples_to_peak == ref[i].peak);
        CHECK(got[i].adc_integral == ref[i].charge);
        CHECK(got[i].adc_peak == ref[i].peak_adc);
    }
}

TEST_CASE("trig tpgchannel streaming in pieces equals whole")
{
    const auto adcs = make_waveform(6000, 999);
    TPGConfig cfg;
    cfg.threshold = 20;
    cfg.pedestal_init = TPGConfig::PedestalInit::fixed;
    cfg.pedestal = 900;

    const int64_t first = 123456;
    TPGChannel whole(cfg);
    std::vector<ThresholdHit> whole_hits;
    whole(adcs, first, whole_hits);
    whole.flush(first + adcs.size(), whole_hits);
    REQUIRE(whole_hits.size() > 10);

    // Cut at awkward places, including inside pulses.
    TPGChannel pieces(cfg);
    std::vector<ThresholdHit> piece_hits;
    const std::vector<size_t> cuts{0, 101, 102, 1000, 2741, 2745, 6000};
    for (size_t ic = 0; ic + 1 < cuts.size(); ++ic) {
        std::vector<int16_t> piece(adcs.begin() + cuts[ic], adcs.begin() + cuts[ic + 1]);
        pieces(piece, first + cuts[ic], piece_hits);
    }
    pieces.flush(first + adcs.size(), piece_hits);

    REQUIRE(piece_hits.size() == whole_hits.size());
    for (size_t i = 0; i < whole_hits.size(); ++i) {
        CHECK(piece_hits[i].sample_start == whole_hits[i].sample_start);
        CHECK(piece_hits[i].samples_over_threshold == whole_hits[i].samples_over_threshold);
        CHECK(piece_hits[i].adc_integral == whole_hits[i].adc_integral);
    }
}

TEST_CASE("trig tpgchannel negative threshold disables and reset forgets")
{
    const auto adcs = make_waveform(2000, 1);
    TPGConfig cfg;
    cfg.threshold = -1;
    TPGChannel off(cfg);
    std::vector<ThresholdHit> hits;
    off(adcs, 0, hits);
    CHECK(hits.empty());

    cfg.threshold = 20;
    TPGChannel tpg(cfg);
    CHECK(tpg.fresh());
    tpg(adcs, 0, hits);
    CHECK(!tpg.fresh());
    tpg.reset();
    CHECK(tpg.fresh());
}

TEST_CASE("trig tpgchannel zero pedestal start converges as in the DAQ")
{
    // The DAQ starts from zero and moves one ADC per accum_limit+1 samples, so
    // a 900 ADC baseline takes about 10k samples to reach.
    TPGConfig cfg;
    cfg.pedestal_init = TPGConfig::PedestalInit::zero;
    cfg.threshold = 1000;  // no hits wanted here
    TPGChannel tpg(cfg);
    std::vector<int16_t> flat(12000, 900);
    std::vector<ThresholdHit> hits;
    tpg(flat, 0, hits);
    CHECK(hits.empty());
    std::vector<int16_t> pulse(20, 900);
    pulse[10] = 2000;
    tpg(pulse, 12000, hits);
    REQUIRE(hits.size() == 1);
    CHECK(hits[0].sample_start == 12010);
    CHECK(hits[0].adc_peak == 1100);
}
