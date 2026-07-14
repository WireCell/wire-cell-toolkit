#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Exceptions.h"

using namespace std;
using namespace WireCell::Waveform;
using namespace WireCell;

TEST_SUITE("issue24") {

    TEST_CASE("median on non-empty waveforms") {
        int nsamples = 10;
        while (nsamples > 0) {
            realseq_t wave(nsamples, 0);
            CHECK_NOTHROW(median(wave));
            --nsamples;
        }
    }

    TEST_CASE("median throws on empty waveform") {
        realseq_t wave;
        CHECK_THROWS_AS(median(wave), ValueError);
    }

    TEST_CASE("percentile throws for out-of-range percentage") {
        realseq_t wave = {6.9, 9.6};

        SUBCASE("negative percentage") {
            CHECK_THROWS_AS(percentile(wave, -0.1), ValueError);
        }

        SUBCASE("percentage over 1") {
            CHECK_THROWS_AS(percentile(wave, 1.1), ValueError);
        }
    }

    TEST_CASE("median correctness") {
        realseq_t wave;
        wave.push_back(6.9);
        wave.push_back(9.6);
        spdlog::debug("median of 2 elements: {}", median(wave));
        CHECK(std::abs(9.6 - median(wave)) < 0.001);

        wave.push_back(0.0);
        spdlog::debug("median of 3 elements: {}", median(wave));
        CHECK(std::abs(6.9 - median(wave)) < 0.001);

        wave.push_back(10.0);
        spdlog::debug("median of 4 elements: {}", median(wave));
        CHECK(std::abs(9.6 - median(wave)) < 0.001);
    }

}  // TEST_SUITE("issue24")
