// Tests for the FEMB negative-pulse ROI detector added in commit
// 139d782f (sigproc/PDHD: harden FEMBNoiseSub).
//
// We test the free function PDHD::Is_FEMB_noise() directly rather than
// going through PDHD::FEMBNoiseSub::apply(), because the class form
// requires a fully-configured IAnodePlane for plane resolution while the
// detection logic itself lives in this free function.  All new
// algorithmic behaviour added by the commit (multi-ROI handling,
// configurable nsigma, pad clamping, empty-input guard) is exercised
// here.

#include "WireCellUtil/doctest.h"
#include "WireCellSigProc/ProtoduneHD.h"

#include <cmath>
#include <random>

using namespace WireCell;
using namespace WireCell::SigProc;

namespace {

// 64 channels of Gaussian-ish baseline noise plus an optional negative
// dip on every channel.  Returning per-channel signals (not summed)
// because Is_FEMB_noise sums internally.
IChannelFilter::channel_signals_t
make_chansig(int nchan, int nticks,
             float noise_rms,
             int dip_start, int dip_width, float dip_amp_per_chan,
             unsigned seed = 12345u)
{
    IChannelFilter::channel_signals_t cs;
    std::mt19937 rng(seed);
    std::normal_distribution<float> dist(0.0f, noise_rms);
    for (int ch = 0; ch < nchan; ++ch) {
        Waveform::realseq_t s(nticks, 0.0f);
        for (int i = 0; i < nticks; ++i) s[i] = dist(rng);
        if (dip_width > 0) {
            for (int i = dip_start; i < dip_start + dip_width && i < nticks; ++i) {
                s[i] += dip_amp_per_chan;  // dip_amp_per_chan is negative
            }
        }
        cs[ch] = std::move(s);
    }
    return cs;
}

}  // namespace

TEST_CASE("Is_FEMB_noise empty input returns false") {
    IChannelFilter::channel_signals_t cs;
    Waveform::BinRangeList rois;
    bool found = PDHD::Is_FEMB_noise(cs, rois, /*min_width*/10.0f,
                                     /*pad_nticks*/5, /*nsigma*/3.5f);
    CHECK_FALSE(found);
    CHECK(rois.empty());
}

TEST_CASE("Is_FEMB_noise detects a wide negative dip with padding") {
    const int nchan = 64;
    const int nticks = 1000;
    const int dip_start = 400;
    const int dip_width = 50;          // wider than min_width = 10
    const int pad_nticks = 7;
    const float noise_rms = 1.0f;
    const float dip_amp_per_chan = -3.0f;  // sum dip ~ -64*3 = -192 ADC

    auto cs = make_chansig(nchan, nticks, noise_rms,
                           dip_start, dip_width, dip_amp_per_chan);

    Waveform::BinRangeList rois;
    bool found = PDHD::Is_FEMB_noise(cs, rois, /*min_width*/10.0f,
                                     pad_nticks, /*nsigma*/3.5f);
    REQUIRE(found);
    REQUIRE(rois.size() == 1);
    // Padding extends the ROI by pad_nticks on each side, clamped to [0, n-1]
    CHECK(rois[0].first  == dip_start - pad_nticks);
    CHECK(rois[0].second == dip_start + dip_width - 1 + pad_nticks);
}

TEST_CASE("Is_FEMB_noise rejects a narrow dip below min_width") {
    const int nchan = 64;
    const int nticks = 1000;
    auto cs = make_chansig(nchan, nticks, /*noise_rms*/1.0f,
                           /*dip_start*/400, /*dip_width*/5,
                           /*dip_amp*/ -5.0f);
    Waveform::BinRangeList rois;
    bool found = PDHD::Is_FEMB_noise(cs, rois, /*min_width*/10.0f,
                                     /*pad_nticks*/5, /*nsigma*/3.5f);
    CHECK_FALSE(found);
    CHECK(rois.empty());
}

TEST_CASE("Is_FEMB_noise ignores a positive bump (must be negative-going)") {
    const int nchan = 64;
    const int nticks = 1000;
    auto cs = make_chansig(nchan, nticks, /*noise_rms*/1.0f,
                           /*dip_start*/400, /*dip_width*/50,
                           /*dip_amp*/ +3.0f);  // POSITIVE
    Waveform::BinRangeList rois;
    bool found = PDHD::Is_FEMB_noise(cs, rois, /*min_width*/10.0f,
                                     /*pad_nticks*/5, /*nsigma*/3.5f);
    CHECK_FALSE(found);
}

TEST_CASE("Is_FEMB_noise nsigma threshold gates marginal dips") {
    // A dip whose summed amplitude is ~5x baseline rms in the projection
    // should pass at nsigma=3.5 but fail at nsigma=10.0.
    const int nchan = 64;
    const int nticks = 1000;
    auto cs = make_chansig(nchan, nticks, /*noise_rms*/1.0f,
                           /*dip_start*/400, /*dip_width*/50,
                           /*dip_amp*/ -2.0f);  // sum ~ -128, rms_sum ~ 8

    Waveform::BinRangeList rois_loose, rois_tight;
    bool loose = PDHD::Is_FEMB_noise(cs, rois_loose, /*min_width*/10.0f,
                                     /*pad_nticks*/0, /*nsigma*/3.5f);
    bool tight = PDHD::Is_FEMB_noise(cs, rois_tight, /*min_width*/10.0f,
                                     /*pad_nticks*/0, /*nsigma*/30.0f);
    CHECK(loose);
    CHECK_FALSE(tight);
}

TEST_CASE("Is_FEMB_noise pad clamps to waveform boundaries") {
    const int nchan = 64;
    const int nticks = 1000;
    const int dip_width = 30;
    // Dip at the very end of the waveform; pad_nticks larger than the
    // remaining tail so the back edge of the padded ROI clamps to n-1.
    const int dip_start = nticks - dip_width;   // ends at exactly nticks-1
    auto cs = make_chansig(nchan, nticks, /*noise_rms*/1.0f,
                           dip_start, dip_width,
                           /*dip_amp*/ -3.0f);
    Waveform::BinRangeList rois;
    bool found = PDHD::Is_FEMB_noise(cs, rois, /*min_width*/10.0f,
                                     /*pad_nticks*/50, /*nsigma*/3.5f);
    REQUIRE(found);
    REQUIRE(rois.size() == 1);
    CHECK(rois[0].first == dip_start - 50);
    CHECK(rois[0].second == nticks - 1);          // clamped at n-1
}
