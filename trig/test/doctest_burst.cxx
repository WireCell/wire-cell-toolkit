// Tests of the burst TC maker core.

#include "WireCellTrig/BurstMaker.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <algorithm>
#include <random>

using namespace WireCell;
using namespace WireCell::Trig;

namespace {

    TriggerActivity make_ta(double t, int detector = 0, size_t ninputs = 5, int algorithm = 1)
    {
        TriggerActivity ta;
        ta.time_start = t;
        ta.time_end = t + 10 * units::us;
        ta.time_peak = ta.time_activity = t + 5 * units::us;
        ta.channel_start = 100;
        ta.channel_end = 104;
        ta.adc_integral = 1000;
        ta.adc_peak = 100;
        ta.detector = detector;
        ta.algorithm = algorithm;
        ta.ninputs = ninputs;
        return ta;
    }

    // Uniform background at "rate" over [0, duration), plus nburst TAs in
    // [tburst, tburst + spread), sorted.  Detectors cycle over ndet.
    std::vector<TriggerActivity> make_tas(unsigned seed, double rate, double duration, size_t nburst,
                                          double tburst, double spread, int ndet = 1)
    {
        std::mt19937 rng(seed);
        std::poisson_distribution<int> pois(rate * duration);
        std::uniform_real_distribution<double> uni(0, 1);
        std::vector<double> times;
        const int nbkg = pois(rng);
        for (int i = 0; i < nbkg; ++i) times.push_back(uni(rng) * duration);
        for (size_t i = 0; i < nburst; ++i) times.push_back(tburst + uni(rng) * spread);
        std::sort(times.begin(), times.end());
        std::vector<TriggerActivity> tas;
        for (size_t i = 0; i < times.size(); ++i) tas.push_back(make_ta(times[i], i % ndet));
        return tas;
    }

    std::vector<TriggerCandidate> run(TCMakerBurst& maker, const std::vector<TriggerActivity>& tas)
    {
        std::vector<TriggerCandidate> tcs;
        for (const auto& ta : tas) maker(ta, tcs);
        maker.flush(std::numeric_limits<double>::infinity(), tcs);
        return tcs;
    }
}  // namespace

TEST_CASE("trig burst matches brute force")
{
    const auto tas = make_tas(1, 20 / units::s, 10 * units::s, 0, 0, 0);
    TCMakerBurst::Config cfg;
    cfg.window = 0.5 * units::s;
    cfg.threshold = 12;
    cfg.holdoff = 0;  // a TC at every TA over threshold
    TCMakerBurst maker(cfg);
    auto tcs = run(maker, tas);

    std::vector<double> expected;
    for (size_t i = 0; i < tas.size(); ++i) {
        size_t count = 0;
        for (size_t j = 0; j <= i; ++j) {
            if (tas[j].time_start > tas[i].time_start - cfg.window) ++count;
        }
        if (count >= cfg.threshold) expected.push_back(tas[i].time_start);
    }
    REQUIRE(!expected.empty());  // mean count is 10, so some fluctuate up
    REQUIRE(tcs.size() == expected.size());
    for (size_t i = 0; i < tcs.size(); ++i) {
        CHECK(tcs[i].time_candidate == expected[i]);
        CHECK(tcs[i].ninputs >= cfg.threshold);
        CHECK(tcs[i].inputs.size() == tcs[i].ninputs);
        CHECK(tcs[i].time_start > tcs[i].time_candidate - cfg.window);
        CHECK(tcs[i].time_end >= tcs[i].time_candidate);
    }
}

TEST_CASE("trig burst finds a burst over background")
{
    const double tburst = 12 * units::s;
    const auto tas = make_tas(2, 10 / units::s, 20 * units::s, 50, tburst, 0.5 * units::s, 3);
    TCMakerBurst::Config cfg;
    cfg.window = 1 * units::s;
    cfg.threshold = 30;  // background mean 10 in the window
    cfg.holdoff = 5 * units::s;
    TCMakerBurst maker(cfg);
    auto tcs = run(maker, tas);
    REQUIRE(tcs.size() == 1);
    const auto& tc = tcs[0];
    CHECK(tc.time_candidate > tburst);
    CHECK(tc.time_candidate < tburst + 0.5 * units::s);
    CHECK(tc.ninputs == 30);
    CHECK(tc.type == 3);
    CHECK(tc.detectors == std::vector<int>{0, 1, 2});
    CHECK(maker.nseen() == tas.size());
    CHECK(maker.nselected() == tas.size());

    // with the default holdoff (the window), TCs are at least a window apart
    cfg.holdoff = -1;
    cfg.threshold = 20;
    TCMakerBurst again(cfg);
    auto more = run(again, tas);
    REQUIRE(more.size() >= 1);
    for (size_t i = 1; i < more.size(); ++i) {
        CHECK(more[i].time_candidate >= more[i - 1].time_candidate + cfg.window);
    }
}

TEST_CASE("trig burst selection")
{
    TCMakerBurst::Config cfg;
    cfg.window = 1 * units::s;
    cfg.threshold = 3;
    cfg.select.algorithms = {1};
    cfg.select.min_inputs = 3;
    cfg.select.max_adc_integral = 5000;
    cfg.select.max_channel_span = 10;
    TCMakerBurst maker(cfg);

    std::vector<TriggerActivity> tas;
    tas.push_back(make_ta(1 * units::ms));                // counts
    tas.push_back(make_ta(2 * units::ms, 0, 5, 11));      // wrong algorithm
    tas.push_back(make_ta(3 * units::ms, 0, 2));          // too few TPs
    auto big = make_ta(4 * units::ms);
    big.adc_integral = 10000;                             // too much charge
    tas.push_back(big);
    auto wide = make_ta(5 * units::ms);
    wide.channel_end = wide.channel_start + 50;           // too wide
    tas.push_back(wide);
    tas.push_back(make_ta(6 * units::ms));                // counts
    auto tcs = run(maker, tas);
    CHECK(tcs.empty());
    CHECK(maker.nselected() == 2);

    std::vector<TriggerCandidate> out;
    maker.reset();
    tas.push_back(make_ta(7 * units::ms));                // the third
    for (const auto& ta : tas) maker(ta, out);
    maker.flush(1 * units::s, out);
    REQUIRE(out.size() == 1);
    CHECK(out[0].time_candidate == 7 * units::ms);
    CHECK(out[0].ninputs == 3);
}

TEST_CASE("trig burst ordering and lateness")
{
    TCMakerBurst::Config cfg;
    cfg.window = 1 * units::s;
    cfg.threshold = 4;
    std::vector<TriggerCandidate> tcs;

    // out of order before the flush is the same as in order
    TCMakerBurst maker(cfg);
    for (double t : {3.0, 1.0, 4.0, 2.0}) maker(make_ta(t * 100 * units::ms), tcs);
    CHECK(tcs.empty());  // nothing before a flush
    maker.flush(250 * units::ms, tcs);
    CHECK(tcs.empty());  // only two processed
    maker.flush(1 * units::s, tcs);
    REQUIRE(tcs.size() == 1);
    CHECK(tcs[0].time_candidate == doctest::Approx(400 * units::ms));
    CHECK(tcs[0].time_start == doctest::Approx(100 * units::ms));

    // late: after the promise, in the window counts, before it is dropped
    TCMakerBurst late(cfg);
    tcs.clear();
    for (double t : {1.0, 1.1, 1.2}) late(make_ta(t * units::s), tcs);
    late.flush(2 * units::s, tcs);
    CHECK(tcs.empty());
    late(make_ta(0.1 * units::s), tcs);  // outside (0.2, 1.2]
    CHECK(late.ndropped() == 1);
    CHECK(tcs.empty());
    late(make_ta(0.9 * units::s), tcs);  // inside
    CHECK(late.nlate() == 2);
    REQUIRE(tcs.size() == 1);
    CHECK(tcs[0].time_candidate == doctest::Approx(1.2 * units::s));
    CHECK(tcs[0].time_start == doctest::Approx(0.9 * units::s));
    CHECK(tcs[0].ninputs == 4);
}
