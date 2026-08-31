#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <algorithm>
#include <complex>

using std::vector;
using namespace WireCell;

TEST_SUITE("waveform") {

TEST_CASE("transform") {
    const int baseline = 1000;
    const int nticks = 100;
    const int period = 10;

    Waveform::realseq_t wf1(nticks), wf2(nticks);

    for (int ind = 0; ind < nticks; ++ind) {
        wf1[ind] = baseline + ind % period;
    }

    transform(wf1.begin(), wf1.end(), wf2.begin(), [](int x) -> int { return x - baseline; });

    for (int ind = 0; ind < nticks; ++ind) {
        CHECK(wf1[ind] - baseline == wf2[ind]);
    }
}

TEST_CASE("mean rms median") {
    Waveform::realseq_t v{1.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 3.0};

    auto us = Waveform::mean_rms(v);
    auto m = Waveform::median(v);

    spdlog::debug("mean={} rms={} median={}", us.first, us.second, m);

    CHECK(us.first == doctest::Approx(2.75));
    CHECK(us.second == doctest::Approx(1.199).epsilon(0.01));
    CHECK(m == doctest::Approx(3.0));
}

TEST_CASE("arithmetic real") {
    using namespace WireCell::Waveform;
    realseq_t v{0.0, 1.0, 2.0};
    auto v2 = v;

    increase(v, 2.0);
    CHECK(v[0] == 2.0);
    CHECK(v[1] == 3.0);
    CHECK(v[2] == 4.0);

    increase(v2, v);
    CHECK(v2[0] == 2.0);
    CHECK(v2[1] == 4.0);
    CHECK(v2[2] == 6.0);

    scale(v, 2.0);
    scale(v2, v);
}

TEST_CASE("arithmetic complex") {
    using namespace WireCell::Waveform;
    compseq_t cv{{1.1, 2.2}, {-3.3, 4.4}, {0, 0}, {1, 0}, {0, 1}, {-1, 0}, {0, -1}};
    auto cv2 = cv;

    increase(cv, complex_t(1.0, 0.0));
    increase(cv2, cv);
    scale(cv, complex_t(1.0, 0.0));
    scale(cv2, cv);

    CHECK(true);  // TODO: add meaningful assertions for complex arithmetic results
}

// Regression test for Waveform::merge(BinRangeList) — sorting input by
// .first does not order .second, so an earlier-but-wider range followed by
// a narrower overlapping one previously shrank the merged range.  The fix
// in util/src/Waveform.cxx takes max(.second).
TEST_CASE("merge binranges") {
    using BR = Waveform::BinRange;

    // Wide range followed by a narrower one whose .first lies inside it but
    // whose .second is *less* than the wide range's .second.  Expected: a
    // single merged range equal to the wide one.
    Waveform::BinRangeList in1{ BR{0, 100}, BR{10, 50} };
    auto out1 = Waveform::merge(in1);
    CHECK(out1.size() == 1);
    CHECK(out1[0].first  == 0);
    CHECK(out1[0].second == 100);

    // Two genuinely overlapping ranges where the second extends past the
    // first.  Expected: union.
    Waveform::BinRangeList in2{ BR{0, 50}, BR{30, 80} };
    auto out2 = Waveform::merge(in2);
    CHECK(out2.size() == 1);
    CHECK(out2[0].first  == 0);
    CHECK(out2[0].second == 80);

    // Disjoint ranges stay separate.
    Waveform::BinRangeList in3{ BR{0, 10}, BR{20, 30} };
    auto out3 = Waveform::merge(in3);
    CHECK(out3.size() == 2);
    CHECK(out3[0].first  == 0);
    CHECK(out3[0].second == 10);
    CHECK(out3[1].first  == 20);
    CHECK(out3[1].second == 30);

    // Unsorted input with the wide-then-narrow pattern hidden by ordering.
    Waveform::BinRangeList in4{ BR{10, 50}, BR{0, 100} };
    auto out4 = Waveform::merge(in4);
    CHECK(out4.size() == 1);
    CHECK(out4[0].first  == 0);
    CHECK(out4[0].second == 100);
}

}  // TEST_SUITE("waveform")
