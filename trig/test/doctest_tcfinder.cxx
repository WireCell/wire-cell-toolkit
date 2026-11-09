// Tests of the TC finder and TA set fanin nodes.

#include "WireCellTrig/DaqConvert.h"
#include "WireCellTrig/SupernovaMakers.h"
#include "WireCellTrig/TASetFanin.h"
#include "WireCellTrig/TCFinder.h"

#include "WireCellAux/SimpleTriggerSets.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <random>

using namespace WireCell;
using namespace WireCell::Trig;

namespace {

    // nsets TA sets of contiguous coverage "span" from t0, with background
    // TAs at "rate" and nburst more in the set "iburst".
    std::vector<ITriggerActivitySet::pointer> ta_sets(unsigned seed, int detector, size_t nsets, double span,
                                                      double rate, size_t iburst = 0, size_t nburst = 0)
    {
        const double t0 = 1 * units::s;
        std::mt19937 rng(seed);
        std::poisson_distribution<int> pois(rate * span);
        std::uniform_real_distribution<double> uni(0, 1);
        std::vector<ITriggerActivitySet::pointer> sets;
        for (size_t is = 0; is < nsets; ++is) {
            const double begin = t0 + is * span;
            size_t n = pois(rng);
            if (is == iburst) n += nburst;
            ITriggerActivitySet::activities_t tas;
            for (size_t i = 0; i < n; ++i) {
                TriggerActivity ta;
                ta.time_start = begin + uni(rng) * span;
                ta.time_end = ta.time_start + 20 * units::us;
                ta.time_peak = ta.time_activity = ta.time_start + 8 * units::us;
                ta.channel_start = 10 * i;
                ta.channel_end = ta.channel_peak = 10 * i + 3;
                ta.adc_integral = 2000;
                ta.adc_peak = 200;
                ta.detector = detector;
                ta.algorithm = 1;
                ta.ninputs = 4;
                ta.inputs.resize(4);
                for (auto& tp : ta.inputs) {
                    tp.channel = ta.channel_start;
                    tp.time_start = ta.time_start;
                }
                tas.push_back(ta);
            }
            std::sort(tas.begin(), tas.end(), [](const auto& a, const auto& b) { return a.time_start < b.time_start; });
            sets.push_back(std::make_shared<Aux::SimpleTriggerActivitySet>(is, detector, begin, begin + span, tas));
        }
        return sets;
    }

    std::vector<ITriggerCandidateSet::pointer> run_tc(TCFinder& tcf, const std::vector<ITriggerActivitySet::pointer>& in)
    {
        ITriggerCandidateFinder::output_queue outq;
        for (auto set : in) tcf(set, outq);
        tcf(nullptr, outq);
        REQUIRE(outq.back() == nullptr);
        std::vector<ITriggerCandidateSet::pointer> sets;
        for (auto set : outq)
            if (set) sets.push_back(set);
        return sets;
    }
}  // namespace

TEST_CASE("trig tc finder burst")
{
    const double span = 100 * units::ms;
    const auto in = ta_sets(1, 4, 50, span, 20 / units::s, 30, 40);

    TCFinder tcf;
    auto cfg = tcf.default_configuration();
    CHECK(cfg["algorithm"].asString() == "Burst");
    cfg["store_inputs"] = false;
    cfg["config"]["window"] = 1 * units::s;  // background mean 20
    cfg["config"]["threshold"] = 45;
    cfg["config"]["holdoff"] = 10 * units::s;
    cfg["config"]["min_inputs"] = 3;
    tcf.configure(cfg);
    auto sets = run_tc(tcf, in);

    REQUIRE(sets.size() == in.size() + 1);
    std::vector<std::pair<size_t, TriggerCandidate>> found;
    for (size_t is = 0; is < sets.size(); ++is) {
        if (is < in.size()) {
            CHECK(sets[is]->time_begin() == in[is]->time_begin());
            CHECK(sets[is]->time_end() == in[is]->time_end());
        }
        CHECK(sets[is]->detector() == 4);
        for (const auto& tc : sets[is]->candidates()) found.emplace_back(is, tc);
    }
    REQUIRE(found.size() == 1);
    const auto& [iset, tc] = found[0];
    // with no lateness, it is made in the set where the threshold is crossed
    CHECK(iset == 30);
    CHECK(tc.time_candidate >= in[30]->time_begin());
    CHECK(tc.time_candidate < in[30]->time_end());
    CHECK(tc.ninputs == 45);
    CHECK(tc.inputs.empty());
    CHECK(tc.detectors == std::vector<int>{4});
    CHECK(tc.type == 3);
}

TEST_CASE("trig tc finder burst lateness delays")
{
    const double span = 100 * units::ms;
    const auto in = ta_sets(2, 0, 20, span, 0, 10, 20);
    TCFinder tcf;
    auto cfg = tcf.default_configuration();
    cfg["lateness"] = 2.5 * span;
    cfg["config"]["threshold"] = 20;
    tcf.configure(cfg);
    auto sets = run_tc(tcf, in);
    size_t iset = sets.size();
    for (size_t is = 0; is < sets.size(); ++is)
        if (!sets[is]->candidates().empty()) iset = is;
    CHECK(iset == 13);  // the TAs of set 10 are processed once 2.5 sets later
    REQUIRE(iset < sets.size());
    CHECK(sets[iset]->candidates()[0].inputs.size() == 20);
    CHECK(sets[iset]->candidates()[0].inputs[0].inputs.size() == 4);
}

TEST_CASE("trig tc finder supernova matches its core")
{
    const auto in = ta_sets(3, 2, 10, 10 * units::ms, 2000 / units::s);
    TCFinder tcf;
    auto cfg = tcf.default_configuration();
    cfg["algorithm"] = "Supernova";
    cfg["config"]["time_window"] = 62500;  // 1 ms
    cfg["config"]["threshold"] = 3;
    cfg["config"]["hit_threshold"] = 2;
    tcf.configure(cfg);
    auto sets = run_tc(tcf, in);

    size_t nfound = 0;
    for (auto set : sets) {
        for (const auto& tc : set->candidates()) {
            ++nfound;
            CHECK(tc.ninputs == 4);
            CHECK(tc.detectors == std::vector<int>{2});
            CHECK(tc.type == 3);
            CHECK(tc.algorithm == 1);
        }
    }

    DaqClock clock;
    clock.tick = 16 * units::ns;
    TCMakerSupernova::Config sncfg;
    sncfg.time_window = 62500;
    TCMakerSupernova core(sncfg);
    std::vector<Daq::TriggerCandidate> direct;
    for (auto set : in)
        for (const auto& ta : set->activities()) core(to_daq(ta, clock, 500 * units::ns), direct);
    CHECK(nfound == direct.size());
    CHECK(nfound > 10);
}

TEST_CASE("trig tc finder rejects bad configuration")
{
    TCFinder tcf;
    auto cfg = tcf.default_configuration();
    cfg["algorithm"] = "Nope";
    CHECK_THROWS(tcf.configure(cfg));
    cfg["algorithm"] = "Burst";
    cfg["config"]["threshold"] = 0;
    CHECK_THROWS(tcf.configure(cfg));
}

TEST_CASE("trig ta set fanin")
{
    const double span = 50 * units::ms;
    const auto a = ta_sets(4, 1, 3, span, 200 / units::s);
    const auto b = ta_sets(5, 2, 3, span, 200 / units::s);
    TASetFanin fanin;
    auto cfg = fanin.default_configuration();
    cfg["multiplicity"] = 2;
    fanin.configure(cfg);
    CHECK(fanin.input_types().size() == 2);

    for (size_t is = 0; is < a.size(); ++is) {
        ITriggerActivitySetFanin::output_pointer out;
        REQUIRE(fanin({a[is], b[is]}, out));
        REQUIRE(out);
        CHECK(out->detector() == -1);
        CHECK(out->time_begin() == a[is]->time_begin());
        CHECK(out->time_end() == a[is]->time_end());
        const auto& tas = out->activities();
        CHECK(tas.size() == a[is]->activities().size() + b[is]->activities().size());
        CHECK(std::is_sorted(tas.begin(), tas.end(),
                             [](const auto& x, const auto& y) { return x.time_start < y.time_start; }));
    }
    // same detector stays
    ITriggerActivitySetFanin::output_pointer out;
    fanin({a[0], a[0]}, out);
    CHECK(out->detector() == 1);
    // mismatched coverage
    CHECK_THROWS(fanin({a[0], b[1]}, out));
    // end of stream
    fanin({nullptr, nullptr}, out);
    CHECK(out == nullptr);
}
