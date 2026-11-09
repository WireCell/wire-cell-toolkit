// Tests of the TP and TA finder nodes, without an anode.

#include "WireCellTrig/TAFinder.h"
#include "WireCellTrig/TPFinder.h"
#include "WireCellTrig/SupernovaMakers.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleTriggerSets.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <random>

using namespace WireCell;
using namespace WireCell::Trig;

namespace {

    const double period = 500 * units::ns;

    // A waveform with baseline 900, a little noise and unipolar pulses.
    std::vector<float> waveform(size_t n, unsigned seed, std::vector<size_t> pulses)
    {
        std::mt19937 rng(seed);
        std::normal_distribution<double> noise(0.0, 2.0);
        std::vector<float> wf(n);
        for (auto& w : wf) w = 900 + std::round(noise(rng));
        for (size_t p : pulses) {
            for (int k = 0; k < 10 && p + k < n; ++k) wf[p + k] += 60 * (k < 5 ? k : 10 - k);
        }
        return wf;
    }

    Configuration tpfinder_config(bool streaming)
    {
        TPFinder probe;
        auto cfg = probe.default_configuration();
        cfg["streaming"] = streaming;
        cfg["tpg"]["threshold"] = 30;
        cfg["tpg"]["pedestal_init"] = "fixed";
        cfg["tpg"]["pedestal"] = 900;
        return cfg;
    }

    IFrame::pointer make_frame(int ident, double time, std::map<int, std::vector<float>> chwf, int tbin = 0)
    {
        ITrace::vector traces;
        for (auto& [ch, wf] : chwf) traces.push_back(std::make_shared<Aux::SimpleTrace>(ch, tbin, wf));
        return std::make_shared<Aux::SimpleFrame>(ident, time, traces, period);
    }

    std::vector<TriggerPrimitive> run_tp(TPFinder& tpf, const std::vector<IFrame::pointer>& frames,
                                         std::vector<ITriggerPrimitiveSet::pointer>& sets)
    {
        std::vector<TriggerPrimitive> all;
        ITriggerPrimitiveFinder::output_queue outq;
        for (auto frame : frames) tpf(frame, outq);
        tpf(nullptr, outq);
        for (auto set : outq) {
            if (!set) continue;
            sets.push_back(set);
            for (const auto& tp : set->primitives()) all.push_back(tp);
        }
        REQUIRE(outq.back() == nullptr);
        return all;
    }
}  // namespace

TEST_CASE("trig finders tp finder per frame")
{
    TPFinder tpf;
    tpf.configure(tpfinder_config(false));
    const double t0 = 1 * units::ms;
    auto frame = make_frame(7, t0, {{100, waveform(1000, 1, {200, 600})}, {101, waveform(1000, 2, {990})}}, 10);
    std::vector<ITriggerPrimitiveSet::pointer> sets;
    auto tps = run_tp(tpf, {frame}, sets);
    REQUIRE(sets.size() == 1);
    CHECK(sets[0]->time_begin() == doctest::Approx(t0 + 10 * period));
    CHECK(sets[0]->time_end() == doctest::Approx(t0 + 1010 * period));
    // two pulses on 100, and the one on 101 cut by the frame end is closed by
    // the per-frame flush
    REQUIRE(tps.size() == 3);
    CHECK(tps[0].channel == 100);
    CHECK(tps[1].channel == 100);
    CHECK(tps[2].channel == 101);
    // pulse samples 60*k: over 30 from k=1 to k=9
    CHECK(tps[0].time_start == doctest::Approx(t0 + (10 + 201) * period));
    CHECK(tps[0].samples_over_threshold == 9);
    CHECK(tps[0].samples_to_peak == 4);
    CHECK(std::abs((int) tps[0].adc_peak - 300) <= 8);  // 5 x 60 plus noise
    CHECK(tps[2].time_start == doctest::Approx(t0 + (10 + 991) * period));
}

TEST_CASE("trig finders tp finder streaming equals whole")
{
    const double t0 = 2 * units::ms;
    const auto wf = waveform(3000, 3, {100, 995, 1500, 2990});

    TPFinder whole;
    whole.configure(tpfinder_config(true));
    std::vector<ITriggerPrimitiveSet::pointer> whole_sets;
    auto whole_tps = run_tp(whole, {make_frame(0, t0, {{5, wf}})}, whole_sets);

    // the same waveform in three contiguous frames, cut inside a pulse
    TPFinder pieces;
    pieces.configure(tpfinder_config(true));
    std::vector<IFrame::pointer> frames;
    const std::vector<size_t> cuts{0, 1000, 2000, 3000};
    for (size_t ic = 0; ic + 1 < cuts.size(); ++ic) {
        std::vector<float> part(wf.begin() + cuts[ic], wf.begin() + cuts[ic + 1]);
        frames.push_back(make_frame(ic, t0 + cuts[ic] * period, {{5, part}}));
    }
    std::vector<ITriggerPrimitiveSet::pointer> piece_sets;
    auto piece_tps = run_tp(pieces, frames, piece_sets);

    REQUIRE(whole_tps.size() == 4);
    REQUIRE(piece_tps.size() == whole_tps.size());
    for (size_t i = 0; i < whole_tps.size(); ++i) {
        CHECK(piece_tps[i].time_start == doctest::Approx(whole_tps[i].time_start));
        CHECK(piece_tps[i].samples_over_threshold == whole_tps[i].samples_over_threshold);
        CHECK(piece_tps[i].adc_integral == whole_tps[i].adc_integral);
    }
    // contiguous coverage, then a last empty-coverage set with the open TP
    REQUIRE(piece_sets.size() == 4);
    CHECK(piece_sets[1]->time_begin() == doctest::Approx(piece_sets[0]->time_end()));
    CHECK(piece_sets[2]->time_begin() == doctest::Approx(piece_sets[1]->time_end()));
    CHECK(piece_sets[3]->time_begin() == piece_sets[3]->time_end());
    CHECK(piece_sets[3]->primitives().size() == 1);
    // the pulse cut at 1000 comes in the set after the one it started in
    CHECK(piece_sets[1]->primitives().front().time_start < piece_sets[1]->time_begin());
}

TEST_CASE("trig finders tp finder rejects another period")
{
    TPFinder tpf;
    tpf.configure(tpfinder_config(false));
    ITrace::vector traces{std::make_shared<Aux::SimpleTrace>(1, 0, std::vector<float>(10, 900))};
    auto frame = std::make_shared<Aux::SimpleFrame>(0, 0, traces, 512 * units::ns);
    ITriggerPrimitiveFinder::output_queue outq;
    CHECK_THROWS(tpf(frame, outq));
}

namespace {
    // Sets of clustered TPs with contiguous coverage of "span" each.
    std::vector<ITriggerPrimitiveSet::pointer> tp_sets(size_t nsets, double span, unsigned seed)
    {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> uni(0, 1);
        std::vector<ITriggerPrimitiveSet::pointer> sets;
        for (size_t is = 0; is < nsets; ++is) {
            const double begin = 1 * units::s + is * span;
            ITriggerPrimitiveSet::primitives_t tps;
            for (int ic = 0; ic < 20; ++ic) {
                const double t = begin + uni(rng) * span * 0.9;
                const int ch = (int) (uni(rng) * 400);
                const int n = 1 + (int) (uni(rng) * 6);
                for (int k = 0; k < n; ++k) {
                    TriggerPrimitive tp;
                    tp.channel = ch + k;
                    tp.time_start = t + k * period;
                    tp.samples_over_threshold = 10;
                    tp.samples_to_peak = 5;
                    tp.adc_peak = 100;
                    tp.adc_integral = 500;
                    tps.push_back(tp);
                }
            }
            std::sort(tps.begin(), tps.end(), [](const auto& a, const auto& b) { return a.time_start < b.time_start; });
            sets.push_back(std::make_shared<Aux::SimpleTriggerPrimitiveSet>(is, 3, begin, begin + span, tps));
        }
        return sets;
    }

    std::vector<ITriggerActivitySet::pointer> run_ta(TAFinder& taf, const std::vector<ITriggerPrimitiveSet::pointer>& in)
    {
        std::vector<ITriggerActivitySet::pointer> sets;
        ITriggerActivityFinder::output_queue outq;
        for (auto set : in) taf(set, outq);
        taf(nullptr, outq);
        REQUIRE(outq.back() == nullptr);
        for (auto set : outq) if (set) sets.push_back(set);
        return sets;
    }
}  // namespace

TEST_CASE("trig finders ta finder supernova matches its core")
{
    const auto in = tp_sets(5, 1 * units::ms, 11);

    TAFinder taf;
    auto cfg = taf.default_configuration();
    cfg["algorithm"] = "Supernova";
    taf.configure(cfg);
    auto sets = run_ta(taf, in);
    REQUIRE(sets.size() == in.size() + 1);
    size_t nfound = 0;
    for (size_t is = 0; is < in.size(); ++is) {
        CHECK(sets[is]->time_begin() == in[is]->time_begin());
        CHECK(sets[is]->time_end() == in[is]->time_end());
        CHECK(sets[is]->detector() == 3);
    }
    for (auto set : sets) nfound += set->activities().size();

    // the core directly, with the same conversion
    DaqClock clock;
    clock.tick = 16 * units::ns;
    TAMakerSupernova::Config sncfg;
    sncfg.ticks_per_sample = period / clock.tick;
    TAMakerSupernova core(sncfg);
    std::vector<Daq::TriggerActivity> direct;
    for (auto set : in)
        for (const auto& tp : set->primitives()) core(to_daq(tp, clock, period), direct);
    core.flush(Daq::invalid_timestamp, direct);
    CHECK(nfound == direct.size());
    CHECK(nfound > 20);
}

TEST_CASE("trig finders ta finder store_inputs and swift windows")
{
    const auto in = tp_sets(6, 1 * units::ms, 12);
    TAFinder taf;
    auto cfg = taf.default_configuration();
    cfg["algorithm"] = "SWIFT";
    cfg["store_inputs"] = false;
    cfg["config"]["inspect_energy_threshold_sadc"] = 0;
    cfg["config"]["accept_energy_threshold_sadc"] = 0;  // accept every non-empty window
    taf.configure(cfg);
    auto sets = run_ta(taf, in);
    REQUIRE(sets.size() == in.size() + 1);

    const double window = 32000 * 16 * units::ns;  // 512 us
    size_t ntas = 0;
    for (size_t is = 0; is < sets.size(); ++is) {
        for (const auto& ta : sets[is]->activities()) {
            ++ntas;
            CHECK(ta.inputs.empty());
            CHECK(ta.ninputs > 0);
            CHECK(ta.time_end - ta.time_start == doctest::Approx(window));
            if (is < in.size()) {
                // closed by the watermark of the set it is emitted in
                CHECK(ta.time_end <= sets[is]->time_end() + 1e-6);
            }
        }
    }
    // windows of 512 us tile the 6 ms: about 12, all with TPs
    CHECK(ntas >= 11);
}

TEST_CASE("trig finders ta finder not streaming")
{
    // the same coverage again and again, as independent frames
    auto in = tp_sets(1, 1 * units::ms, 13);
    std::vector<ITriggerPrimitiveSet::pointer> frames;
    for (int ind = 0; ind < 3; ++ind) {
        frames.push_back(std::make_shared<Aux::SimpleTriggerPrimitiveSet>(
            ind, 3, in[0]->time_begin(), in[0]->time_end(), in[0]->primitives()));
    }
    TAFinder taf;
    auto cfg = taf.default_configuration();
    cfg["streaming"] = false;
    taf.configure(cfg);
    auto sets = run_ta(taf, frames);
    REQUIRE(sets.size() == frames.size());  // no end of stream set
    for (const auto& set : sets) {
        // everything of a frame is flushed with it, and each is the same
        CHECK(set->activities().size() == sets[0]->activities().size());
        CHECK(set->activities().size() > 0);
    }
}
