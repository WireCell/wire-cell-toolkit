// Tests of the TA/TC maker ports.  The references transcribe DUNE-DAQ
// triggeralgs develop (450867e): TAMakerSupernovaAlgorithm,
// TCMakerSupernovaAlgorithm and TriggerActivityMakerSWIFT, with TLOG and
// registration removed.

#include "WireCellTrig/SWIFTMaker.h"
#include "WireCellTrig/SupernovaMakers.h"

#include "WireCellUtil/doctest.h"

#include <algorithm>
#include <cstdlib>
#include <random>
#include <vector>

using namespace WireCell::Trig;
using namespace WireCell::Trig::Daq;

namespace {

    // A time-ordered stream of TPs: noise-like singles plus compact clusters.
    std::vector<TriggerPrimitive> make_tps(size_t nclusters, unsigned seed, timestamp_t t0 = 1'000'000)
    {
        std::mt19937 rng(seed);
        std::uniform_int_distribution<int> chan(0, 479);
        std::uniform_int_distribution<int> gap(0, 20000);
        std::uniform_int_distribution<int> size(1, 8);
        std::uniform_int_distribution<int> sot(1, 30);
        std::uniform_int_distribution<int> peak(20, 400);
        std::vector<TriggerPrimitive> tps;
        timestamp_t t = t0;
        for (size_t ic = 0; ic < nclusters; ++ic) {
            t += gap(rng);
            const int c0 = chan(rng);
            const int n = size(rng);
            for (int k = 0; k < n; ++k) {
                TriggerPrimitive tp;
                tp.detid = 3;
                tp.channel = std::min(479, c0 + k);
                tp.time_start = t + 40 * k;
                tp.samples_over_threshold = sot(rng);
                tp.samples_to_peak = tp.samples_over_threshold / 2;
                tp.adc_peak = peak(rng);
                tp.adc_integral = tp.adc_peak * tp.samples_over_threshold / 2;
                tps.push_back(tp);
            }
        }
        std::stable_sort(tps.begin(), tps.end(),
                         [](const auto& a, const auto& b) { return a.time_start < b.time_start; });
        return tps;
    }

    // ---- transcription of TAMakerSupernovaAlgorithm ----
    struct RefSupernovaTA {
        timestamp_diff_t m_time_tolerance = 250;
        channel_t m_channel_tolerance = 2;
        timestamp_t m_time_start = 0, m_time_end = 0, m_time_peak = 0, m_time_activity = 0;
        channel_diff_t m_channel_start = 0, m_channel_end = 0, m_channel_peak = 0;
        uint64_t m_adc_integral = 0;
        uint16_t m_adc_peak = 0;
        detid_t m_detid = 0;
        std::vector<TriggerPrimitive> m_tp_list;

        bool is_time_consistent(const TriggerPrimitive& input_tp) const
        {
            timestamp_t tend = input_tp.time_start + input_tp.samples_over_threshold * 32;
            bool is_close_to_edge =
                (m_time_tolerance > std::abs(timestamp_diff_t(input_tp.time_start) - timestamp_diff_t(m_time_end))) ||
                m_time_tolerance > std::abs(timestamp_diff_t(input_tp.time_start) - timestamp_diff_t(m_time_start)) ||
                m_time_tolerance > std::abs(timestamp_diff_t(tend) - timestamp_diff_t(m_time_end)) ||
                m_time_tolerance > std::abs(timestamp_diff_t(tend) - timestamp_diff_t(m_time_start));
            bool is_in_between_edge = ((tend > m_time_start && tend < m_time_end) ||
                                       (input_tp.time_start > m_time_start && input_tp.time_start < m_time_end));
            return is_in_between_edge || is_close_to_edge;
        }
        bool is_channel_consistent(const TriggerPrimitive& input_tp) const
        {
            bool is_close_to_edge =
                ((channel_diff_t) m_channel_tolerance > std::abs(channel_diff_t(input_tp.channel) - m_channel_end) ||
                 (channel_diff_t) m_channel_tolerance > std::abs(channel_diff_t(input_tp.channel) - m_channel_start));
            bool is_in_between_edge = ((channel_diff_t) input_tp.channel > m_channel_start &&
                                       (channel_diff_t) input_tp.channel < m_channel_end);
            return is_in_between_edge || is_close_to_edge;
        }
        TriggerActivity make() const
        {
            TriggerActivity ta;
            ta.time_start = m_time_start;
            ta.time_end = m_time_end;
            ta.time_peak = m_time_peak;
            ta.time_activity = m_time_activity;
            ta.channel_start = m_channel_start;
            ta.channel_end = m_channel_end;
            ta.channel_peak = m_channel_peak;
            ta.adc_integral = m_adc_integral;
            ta.adc_peak = m_adc_peak;
            ta.detid = m_detid;
            ta.type = TriggerActivityData::Type::kTPC;
            ta.algorithm = TriggerActivityData::Algorithm::kSupernova;
            ta.inputs = m_tp_list;
            return ta;
        }
        void start(const TriggerPrimitive& input_tp, timestamp_t tend)
        {
            m_tp_list.erase(m_tp_list.begin(), m_tp_list.end());
            m_tp_list.push_back(input_tp);
            m_time_start = input_tp.time_start;
            m_time_end = tend;
            m_time_peak = input_tp.samples_to_peak * 32 + input_tp.time_start;
            m_channel_start = input_tp.channel;
            m_channel_end = input_tp.channel;
            m_channel_peak = input_tp.channel;
            m_adc_integral = input_tp.adc_integral;
            m_adc_peak = input_tp.adc_peak;
            m_detid = input_tp.detid;
        }
        void process(const TriggerPrimitive& input_tp, std::vector<TriggerActivity>& output_ta)
        {
            timestamp_t tend = input_tp.time_start + input_tp.samples_over_threshold * 32;
            if (m_time_start == 0) {
                start(input_tp, tend);
                return;
            }
            bool time_ok = is_time_consistent(input_tp);
            bool channel_ok = is_channel_consistent(input_tp);
            if (!time_ok && !channel_ok) {
                output_ta.push_back(make());
                start(input_tp, tend);
                return;
            }
            if (input_tp.time_start < m_time_start) m_time_start = input_tp.time_start;
            if (tend > m_time_end) m_time_end = tend;
            if (input_tp.adc_peak > m_adc_peak) {
                m_time_peak = input_tp.samples_to_peak * 32 + input_tp.time_start;
                m_adc_peak = input_tp.adc_peak;
                m_channel_peak = input_tp.channel;
            }
            if ((channel_diff_t) input_tp.channel > m_channel_end) m_channel_end = input_tp.channel;
            if ((channel_diff_t) input_tp.channel < m_channel_start) m_channel_start = input_tp.channel;
            m_tp_list.push_back(input_tp);
            m_adc_integral += input_tp.adc_integral;
            m_detid |= input_tp.detid;
        }
    };

    // ---- transcription of TCMakerSupernovaAlgorithm ----
    struct RefSupernovaTC {
        std::vector<TriggerActivityData> m_activity;
        int64_t m_time_window = 500'000'000;
        uint16_t m_threshold = 3;
        uint16_t m_hit_threshold = 2;
        void process(const TriggerActivity& activity, std::vector<TriggerCandidate>& cand)
        {
            timestamp_t time = activity.time_start;
            timestamp_diff_t how_far = time - m_time_window;
            auto end = std::remove_if(m_activity.begin(), m_activity.end(), [how_far](auto& c) -> bool {
                return (static_cast<timestamp_diff_t>(c.time_start) < how_far);
            });
            m_activity.erase(end, m_activity.end());
            if (activity.inputs.size() > m_hit_threshold)
                m_activity.push_back(static_cast<const TriggerActivityData&>(activity));
            if (m_activity.size() > m_threshold) {
                TriggerCandidate tc;
                tc.time_start = time - 500'000'000;
                tc.time_end = activity.time_end;
                tc.time_candidate = time;
                tc.detid = whole_detector;
                tc.type = TriggerCandidateData::Type::kSupernova;
                tc.algorithm = TriggerCandidateData::Algorithm::kSupernova;
                tc.inputs = m_activity;
                m_activity.clear();
                cand.push_back(tc);
            }
        }
    };

    void check_same_ta(const TriggerActivity& a, const TriggerActivity& b)
    {
        CHECK(a.time_start == b.time_start);
        CHECK(a.time_end == b.time_end);
        CHECK(a.time_peak == b.time_peak);
        CHECK(a.time_activity == b.time_activity);
        CHECK(a.channel_start == b.channel_start);
        CHECK(a.channel_end == b.channel_end);
        CHECK(a.channel_peak == b.channel_peak);
        CHECK(a.adc_integral == b.adc_integral);
        CHECK(a.adc_peak == b.adc_peak);
        CHECK(a.detid == b.detid);
        CHECK(a.algorithm == b.algorithm);
        CHECK(a.inputs.size() == b.inputs.size());
    }
}  // namespace

TEST_CASE("trig makers supernova TA and TC match triggeralgs")
{
    // TC time_start would wrap in triggeralgs before 5e8 ticks; start later.
    const auto tps = make_tps(3000, 42, 600'000'000);
    RefSupernovaTA ref;
    TAMakerSupernova port;
    std::vector<TriggerActivity> ref_tas, got_tas;
    for (const auto& tp : tps) {
        ref.process(tp, ref_tas);
        port(tp, got_tas);
    }
    REQUIRE(ref_tas.size() > 100);
    REQUIRE(got_tas.size() == ref_tas.size());
    for (size_t i = 0; i < got_tas.size(); ++i) check_same_ta(got_tas[i], ref_tas[i]);

    RefSupernovaTC ref_tc;
    TCMaker::Config tccfg;
    tccfg.tc_type = TriggerCandidateData::Type::kSupernova;
    TCMakerSupernova port_tc;
    port_tc.set_common(tccfg);
    std::vector<TriggerCandidate> ref_tcs, got_tcs;
    for (const auto& ta : ref_tas) {
        ref_tc.process(ta, ref_tcs);
        port_tc(ta, got_tcs);
    }
    REQUIRE(ref_tcs.size() > 0);
    REQUIRE(got_tcs.size() == ref_tcs.size());
    for (size_t i = 0; i < got_tcs.size(); ++i) {
        CHECK(got_tcs[i].time_start == ref_tcs[i].time_start);
        CHECK(got_tcs[i].time_end == ref_tcs[i].time_end);
        CHECK(got_tcs[i].time_candidate == ref_tcs[i].time_candidate);
        CHECK(got_tcs[i].detid == ref_tcs[i].detid);
        CHECK(got_tcs[i].type == ref_tcs[i].type);
        CHECK(got_tcs[i].inputs.size() == ref_tcs[i].inputs.size());
    }
}

TEST_CASE("trig makers supernova deviations")
{
    TriggerPrimitive a, b;
    a.channel = 10;
    a.time_start = 0;  // triggeralgs would take this as "no cluster yet"
    a.samples_over_threshold = 4;
    a.samples_to_peak = 2;
    b = a;
    b.channel = 300;  // far in channel
    b.time_start = 100;  // near in time

    // Either-or (DAQ) merges b, both=true splits.
    {
        TAMakerSupernova either;
        std::vector<TriggerActivity> tas;
        either(a, tas);
        either(b, tas);
        CHECK(tas.empty());
        either.flush(invalid_timestamp, tas);
        REQUIRE(tas.size() == 1);
        CHECK(tas[0].inputs.size() == 2);
        // flush reset: a second flush gives nothing
        either.flush(invalid_timestamp, tas);
        CHECK(tas.size() == 1);
    }
    {
        TAMakerSupernova::Config cfg;
        cfg.both = true;
        TAMakerSupernova both(cfg);
        std::vector<TriggerActivity> tas;
        both(a, tas);
        both(b, tas);
        REQUIRE(tas.size() == 1);
        CHECK(tas[0].time_start == 0);
        CHECK(tas[0].inputs.size() == 1);
    }
    {
        // flush(until) closes a both-mode cluster only once no TP can join.
        TAMakerSupernova::Config cfg;
        cfg.both = true;
        TAMakerSupernova both(cfg);
        TriggerPrimitive c = a;
        c.time_start = 1000;  // time_end = 1000 + 4*32 = 1128
        std::vector<TriggerActivity> tas;
        both(c, tas);
        both.flush(1128 + 250 - 1, tas);
        CHECK(tas.empty());
        both.flush(1128 + 250, tas);
        CHECK(tas.size() == 1);
        // either-or mode only closes at end of stream
        TAMakerSupernova either;
        tas.clear();
        either(c, tas);
        either.flush(1'000'000'000, tas);
        CHECK(tas.empty());
        either.flush(invalid_timestamp, tas);
        CHECK(tas.size() == 1);
    }
    {
        // Flushing a fresh maker gives no empty TA.
        TAMakerSupernova fresh;
        std::vector<TriggerActivity> tas;
        fresh.flush(invalid_timestamp, tas);
        CHECK(tas.empty());
    }
    {
        // The TC start is clamped at zero.
        TCMakerSupernova::Config cfg;
        cfg.threshold = 0;
        cfg.hit_threshold = 0;
        TCMakerSupernova tcm(cfg);
        TriggerActivity ta;
        ta.time_start = 1000;
        ta.time_end = 2000;
        ta.inputs.resize(1);
        std::vector<TriggerCandidate> tcs;
        tcm(ta, tcs);
        REQUIRE(tcs.size() == 1);
        CHECK(tcs[0].time_start == 0);
    }
}

TEST_CASE("trig makers base prescale and store_inputs")
{
    const auto tps = make_tps(500, 7, 600'000'000);
    TAMakerSupernova all;
    std::vector<TriggerActivity> all_tas;
    for (const auto& tp : tps) all(tp, all_tas);
    REQUIRE(all_tas.size() > 30);

    TAMaker::Config cfg;
    cfg.prescale = 3;
    cfg.store_inputs = false;
    TAMakerSupernova some;
    some.set_common(cfg);
    std::vector<TriggerActivity> some_tas;  // accumulated over calls on purpose
    for (const auto& tp : tps) some(tp, some_tas);
    REQUIRE(some_tas.size() == all_tas.size() / 3);
    for (size_t i = 0; i < some_tas.size(); ++i) {
        CHECK(some_tas[i].time_start == all_tas[3 * i + 2].time_start);
        CHECK(some_tas[i].inputs.empty());
    }

    // max_samples_over_threshold drops long TPs before the algorithm.
    TAMaker::Config cut;
    cut.max_samples_over_threshold = 10;
    TAMakerSupernova cutter;
    cutter.set_common(cut);
    std::vector<TriggerActivity> cut_tas;
    for (const auto& tp : tps) cutter(tp, cut_tas);
    cutter.flush(invalid_timestamp, cut_tas);
    for (const auto& ta : cut_tas) {
        for (const auto& tp : ta.inputs) CHECK(tp.samples_over_threshold <= 10);
    }
}

namespace {
    // ---- transcription of TriggerActivityMakerSWIFT (develop) ----
    struct RefSWIFT {
        uint64_t m_window_length = 32000, m_inspect = 15000, m_accept = 55000;
        uint16_t m_min_adc_peak = 80, m_min_sot = 8;
        float m_cm_per_tick = 0.016 * 0.16, m_wire_pitch = 0.48;
        int m_db_min_samples = 2;
        float m_db_eps = 2;
        uint64_t m_cut = 22000;
        bool m_initialised = false;
        uint64_t m_window_start = 0, m_window_energy = 0;
        TriggerActivity m_current_ta;

        bool preprocess(const TriggerPrimitive& tp)
        {
            if ((tp.adc_peak < m_min_adc_peak) && (tp.samples_over_threshold < m_min_sot)) return false;
            return true;
        }
        void reset_window_state(uint64_t s)
        {
            m_window_start = s;
            m_window_energy = 0;
            m_current_ta = TriggerActivity();
            m_current_ta.time_start = m_window_start;
        }
        void process(const TriggerPrimitive& tp, std::vector<TriggerActivity>& out)
        {
            if (!preprocess(tp)) return;
            const uint64_t ws = (tp.time_start / m_window_length) * m_window_length;
            if (!m_initialised) {
                m_current_ta = TriggerActivity();
                reset_window_state(ws);
                m_initialised = true;
            }
            if (ws < m_window_start) return;
            if (ws > m_window_start) {
                close_window(out);
                reset_window_state(ws);
            }
            m_current_ta.inputs.push_back(tp);
            m_window_energy += tp.adc_integral;
        }
        void close_window(std::vector<TriggerActivity>& out)
        {
            if (m_current_ta.inputs.empty()) return;
            int decision;
            if (m_window_energy >= m_accept) decision = 2;
            else if (m_window_energy >= m_inspect) decision = 1;
            else return;
            if (decision == 1) {
                if (cluster(m_current_ta.inputs, m_db_eps, m_db_min_samples) <= m_cut) return;
            }
            set_ta_attributes();
            out.push_back(m_current_ta);
        }
        uint64_t cluster(const std::vector<TriggerPrimitive>& tps, float eps, int min_samples)
        {
            struct Point { float z; float t; uint64_t adc; };
            std::vector<Point> points;
            for (const auto& tp : tps)
                points.push_back({tp.channel * m_wire_pitch, (tp.time_start - m_window_start) * m_cm_per_tick, tp.adc_integral});
            const int N = points.size();
            int cluster_id = 0;
            std::vector<int> labels(N, -1);
            std::vector<uint8_t> visited(N, 0);
            float eps2 = eps * eps;
            for (int i = 0; i < N; ++i) {
                if (visited[i]) continue;
                visited[i] = 1;
                std::vector<int> neigh;
                for (int j = 0; j < N; ++j) {
                    float dz = points[j].z - points[i].z, dt = points[j].t - points[i].t;
                    if (dz * dz + dt * dt <= eps2) neigh.push_back(j);
                }
                if (neigh.size() < static_cast<size_t>(min_samples)) { labels[i] = -1; continue; }
                labels[i] = cluster_id;
                size_t k = 0;
                while (k < neigh.size()) {
                    int j = neigh[k];
                    if (!visited[j]) {
                        visited[j] = 1;
                        std::vector<int> neigh2;
                        for (int m = 0; m < N; ++m) {
                            float dz = points[m].z - points[j].z, dt = points[m].t - points[j].t;
                            if (dz * dz + dt * dt <= eps2) neigh2.push_back(m);
                        }
                        if (neigh2.size() >= static_cast<size_t>(min_samples))
                            for (int m : neigh2)
                                if (std::find(neigh.begin(), neigh.end(), m) == neigh.end()) neigh.push_back(m);
                    }
                    if (labels[j] == -1) labels[j] = cluster_id;
                    ++k;
                }
                ++cluster_id;
            }
            std::vector<uint64_t> sums(cluster_id, 0);
            for (int i = 0; i < N; ++i) if (labels[i] >= 0) sums[labels[i]] += points[i].adc;
            if (sums.empty()) return 0;
            return *std::max_element(sums.begin(), sums.end());
        }
        void set_ta_attributes()
        {
            m_current_ta.time_start = m_window_start;
            m_current_ta.time_end = m_window_start + m_window_length;
            m_current_ta.adc_integral = m_window_energy;
            const auto& first_tp = m_current_ta.inputs.front();
            m_current_ta.detid = first_tp.detid;
            m_current_ta.type = TriggerActivityData::Type::kTPC;
            m_current_ta.algorithm = TriggerActivityData::Algorithm::kSWIFT;
            channel_t min_ch = first_tp.channel, max_ch = first_tp.channel;
            m_current_ta.adc_peak = 0;
            for (const auto& tp : m_current_ta.inputs) {
                if (tp.channel < min_ch) min_ch = tp.channel;
                if (tp.channel > max_ch) max_ch = tp.channel;
                if (tp.adc_peak > m_current_ta.adc_peak) {
                    m_current_ta.adc_peak = tp.adc_peak;
                    m_current_ta.channel_peak = tp.channel;
                    m_current_ta.time_peak = tp.time_start;
                }
            }
            m_current_ta.channel_start = min_ch;
            m_current_ta.channel_end = max_ch;
            m_current_ta.time_activity = m_current_ta.time_peak;
        }
    };
}  // namespace

TEST_CASE("trig makers swift matches triggeralgs")
{
    // Lower thresholds so all three window decisions occur.
    auto tps = make_tps(4000, 2024);
    RefSWIFT ref;
    ref.m_inspect = 2000;
    ref.m_accept = 6000;
    ref.m_cut = 1500;
    TAMakerSWIFT::Config cfg;
    cfg.inspect_energy_threshold = 2000;
    cfg.accept_energy_threshold = 6000;
    cfg.cluster_energy_cut = 1500;
    TAMakerSWIFT port(cfg);

    std::vector<TriggerActivity> ref_tas, got_tas;
    for (const auto& tp : tps) {
        ref.process(tp, ref_tas);
        port(tp, got_tas);
    }
    REQUIRE(ref_tas.size() > 50);
    REQUIRE(got_tas.size() == ref_tas.size());
    for (size_t i = 0; i < got_tas.size(); ++i) check_same_ta(got_tas[i], ref_tas[i]);

    // Some accepted windows went through inspection.
    size_t inspected = 0;
    for (const auto& ta : got_tas) {
        if (ta.adc_integral < 6000) ++inspected;
    }
    CHECK(inspected > 0);
}

TEST_CASE("trig makers swift flush semantics")
{
    TAMakerSWIFT::Config cfg;
    cfg.inspect_energy_threshold = 0;
    cfg.accept_energy_threshold = 0;  // accept every non-empty window
    TAMakerSWIFT swift(cfg);
    TriggerPrimitive tp;
    tp.channel = 5;
    tp.time_start = 32000 * 10 + 100;  // in window [320000, 352000)
    tp.samples_over_threshold = 20;
    tp.adc_peak = 100;
    tp.adc_integral = 1000;
    std::vector<TriggerActivity> tas;
    swift(tp, tas);
    swift.flush(352000 - 1, tas);  // window not yet complete
    CHECK(tas.empty());
    swift.flush(352000, tas);
    REQUIRE(tas.size() == 1);
    CHECK(tas[0].time_start == 320000);
    CHECK(tas[0].time_end == 352000);
    swift.flush(10'000'000, tas);  // nothing more
    CHECK(tas.size() == 1);
    // A TP in the flushed window is late and dropped.
    swift(tp, tas);
    swift.flush(invalid_timestamp, tas);
    CHECK(tas.size() == 1);

    // Preprocess keeps a TP if either peak or length passes.
    TAMakerSWIFT keep(cfg);
    TriggerPrimitive weak = tp;
    weak.adc_peak = 10;
    weak.samples_over_threshold = 2;
    std::vector<TriggerActivity> none;
    keep(weak, none);
    keep.flush(invalid_timestamp, none);
    CHECK(none.empty());
    weak.samples_over_threshold = 8;
    keep(weak, none);
    keep.flush(invalid_timestamp, none);
    CHECK(none.size() == 1);
}
