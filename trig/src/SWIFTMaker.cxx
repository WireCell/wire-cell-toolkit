#include "WireCellTrig/SWIFTMaker.h"

#include <algorithm>

using namespace WireCell::Trig;
using namespace WireCell::Trig::Daq;

TAMakerSWIFT::TAMakerSWIFT()
  : m_cfg()
{
}

TAMakerSWIFT::TAMakerSWIFT(const Config& cfg)
  : m_cfg(cfg)
{
}

bool TAMakerSWIFT::preprocess(const input_type& tp)
{
    // As triggeralgs, this replaces the base max_samples_over_threshold cut.
    return !(tp.adc_peak < m_cfg.min_adc_peak && tp.samples_over_threshold < m_cfg.min_samples_over_threshold);
}

void TAMakerSWIFT::reset_window(timestamp_t window_start)
{
    m_window_start = window_start;
    m_window_energy = 0;
    m_ta = TriggerActivity{};
    m_ta.time_start = window_start;
}

void TAMakerSWIFT::process(const input_type& tp, output_vector& out)
{
    const timestamp_t tp_window_start = (tp.time_start / m_cfg.window_length) * m_cfg.window_length;
    if (!m_initialised) {
        reset_window(tp_window_start);
        m_initialised = true;
    }
    if (tp_window_start < m_window_start) {
        return;  // late TP, its window is closed
    }
    if (tp_window_start > m_window_start) {
        close_window(out);
        reset_window(tp_window_start);
    }
    m_ta.inputs.push_back(tp);
    m_window_energy += tp.adc_integral;
}

void TAMakerSWIFT::do_flush(timestamp_t until, output_vector& out)
{
    if (!m_initialised) return;
    const timestamp_t window_end = m_window_start + m_cfg.window_length;
    if (until < window_end) return;
    close_window(out);
    reset_window(window_end);
}

void TAMakerSWIFT::close_window(output_vector& out)
{
    if (m_ta.inputs.empty()) return;
    if (m_window_energy < m_cfg.inspect_energy_threshold) return;  // reject
    if (m_window_energy < m_cfg.accept_energy_threshold) {         // inspect
        if (dominant_cluster_energy(m_ta.inputs, m_window_start) <= m_cfg.cluster_energy_cut) return;
    }
    set_ta_attributes();
    out.push_back(m_ta);
}

uint64_t TAMakerSWIFT::dominant_cluster_energy(const std::vector<input_type>& tps, timestamp_t window_start) const
{
    struct Point {
        float z, t;
        uint64_t adc;
    };
    std::vector<Point> points;
    points.reserve(tps.size());
    for (const auto& tp : tps) {
        points.push_back({tp.channel * m_cfg.wire_pitch, (tp.time_start - window_start) * m_cfg.cm_per_tick,
                          tp.adc_integral});
    }
    const int npoints = points.size();
    const float eps2 = m_cfg.epsilon * m_cfg.epsilon;
    auto neighbors = [&](int ind) {
        std::vector<int> ret;
        for (int other = 0; other < npoints; ++other) {
            const float dz = points[other].z - points[ind].z;
            const float dt = points[other].t - points[ind].t;
            if (dz * dz + dt * dt <= eps2) ret.push_back(other);
        }
        return ret;
    };

    int cluster_id = 0;
    std::vector<int> labels(npoints, -1);
    std::vector<uint8_t> visited(npoints, 0);
    for (int ind = 0; ind < npoints; ++ind) {
        if (visited[ind]) continue;
        visited[ind] = 1;
        auto neigh = neighbors(ind);
        if (neigh.size() < static_cast<size_t>(m_cfg.min_samples)) {
            labels[ind] = -1;
            continue;
        }
        labels[ind] = cluster_id;
        for (size_t k = 0; k < neigh.size(); ++k) {
            const int jnd = neigh[k];
            if (!visited[jnd]) {
                visited[jnd] = 1;
                auto neigh2 = neighbors(jnd);
                if (neigh2.size() >= static_cast<size_t>(m_cfg.min_samples)) {
                    for (int mnd : neigh2) {
                        if (std::find(neigh.begin(), neigh.end(), mnd) == neigh.end()) neigh.push_back(mnd);
                    }
                }
            }
            if (labels[jnd] == -1) labels[jnd] = cluster_id;
        }
        ++cluster_id;
    }

    std::vector<uint64_t> sums(cluster_id, 0);
    for (int ind = 0; ind < npoints; ++ind) {
        if (labels[ind] >= 0) sums[labels[ind]] += points[ind].adc;
    }
    if (sums.empty()) return 0;
    return *std::max_element(sums.begin(), sums.end());
}

void TAMakerSWIFT::set_ta_attributes()
{
    m_ta.time_start = m_window_start;
    m_ta.time_end = m_window_start + m_cfg.window_length;
    m_ta.adc_integral = m_window_energy;
    const auto& first = m_ta.inputs.front();
    m_ta.detid = first.detid;
    m_ta.type = TriggerActivityData::Type::kTPC;
    m_ta.algorithm = TriggerActivityData::Algorithm::kSWIFT;
    channel_t min_ch = first.channel, max_ch = first.channel;
    m_ta.adc_peak = 0;
    for (const auto& tp : m_ta.inputs) {
        min_ch = std::min(min_ch, tp.channel);
        max_ch = std::max(max_ch, tp.channel);
        if (tp.adc_peak > m_ta.adc_peak) {
            m_ta.adc_peak = tp.adc_peak;
            m_ta.channel_peak = tp.channel;
            m_ta.time_peak = tp.time_start;
        }
    }
    m_ta.channel_start = min_ch;
    m_ta.channel_end = max_ch;
    m_ta.time_activity = m_ta.time_peak;
}
