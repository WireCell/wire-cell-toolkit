#include "WireCellTrig/TPGChannel.h"

#include <algorithm>

using namespace WireCell::Trig;

TPGChannel::TPGChannel(const TPGConfig& cfg)
  : m_cfg(cfg)
{
    reset();
}

void TPGChannel::reset()
{
    m_fresh = true;
    m_pedestal = FrugalPedestal{};
    m_pedestal.accum_limit = m_cfg.accum_limit;
    m_rs = RunningSum{};
    m_rs.r_value = m_cfg.r_value;
    m_rs.scale_factor = m_cfg.scale_factor;
    m_rs.saturate_low = m_cfg.saturate_low;
    m_ars = AbsRunningSum{};
    m_ars.r_value = m_cfg.r_value;
    m_ars.scale_factor = m_cfg.scale_factor;
    m_ars.saturate_low = m_cfg.saturate_low;
    m_filter_pedestal = FrugalPedestal{};
    m_filter_pedestal.accum_limit = m_cfg.accum_limit;
    m_hitter = ThresholdHitter{};
    m_hitter.threshold = m_cfg.threshold;
    m_hitter.min_samples_over_threshold = m_cfg.min_samples_over_threshold;
}

void TPGChannel::initialize(const std::vector<int16_t>& samples)
{
    switch (m_cfg.pedestal_init) {
    case TPGConfig::PedestalInit::zero:
        m_pedestal.pedestal = 0;
        break;
    case TPGConfig::PedestalInit::mode:
        m_pedestal.pedestal = mode_pedestal(samples);
        break;
    case TPGConfig::PedestalInit::fixed:
        m_pedestal.pedestal = m_cfg.pedestal;
        break;
    }
    m_fresh = false;
}

int16_t TPGChannel::filter(int16_t sample)
{
    sample -= m_pedestal(sample);
    switch (m_cfg.filter) {
    case TPGConfig::Filter::none:
        break;
    case TPGConfig::Filter::running_sum:
        sample = m_rs(sample);
        break;
    case TPGConfig::Filter::abs_running_sum:
        sample = m_ars(sample);
        break;
    }
    if (m_cfg.filter_pedestal) {
        sample -= m_filter_pedestal(sample);
    }
    return sample;
}

void TPGChannel::settle(const std::vector<int16_t>& samples, size_t nsamples)
{
    if (m_fresh) initialize(samples);
    nsamples = std::min(nsamples, samples.size());
    for (size_t ind = 0; ind < nsamples; ++ind) {
        filter(samples[ind]);
    }
}

void TPGChannel::operator()(const std::vector<int16_t>& samples, int64_t first, std::vector<ThresholdHit>& hits)
{
    if (m_cfg.threshold < 0) return;
    if (m_fresh) initialize(samples);
    ThresholdHit hit;
    const size_t nsamples = samples.size();
    for (size_t ind = 0; ind < nsamples; ++ind) {
        const int16_t sample = filter(samples[ind]);
        if (m_hitter(first + (int64_t) ind, sample, hit)) {
            hits.push_back(hit);
        }
    }
}

void TPGChannel::flush(int64_t next, std::vector<ThresholdHit>& hits)
{
    ThresholdHit hit;
    if (m_hitter.flush(next, hit)) {
        hits.push_back(hit);
    }
}
