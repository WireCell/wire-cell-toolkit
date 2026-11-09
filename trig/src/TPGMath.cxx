#include "WireCellTrig/TPGMath.h"

#include <cstdlib>
#include <limits>
#include <map>

using namespace WireCell::Trig;

int16_t WireCell::Trig::avx2_divide(int16_t a, int16_t b)
{
    int16_t vb = (1 << 15) / b;   //  1 / b * 2^15
    int32_t mulhrs = a * vb;      //  a / b * 2^15
    mulhrs = (mulhrs >> 14) + 1;  // ~ a / b * 2 + 1
    mulhrs = mulhrs >> 1;         // ~ a / b, the +1 gives the odd rounding
    return (int16_t) mulhrs;
}

int16_t WireCell::Trig::mode_pedestal(const std::vector<int16_t>& adcs)
{
    std::map<int16_t, size_t> counts;
    for (const auto& adc : adcs) {
        counts[adc] += 1;
    }
    int16_t mode = 0;
    size_t max_counts = 0;
    for (const auto& [value, count] : counts) {
        if (count > max_counts) {
            mode = value;
            max_counts = count;
        }
    }
    return mode;
}

int16_t FrugalPedestal::operator()(int16_t sample)
{
    if (sample > pedestal) ++accum;
    if (sample < pedestal) --accum;
    if (accum > accum_limit) {
        ++pedestal;
        accum = 0;
    }
    if (accum < -1 * accum_limit) {
        --pedestal;
        accum = 0;
    }
    return pedestal;
}

int16_t RunningSum::operator()(int16_t sample)
{
    const int32_t next = r_value * avx2_divide(sum, 10) + avx2_divide(sample, 10) * scale_factor;
    if (next > std::numeric_limits<int16_t>::max()) {
        sum = std::numeric_limits<int16_t>::max();
    }
    else if (saturate_low && next < std::numeric_limits<int16_t>::min()) {
        sum = std::numeric_limits<int16_t>::min();
    }
    else {
        sum = (int16_t) next;  // wraps below INT16_MIN, as tpglibs does
    }
    return sum;
}

int16_t AbsRunningSum::operator()(int16_t sample)
{
    return RunningSum::operator()((int16_t) std::abs(sample));
}

bool ThresholdHitter::operator()(int64_t sample_count, int16_t sample, ThresholdHit& hit)
{
    const bool is_over = sample > threshold;
    if (is_over) {
        m_charge += sample;
        if (sample > m_peak_adc) {
            m_peak_adc = (uint16_t) sample;
            m_peak_time = m_tover;
        }
        ++m_tover;
    }
    bool ended = false;
    if (m_over && !is_over) {
        ended = flush(sample_count, hit);
    }
    m_over = is_over;
    return ended;
}

bool ThresholdHitter::flush(int64_t sample_count, ThresholdHit& hit)
{
    if (m_tover == 0 || m_tover < min_samples_over_threshold) {
        reset();
        return false;
    }
    hit.sample_start = sample_count - m_tover;
    hit.samples_over_threshold = m_tover;
    hit.samples_to_peak = m_peak_time;
    hit.adc_integral = m_charge;
    hit.adc_peak = m_peak_adc;
    reset();
    return true;
}

void ThresholdHitter::reset()
{
    m_over = false;
    m_tover = 0;
    m_peak_time = 0;
    m_peak_adc = 0;
    m_charge = 0;
}
