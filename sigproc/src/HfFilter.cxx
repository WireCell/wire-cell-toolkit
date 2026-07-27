#include "WireCellSigProc/HfFilter.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(HfFilter, WireCell::SigProc::HfFilter, WireCell::IFilterWaveform, WireCell::IConfigurable)

using namespace WireCell;

SigProc::HfFilter::HfFilter(double max_freq, double sigma, double power, bool flag)
  : m_max_freq(max_freq)
  , m_sigma(sigma)
  , m_power(power)
  , m_flag(flag)
{
}

SigProc::HfFilter::~HfFilter() {}

WireCell::Configuration SigProc::HfFilter::default_configuration() const
{
    Configuration cfg;
    cfg["sigma"] = m_sigma;
    cfg["power"] = m_power;
    cfg["flag"] = m_flag;
    cfg["max_freq"] = m_max_freq;
    cfg["use_negative_freqs"] = m_use_negative_freqs;
    return cfg;
}

void SigProc::HfFilter::configure(const WireCell::Configuration& cfg)
{
    m_sigma = get(cfg, "sigma", m_sigma);
    m_power = get(cfg, "power", m_power);
    m_flag = get(cfg, "flag", m_flag);
    m_use_negative_freqs = get(cfg, "use_negative_freqs", m_use_negative_freqs);
    m_max_freq = get(cfg, "max_freq", m_max_freq);
}

const Waveform::realseq_t SigProc::HfFilter::filter_waveform(int nbins) const
{
    Waveform::realseq_t m_wfs(nbins);

    Response::HfFilter hf_filter(m_sigma, m_power, m_flag);

    double freq_factor = (m_use_negative_freqs ? 2. : 1.);

    for (size_t i = 0; i != m_wfs.size(); i++) {
        double freq = i * 1.0 / int(m_wfs.size()) * freq_factor * m_max_freq;
        if ((freq > m_max_freq) && m_use_negative_freqs) {
            freq = freq - 2 * m_max_freq;
        }
        else if ((freq > m_max_freq) && !m_use_negative_freqs) {
            THROW(ValueError() << errmsg{"HfFilter::filter_waveform -- Passed max freq when not using negative freqs"});

        }
        m_wfs.at(i) = hf_filter(fabs(freq));

    }
    // std::cout << m_wfs.size() << std::endl;

    return m_wfs;
}
