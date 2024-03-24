#include "WireCellGen/RCResponse.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Response.h"

WIRECELL_FACTORY(RCResponse, WireCell::Gen::RCResponse, WireCell::IWaveform, WireCell::IConfigurable)

using namespace WireCell;

Gen::RCResponse::RCResponse(int nticks, double t0, double width, double tick)
{
    m_cfg["start"] = t0;
    m_cfg["tick"] = tick;
    m_cfg["nticks"] = nticks;
    m_cfg["width"] = width;
}
WireCell::Configuration Gen::RCResponse::default_configuration() const { return m_cfg; }
void Gen::RCResponse::configure(const WireCell::Configuration& cfg)
{
    m_cfg = cfg;

    // fixme: why give SimpleRC tick twice?  Once in ctor and once in
    // generate()?
    m_rc = std::make_unique<Response::SimpleRC>(m_cfg["width"].asDouble(), waveform_period(), waveform_start());

    const int nbins = m_cfg["nticks"].asInt();
    const double t0 = waveform_start();
    const double tick = waveform_period();
    Binning bins(nbins, t0, t0 + nbins * tick);
    m_wave = m_rc->generate(bins);
}

double Gen::RCResponse::waveform_start() const { return m_cfg["start"].asDouble(); }

double Gen::RCResponse::waveform_period() const { return m_cfg["tick"].asDouble(); }

const IWaveform::sequence_type& Gen::RCResponse::waveform_samples() const { return m_wave; }

IWaveform::sequence_type Gen::RCResponse::waveform_samples(const WireCell::Binning& tbins) const
{
    sequence_type rebinned_wave = m_rc->generate(tbins);
    return rebinned_wave;
}
