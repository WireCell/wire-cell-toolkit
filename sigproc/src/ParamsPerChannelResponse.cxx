#include "WireCellSigProc/ParamsPerChannelResponse.h"

#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Response.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(ParamsPerChannelResponse, WireCell::SigProc::ParamsPerChannelResponse, WireCell::IChannelResponse,
                 WireCell::IConfigurable)

using namespace WireCell;

SigProc::ParamsPerChannelResponse::ParamsPerChannelResponse(const char* filename)
  : m_filename(filename)
{
}

SigProc::ParamsPerChannelResponse::~ParamsPerChannelResponse() {}

WireCell::Configuration SigProc::ParamsPerChannelResponse::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_filename;
    return cfg;
}

void SigProc::ParamsPerChannelResponse::configure(const WireCell::Configuration& cfg)
{
    m_filename = get(cfg, "filename", m_filename);
    if (m_filename.empty()) {
        THROW(ValueError() << errmsg{"must supply a ParamsPerChannelResponse filename"});
    }

    auto top = Persist::load(m_filename);
    const double tick = top["tick"].asFloat();
    const double t0 = top["t0"].asFloat();
    const int nsamp = top["nsamp"].asInt();
    if (!m_bins.nbins()) {  // first time
        m_bins = Binning(nsamp, t0, t0 + nsamp * tick);
    }

    auto jchannels = top["channel_info"];
    if (jchannels.isNull()) {
        THROW(ValueError() << errmsg{"no channels given in file " + m_filename});
    }

    for (auto jchinfo : jchannels) {
        const double gain = jchinfo["gain"].asFloat();
        const double shaping = jchinfo["shaping"].asFloat();
        const double k3 = jchinfo["k3"].asFloat();
        const double k4 = jchinfo["k4"].asFloat();
        const double k5 = jchinfo["k5"].asFloat();
        const double k6 = jchinfo["k6"].asFloat();
        Response::ParamsColdElec ce(gain, shaping, k3, k4, k5, k6); // 14mV/fC, 2.0us, 0.1, 0.1, 0.0, 0.0
        Waveform::realseq_t resp = ce.generate(m_bins);

        if (jchinfo["channels"].isArray()) {
          for (auto& ch: jchinfo["channels"]){
            m_cr[ch.asInt()] = resp;
            // std::cout << "[ParamsPerChannelResponse] channel: " << ch.asInt()
            //           << " gain: " << gain
            //           << " shaping: " << shaping
            //           << " k3: " << k3
            //           << " k4: " << k4
            //           << " k5: " << k5
            //           << " k6: " << k6
            //           << std::endl;
          }
        }
    }
}

const Waveform::realseq_t& SigProc::ParamsPerChannelResponse::channel_response(int channel_ident) const
{
    const auto& it = m_cr.find(channel_ident);
    if (it == m_cr.end()) {
        THROW(KeyError() << errmsg{String::format("no response for channel %d", channel_ident)});
    }
    return it->second;
}

Binning SigProc::ParamsPerChannelResponse::channel_response_binning() const { return m_bins; }
