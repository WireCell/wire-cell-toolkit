#include "WireCellSigProc/FilterResponse.h"

#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Response.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(FilterResponse, WireCell::SigProc::FilterResponse, WireCell::IChannelResponse,
                 WireCell::IConfigurable)

using namespace WireCell;

SigProc::FilterResponse::FilterResponse(const char* filename, const int planeid)
  : m_filename(filename), m_planeid(planeid)
{
}

SigProc::FilterResponse::~FilterResponse() {}

WireCell::Configuration SigProc::FilterResponse::default_configuration() const
{
    Configuration cfg;
    cfg["filename"] = m_filename;
    cfg["planeid"] = m_planeid;
    return cfg;
}

void SigProc::FilterResponse::configure(const WireCell::Configuration& cfg)
{

    m_filename = get(cfg, "filename", m_filename);
    if (m_filename.empty()) {
        THROW(ValueError() << errmsg{"must supply a FilterResponse filename"});
    }
    m_planeid = get(cfg, "planeid", m_planeid);

    auto top = Persist::load(m_filename);
    // const int nwires = top["nwires"].asInt();
    const int nticks = top["nticks"].asInt();
    if (!m_bins.nbins()) {  // first time
        m_bins = Binning(nticks, 0, nticks); // tick range for filter
    }
    auto jflts = top["filters"];
    if (jflts.isNull()) {
        THROW(ValueError() << errmsg{"no filters given in file " + m_filename});
    }

    for (auto jflt : jflts) {
        const int plane = jflt["plane"].asInt();
        if (plane != m_planeid) continue;

        const int wire = jflt["wire"].asInt();
        Waveform::realseq_t resp(nticks, 0.0);
        std::transform(jflt["values"].begin(), jflt["values"].end(), resp.begin(),
                       [](const auto& v) { return v.asFloat(); });
        m_cr[wire] = resp;
    }
}

const Waveform::realseq_t& SigProc::FilterResponse::channel_response(int channel_ident) const
{
    const auto& it = m_cr.find(channel_ident);
    if (it == m_cr.end()) {
        THROW(KeyError() << errmsg{String::format("no response for channel %d", channel_ident)});
    }
    return it->second;
}

Binning SigProc::FilterResponse::channel_response_binning() const { return m_bins; }
