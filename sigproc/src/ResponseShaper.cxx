/**
 * A generic filter to correct reponse functions 
 * Note: the purpose of this filter is not to replace the perchannel response
 * correction in the OmnibusSigProc, instead it's implemented for validation purpose
 */

#include "WireCellSigProc/ResponseShaper.h"
#include "WireCellSigProc/Derivations.h"

#include "WireCellAux/DftTools.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Response.h"
#include "WireCellIface/IChannelResponse.h"
#include "WireCellIface/IWaveform.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <set>

WIRECELL_FACTORY(perChannelShaper, WireCell::SigProc::ResponseShaper::OneChannelResponse, WireCell::IChannelFilter,
                 WireCell::IConfigurable)

using namespace WireCell::SigProc;
using WireCell::Aux::DftTools::fwd_r2c;
using WireCell::Aux::DftTools::inv_c2r;


ResponseShaper::OneChannelResponse::OneChannelResponse(const std::string& anode, const std::string& noisedb)
  : m_anode_tn(anode)
  , m_noisedb_tn(noisedb)
{
}
ResponseShaper::OneChannelResponse::~OneChannelResponse() {}

void ResponseShaper::OneChannelResponse::configure(const WireCell::Configuration& cfg)
{
    m_elecresponse_tn = get(cfg, "elecresponse", m_elecresponse_tn); 
    m_anode_tn = get(cfg, "anode", m_anode_tn);
    m_anode = Factory::find_tn<IAnodePlane>(m_anode_tn);
    m_noisedb_tn = get(cfg, "noisedb", m_noisedb_tn);
    m_noisedb = Factory::find_tn<IChannelNoiseDatabase>(m_noisedb_tn);

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);
}
WireCell::Configuration ResponseShaper::OneChannelResponse::default_configuration() const
{
    Configuration cfg;
    cfg["elecresponse"] = m_elecresponse_tn;
    cfg["anode"] = m_anode_tn;
    cfg["noisedb"] = m_noisedb_tn;
    cfg["dft"] = "FftwDFT";     // type-name for the DFT to use
    return cfg;
}

WireCell::Waveform::ChannelMaskMap ResponseShaper::OneChannelResponse::apply(int ch, signal_t& signal) const
{
    WireCell::Waveform::ChannelMaskMap ret;

    // correct rc undershoot
    auto spectrum = fwd_r2c(m_dft, signal);

    // now apply the ch-by-ch response ...
    auto cr = Factory::find_tn<IChannelResponse>("ParamsPerChannelResponse");
    auto cr_bins = cr->channel_response_binning();
    // if (cr_bins.binsize() != m_period) { // FIXME
    //     THROW(ValueError() << errmsg{"ProtoduneHD::OneChannelResponse channel response size mismatch"});
    // }
    auto period = cr_bins.binsize();
    WireCell::Binning tbins(signal.size(), cr_bins.min(), cr_bins.min() + signal.size() * period);

    // Option 1: declare response function explicitly
    // Response::ColdElec ce(14 * units::mV / units::fC, 2.2 * units::us);
    // auto ewave = ce.generate(tbins);
    
    // Option 2: retrieve response function from configuration
    auto elecresponse = Factory::find_tn<IWaveform>(m_elecresponse_tn);
    auto ewave = elecresponse->waveform_samples(tbins);

    const WireCell::Waveform::compseq_t elec = fwd_r2c(m_dft, ewave);

    Waveform::realseq_t tch_resp = cr->channel_response(ch);
    tch_resp.resize(signal.size(), 0);
    const WireCell::Waveform::compseq_t ch_elec = fwd_r2c(m_dft, tch_resp);

    for (unsigned int ind = 0; ind != elec.size(); ind ++) {
        const auto four = ch_elec.at(ind);
        if (std::abs(four) != 0) {
            spectrum.at(ind) *= elec.at(ind) / four;
        }
        else {
            spectrum.at(ind) = 0;
        }
    }


    // remove the DC component
    spectrum.front() = 0;
    signal = inv_c2r(m_dft, spectrum);

    // Now calculate the baseline ...
    std::pair<double, double> temp = WireCell::Waveform::mean_rms(signal);
    auto temp_signal = signal;
    for (size_t i = 0; i != temp_signal.size(); i++) {
        if (fabs(temp_signal.at(i) - temp.first) > 6 * temp.second) {
            temp_signal.at(i) = temp.first;
        }
    }
    float baseline = WireCell::Waveform::median_binned(temp_signal);
    // correct baseline
    WireCell::Waveform::increase(signal, baseline * (-1));


    return ret;
}

WireCell::Waveform::ChannelMaskMap ResponseShaper::OneChannelResponse::apply(channel_signals_t& chansig) const
{
    return WireCell::Waveform::ChannelMaskMap();
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
