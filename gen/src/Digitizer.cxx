#include "WireCellGen/Digitizer.h"

#include "WireCellIface/IWireSelectors.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"

#include "WireCellAux/FrameTools.h"

#include "WireCellUtil/Testing.h"
#include "WireCellUtil/NamedFactory.h"
#include <sstream>

WIRECELL_FACTORY(Digitizer, WireCell::Gen::Digitizer,
                 WireCell::INamed,
                 WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

Gen::Digitizer::Digitizer()
  : Aux::Logger("Digitizer", "gen")
{
}

Gen::Digitizer::~Digitizer() {}

WireCell::Configuration Gen::Digitizer::default_configuration() const
{
    Configuration cfg;
    put(cfg, "anode", "AnodePlane");

    put(cfg, "resolution", m_resolution);
    put(cfg, "gain", m_gain);
    Configuration fs(Json::arrayValue);  // fixme: sure would be nice if we had some templated sugar for this
    for (int ind = 0; ind < 2; ++ind) {
        fs[ind] = m_fullscale[ind];
    }
    cfg["fullscale"] = fs;

    Configuration bl(Json::arrayValue);
    for (int ind = 0; ind < 3; ++ind) {
        bl[ind] = m_baselines[ind];
    }
    cfg["baselines"] = bl;
    cfg["round"] = m_round;

    cfg["frame_tag"] = m_frame_tag;
    return cfg;
}

void Gen::Digitizer::configure(const Configuration& cfg)
{
    auto anode_tn = get<string>(cfg, "anode", "AnodePlane");
    m_anode = Factory::find_tn<IAnodePlane>(anode_tn);

    m_resolution = get(cfg, "resolution", m_resolution);
    m_gain = get(cfg, "gain", m_gain);
    m_fullscale = get(cfg, "fullscale", m_fullscale);
    m_baselines = get(cfg, "baselines", m_baselines);
    m_round = get(cfg, "round", m_round);
    m_frame_tag = get(cfg, "frame_tag", m_frame_tag);

    std::stringstream ss;
    ss << "tag=\"" << m_frame_tag << "\", "
       << "resolution=" << m_resolution << " bits, "
       << "maxvalue=" << (1 << m_resolution) << " counts, "
       << "gain=" << m_gain << ", "
       << "fullscale=[" << m_fullscale[0] / units::mV << "," << m_fullscale[1] / units::mV << "] mV, "
       << "baselines=[" << m_baselines[0] / units::mV << "," << m_baselines[1] / units::mV << ","
       << m_baselines[2] / units::mV << "] mV "
       << "round=" << m_round;
    log->debug(ss.str());
}

double Gen::Digitizer::digitize(double voltage)
{
    const int adcmaxval = (1 << m_resolution) - 1;

    if (voltage <= m_fullscale[0]) {
        return 0;
    }
    if (voltage >= m_fullscale[1]) {
        return adcmaxval;
    }
    const double relvoltage = (voltage - m_fullscale[0]) / (m_fullscale[1] - m_fullscale[0]);
    const double fp_adc = relvoltage * adcmaxval;
    if (m_round) {
        return round(fp_adc);
    }
    return floor(fp_adc);
}

bool Gen::Digitizer::operator()(const input_pointer& vframe, output_pointer& adcframe)
{
    if (!vframe) {  // EOS
        log->debug("see EOS at call={}", m_count);
        adcframe = nullptr;
        ++m_count;
        return true;
    }

    // fixme: maybe make this honor a tag
    // auto vtraces = Aux::untagged_traces(vframe);
    const auto& vtraces = *vframe->traces();
    if (vtraces.empty()) {
        // Empty is not EOS, must output a corresponding empty frame.
        adcframe = std::make_shared<SimpleFrame>(vframe->ident(), vframe->time(), vframe->tick());

        log->warn("no traces in at call={}: {}", m_count, Aux::taginfo(adcframe));
        ++m_count;
        return true;
    }

    // Get extent in channel and tbin
    auto channels = Aux::channels(vtraces);
    std::sort(channels.begin(), channels.end());
    auto chbeg = channels.begin();
    auto chend = std::unique(chbeg, channels.end());
    auto tbinmm = Aux::tbin_range(vtraces);

    const size_t ncols = tbinmm.second - tbinmm.first;
    const size_t nrows = std::distance(chbeg, chend);

    // make a dense array working space.  a row is one trace.  a
    // column is one tick.
    Array::array_xxf arr = Array::array_xxf::Zero(nrows, ncols);
    Aux::fill(arr, vtraces, channels.begin(), chend, tbinmm.first);

    ITrace::vector adctraces(nrows);

    double totadc = 0;

    for (size_t irow = 0; irow < nrows; ++irow) {
        int ch = channels[irow];
        WirePlaneId wpid = m_anode->resolve(ch);
        if (!wpid.valid()) {
            log->warn("got invalid WPID for channel {}: {}, skipping", ch, wpid);
            continue;
        }
        const float baseline = m_baselines[wpid.index()];

        ITrace::ChargeSequence adcwave(ncols);
        for (size_t icol = 0; icol < ncols; ++icol) {
            double voltage = m_gain * arr(irow, icol) + baseline;
            const float adcf = digitize(voltage);
            adcwave[icol] = adcf;
            totadc += adcf;
        }
        adctraces[irow] = make_shared<SimpleTrace>(ch, tbinmm.first, adcwave);
    }
    auto sframe = make_shared<SimpleFrame>(vframe->ident(), vframe->time(), adctraces, vframe->tick(), vframe->masks());
    if (!m_frame_tag.empty()) {
        sframe->tag_frame(m_frame_tag);
    }
    adcframe = sframe;

    log->debug("call={} traces={} frame={} totadc={} outtag=\"{}\"",
               m_count,
               adctraces.size(), vframe->ident(), totadc, m_frame_tag);

    log->debug("input : {}", Aux::taginfo(vframe));
    log->debug("output: {}", Aux::taginfo(adcframe));

    ++m_count;
    return true;
}
