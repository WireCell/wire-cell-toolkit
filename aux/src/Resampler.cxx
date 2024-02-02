#include "WireCellAux/Resampler.h"
#include "WireCellUtil/LMN.h"
#include "WireCellUtil/Units.h"

#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/DftTools.h"

#include "WireCellUtil/NamedFactory.h"

#include <map>
#include <unordered_set>

WIRECELL_FACTORY(Resampler, WireCell::Aux::Resampler,
                 WireCell::INamed,
                 WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;
using WireCell::Aux::DftTools::fwd_r2c;
using WireCell::Aux::DftTools::inv_c2r;
using namespace WireCell::Aux;
using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

Aux::Resampler::Resampler()
    : Aux::Logger("Resampler", "aux")
{
}

Aux::Resampler::~Resampler() {}

WireCell::Configuration Aux::Resampler::default_configuration() const
{
    Configuration cfg;

    cfg["tags"] = Json::arrayValue;
    cfg["period"] = m_period;
    cfg["nticks"] = m_nticks;
    cfg["dft"] = "FftwDFT";

    // fixme: padding strategy string may be added in future

    return cfg;
}

void Aux::Resampler::configure(const WireCell::Configuration& cfg)
{
    m_period = get(cfg, "period", m_period);
    m_nticks = get(cfg, "nticks", m_nticks);
    auto dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);
}


bool Aux::Resampler::operator()(const input_pointer& inframe, output_pointer& outframe)
{
    outframe = nullptr;
    if (!inframe) {
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    const double Ts = inframe->tick();
    const double Tr = m_period;
    const size_t Nrat = LMN::rational(Ts, Tr);

    std::unordered_map< std::string, IFrame::trace_list_t> tag_indicies;

    ITrace::vector out_traces;
    for (const auto& trace : *inframe->traces()) {

        if (trace->tbin()) {
            raise<ValueError>("currently no support for nonzero tbin (fixme)");
        }

        auto wave = trace->charge();

        const size_t Ns_orig = wave.size();
        const size_t Ns_pad = LMN::nbigger(Ns_orig, Nrat);
        const double duration = Ns_pad * Ts;
        const size_t Nr = duration / Tr; // check for error?

        wave = LMN::resize(wave, Ns_pad);
        auto spec = fwd_r2c(m_dft, wave);
        spec = LMN::resample(spec, Nr);
        wave = inv_c2r(m_dft, spec);

        // Interpolation interpretation.
        double norm = wave.size() / (double)Ns_orig; 
        Waveform::scale(wave, norm);

        out_traces.push_back(std::make_shared<SimpleTrace>(trace->channel(), 0, wave));
    }

    auto sf = std::make_shared<SimpleFrame>(inframe->ident(), inframe->time(), out_traces, Tr);

    for (const auto& frame_tag : inframe->frame_tags()) {
        sf->tag_frame(frame_tag);
    }
    for (const auto& trace_tag : inframe->trace_tags()) {
        auto tl = inframe->tagged_traces(trace_tag);
        auto ts = inframe->trace_summary(trace_tag);
        sf->tag_traces(trace_tag, tl, ts);
    }
    outframe = sf;
    ++m_count;
    return true;
}
