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

    cfg["period"] = m_period;
    cfg["dft"] = "FftwDFT";
    cfg["time_pad"] = m_time_pad;
    cfg["time_sizing"] = m_time_sizing;
    return cfg;
}

void Aux::Resampler::configure(const WireCell::Configuration& cfg)
{
    m_period = get(cfg, "period", m_period);
    m_time_pad = get(cfg, "time_pad", m_time_pad);
    m_time_sizing = get(cfg, "time_sizing", m_time_sizing);
    auto dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);
    log->debug("resample period={}, time_pad={} time_sizing={}",
               m_period, m_time_pad, m_time_sizing);
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

    if (Ts == Tr) {
        log->warn("same periods: Ts={} Tr={} at call={} is probably not what you want", Ts, Tr, m_count);
        outframe = inframe;
        ++m_count;
        return true;
    }

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
        if (Ns_pad > Ns_orig) {
            const float first = wave[0];
            const float last = wave[Ns_orig-1];

            auto start = wave.begin()+Ns_orig;

            if (m_time_pad == "linear") {
                LMN::fill_linear(start, wave.end(), last, first);
            }
            // else if (m_time_pad == "cosine") {
            //     LMN::fill_cosine(start, wave.end(), last, first);
            // }
            else {
                float pad = 0;      // "zero"
                if (m_time_pad == "first") {
                    pad = first;
                }
                else if (m_time_pad == "last") {
                    pad = last;
                }
                LMN::fill_constant(start, wave.end(), pad);
            }
        }

        auto spec = fwd_r2c(m_dft, wave);
        spec = LMN::resample(spec, Nr);
        wave = inv_c2r(m_dft, spec);

        // Interpolation interpretation.
        double norm = wave.size() / (double)Ns_pad; 
        Waveform::scale(wave, norm);

        if (m_time_sizing == "duration") {
            const size_t Nr_unpadded = (int)(Ns_orig * (Ts / Tr));
            wave = LMN::resize(wave, Nr_unpadded);
        }
        // else: "padded" is as-is

        if (out_traces.empty()) {
            log->debug("first ch={} Ts={} Ns={} Ns_pad={} Nrat={} Tr={} Nr={} Nout={} padding:{}",
                       trace->channel(), Ts, Ns_orig, Ns_pad, Nrat, Tr, Nr, wave.size(), m_time_pad);
        }
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
    log->debug("resample {} traces at call={}", out_traces.size(), m_count);

    ++m_count;
    return true;
}
