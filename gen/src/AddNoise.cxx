#include "WireCellGen/AddNoise.h"

#include "WireCellAux/RandTools.h"
#include "WireCellAux/NoiseTools.h"

#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleFrame.h"


#include <unordered_map>

WIRECELL_FACTORY(IncoherentAddNoise, WireCell::Gen::IncoherentAddNoise,
                 WireCell::INamed, WireCell::IFrameFilter, WireCell::IConfigurable)
WIRECELL_FACTORY(CoherentAddNoise, WireCell::Gen::CoherentAddNoise,
                 WireCell::INamed, WireCell::IFrameFilter, WireCell::IConfigurable)

// Historical aliases 
WIRECELL_FACTORY(AddNoise, WireCell::Gen::IncoherentAddNoise,
                 WireCell::INamed, WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace std;
using namespace WireCell;

using WireCell::Aux::SimpleTrace;
using WireCell::Aux::SimpleFrame;

using namespace WireCell::Aux::RandTools;
using namespace WireCell::Aux::NoiseTools;

Gen::IncoherentAddNoise::IncoherentAddNoise()
    : Gen::NoiseBaseT<IChannelSpectrum>("IncoherentAddNoise")
{
}
Gen::IncoherentAddNoise::~IncoherentAddNoise() {}

Gen::CoherentAddNoise::CoherentAddNoise()
    : Gen::NoiseBaseT<IGroupSpectrum>("CoherentAddNoise")
{
}
Gen::CoherentAddNoise::~CoherentAddNoise() {}


bool Gen::IncoherentAddNoise::operator()(const input_pointer& inframe, output_pointer& outframe)
{
    if (!inframe) {
        outframe = nullptr;
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    // We "recycle" the randoms in order to speed up generation.  Each
    // array will start from a random location in the recycled buffer
    // and a few percent will be "freshened".  This results in a small
    // amount of coherency between nearby channels.
    const double BUG=m_bug202; // 0.0 is correct, passing rec. perc. 0.04 opened issue 202
    auto rn = Normals::make_recycling(m_rng, 2*m_nsamples, BUG, 1, 2*m_rep_percent);
    GeneratorN rwgen(m_dft, rn);

    // Limit number of warnings below
    static bool warned = false;

    // Make waveforms of size nsample from each model, adding only
    // ncharge of their element to the trace charge.  This
    // full-nsample followed by ncharge-truncation may CPU-wasteful in
    // the sparse traces case.

    const float sqrt2opi = sqrt(2.0/3.141592);

    ITrace::vector outtraces;
    for (const auto& intrace : *inframe->traces()) {

        const int chid = intrace->channel();
        auto charge = intrace->charge(); // copies
        const size_t ncharge = charge.size();

        for (auto& [mtn, model] : m_models) {
            const auto& spec = model->channel_spectrum(chid);

            if (spec.empty()) {
                continue;       // channel not in model
            }
            // The model spec size may differ than expected nsamples.
            // We could interpolate to correct for that which would
            // slow things down.  Better to correct the model(s) code
            // and configuration.
            if (not warned and spec.size() != m_nsamples) {
                log->warn("model {} produced {} samples instead of expected {}, future warnings muted",
                          mtn, spec.size(), m_nsamples);
                warned = true;
            }

            real_vector_t sigmas(m_nsamples);
            for (size_t ind=0; ind<m_nsamples; ++ind) {
                sigmas[ind] = spec[ind]*sqrt2opi;
            }
            auto wave = rwgen.wave(sigmas);
            wave.resize(ncharge);
            Waveform::increase(charge, wave);
        }

        auto trace = make_shared<SimpleTrace>(chid, intrace->tbin(), charge);
        outtraces.push_back(trace);
    }
    outframe = make_shared<SimpleFrame>(inframe->ident(), inframe->time(), outtraces, inframe->tick());
    log->debug("call={} frame={} {} traces",
               m_count, inframe->ident(), outtraces.size());
    ++m_count;
    return true;
}

bool Gen::CoherentAddNoise::operator()(const input_pointer& inframe, output_pointer& outframe)
{
    if (!inframe) {
        outframe = nullptr;
        log->debug("EOS at call={}", m_count);
        ++m_count;
        return true;
    }

    // We "recycle" the randoms in order to speed up generation.  Each
    // array will start from a random location in the recycled buffer
    // and a few percent will be "freshened".  This results in a small
    // amount of coherency between nearby channels.
    const double BUG=m_bug202; // 0.0 is correct, passing rec. perc. 0.04 opened issue 202
    auto rn = Normals::make_recycling(m_rng, 2*m_nsamples, BUG, 1, 2*m_rep_percent);
    GeneratorN rwgen(m_dft, rn);

    // Look up the generated wave for a group.
    using group_wave_lu = std::unordered_map<int, real_vector_t>;
    // Models may not be coherent across their groups so we have a LU
    // per model.
    std::unordered_map<std::string, group_wave_lu> model_group_waves;

    // Limit number of warnings below
    static bool warned = false;

    const float sqrt2opi = sqrt(2.0/3.141592);

    // Make waveforms of size nsample from each model, adding only
    // ncharge of their element to the trace charge.  This
    // full-nsample followed by ncharge-truncation may CPU-wasteful in
    // the sparse traces case.

    ITrace::vector outtraces;
    for (const auto& intrace : *inframe->traces()) {

        const int chid = intrace->channel();
        auto charge = intrace->charge(); // copies
        const size_t ncharge = charge.size();

        for (auto& [mtn, model] : m_models) {
            auto& gwlu = model_group_waves[mtn];
            int grpid = model->groupid(chid);
            if (gwlu.find(grpid) == gwlu.end()) {
                const auto& spec = model->group_spectrum(grpid);
                
                if (spec.empty()) {
                    continue;       // channel not in model
                }

                // The model spec size may differ than expected nsamples.
                // We could interpolate to correct for that which would
                // slow things down.  Better to correct the model(s) code
                // and configuration.
                if (not warned and spec.size() != m_nsamples) {
                    log->warn("model {} produced {} samples instead of expected {}, future warnings muted",
                              mtn, spec.size(), m_nsamples);
                    warned = true;
                }
                real_vector_t sigmas(m_nsamples);
                for (size_t ind=0; ind<m_nsamples; ++ind) {
                    sigmas[ind] = spec[ind]*sqrt2opi;
                }
                auto wave = rwgen.wave(sigmas);
                wave.resize(ncharge);
                gwlu[grpid] = wave;
            }
            Waveform::increase(charge, gwlu[grpid]);
        }

        auto trace = make_shared<SimpleTrace>(chid, intrace->tbin(), charge);
        outtraces.push_back(trace);
    }
    outframe = make_shared<SimpleFrame>(inframe->ident(), inframe->time(), outtraces, inframe->tick());
    log->debug("call={} frame={} {} traces",
               m_count, inframe->ident(), outtraces.size());
    ++m_count;
    return true;

}

