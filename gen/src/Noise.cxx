#include "WireCellGen/Noise.h"

#include "WireCellAux/DftTools.h"
#include "WireCellAux/RandTools.h"

#include "WireCellUtil/NamedFactory.h"

#include <unordered_set>


using namespace std;
using namespace WireCell;
using namespace WireCell::Aux::RandTools;

Gen::NoiseBase::NoiseBase(const std::string& tn)
    : Aux::Logger(tn, "gen")
{
}
Gen::NoiseBase::~NoiseBase() 
{
}

WireCell::Configuration Gen::NoiseBase::default_configuration() const
{
    Configuration cfg;
    // The name of the noise model.  Or, if an array is given, the
    // names of multiple noise models.  The plural spelling "models"
    // is also checked for an array of model names.
    cfg["model"] = "";
    // The IRandom
    cfg["rng"] = "Random";
    // The IDFT
    cfg["dft"] = "FftwDFT";
    // The size (number of ticks) of the noise waveforms
    cfg["nsamples"] = m_nsamples;
    // Percentage of fresh randomness added in using the recyled
    // random normal-Gaussian distribution.  Note, to match prior
    // interpreation of "percent" we apply twice this amount of
    // refreshing.  0.02 remains a reasonable choice.
    cfg["replacement_percentage"] = m_rep_percent;
    return cfg;
}

void Gen::NoiseBase::configure(const WireCell::Configuration& cfg)
{
    std::string rng_tn = get<std::string>(cfg, "rng", "Random");
    m_rng = Factory::find_tn<IRandom>(rng_tn);

    std::string dft_tn = get<std::string>(cfg, "dft", "FftwDFT");
    m_dft = Factory::find_tn<IDFT>(dft_tn);

    // Specify a single model name or an array of model names.  We
    // also accept an array given as a pluralized "models".
    m_model_tns.clear();
    auto jmodel = cfg["model"];
    if (jmodel.isString()) {
        m_model_tns.insert(jmodel.asString());
    }
    else {
        for (const auto& jmod : jmodel) {
            m_model_tns.insert(jmod.asString());
        }
    }
    if (cfg["models"].isArray()) {
        for (const auto& jmod : cfg["models"]) {
            m_model_tns.insert(jmod.asString());
        }
    }
    for (const auto& mtn : m_model_tns) {
        log->debug("using noise model: \"{}\"", mtn);
    }

    m_nsamples = get<int>(cfg, "nsamples", m_nsamples);
    m_rep_percent = get<double>(cfg, "replacement_percentage", m_rep_percent);
    m_bug202 = get<double>(cfg, "bug202", 0.0);
    if (m_bug202>0) {
        log->warn("BUG 202 IS ACTIVATED: {}", m_bug202);
    }
}

