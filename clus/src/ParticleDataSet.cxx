#include "WireCellClus/ParticleDataSet.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(ParticleDataSet, WireCell::CLUS::ParticleDataSet, WireCell::IConfigurable)

using namespace WireCell;

CLUS::ParticleDataSet::ParticleDataSet() {}
CLUS::ParticleDataSet::~ParticleDataSet() {}

void CLUS::ParticleDataSet::configure(const WireCell::Configuration& config) {
    // Configure dE/dx functions
    if (config.isMember("dedx_functions")) {
        const auto& dedx_config = config["dedx_functions"];
        for (const auto& particle : dedx_config.getMemberNames()) {
            const auto& func_name = dedx_config[particle].asString();
            auto func = Factory::find_tn<IScalarFunction>(func_name);
            if (func) {
                m_dedx_functions[particle] = func;
            }
        }
    }
    
    // Configure range functions  
    if (config.isMember("range_functions")) {
        const auto& range_config = config["range_functions"];
        for (const auto& particle : range_config.getMemberNames()) {
            const auto& func_name = range_config[particle].asString();
            auto func = Factory::find_tn<IScalarFunction>(func_name);
            if (func) {
                m_range_functions[particle] = func;
            }
        }
    }
}

WireCell::Configuration CLUS::ParticleDataSet::default_configuration() const {
    Configuration cfg;
    cfg["dedx_functions"] = Json::Value(Json::objectValue);
    cfg["range_functions"] = Json::Value(Json::objectValue);
    return cfg;
}

IScalarFunction::pointer CLUS::ParticleDataSet::get_dEdx_function(const std::string& particle) const {
    auto it = m_dedx_functions.find(particle);
    return (it != m_dedx_functions.end()) ? it->second : nullptr;
}

IScalarFunction::pointer CLUS::ParticleDataSet::get_range_function(const std::string& particle) const {
    auto it = m_range_functions.find(particle);
    return (it != m_range_functions.end()) ? it->second : nullptr;
}

std::vector<std::string> CLUS::ParticleDataSet::get_particles() const {
    std::vector<std::string> particles;
    for (const auto& pair : m_dedx_functions) {
        particles.push_back(pair.first);
    }
    return particles;
}