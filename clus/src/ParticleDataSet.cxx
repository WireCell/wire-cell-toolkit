#include "WireCellClus/ParticleDataSet.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"

WIRECELL_FACTORY(ParticleDataSet, WireCell::Clus::ParticleDataSet, WireCell::IConfigurable)

using namespace WireCell;

Clus::ParticleDataSet::ParticleDataSet() {}
Clus::ParticleDataSet::~ParticleDataSet() {}

void Clus::ParticleDataSet::configure(const WireCell::Configuration& config) {
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

WireCell::Configuration Clus::ParticleDataSet::default_configuration() const {
    Configuration cfg;
    cfg["dedx_functions"] = Json::Value(Json::objectValue);
    cfg["range_functions"] = Json::Value(Json::objectValue);
    return cfg;
}

IScalarFunction::pointer Clus::ParticleDataSet::get_dEdx_function(const std::string& particle) const {
    auto it = m_dedx_functions.find(particle);
    return (it != m_dedx_functions.end()) ? it->second : nullptr;
}

IScalarFunction::pointer Clus::ParticleDataSet::get_range_function(const std::string& particle) const {
    auto it = m_range_functions.find(particle);
    return (it != m_range_functions.end()) ? it->second : nullptr;
}

std::vector<std::string> Clus::ParticleDataSet::get_particles() const {
    std::vector<std::string> particles;
    for (const auto& pair : m_dedx_functions) {
        particles.push_back(pair.first);
    }
    return particles;
}

double Clus::ParticleDataSet::get_particle_mass(int pdg_code) const {
    // Particle Data Group (PDG) codes and their masses in MeV/c^2
    static const std::map<int, double> pdg_mass_map = {
        {11, 0.5109989461},    // electron
        {-11, 0.5109989461},   // positron
        {13, 105.6583745},     // muon
        {-13, 105.6583745},    // anti-muon
        {211, 139.57039},      // charged pion
        {-211, 139.57039},     // charged pion
        {321, 493.677},        // charged kaon
        {-321, 493.677},       // charged kaon
        {2212, 938.2720813},   // proton
        {-2212, 938.2720813}   // anti-proton
    };
    
    auto it = pdg_mass_map.find(pdg_code);
    if (it != pdg_mass_map.end()) {
        return it->second * units::MeV; // Convert to internal units (MeV/c^2)
    } else {
        return 0.0; // Unknown particle type
    }

}
