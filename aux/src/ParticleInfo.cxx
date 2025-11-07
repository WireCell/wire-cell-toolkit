#include "WireCellAux/ParticleInfo.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Units.h"

#include <cmath>
#include <stdexcept>

using namespace WireCell;
using namespace WireCell::Aux;

// Default constructor
ParticleInfo::ParticleInfo()
    : m_pdg_code(0)
    , m_mass(0.0)
    , m_name("unknown")
    , m_kinetic_energy(0.0)
    , m_four_momentum(0.0, 0.0, 0.0, 0.0)
    , m_id(0)
    , m_charge(0.0)
{
}

// Comprehensive constructor
ParticleInfo::ParticleInfo(int pdg_code, 
                          double mass,
                          const std::string& name,
                          double kinetic_energy,
                          const WireCell::Point& momentum_3vec,
                          int id,
                          double charge)
    : m_pdg_code(pdg_code)
    , m_mass(mass)
    , m_name(name)
    , m_kinetic_energy(kinetic_energy)
    , m_id(id)
    , m_charge(charge)
{
    // Calculate total energy from kinetic energy and mass
    double total_energy = m_kinetic_energy + m_mass;
    m_four_momentum.set(total_energy, momentum_3vec.x(), momentum_3vec.y(), momentum_3vec.z());
    validate_inputs();
}

// Constructor from 4-momentum
ParticleInfo::ParticleInfo(int pdg_code,
                          double mass,
                          const std::string& name,
                          const WireCell::D4Vector<double>& four_momentum,
                          int id,
                          double charge)
    : m_pdg_code(pdg_code)
    , m_mass(mass)
    , m_name(name)
    , m_four_momentum(four_momentum)
    , m_id(id)
    , m_charge(charge)
{
    // Calculate kinetic energy from total energy and mass
    m_kinetic_energy = m_four_momentum.e() - m_mass;
    validate_inputs();
}

bool ParticleInfo::is_stable() const {
    // Common stable particles by PDG code
    static const std::set<int> stable_pdgs = {
        11, -11,    // electron, positron
        12, -12,    // electron neutrino, anti-electron neutrino
        13, -13,    // muon, anti-muon (long-lived)
        14, -14,    // muon neutrino, anti-muon neutrino
        16, -16,    // tau neutrino, anti-tau neutrino
        22,         // photon
        2112,       // neutron (relatively stable)
        2212,       // proton
        -2212       // anti-proton
    };
    return stable_pdgs.find(std::abs(m_pdg_code)) != stable_pdgs.end();
}

// Setters with kinematic updates
void ParticleInfo::set_momentum(const WireCell::Point& momentum_3vec) {
    m_four_momentum.set(m_four_momentum.e(), momentum_3vec.x(), momentum_3vec.y(), momentum_3vec.z());
    update_kinematics();
}

void ParticleInfo::set_four_momentum(const WireCell::D4Vector<double>& four_momentum) {
    m_four_momentum = four_momentum;
    m_kinetic_energy = m_four_momentum.e() - m_mass;
}

void ParticleInfo::set_kinetic_energy(double ke) {
    m_kinetic_energy = ke;
    double total_energy = m_kinetic_energy + m_mass;
    
    // Recalculate momentum magnitude if needed
    double p_mag = std::sqrt(total_energy * total_energy - m_mass * m_mass);
    
    if (momentum_magnitude() > 0) {
        // Scale existing momentum direction
        WireCell::Point current_momentum = momentum();
        double current_p_mag = momentum_magnitude();
        double scale = p_mag / current_p_mag;
        m_four_momentum.set(total_energy, 
                           current_momentum.x() * scale,
                           current_momentum.y() * scale, 
                           current_momentum.z() * scale);
    } else {
        // Set momentum along z-axis if no direction exists
        m_four_momentum.set(total_energy, 0, 0, p_mag);
    }
}

// Private helper methods
void ParticleInfo::update_kinematics() {
    double p_mag = momentum_magnitude();
    double total_energy = std::sqrt(p_mag * p_mag + m_mass * m_mass);
    m_four_momentum.e(total_energy);
    m_kinetic_energy = total_energy - m_mass;
}

void ParticleInfo::validate_inputs() {
    if (m_mass < 0.0) {
        raise<ValueError>("ParticleInfo: mass cannot be negative");
    }
    if (m_four_momentum.e() < m_mass) {
        raise<ValueError>("ParticleInfo: total energy cannot be less than rest mass");
    }
    if (m_kinetic_energy < 0.0) {
        raise<ValueError>("ParticleInfo: kinetic energy cannot be negative");
    }
    
    // Check energy-momentum relation using D4Vector's mass calculation
    double calculated_mass = m_four_momentum.mass();
    if (std::abs(calculated_mass - m_mass) > 1e-6 * m_mass) {
        raise<ValueError>("ParticleInfo: energy-momentum relation violated");
    }
}

// Static utility methods for PDG lookups
std::string ParticleInfo::pdg_to_name(int pdg_code) {
    const auto& name_map = get_pdg_name_map();
    auto it = name_map.find(std::abs(pdg_code));
    if (it != name_map.end()) {
        return pdg_code < 0 ? "anti-" + it->second : it->second;
    }
    return "unknown";
}

double ParticleInfo::pdg_to_mass(int pdg_code) {
    const auto& mass_map = get_pdg_mass_map();
    auto it = mass_map.find(std::abs(pdg_code));
    return it != mass_map.end() ? it->second : 0.0;
}

double ParticleInfo::pdg_to_charge(int pdg_code) {
    const auto& charge_map = get_pdg_charge_map();
    auto it = charge_map.find(pdg_code);  // Don't use abs() here as charge depends on sign
    return it != charge_map.end() ? it->second : 0.0;
}

ParticleInfo ParticleInfo::from_pdg(int pdg_code, 
                                   const WireCell::Point& momentum_3vec,
                                   int id) {
    return ParticleInfo(pdg_code,
                       pdg_to_mass(pdg_code),
                       pdg_to_name(pdg_code),
                       0.0,  // KE will be calculated
                       momentum_3vec,
                       id,
                       pdg_to_charge(pdg_code));
}

// Static data for PDG lookups
const std::map<int, std::string>& ParticleInfo::get_pdg_name_map() {
    static const std::map<int, std::string> pdg_names = {
        {11, "electron"},
        {12, "electron_neutrino"},
        {13, "muon"},
        {14, "muon_neutrino"},
        {15, "tau"},
        {16, "tau_neutrino"},
        {22, "photon"},
        {111, "pi0"},
        {211, "pi_plus"},
        {311, "K0"},
        {321, "K_plus"},
        {2112, "neutron"},
        {2212, "proton"},
        {3122, "lambda"},
        {3222, "sigma_plus"},
        {3112, "sigma_minus"},
        {3322, "xi_minus"},
        {3334, "omega_minus"}
    };
    return pdg_names;
}

const std::map<int, double>& ParticleInfo::get_pdg_mass_map() {
    static const std::map<int, double> pdg_masses = {
        {11, 0.511 * units::MeV},     // electron
        {12, 0.0},                    // electron neutrino (massless)
        {13, 105.658 * units::MeV},   // muon
        {14, 0.0},                    // muon neutrino (massless)
        {15, 1776.86 * units::MeV},   // tau
        {16, 0.0},                    // tau neutrino (massless)
        {22, 0.0},                    // photon
        {111, 134.977 * units::MeV},  // pi0
        {211, 139.570 * units::MeV},  // pi_plus
        {311, 497.648 * units::MeV},  // K0
        {321, 493.677 * units::MeV},  // K_plus
        {2112, 939.565 * units::MeV}, // neutron
        {2212, 938.272 * units::MeV}, // proton
        {3122, 1115.683 * units::MeV}, // lambda
        {3222, 1189.37 * units::MeV},  // sigma_plus
        {3112, 1197.449 * units::MeV}, // sigma_minus
        {3322, 1321.71 * units::MeV},  // xi_minus
        {3334, 1672.45 * units::MeV}   // omega_minus
    };
    return pdg_masses;
}

const std::map<int, double>& ParticleInfo::get_pdg_charge_map() {
    static const std::map<int, double> pdg_charges = {
        {11, -1.0},   // electron
        {-11, 1.0},   // positron
        {12, 0.0},    // electron neutrino
        {-12, 0.0},   // anti-electron neutrino
        {13, -1.0},   // muon
        {-13, 1.0},   // anti-muon
        {14, 0.0},    // muon neutrino
        {-14, 0.0},   // anti-muon neutrino
        {15, -1.0},   // tau
        {-15, 1.0},   // anti-tau
        {16, 0.0},    // tau neutrino
        {-16, 0.0},   // anti-tau neutrino
        {22, 0.0},    // photon
        {111, 0.0},   // pi0
        {211, 1.0},   // pi_plus
        {-211, -1.0}, // pi_minus
        {311, 0.0},   // K0
        {321, 1.0},   // K_plus
        {-321, -1.0}, // K_minus
        {2112, 0.0},  // neutron
        {2212, 1.0},  // proton
        {-2212, -1.0}, // anti-proton
        {3122, 0.0},  // lambda
        {3222, 1.0},  // sigma_plus
        {3112, -1.0}, // sigma_minus
        {3322, -1.0}, // xi_minus
        {3334, -1.0}  // omega_minus
    };
    return pdg_charges;
}

// Helper functions for common particle types
namespace WireCell::Aux::ParticleHelpers {

    ParticleInfo electron(const WireCell::Point& momentum, int id) {
        return ParticleInfo::from_pdg(11, momentum, id);
    }

    ParticleInfo muon(const WireCell::Point& momentum, int charge_sign, int id) {
        return ParticleInfo::from_pdg(charge_sign < 0 ? 13 : -13, momentum, id);
    }

    ParticleInfo pion_charged(const WireCell::Point& momentum, int charge_sign, int id) {
        return ParticleInfo::from_pdg(charge_sign > 0 ? 211 : -211, momentum, id);
    }

    ParticleInfo proton(const WireCell::Point& momentum, int id) {
        return ParticleInfo::from_pdg(2212, momentum, id);
    }

    ParticleInfo photon(const WireCell::Point& momentum, int id) {
        return ParticleInfo::from_pdg(22, momentum, id);
    }

}