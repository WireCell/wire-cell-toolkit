#ifndef WIRECELLAUX_PARTICLEINFO
#define WIRECELLAUX_PARTICLEINFO

#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/D4Vector.h"
#include <string>
#include <map>

namespace WireCell::Aux {

    /**
     * @brief A comprehensive particle information storage class
     * 
     * This class stores all relevant particle physics information including
     * PDG codes, 4-momentum, mass, name, and kinetic energy. It follows
     * Wire-Cell toolkit conventions and integrates well with existing
     * deposition and clustering infrastructure.
     */
    class ParticleInfo {
    public:
        // Default constructor - creates an invalid/empty particle
        ParticleInfo();
        
        // Comprehensive constructor
        ParticleInfo(int pdg_code, 
                    double mass,
                    const std::string& name,
                    double kinetic_energy,
                    const WireCell::Point& momentum_3vec,
                    int id = 0,
                    double charge = 0.0);
                    
        // Constructor from 4-momentum
        ParticleInfo(int pdg_code,
                    double mass, 
                    const std::string& name,
                    const WireCell::D4Vector<double>& four_momentum,
                    int id = 0,
                    double charge = 0.0);

        virtual ~ParticleInfo() = default;

        // Core particle properties
        int pdg() const { return m_pdg_code; }
        double mass() const { return m_mass; }
        const std::string& name() const { return m_name; }
        double kinetic_energy() const { return m_kinetic_energy; }
        int id() const { return m_id; }
        double charge() const { return m_charge; }

        // 4-momentum access
        double energy() const { return m_four_momentum.e(); }
        WireCell::Point momentum() const { return WireCell::Point(m_four_momentum.px(), m_four_momentum.py(), m_four_momentum.pz()); }
        const WireCell::D4Vector<double>& four_momentum() const { return m_four_momentum; }
        
        // Derived quantities
        double momentum_magnitude() const { return m_four_momentum.p(); }
        double beta() const { return m_four_momentum.beta(); }           // v/c
        double gamma() const { return m_four_momentum.gamma(); }          // Lorentz factor
        double rapidity() const { return m_four_momentum.rapidity(); }       // 0.5 * ln((E+pz)/(E-pz))
        
        // Utility methods
        bool is_valid() const { return m_pdg_code != 0; }
        bool is_charged() const { return std::abs(m_charge) > 1e-6; }
        bool is_stable() const;        // Based on PDG code
        
        // Setters (for cases where you need to modify after construction)
        void set_pdg(int pdg) { m_pdg_code = pdg; }
        void set_mass(double mass) { m_mass = mass; update_kinematics(); }
        void set_name(const std::string& name) { m_name = name; }
        void set_momentum(const WireCell::Point& momentum_3vec);
        void set_four_momentum(const WireCell::D4Vector<double>& four_momentum);
        void set_kinetic_energy(double ke);
        void set_id(int id) { m_id = id; }
        void set_charge(double charge) { m_charge = charge; }

        // Static utility methods for PDG lookups
        static std::string pdg_to_name(int pdg_code);
        static double pdg_to_mass(int pdg_code);
        static double pdg_to_charge(int pdg_code);
        
        // Factory method to create from PDG code
        static ParticleInfo from_pdg(int pdg_code, 
                                   const WireCell::Point& momentum_3vec,
                                   int id = 0);

    private:
        // Core data members
        int m_pdg_code;
        double m_mass;
        std::string m_name;
        double m_kinetic_energy;
        WireCell::D4Vector<double> m_four_momentum;  // (E, px, py, pz)
        int m_id;
        double m_charge;

        // Internal helper methods
        void update_kinematics();     // Recalculate energy/momentum relationships
        void validate_inputs();       // Check for physical consistency
        
        // Static data for PDG lookups
        static const std::map<int, std::string>& get_pdg_name_map();
        static const std::map<int, double>& get_pdg_mass_map();
        static const std::map<int, double>& get_pdg_charge_map();
    };

    // Convenience typedefs
    using ParticleInfoPtr = std::shared_ptr<ParticleInfo>;
    using ParticleInfoVector = std::vector<ParticleInfo>;
    using ParticleInfoSelection = std::vector<ParticleInfo*>;

    // Helper functions for common particle types
    namespace ParticleHelpers {
        ParticleInfo electron(const WireCell::Point& momentum, int id = 0);
        ParticleInfo muon(const WireCell::Point& momentum, int charge_sign = 1, int id = 0);
        ParticleInfo pion_charged(const WireCell::Point& momentum, int charge_sign = 1, int id = 0);
        ParticleInfo proton(const WireCell::Point& momentum, int id = 0);
        ParticleInfo photon(const WireCell::Point& momentum, int id = 0);
    }

} // namespace WireCell::Aux

#endif // WIRECELLAUX_PARTICLEINFO