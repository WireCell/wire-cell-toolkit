#ifndef WIRECELLCLUS_PARTICLEDATASET
#define WIRECELLCLUS_PARTICLEDATASET

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IScalarFunction.h"
#include "WireCellUtil/NamedFactory.h"
#include <map>
#include <string>

namespace WireCell {
    namespace Clus {
        
        class ParticleDataSet : public IConfigurable {
        public:
            typedef std::shared_ptr<ParticleDataSet> pointer;
            
            ParticleDataSet();
            virtual ~ParticleDataSet();
            
            virtual void configure(const WireCell::Configuration& config);
            virtual WireCell::Configuration default_configuration() const;
            
            // Access functions
            IScalarFunction::pointer get_dEdx_function(const std::string& particle) const;
            IScalarFunction::pointer get_range_function(const std::string& particle) const;
            double get_particle_mass(int pdg_code) const;
            std::string pdg_to_name(int pdg_code) const;
            
            // Get available particles
            std::vector<std::string> get_particles() const;
            
        private:
            std::map<std::string, IScalarFunction::pointer> m_dedx_functions;
            std::map<std::string, IScalarFunction::pointer> m_range_functions;
        };
        
    }
}

#endif