#ifndef WIRECELLAUX_SIMPLEPRIMARYPARTICLE
#define WIRECELLAUX_SIMPLEPRIMARYPARTICLE

#include "WireCellIface/IPrimaryParticle.h"

namespace WireCell::Aux {

    // A primary particle that simply holds all the data it presents.
    class SimplePrimaryParticle : public WireCell::IPrimaryParticle {
       public:
        SimplePrimaryParticle(int pdg,
                              const WireCell::Vector& momentum,
                              double energy,
                              int id = -1,
                              const WireCell::Configuration& metadata = WireCell::Configuration());
        virtual ~SimplePrimaryParticle();

        virtual int pdg() const;
        virtual WireCell::Vector momentum() const;
        virtual double energy() const;
        virtual int id() const;
        virtual WireCell::Configuration metadata() const;

       private:
        int m_pdg;
        WireCell::Vector m_momentum;
        double m_energy;
        int m_id;
        WireCell::Configuration m_metadata;
    };

}  // namespace WireCell::Aux

#endif
