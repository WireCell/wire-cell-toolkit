#ifndef WIRECELLAUX_SIMPLEPRIMARYVERTEX
#define WIRECELLAUX_SIMPLEPRIMARYVERTEX

#include "WireCellIface/IPrimaryVertex.h"

namespace WireCell::Aux {

    // A primary vertex that simply holds all the data it presents.
    class SimplePrimaryVertex : public WireCell::IPrimaryVertex {
       public:
        SimplePrimaryVertex(const WireCell::Point& position,
                            double time,
                            const WireCell::IPrimaryParticle::vector& particles = {},
                            const WireCell::Configuration& metadata = WireCell::Configuration());
        virtual ~SimplePrimaryVertex();

        virtual WireCell::Point position() const;
        virtual double time() const;
        virtual WireCell::IPrimaryParticle::shared_vector particles() const;
        virtual WireCell::Configuration metadata() const;

       private:
        WireCell::Point m_position;
        double m_time;
        WireCell::IPrimaryParticle::shared_vector m_particles;
        WireCell::Configuration m_metadata;
    };

}  // namespace WireCell::Aux

#endif
