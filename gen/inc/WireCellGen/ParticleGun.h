#ifndef WIRECELLGEN_PARTICLEGUN
#define WIRECELLGEN_PARTICLEGUN

#include "WireCellIface/IPrimaryVertexSetSource.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IRandom.h"
#include "WireCellAux/Logger.h"
#include "WireCellUtil/Point.h"

#include <string>

namespace WireCell::Gen {

    /** ParticleGun: a source of single-particle primary-vertex events.
     *
     *  Each emitted IPrimaryVertexSet holds one vertex with one particle.  The
     *  particle type, energy, vertex position and direction are configurable,
     *  and each of energy/position may be a fixed scalar or sampled uniformly
     *  from bounds.  Direction is given as a central pair of angles plus an
     *  optional cone half-angle; when the half-angle is > 0 the direction is
     *  sampled uniformly in solid angle (uniform in cos-theta) within the cone.
     *
     *  Configuration:
     *   - particle    : PDG code (int) or canonical name ("electron", "muon",
     *                   "antimuon", ...).
     *   - energy      : scalar, or {min,max} / [min,max] for a uniform sample.
     *   - energy_type : "kinetic" (default), "total", or "momentum" -- how the
     *                   energy value is interpreted.
     *   - position    : [x,y,z]; each component a scalar or {min,max}/[min,max].
     *   - direction   : { theta, phi, max_angle } or { dir:[x,y,z], max_angle }.
     *                   theta/phi (spherical about +z) or the vector "dir" give
     *                   the central direction; max_angle (default 0) is the cone
     *                   half-angle for uniform-in-cos-theta sampling about it.
     *   - count       : number of events to emit (default 1).
     *   - event0      : id of the first event (default 1; must be > 0 since a
     *                   downstream tracker may seed its RNG from it).
     *   - rng         : IRandom component name (default "Random"); only needed
     *                   if any distribution is configured.
     */
    class ParticleGun : public Aux::Logger, public IPrimaryVertexSetSource, public IConfigurable {
       public:
        ParticleGun();
        virtual ~ParticleGun();

        // ISourceNode
        virtual bool operator()(output_pointer& out);

        // IConfigurable
        virtual void configure(const WireCell::Configuration& config);
        virtual WireCell::Configuration default_configuration() const;

       private:
        // A scalar-or-uniform-range value.
        struct Spec {
            bool range{false};
            double lo{0.0};
            double hi{0.0};
        };
        double sample(const Spec& s) const;                 // draw a value
        WireCell::Vector direction() const;                 // draw a unit direction

        std::string m_rng_tn{"Random"};
        IRandom::pointer m_rng;
        bool m_needs_rng{false};

        int m_pdg{13};
        double m_mass{0.0};
        std::string m_name{"muon"};

        Spec m_energy;
        std::string m_energy_type{"kinetic"};
        Spec m_pos[3];

        WireCell::Vector m_dir0{0.0, 0.0, 1.0};  // central unit direction
        double m_max_angle{0.0};                 // cone half-angle (0 => fixed direction)

        int m_number{1};           // events to emit
        int m_event0{1};           // first event id
        int m_emitted{0};
        bool m_eos{false};
    };

}  // namespace WireCell::Gen

#endif
