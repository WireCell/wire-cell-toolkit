#ifndef WIRECELL_IPRIMARYPARTICLE
#define WIRECELL_IPRIMARYPARTICLE

#include "WireCellIface/IData.h"
#include "WireCellUtil/Point.h"          // Vector
#include "WireCellUtil/Configuration.h"  // Configuration (rarely-used attrs)

namespace WireCell {

    /** An interface to a primary particle: one particle handed to a
     *  tracking simulation (eg Geant4) to be tracked, as emitted from an
     *  IPrimaryVertex.
     *
     *  The interface models the commonly-used core (PDG code, momentum,
     *  energy and a track id linking to the truth trajectory record) as
     *  first-class accessors, striking a balance between edep-sim's
     *  TG4PrimaryParticle (PDG + 4-momentum only) and Geant4's much richer
     *  G4PrimaryParticle.  The less-common Geant4 attributes are carried,
     *  when a producer provides them, in the metadata() Configuration:
     *
     *    "name"         (string)  particle name
     *    "mass"         (double)  rest mass, if not implied by pdg()
     *    "charge"       (double)  charge, if not implied by pdg()
     *    "polarization" ([x,y,z]) polarization vector
     *    "proper_time"  (double)  pre-tracking proper time
     *    "weight"       (double)  statistical weight
     *
     *  (Pre-assigned decay daughters are not modeled here yet.)
     *
     *  Values are expressed in the WCT system of units.
     */
    class IPrimaryParticle : public IData<IPrimaryParticle> {
       public:
        virtual ~IPrimaryParticle();

        /// The PDG particle code.
        virtual int pdg() const = 0;

        /// The initial 3-momentum (px, py, pz).  Returned by value so a facade
        /// implementation may synthesize it from its backing store without
        /// having to hold a WireCell::Vector member.
        virtual Vector momentum() const = 0;

        /// The total (kinetic + rest) energy.  The rest mass is
        /// sqrt(energy^2 - |momentum|^2).
        virtual double energy() const = 0;

        /// The track id linking this primary to a trajectory in the truth
        /// record.  Negative if the particle is not tracked or the producer
        /// provides no association.
        virtual int id() const = 0;

        /// Less-common producer-provided attributes (see class docs).  An
        /// empty (null) Configuration if the producer provides none.
        virtual Configuration metadata() const { return Configuration(); }
    };

}  // namespace WireCell

#endif
