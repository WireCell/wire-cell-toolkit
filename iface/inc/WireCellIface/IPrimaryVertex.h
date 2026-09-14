#ifndef WIRECELL_IPRIMARYVERTEX
#define WIRECELL_IPRIMARYVERTEX

#include "WireCellIface/IData.h"
#include "WireCellIface/IPrimaryParticle.h"
#include "WireCellUtil/Point.h"          // Point
#include "WireCellUtil/Configuration.h"  // Configuration (rarely-used attrs)

namespace WireCell {

    /** An interface to a primary vertex: a point in space and time from
     *  which one or more IPrimaryParticle are emitted into a tracking
     *  simulation.
     *
     *  As with IPrimaryParticle, the interface models the common core
     *  (position, time and the emitted particles) directly and defers the
     *  less-common, generator-provenance attributes to metadata().  This
     *  strikes a balance between edep-sim's TG4PrimaryVertex and Geant4's
     *  G4PrimaryVertex.  Recognized metadata() keys (all optional):
     *
     *    "generator"          (string)  name of the generator
     *    "reaction"           (string)  reaction that made the vertex
     *    "filename"           (string)  input kinematics file
     *    "interaction_number" (int)     index in the kinematics file
     *    "cross_section"      (double)  cross section
     *    "diff_cross_section" (double)  differential cross section
     *    "weight"             (double)  interaction weight
     *    "probability"        (double)  overall interaction probability
     *
     *  (edep-sim's nested "informational" vertices are not modeled here yet.)
     *
     *  Values are expressed in the WCT system of units.
     */
    class IPrimaryVertex : public IData<IPrimaryVertex> {
       public:
        virtual ~IPrimaryVertex();

        /// The position of the vertex.  Returned by value so a facade
        /// implementation may synthesize it from its backing store without
        /// having to hold a WireCell::Point member.
        virtual Point position() const = 0;

        /// The time of the vertex.
        virtual double time() const = 0;

        /// The primary particles emitted from this vertex.
        virtual IPrimaryParticle::shared_vector particles() const = 0;

        /// Less-common producer-provided attributes (see class docs).  An
        /// empty (null) Configuration if the producer provides none.
        virtual Configuration metadata() const { return Configuration(); }
    };

}  // namespace WireCell

#endif
