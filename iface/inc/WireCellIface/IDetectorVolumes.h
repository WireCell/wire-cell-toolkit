/** A view of the basic spatial composition of an entire detector.
    
    The information returned by implementations of this interface MUST adhere to
    the "wire schema" and its convention.  See the util/docs/wire-schema.org.

 */
#ifndef WIRECELL_IDETECTORVOLUMES
#define WIRECELL_IDETECTORVOLUMES

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/Point.h"
#include "WireCellIface/WirePlaneId.h"


namespace WireCell {

    class IDetectorVolumes : public IComponent<IDetectorVolumes> {
    public:
        virtual ~IDetectorVolumes() {}

        /// Return the the wpid for the detector unit that contains the point in
        /// its sensitive volume.  The wpid is evaluates to false if the point
        /// is not contained in any sensitive volume.  The Point must be
        /// provided in the global spatial coordinate system (that which is
        /// defined by "wires").
        virtual WirePlaneId contained_by(const Point& point) const = 0;

        /// Return the sign (+/-1) of the direction along the global X
        /// coordinate direction in which the sensitive face of the detector
        /// volume is "looking".  Note, this is always opposite of the nominal
        /// direction of electron drift toward that face.  If wpid is illegal or
        /// unknown, 0 is returned.
        virtual int face_dirx(WirePlaneId wpid) const = 0;

        /// Return a unit vector along the "wire" (or strip) direction for the
        /// plane.  If wpid is illegal or unknown, a zero vector is returned.
        virtual Vector wire_direction(WirePlaneId wpid) const = 0;

        /// Return a vector with magnitude equal to the pitch between "wires"
        /// (or strips) and with a direction that is perpendicular to the "wire"
        /// direction and along increasing pitch. The "layer" component of the
        /// wpid must be well determined.  If wpid is illegal or unknown, a zero
        /// vector is returned.
        virtual Vector pitch_vector(WirePlaneId wpid) const = 0;

    };

}
#endif
