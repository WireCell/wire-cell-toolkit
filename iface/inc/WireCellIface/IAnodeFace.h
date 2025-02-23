/**
   Provides information about an "anode face" which consists of a
   number of parallel wire planes, each consisting of a number of
   parallel wires.

   Information includes:

   - Internal geometry of the wires in their planes and how they
   relate to a coordinate system of the larger enclosing volume.

   - field response information

   - wire/channel connectivity and numbering


 */

#ifndef WIRECELLIFACES_IANODEFACE
#define WIRECELLIFACES_IANODEFACE

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/BoundingBox.h"
#include "WireCellUtil/RayGrid.h"
#include "WireCellIface/IWirePlane.h"

namespace WireCell {

    class IAnodeFace : public IComponent<IAnodeFace> {
       public:
        virtual ~IAnodeFace();

        /// Return an application-specific, externally defined identifier number
        /// for this face.  WCT should NOT interpret this number.  If you find a
        /// place where it does, please file a bug.
        virtual int ident() const = 0;

        /// This returns an index of the face in the context of the anode.  The
        /// order is typically determined by user configuration.  There is NO
        /// EXPECTATION for this number to relate to face direction.  If you
        /// find a place where this is not true, please file a bug.  See dirx().
        virtual int which() const = 0;

        /// Return +1 when the face-normal vector points in the positive
        /// global-X direction and -1 when pointing in the negative global-X
        /// direction.  Note, in all cases, the face-normal points opposite to
        /// the drift direction.
        virtual int dirx() const = 0;

        /// Return the anode plane ident number in which this face
        /// resides.
        virtual int anode() const = 0;

        /// Return the number of wire planes in the given side
        virtual int nplanes() const = 0;

        /// Return the wire plane with the given ident or nullptr if unknown.
        virtual IWirePlane::pointer plane(int ident) const = 0;

        /// Return all wires planes
        virtual IWirePlane::vector planes() const = 0;

        /// Return a bounding box containing the volume to which this
        /// face is sensitive.  BB may be zero volume.
        virtual BoundingBox sensitive() const = 0;

        /// Return a RayGrid::Coordinates corresponding to this face.
        virtual const RayGrid::Coordinates& raygrid() const = 0;
    };
}  // namespace WireCell

#endif
