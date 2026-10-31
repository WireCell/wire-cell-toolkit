/** A drift-planes provider gives the locations along the drift (X) axis of
 * the three conceptual planes that bound and partition each sensitive anode
 * face:
 *
 * - anode :: bounds the drift region on the side of the wires.
 * - response :: where the field response begins; depos are drifted to it.
 * - cathode :: bounds the drift region opposite the wires.
 *
 * Without a provider these planes are given directly in the configuration of
 * AnodePlane ("faces") and Drifter ("xregions").  A provider lets both take
 * them from one source, such as the geometry of a host application.  Values
 * are returned in exactly those configuration schemas.
 */

#ifndef WIRECELLIFACE_IDRIFTPLANES
#define WIRECELLIFACE_IDRIFTPLANES

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/Configuration.h"

namespace WireCell {

    class IDriftPlanes : public IComponent<IDriftPlanes> {
       public:
        virtual ~IDriftPlanes();

        /// Return the planes of the faces of the anode with the given ident
        /// in the schema of AnodePlane's "faces" parameter: an array indexed
        /// by face, each element null (insensitive face) or an object with
        /// "anode", "response" and "cathode" attributes.
        virtual Configuration faces(int anode_ident) const = 0;

        /// Return the distinct regions in the schema of Drifter's "xregions"
        /// parameter: an array of objects, each as an element of faces().
        virtual Configuration xregions() const = 0;
    };

}  // namespace WireCell

#endif
