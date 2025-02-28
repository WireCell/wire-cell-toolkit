/** Determine if a point is inside a closed 3D region */
#ifndef WIRECELL_IFIDUCIAL
#define WIRECELL_IFIDUCIAL

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/Point.h"


namespace WireCell {

    class IFiducial : public IComponent<IFiducial> {
    public:
        virtual ~IFiducial() {}

        /// Return true if the point is inside the region, else false.  The
        /// region is considered closed.  That is, a point on the region
        /// boundary is inside.
        virtual bool contained(const Point& point) const = 0;
    };

}
#endif
