#ifndef WIRECELLIFACE_ISCEFIELD
#define WIRECELLIFACE_ISCEFIELD

#include "WireCellUtil/IComponent.h"

namespace WireCell {

    /** Per-TPC backward (reco -> true) Space-Charge-Effect displacement field.
     *
     *  Implementations expose, for a given APA/TPC index and 3D position,
     *  the SIGNED backward displacement to be added directly to a t0-corrected
     *  coordinate to recover the true coordinate.  Out-of-domain queries are
     *  clamped internally by the implementation.  Implementations that carry no
     *  map should return 0 (no-op safe).
     *
     *  Units: WCT internal (mm) in, (mm) out.  Implementations handle any
     *  conversion required by their underlying map.
     *
     *  displacement_x is mandatory.  displacement_y / displacement_z are
     *  optional: implementations whose underlying map lacks transverse
     *  components inherit the default no-op (return 0).
     */
    class ISCEField : public IComponent<ISCEField> {
       public:
        virtual ~ISCEField() = default;

        /// Backward X-displacement at (x,y,z) for APA index `apa`.
        /// SBND convention: apa=0 East, apa=1 West.  Units: mm.
        virtual double displacement_x(int apa, double x, double y, double z) const = 0;

        /// Backward Y-displacement at (x,y,z).  Units: mm.  Optional (default no-op).
        virtual double displacement_y(int /*apa*/, double /*x*/, double /*y*/, double /*z*/) const { return 0.0; }

        /// Backward Z-displacement at (x,y,z).  Units: mm.  Optional (default no-op).
        virtual double displacement_z(int /*apa*/, double /*x*/, double /*y*/, double /*z*/) const { return 0.0; }
    };

}  // namespace WireCell

#endif
