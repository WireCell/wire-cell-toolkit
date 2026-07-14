#ifndef WIRECELLIFACE_IWIREPARAMETERS
#define WIRECELLIFACE_IWIREPARAMETERS

#include "WireCellIface/WirePlaneId.h"

#include "WireCellUtil/IComponent.h"
#include "WireCellUtil/Point.h"

namespace WireCell {

    // fixme: move this blah blah into the docs.

    /** Interface by which parameters describing wire planes can be
     * accessed.
     *
     * Assumptions and conventions:
     *
     * - Each plane is parallel to the others and are perpendicular to
     * the X-axis.  In decreasing X value they are labeled: U, V and
     * W.
     *
     * - The wires in a plane run parallel to each other, this is
     * called the plane's direction which is also labeled U, V and W.
     *
     * - The angle of the plane is measured from the Y-axis and is
     * signed following the right-hand-rule. (Y x Wire = X).
     *
     * - The perpendicular distances ("pitch") between adjacent wires
     * is regular.  The pitch direction points toward increasing wire
     * index and is take such that (Wire x Pitch) = X
     *
     * The above is sufficient to define the wire planes, but further
     * convention helps to conceptualize things.
     *
     * U wires:
     * - angle: [0, +90]
     * - points: (+Y,+Z)
     * - pitch:  (-Y,+Z)
     * - starts: (-Y,-Z)
     *
     * V wires:
     * - angle: [-90, 0]
     * - points: (+Y,-Z)
     * - pitch:  (-Y,-Z)
     * - starts: (+Y,-Z)
     *
     * W wires:
     * - angle: 0
     * - points: +Y
     * - pitch:  +Z
     * - starts: -Z
     *
     * Given this, the wires planes are fully specified by four rays.
     *
     * \param bounds is a WireCell::Ray with its tail at the
     * negative-most corner and its head at the positive-most
     * corner of a rectangular bounding box.
     *
     * \param U is a WireCell::Ray with direction giving the plane's
     * pitch and with endpoints on wires within the bounding box.
     *
     * \param V is same but same for V wire plane.
     *
     * \param W is same but for for W wire plane.
     *
     */

    // fixme: this doesn't need to be so grand as to be an IComponent,
    // it is just to give an initial test vehicle.
    class IWireParameters : virtual public IComponent<IWireParameters> {
       public:
        virtual ~IWireParameters();

        /** Provide access to the rays which were used to define the wires. */

        /// Defines a box that bounds the set of wires.
        virtual const Ray& bounds() const = 0;

        /// A ray going from the center of the first U wire to the
        /// second and perpendicular to both.
        virtual const Ray& pitchU() const = 0;

        /// A ray going from the center of the first V wire to the
        /// second and perpendicular to both.
        virtual const Ray& pitchV() const = 0;

        /// A ray going from the center of the first W wire to the
        /// second and perpendicular to both.
        virtual const Ray& pitchW() const = 0;

        // helper to return the pitch based on dynamic layer value
        virtual const Ray& pitch(WireCell::WirePlaneLayer_t layer) const
        {
            static Ray bogus;
            switch (layer) {
            case kUlayer:
                return pitchU();
            case kVlayer:
                return pitchV();
            case kWlayer:
                return pitchW();
            default:
                return bogus;
            }
            return bogus;
        }
    };

}  // namespace WireCell

#endif
