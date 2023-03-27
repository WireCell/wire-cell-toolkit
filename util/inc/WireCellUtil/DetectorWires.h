/**
   Provide wire schema generators for known detectors
 */

#ifndef WIRECELL_DETECTORWIRES
#define WIRECELL_DETECTORWIRES

#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

namespace WireCell::DetectorWires
{
    constexpr double pi = 3.14159265;

    struct Plane {
        double pitch;           // pitch magnitude
        double pangle;          // pitch angle w.r.t. z axis in radians
        Point center0;          // center of wire zero
        Ray bounds;             // bounding box
    };

    const std::vector<Plane> uboone = {
        {3*units::mm, -pi/3,
         Point(0, 1173.0155, 2.919594)*units::mm,
         Ray(Point( 0, -1155.1, 0.352608)*units::mm,
             Point( 0,  1174.5, 10369.6)*units::mm)},
        {3*units::mm, +pi/3, 
         Point(-3, -1153.615, 2.919594)*units::mm,
         Ray(Point(-3, -1155.1, 0.352608)*units::mm,
             Point(-3,  1174.5, 10369.6)*units::mm)},
        {3*units::mm, 0.0,
         Point(-6, 9.7, 2.5)*units::mm,
         Ray(Point(-6, -1155.3, 2.5)*units::mm,
             Point(-6,  1174.7, 10367.5)*units::mm)},
    };


    // Fill one face of wires given pitch magnitude, pitch angles
    // w.r.t Z axis and wire0 centers.
    size_t plane(WireSchema::StoreDB& storedb,
                 WireSchema::Plane& plane,
                 const Plane& param);

    
    
    
    // Set channels of one planes wires by counting wires in plane.
    void flat_channels(WireSchema::StoreDB& storedb, WireSchema::Plane& plane, int ch0=0);

    // Set channels by counting their wires in plane in hierarchical
    // order.
    void flat_channels(WireSchema::StoreDB& storedb, int ch0=0);

}

#endif
