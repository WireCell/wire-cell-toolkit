#ifndef WIRECELL_IFACE_IYZMAP
#define WIRECELL_IFACE_IYZMAP

#include "WireCellUtil/IComponent.h"

#include <string>

namespace WireCell {
    class IYZMap : public IComponent<IYZMap> {
    public:
        virtual ~IYZMap();

        // Return map value for the given anode, plane, and detector (y, z)
        // coordinates in WireCell units.  Handles binning and boundary
        // clamping internally.
        virtual double value(const std::string& anode_name,
                             int plane,
                             double y, double z) const = 0;
    };
}

#endif
