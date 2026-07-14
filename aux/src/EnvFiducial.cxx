#include "WireCellIface/IFiducial.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/BoundingBox.h"

#include <vector>
//#include <iostream>             // only debug

// Implementation is totally local to this comp. unit so no need for namespacing.
class EnvFiducial;

WIRECELL_FACTORY(EnvFiducial, EnvFiducial,
                 WireCell::IFiducial, WireCell::IConfigurable)


using namespace WireCell;

/// The EnvFiducial is a simple bounding box that contains all anode sensitive
/// volumes.  A point which is inside the gaps between these volumes is inside
/// this fiducial.  If you want to exclude the gaps, use the IFiducial interface
/// of the DetectorVolumes component.
class EnvFiducial : public IFiducial, public IConfigurable {

    BoundingBox m_bb;

public:

    EnvFiducial() {}
    virtual ~EnvFiducial() {}

    virtual Configuration default_configuration() const {
        Configuration cfg;
        // A list of IAnodePlane "type:name" identifiers to match up to wpids. 
        cfg["anodes"] = Json::arrayValue;
        return cfg;
    }

    virtual void configure(const Configuration& cfg) {
        const auto& janodes = cfg["anodes"];

        if (janodes.empty()) {
            raise<ValueError>("EnvFiducial 'anodes' list is empty, check your configuration");
        }                              

        m_bb = BoundingBox();
        for (const auto& janode : cfg["anodes"]) {
            const std::string anode_tn = janode.asString();
            auto ianode = Factory::find_tn<IAnodePlane>(anode_tn);
            for (auto iface : ianode->faces()) {
                m_bb(iface->sensitive().bounds());
            }
        }
        // std::cerr << "EnvFiducial bounds: " << m_bb.bounds() << "\n";
    }

    // IFiducial
    virtual bool contained(const Point& point) const {
        return m_bb.inside(point);
    }
};
