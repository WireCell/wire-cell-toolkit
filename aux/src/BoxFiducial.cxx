#include "WireCellIface/IFiducial.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/BoundingBox.h"

#include <vector>
// #include <iostream>             // only debug

// Implementation is totally local to this comp. unit so no need for namespacing.
class BoxFiducial;

WIRECELL_FACTORY(BoxFiducial, BoxFiducial,
                 WireCell::IFiducial, WireCell::IConfigurable)


using namespace WireCell;

class BoxFiducial : public IFiducial, public IConfigurable {

    BoundingBox m_bb;

public:

    BoxFiducial() {}
    virtual ~BoxFiducial() {}

    virtual Configuration default_configuration() const {
        Configuration cfg;
        // This is a standard "Ray" form:
        //   {tail:{x:1,y:2,z:3}, head:{x:10,y:20:z:30}}
        cfg["bounds"] = Json::objectValue;
        return cfg;
    }

    virtual void configure(const Configuration& cfg) {
        Ray ray = get<Ray>(cfg, "bounds");
        m_bb = BoundingBox(ray);
    }


    // IFiducial
    virtual bool contained(const Point& point) const {
        return m_bb.inside(point);
    }
};
