#include "WireCellUtil/DetectorWires.h"
using namespace WireCell;
using namespace WireCell::WireSchema;

size_t DetectorWires::plane(WireSchema::StoreDB& storedb, WireSchema::Plane& plane,
                            const Plane& param)
{
    const Vector pdir(0, sin(param.pangle), cos(param.pangle));
    const Vector pvec = param.pitch * pdir;
    const Ray pray(param.center0, param.center0 + pvec);
    return generate(storedb, plane, pray, param.bounds);
}

template<typename T>
T& element(std::vector<T>& arr, int ind)
{
    return arr[ind];
}

void DetectorWires::flat_channels(WireSchema::StoreDB& storedb, WireSchema::Plane& plane, int ch0)
{
    for (int iwire : plane.wires) {
        auto& wire = element(storedb.wires, iwire);
        wire.channel = ch0;
        ++ch0;
    }
}


void DetectorWires::flat_channels(WireSchema::StoreDB& storedb, int ch0)
{
    for (auto& detector : storedb.detectors) {
        for (int ianode : detector.anodes) {
            auto& anode = element(storedb.anodes, ianode);
            for (int iface : anode.faces) {
                auto& face = element(storedb.faces, iface);
                for (int iplane : face.planes) {
                    auto& plane = element(storedb.planes, iplane);
                    flat_channels(storedb, plane, ch0);
                    ch0 += plane.wires.size();
                }
            }
        }
    }
}

