#include "WireCellAux/PlaneTools.h"

using namespace WireCell;

IChannel::vector Aux::plane_channels(IAnodePlane::pointer anode,
                                     int wire_plane_index)
{
    IChannel::vector ret;    
    for (auto face : anode->faces()) {
        if (!face) {   // A null face means one sided AnodePlane.
            continue;  // Can be "back" or "front" face.
        }
        for (auto plane : face->planes()) {
            if (wire_plane_index != plane->planeid().index()) {
                continue;
            }
            // These IChannel vectors are ordered in same order as wire-in-plane.
            const auto& ichans = plane->channels();
            // Append
            ret.reserve(ret.size() + ichans.size());
            ret.insert(ret.end(), ichans.begin(), ichans.end());
        }
    }
    return ret;
}


// Get wire information for a specific plane from IAnodeFace
Aux::WirePlaneInfo Aux::get_wire_plane_info(IAnodeFace::pointer face, WirePlaneLayer_t layer) {
    WirePlaneInfo info = {0, 0, 0};
    
    // Get the wire plane for this layer
    auto planes = face->planes();
    IWirePlane::pointer target_plane = nullptr;
    for (auto plane : planes) {
        if (plane->planeid().layer() == layer) {
            target_plane = plane;
            break;
        }
    }
    
    if (!target_plane) {
        return info;
    }

    // Get all wires in this plane
    const auto& wires = target_plane->wires();
    if (wires.empty()) {
        return info;
    }

    // Find min and max wire indices
    info.start_index = wires.front()->index();
    info.end_index = wires.back()->index();
    info.total_wires = wires.size();

    return info;
}

