#include "WireCellAux/PlaneTools.h"

using namespace WireCell;

IChannel::vector Aux::plane_channels(IAnodePlane::pointer anode,
                                     int wire_plane_index)
{
    IChannel::vector ret;    
    // std::cout << "plane_channels called for plane index " << wire_plane_index << std::endl;
    for (auto face : anode->faces()) {
        if (!face) {   // A null face means one sided AnodePlane.
            continue;  // Can be "back" or "front" face.
        }
        // int face_num = face->which();
        // std::cout << "  Checking face " << face_num << std::endl;

        for (auto plane : face->planes()) {
            if (wire_plane_index != plane->planeid().index()) {
                continue;
            }
            // auto plane_wpid = plane->planeid();
            // std::cout << "    Found matching plane: " << plane_wpid.name() 
            //          << " (face=" << plane_wpid.face() << ", layer=" << plane_wpid.layer() << ")" << std::endl;
            
            const auto& ichans = plane->channels();
            // std::cout << "    Adding " << ichans.size() << " channels from this plane" << std::endl;
            
            // Debug: Check what face the channels actually belong to
            // if (wire_plane_index == 1 && ichans.size() > 0) { // W-layer
            //     std::cout << "    Sample channel check:" << std::endl;
            //     auto sample_channel = ichans[0];
            //     auto sample_wires = anode->wires(sample_channel->ident());
            //     for (const auto& wire : sample_wires) {
            //         auto wire_wpid = wire->planeid();
            //         std::cout << "      Channel " << sample_channel->ident() 
            //                  << " -> Wire face=" << wire_wpid.face() << std::endl;
            //     }
            // }            
            // These IChannel vectors are ordered in same order as wire-in-plane.
            // const auto& ichans = plane->channels();
            // Append
            ret.reserve(ret.size() + ichans.size());
            ret.insert(ret.end(), ichans.begin(), ichans.end());
        }
    }
    //     std::cerr << "PlaneTools: ret contains " << ret.size() << " elements: ";
    // for (const auto& val : ret) {
    //     std::cerr << val->ident() << " ";
    // }
    // std::cerr << std::endl;
    return ret;
}

IChannel::vector Aux::plane_channels_for_face(IAnodePlane::pointer anode,
                                               int wire_plane_index,
                                               int target_face)
{
    IChannel::vector ret;    

    for (auto face : anode->faces()) {
        if (!face || face->which() != target_face) continue;
              int face_num = face->which(); 
            std::cout << "  Checking face " << face_num << std::endl;

        for (auto plane : face->planes()) {
            if (wire_plane_index != plane->planeid().index()) {
                continue;
            }
            auto plane_wpid = plane->planeid();
            std::cout << "    Found matching plane: " << plane_wpid.name() 
                     << " (face=" << plane_wpid.face() << ", layer=" << plane_wpid.layer() << ")" << std::endl;
            
            const auto& ichans = plane->channels();
            std::cout << "    Adding " << ichans.size() << " channels from this plane" << std::endl;
            
            // Debug: Check what face the channels actually belong to
            if (wire_plane_index == 0 && ichans.size() > 0) { // W-layer

                std::cout << "    Sample channel check:" << std::endl;
                auto sample_channel = ichans[0];
                auto sample_wires = anode->wires(sample_channel->ident());
                for (const auto& wire : sample_wires) {
                    auto wire_wpid = wire->planeid();
                    std::cout << "      Channel " << sample_channel->ident() 
                             << " -> Wire face=" << wire_wpid.face() << std::endl;
                }
            }             
            // const auto& ichans = plane->channels();
            ret.reserve(ret.size() + ichans.size());
            ret.insert(ret.end(), ichans.begin(), ichans.end());
        }
    }
    // std::sort(ret.begin(), ret.end(), [](const auto& a, const auto& b) {
    //     return a->ident() < b->ident();
    // });
    // std::cerr << "PlaneTools: ret contains " << ret.size() << " elements: ";
    // for (const auto& val : ret) {
    //     std::cerr << val->ident() << " ";
    // }
    // std::cerr << std::endl;
    // exit(0); // Debugging exit point
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

