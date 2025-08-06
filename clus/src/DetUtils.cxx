#include "WireCellClus/DetUtils.h"

namespace WireCell::Clus {

    std::set<int> apa_idents(IDetectorVolumes::pointer dv)
    {
        std::set<int> apas;
        for (auto& [wpid_ident, iface] : dv->wpident_faces()) {
            const WirePlaneId wpid(wpid_ident);
            apas.insert(wpid.apa());
        }
        return apas;
    }

    wpid_faceparams_map face_parameters(IDetectorVolumes::pointer dv)
    {
        wpid_faceparams_map ret;
        for (auto& [wpid_ident, iface] : dv->wpident_faces()) {
            const WirePlaneId wpid(wpid_ident);
            int apa = wpid.apa();
            int face = wpid.face();

            // Create wpids for all three planes with this APA and face
            WirePlaneId wpid_u(kUlayer, face, apa);
            WirePlaneId wpid_v(kVlayer, face, apa);
            WirePlaneId wpid_w(kWlayer, face, apa);
     
            // Get drift direction based on face orientation
            int face_dirx = dv->face_dirx(wpid_u);
        
            // Get wire directions for all planes
            Vector wire_dir_u = dv->wire_direction(wpid_u);
            Vector wire_dir_v = dv->wire_direction(wpid_v);
            Vector wire_dir_w = dv->wire_direction(wpid_w);

            // Calculate angles
            double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
            double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
            double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());

            ret[wpid] = {Vector(face_dirx,0,0), angle_u, angle_v, angle_w};
        }
        return ret;
    }
    
    std::shared_ptr<Facade::DynamicPointCloud> make_dynamicpointcloud(IDetectorVolumes::pointer dv)
    {
        return std::make_shared<Facade::DynamicPointCloud>(face_parameters(dv));
    }

}
