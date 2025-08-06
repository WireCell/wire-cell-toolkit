#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/Facade.h"
#include "WireCellClus/Facade_Cluster.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

namespace WireCell::Clus::Facade {

// Function to compute wire plane parameters
void compute_wireplane_params(
    const std::set<WirePlaneId>& wpids,
    const IDetectorVolumes::pointer dv,
    std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>>& wpid_params,
    std::map<WirePlaneId, std::pair<geo_point_t, double>>& wpid_U_dir,
    std::map<WirePlaneId, std::pair<geo_point_t, double>>& wpid_V_dir,
    std::map<WirePlaneId, std::pair<geo_point_t, double>>& wpid_W_dir,
    std::set<int>& apas)
{
    for (const auto& wpid : wpids) {
        int apa = wpid.apa();
        int face = wpid.face();
        apas.insert(apa);

        // Create wpids for all three planes with this APA and face
        WirePlaneId wpid_u(kUlayer, face, apa);
        WirePlaneId wpid_v(kVlayer, face, apa);
        WirePlaneId wpid_w(kWlayer, face, apa);

        // Get drift direction based on face orientation
        int face_dirx = dv->face_dirx(wpid_u);
        geo_point_t drift_dir(face_dirx, 0, 0);

        // Get wire directions for all planes
        Vector wire_dir_u = dv->wire_direction(wpid_u);
        Vector wire_dir_v = dv->wire_direction(wpid_v);
        Vector wire_dir_w = dv->wire_direction(wpid_w);

        // Calculate angles
        double angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
        double angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
        double angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());

        wpid_params[wpid] = std::make_tuple(drift_dir, angle_u, angle_v, angle_w);
        wpid_U_dir[wpid] = std::make_pair(geo_point_t(0, std::cos(angle_u), std::sin(angle_u)), angle_u);
        wpid_V_dir[wpid] = std::make_pair(geo_point_t(0, std::cos(angle_v), std::sin(angle_v)), angle_v);
        wpid_W_dir[wpid] = std::make_pair(geo_point_t(0, std::cos(angle_w), std::sin(angle_w)), angle_w);
    }
}

}  // namespace WireCell::Clus::Facade