#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellAux/PlaneTools.h"
#include <boost/container_hash/hash.hpp>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud;
// using WireCell::PointCloud::Dataset;
using namespace WireCell::PointCloud::Tree;  // for "Points" node value type
// using WireCell::PointCloud::Tree::named_pointclouds_t;

#include "WireCellUtil/Logging.h"
using spdlog::debug;

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

std::ostream& Facade::operator<<(std::ostream& os, const Facade::Grouping& grouping)
{
    os << "<Grouping [" << (void*) &grouping << "]:" << " nclusters=" << grouping.nchildren() << ">";
    return os;
}

std::string Facade::dump(const Facade::Grouping& grouping, int level)
{
    std::stringstream ss;

    ss << grouping;
    if (level == 0) {
        return ss.str();
    }
    ss << "\n";
    size_t nc = 0;
    for (const auto* cluster : grouping.children()) {
        ss << nc++ << "\t" << *cluster << "\n";
        if (level == 1) {
            continue;
        }
        size_t nb = 0;
        for (const auto* blob : cluster->children()) {
            ss << nb++ << "\t\t" << *blob << "\n";
        }
    }
    return ss.str();
}


static std::tuple<int, int, int> parse_dead_winds(const std::string& ds_name) {
    int apa, face;
    char plane;
    // Use sscanf to extract the numbers and the plane letter.
    // The format string must match the structure of ds_name.
    if (std::sscanf(ds_name.c_str(), "dead_winds_a%df%dp%c", &apa, &face, &plane) != 3) {
        throw std::runtime_error("Failed to parse string: " + ds_name);
    }
    // Convert the plane letter to an index.
    int plane_index = -1;
    switch (plane) {
        case 'U': plane_index = 0; break;
        case 'V': plane_index = 1; break;
        case 'W': plane_index = 2; break;
        default: 
            throw std::runtime_error("Unexpected plane letter in: " + ds_name);
    }
    return std::make_tuple(apa, face, plane_index);
}

void Grouping::on_construct(node_type* node)
{
    this->NaryTree::Facade<points_t>::on_construct(node);
    const auto& lpcs = m_node->value.local_pcs();
    for (const auto& [name, pc_dead_winds] : lpcs) {
        if (name.find("dead_winds") != std::string::npos) {
            const auto& xbeg = pc_dead_winds.get("xbeg")->elements<float_t>();
            const auto& xend = pc_dead_winds.get("xend")->elements<float_t>();
            const auto& wind = pc_dead_winds.get("wind")->elements<int_t>();
            auto [apa, face, plane] = parse_dead_winds(name);
            for (size_t i = 0; i < xbeg.size(); ++i) {
                m_dead_winds[apa][face][plane][wind[i]] = {xbeg[i], xend[i]};
            }
        }
    }

    // for (const int face : faces) {
    //     for (const int plane : planes) {
    //         const std::string ds_name = String::format("dead_winds_f%dp%d", face, plane);
    //         if (lpcs.find(ds_name) == lpcs.end()) continue;
    //         const auto& pc_dead_winds = lpcs.at(ds_name);
    //         const auto& xbeg = pc_dead_winds.get("xbeg")->elements<float_t>();
    //         const auto& xend = pc_dead_winds.get("xend")->elements<float_t>();
    //         const auto& wind = pc_dead_winds.get("wind")->elements<int_t>();
    //         for (size_t i = 0; i < xbeg.size(); ++i) {
    //             m_dead_winds[face][plane][wind[i]] = {xbeg[i], xend[i]};
    //         }
    //     }
    // }
}


void Facade::Grouping::from(const Grouping& other)
{
    m_anodes = other.m_anodes;
    m_dv = other.m_dv;
}

void Facade::Grouping::enumerate_idents(int id)
{
    for (auto* cluster : children()) {
        cluster->set_cluster_id(id++);
    }
}


void Grouping::fill_cache(GroupingCache& gc) const
{
    {
        // In pre-cached code this was Grouping::fill_proj_centers_pitch_mags() const

        const int ndummy_layers = 2;
        if (m_anodes.size()==0) {
            raise<ValueError>("anode is null");
        }
        for (const auto& [ident, anode] : m_anodes) {
            for (const auto& face : anode->faces()) {
                const auto& coords = face->raygrid();
                // skip dummy layers so the vector matches 0, 1, 2 plane order
                for (int layer=ndummy_layers; layer<coords.nlayers(); ++layer) {
                    const auto& pitch_dir = coords.pitch_dirs()[layer];
                    const auto& center = coords.centers()[layer];
                    double proj_center = center.dot(pitch_dir);
                    gc.proj_centers[anode->ident()][face->ident()][layer - ndummy_layers] = proj_center;
                    gc.pitch_mags[anode->ident()][face->ident()][layer - ndummy_layers] = coords.pitch_mags()[layer];
                }
            }
        }
    }

    {
        for (size_t iclus = 0; iclus != children().size(); iclus++) {
            const Cluster* cluster = children().at(iclus);
            const auto& wpids = cluster->wpids_blob();
            gc.cluster_wpids.insert(wpids.begin(), wpids.end());
        }
        // for (const auto wpid : gc.cluster_wpids) {
        //     std::cout << "Grouping::fill_cache wpid: " << wpid.name() << std::endl;
        // }
    }

    // fill cache related to the detector volume
    // for (const auto wpid : gc.dv_wpids) {
    //     std::cout << "DEBUG Grouping::fill_cache wpid: " << wpid.name() << std::endl;
    // }
    fill_dv_cache(gc);
    // for (const auto wpid : gc.dv_wpids) {
    //     std::cout << "DEBUG Grouping::fill_cache wpid: " << wpid.name() << std::endl;
    // }

    // fill plane_channels ...
    fill_plane_channels_cache(gc);
}

void Grouping::fill_plane_channels_cache(GroupingCache& gc) const
{
    for (auto wpid : gc.dv_wpids) {
        int face = wpid.face();
        int apa = wpid.apa();

        // Create wpids for all three planes with the same APA and face
        WirePlaneId wpid_u(kUlayer, face, apa);
        WirePlaneId wpid_v(kVlayer, face, apa);
        WirePlaneId wpid_w(kWlayer, face, apa);
        auto& anode = get_anode(apa);

        gc.map_plane_channels[apa][face][kUlayer] = Aux::plane_channels(anode, wpid_u.index());
        gc.map_plane_channels[apa][face][kVlayer] = Aux::plane_channels(anode, wpid_v.index());
        gc.map_plane_channels[apa][face][kWlayer] = Aux::plane_channels(anode, wpid_w.index());

        // std::cout << "Test: " << apa << " " << face << " " << wpid_u.index() << " " << gc.map_plane_channels[apa][face][kUlayer].size() << " " << wpid_v.index() << " " << gc.map_plane_channels[apa][face][kVlayer].size() <<  " " << wpid_w.index() << " " << gc.map_plane_channels[apa][face][kWlayer].size() << " " << wpid_u.layer() << " " << wpid_v.layer() << " " << wpid_w.layer() << std::endl;
    }
}

void Grouping::fill_dv_cache(GroupingCache& gc) const
{
    if (m_dv != nullptr) {
        for (auto& [wpid_ident, iface] : m_dv->wpident_faces()) {
            const WirePlaneId wpid(wpid_ident);
            // std::cout << "DEBUG Grouping::fill_dv_cache wpid: " << wpid.name() << std::endl;
            gc.dv_wpids.insert(wpid);
            // std::cout << "DEBUG Grouping::fill_dv_cache gc.dv_wpids.size() " << gc.dv_wpids.size() << std::endl;
            int face = wpid.face();
            int apa = wpid.apa();
            // int plane = wpid.index();
            // std::cout << "Test: " << apa << " " << face << " " << plane << " " << kAllLayers << " " << m_dv << std::endl;
            WirePlaneId wpid_all(kAllLayers, face, apa);
            gc.map_time_offset[apa][face] = m_dv->metadata(wpid_all)["time_offset"].asDouble();
            gc.map_drift_speed[apa][face] = m_dv->metadata(wpid_all)["drift_speed"].asDouble();
            gc.map_tick[apa][face] = m_dv->metadata(wpid_all)["tick"].asDouble();

            // Create wpids for all three planes with the same APA and face
            WirePlaneId wpid_u(kUlayer, face, apa);
            WirePlaneId wpid_v(kVlayer, face, apa);
            WirePlaneId wpid_w(kWlayer, face, apa);
            
            // Get wire directions for all planes
            Vector wire_dir_u = m_dv->wire_direction(wpid_u);
            Vector wire_dir_v = m_dv->wire_direction(wpid_v);
            Vector wire_dir_w = m_dv->wire_direction(wpid_w);
            
            // Calculate angles
            gc.map_wire_angles[apa][face][0] = std::atan2(wire_dir_u.z(), wire_dir_u.y());
            gc.map_wire_angles[apa][face][1] = std::atan2(wire_dir_v.z(), wire_dir_v.y());
            gc.map_wire_angles[apa][face][2] = std::atan2(wire_dir_w.z(), wire_dir_w.y());

            gc.map_drift_dir[apa][face] = m_dv->face_dirx(wpid);

            gc.map_nticks_per_slice[apa][face] = m_dv->metadata(wpid_all)["nticks_live_slice"].asInt();

            // std::cout << "Test: " << gc.map_time_offset[apa][face] << " " << gc.map_drift_speed[apa][face] << " " << gc.map_tick[apa][face] << " " << gc.map_drift_dir[apa][face]  << std::endl;
        }
        // for (auto wpid : gc.dv_wpids) {
        //     std::cout << "DEBUG Grouping::fill_dv_cache gc.dv_wpids wpid: " << wpid.name() << std::endl;
        // }
        // double time_offset = m_dv->metadata(wpid_all)["time_offset"].asDouble(); 
        // std::map<int, std::map<int, std::map<string, double> > > map_time_offset;
        // std::map<int, std::map<int, std::map<string, double> > > map_drift_speed;
        // std::map<int, std::map<int, std::map<string, double> > > map_tick;
    }
}




void Grouping::set_anodes(const std::vector<IAnodePlane::pointer>& anodes) {
    for (auto anode : anodes) {
        m_anodes[anode->ident()] = anode;
    }
}

const IAnodePlane::pointer Grouping::get_anode(const int ident) const {
    if (m_anodes.find(ident) == m_anodes.end()) {
        raise<ValueError>("anode %d not found", ident);
    }
    return m_anodes.at(ident);
}

size_t Grouping::hash() const
{
    std::size_t h = 0;
    for (auto wpid : cache().dv_wpids) {
        boost::hash_combine(h, wpid.ident());
    }
    auto clusters = children();  // copy vector
    // sort_clusters(clusters);
    for (const Cluster* cluster : clusters) {
        boost::hash_combine(h, cluster->hash());
    }
    return h;
}

const Grouping::kd2d_t& Grouping::kd2d(const int apa, const int face, const int pind) const
{
    std::vector<std::string> plane_names = {"U", "V", "W"};
    const auto sname = String::format("ctpc_a%df%dp%d",apa, face, plane_names[pind]);
    // const auto sname = String::format("ctpc_f%dp%d", face, pind);
    Tree::Scope scope = {sname, {"x", "y"}, 1};
    const auto& sv = m_node->value.scoped_view(scope);
    // std::cout << "sname: " << sname << " npoints: " << sv.kd().npoints() << std::endl;
    return sv.kd();
}


bool Grouping::is_good_point(const geo_point_t& point, const int apa, const int face, double radius, int ch_range, int allowed_bad) const {
    const int nplanes = 3;
    int matched_planes = 0;
    for (int pind = 0; pind < nplanes; ++pind) {
        if (get_closest_points(point, radius, apa, face, pind).size() > 0) {
            matched_planes++;
        } else if (get_closest_dead_chs(point, ch_range, apa, face, pind)) {
            matched_planes++;
        }
    }
    // std::cout << "matched_planes: " << matched_planes << std::endl;
    if (matched_planes >= nplanes - allowed_bad) {
        return true;
    }
    return false;
}

bool Grouping::is_good_point_wc(const geo_point_t& point, const int apa, const int face, double radius, int ch_range, int allowed_bad) const 
{
    const int nplanes = 3;
    int matched_planes = 0;
    
    // Loop through U,V,W planes
    for (int pind = 0; pind < nplanes; pind++) {
        int weight = (pind == 2) ? 2 : 1; // W plane counts double
        if (get_closest_points(point, radius, apa, face, pind).size() > 0) {
            matched_planes += weight;
        }
        else if (get_closest_dead_chs(point, ch_range, apa, face, pind)) {
            matched_planes += weight;
        }
    }

    return matched_planes >= 4 - allowed_bad;
}

std::vector<int> Grouping::test_good_point(const geo_point_t& point, const int apa, const int face, 
    double radius, int ch_range) const 
{
    std::vector<int> num_planes(6, 0);  // Initialize with 6 zeros
    // std::cout << "abc: " << point << " " << radius << " " << ch_range << std::endl;
    // Check each plane (0,1,2)
    for (int pind = 0; pind < 3; ++pind) {
        // Get closest points for this plane
        const auto closest_pts = get_closest_points(point, radius, apa, face, pind);
        
        if (closest_pts.size() > 0) {
            // Has hits in this plane
            num_planes[pind]++;
        }
        else {
            // No hits, check if it's in dead region
            if (get_closest_dead_chs(point, ch_range, apa, face, pind)) {
                num_planes[pind + 3]++;
            }
        }
        // std::cout << closest_pts.size() << " " << get_closest_dead_chs(point, ch_range, face, pind) << " " << num_planes[pind] << " " << num_planes[pind+3] << std::endl;
    }
    
    return num_planes;
}

double Facade::Grouping::get_ave_3d_charge(const geo_point_t& point, const int apa, const int face,  const double radius) const {
    double charge = 0;
    int ncount = 0;
    const int nplanes = 3;
    // Check all three planes
    for (int pind = 0; pind < nplanes; ++pind) {
        if (!get_closest_dead_chs(point, 1, apa, face, pind)) {
            charge += get_ave_charge(point, apa, face, pind, radius);
            ncount++;
        }
    }

    if (ncount != 0) {
        charge /= ncount;
    }
    return charge;
}

double Facade::Grouping::get_ave_charge(const geo_point_t& point, const int apa, const int face, const int pind, const double radius) const {
    double sum_charge = 0;
    double ncount = 0;

    // Get closest points within radius
    auto nearby_points = get_closest_points(point, radius, apa, face, pind);

    // Access the charge information from ctpc dataset
    std::vector<std::string> plane_names = {"U", "V", "W"};
    const std::string ds_name = String::format("ctpc_a%df%dp%d",apa, face, plane_names[pind]);
    // const std::string ds_name = String::format("ctpc_f%dp%d", face, pind);
    const auto& local_pcs = m_node->value.local_pcs();
    
    if (local_pcs.find(ds_name) == local_pcs.end()) {
        return 0.0;
    }

    const auto& ds = local_pcs.at(ds_name);
    const auto& charges = ds.get("charge")->elements<float_t>();
    
    // Sum charges for nearby points
    for (const auto& [ind, dist] : nearby_points) {
        sum_charge += charges[ind];
        ncount++;
    }

    if (ncount != 0) {
        sum_charge /= ncount;
    }
    return sum_charge;
}



Grouping::kd_results_t Grouping::get_closest_points(const geo_point_t& point, const double radius, const int apa, const int face,
                                                    int pind) const
{
    double x = point[0];
    const auto [angle_u,angle_v,angle_w] = wire_angles(apa, face);
    std::vector<double> angles = {angle_u, angle_v, angle_w};
    double y = cos(angles[pind]) * point[2] - sin(angles[pind]) * point[1];
    const auto& skd = kd2d(apa, face, pind);
    return skd.radius<std::vector<double>>(radius * radius, {x, y});
}

bool Grouping::get_closest_dead_chs(const geo_point_t& point, const int ch_range, const int apa, const int face, int pind) const {
    const auto [tind, wind] = convert_3Dpoint_time_ch(point, apa, face, pind);
    const auto& ch2xrange = get_dead_winds(apa, face, pind);
    for (int ch = wind - ch_range; ch <= wind + ch_range; ++ch) {
        if (ch2xrange.find(ch) ==  ch2xrange.end()) continue;
        const auto [xmin, xmax] = ch2xrange.at(ch);
        if (point[0] >= xmin && point[0] <= xmax) {
            // std::cout << "ch " << ch << " x " << point[0] << " xmin " << xmin << " xmax " << xmax << std::endl;
            return true;
        }
    }
    return false;
}

std::tuple<int, int> Grouping::convert_3Dpoint_time_ch(const geo_point_t& point, const int apa, const int face, const int pind) const {
    if (m_anodes.size()==0) {
        raise<ValueError>("Anode is null");
    }
    const auto iface = m_anodes.at(apa)->faces()[face];
    if (iface == nullptr) {
        raise<ValueError>("anode %d has no face %d", m_anodes.at(apa)->ident(), face);
    }

    const auto [angle_u,angle_v,angle_w] = wire_angles(apa, face);
    std::vector<double> angles = {angle_u, angle_v, angle_w};
    const double angle = angles[pind];
    const double pitch = pitch_mags().at(apa).at(face).at(pind);
    const double center = proj_centers().at(apa).at(face).at(pind);

    // std::cout << "Test: " << pitch/units::cm << " " << center/units::cm << std::endl;

    const int wind = point2wind(point, angle, pitch, center);

    // const auto params = get_params();
    double time_offset = cache().map_time_offset.at(apa).at(face);
    double drift_speed = cache().map_drift_speed.at(apa).at(face);
    double tick = cache().map_tick.at(apa).at(face);

    //std::cout << "Test: " << params.time_offset/units::us << " " << params.drift_speed/(units::mm/units::us) << " " << point[0] << std::endl;

    const double time = drift2time(iface, time_offset, drift_speed, point[0]);
    const int tind = std::round(time / tick);

    return {tind, wind};
}

std::pair<double,double> Grouping::convert_time_ch_2Dpoint(const int timeslice, const int channel, const int apa, const int face, const int plane) const 
{
    if (m_anodes.size() == 0) {
        raise<ValueError>("Anode is null");
    }
    const auto iface = m_anodes.at(apa)->faces()[face];
    if (iface == nullptr) {
        raise<ValueError>("anode %d has no face %d", m_anodes.at(apa)->ident(), face);
    }
    const int nplanes = 3;
    // const auto params = get_params();
    const auto& pitch_mags = this->pitch_mags();
    const auto& proj_centers = this->proj_centers();
    
    double time_offset = cache().map_time_offset.at(apa).at(face);
    double drift_speed = cache().map_drift_speed.at(apa).at(face);
    double tick = cache().map_tick.at(apa).at(face);

    // Convert time to x position
    const double x = time2drift(iface, time_offset, drift_speed, timeslice * tick);
    
    // Get y position based on channel and plane
    double y;
    if (plane >= 0 && plane < nplanes) {
        const double pitch = pitch_mags.at(apa).at(face).at(plane);
        const double center = proj_centers.at(apa).at(face).at(plane);
        y = pitch * (channel+0.5) + center;
    }
    else {
        raise<ValueError>("invalid plane index %d", plane);
    }

    return std::make_pair(x, y);
}


size_t Grouping::get_num_points(const int apa, const int face, const int pind) const {
    std::vector<std::string> plane_names = {"U", "V", "W"};
    const auto sname = String::format("ctpc_a%df%dp%d",apa, face, plane_names[pind]);
    // const auto sname = String::format("ctpc_f%dp%d", face, pind);
    Tree::Scope scope = {sname, {"x", "y"}, 1};
    const auto& sv = m_node->value.scoped_view(scope);
    return sv.npoints();
}

std::vector<std::pair<int, int>> Facade::Grouping::get_overlap_dead_chs(const int min_time, const int max_time,
    const int min_ch, const int max_ch, const int apa, const int face, const int pind, const bool flag_ignore_time) const 
{
    if (m_anodes.size() == 0) {
        raise<ValueError>("anode is null");
    }
    // const auto& params = get_params();
    double time_offset = cache().map_time_offset.at(apa).at(face);
    double drift_speed = cache().map_drift_speed.at(apa).at(face);
    // double tick = cache().map_tick.at(apa).at(face);
    
    // Convert time to position
    double min_xpos = time2drift(m_anodes.at(apa)->faces()[face], time_offset, drift_speed, min_time);
    double max_xpos = time2drift(m_anodes.at(apa)->faces()[face], time_offset, drift_speed, max_time);

    // Ensure min_xpos is actually less than max_xpos
    if (min_xpos > max_xpos) {
        std::swap(min_xpos, max_xpos);
    }

    std::set<int> dead_chs;
    const auto& dead_winds = get_dead_winds(apa, face, pind);

    // Find overlapping dead channels
    for (const auto& [wind, xrange] : dead_winds) {
        const int temp_ch = wind;
        const double temp_min_xpos = xrange.first;
        const double temp_max_xpos = xrange.second;

        if (flag_ignore_time) {
            if (temp_ch >= min_ch && temp_ch <= max_ch) {
                dead_chs.insert(temp_ch);
            }
        }
        else {
            if (temp_ch >= min_ch && temp_ch <= max_ch &&
                max_xpos >= temp_min_xpos && min_xpos <= temp_max_xpos) {
                dead_chs.insert(temp_ch);
            }
        }
    }

    // Convert set of channels to ranges
    std::vector<std::pair<int, int>> dead_ch_range;
    for (const auto ch : dead_chs) {
        if (dead_ch_range.empty()) {
            dead_ch_range.push_back(std::make_pair(ch, ch));
        }
        else {
            if (ch - dead_ch_range.back().second == 1) {
                dead_ch_range.back().second = ch;
            }
            else {
                dead_ch_range.push_back(std::make_pair(ch, ch));
            }
        }
    }

    return dead_ch_range;
}

std::map<int, std::pair<int, int>> Facade::Grouping::get_all_dead_chs(const int apa, const int face, const int pind, int expand) const
{
    std::map<int, std::pair<int, int>> results;
    
    const auto& dead_winds = get_dead_winds(apa, face, pind);
    
    double time_offset = cache().map_time_offset.at(apa).at(face);
    double drift_speed = cache().map_drift_speed.at(apa).at(face);

    // Add entries for this face/plane's dead channels
    for (const auto& [wind, xrange] : dead_winds) {
        int temp_ch = wind;
        
        // Convert position range to time ticks using drift parameters
        int min_time = std::round(drift2time(m_anodes.at(apa)->faces()[face], 
                                           time_offset,
                                           drift_speed, 
                                           xrange.first)) - expand;
        int max_time = std::round(drift2time(m_anodes.at(apa)->faces()[face],
                                           time_offset,
                                           drift_speed,
                                           xrange.second)) + expand;
        
        results[temp_ch] = std::make_pair(std::min(min_time, max_time), std::max(min_time, max_time));
    }
    
    return results;
}

std::map<std::pair<int,int>, std::pair<double,double>> Facade::Grouping::get_overlap_good_ch_charge(
    int min_time, int max_time, int min_ch, int max_ch, const int apa, 
    const int face, const int pind) const 
{
    std::map<std::pair<int,int>, std::pair<double,double>> map_time_ch_charge;
    
    // Get the point cloud for this face/plane
    std::vector<std::string> plane_names = {"U", "V", "W"};
    const std::string ds_name = String::format("ctpc_a%df%dp%d",apa, face, plane_names[pind]);
    // const std::string ds_name = String::format("ctpc_f%dp%d", face, pind);
    if (m_node->value.local_pcs().find(ds_name) == m_node->value.local_pcs().end()) {
        return map_time_ch_charge; // Return empty if dataset not found
    }
    
    const auto& ctpc = m_node->value.local_pcs().at(ds_name);
    const auto& slice_index = ctpc.get("slice_index")->elements<int_t>();
    const auto& wind = ctpc.get("wind")->elements<int_t>();
    const auto& charge = ctpc.get("charge")->elements<float_t>();
    const auto& charge_err = ctpc.get("charge_err")->elements<float_t>();

    // Fill the map for points within the specified window
    for (size_t i = 0; i < slice_index.size(); ++i) {
        if (slice_index[i] >= min_time && slice_index[i] <= max_time &&
            wind[i] >= min_ch && wind[i] <= max_ch) {
            map_time_ch_charge[std::make_pair(slice_index[i], wind[i])] = 
                std::make_pair(charge[i], charge_err[i]);
        }
    }

    return map_time_ch_charge;
}



void Grouping::clear_cache() const
{
    this->Mixin<Grouping, GroupingCache>::clear_cache();


    // This is utterly broken.  #381.
    m_dead_winds.clear(); 

}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
