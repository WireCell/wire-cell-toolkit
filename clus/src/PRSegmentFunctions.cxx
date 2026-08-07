#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/PRTrajectoryView.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/KSTest.h"
#include "WireCellUtil/Logging.h"
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <numeric>
#include <algorithm>

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");

namespace WireCell::Clus::PR {
    void create_segment_point_cloud(SegmentPtr segment,
                                const std::vector<geo_point_t>& path_points,
                                const IDetectorVolumes::pointer& dv,
                                const std::string& cloud_name,
                                const std::vector<size_t>& global_indices)
    {
        if (!segment || !segment->cluster()) {
            raise<RuntimeError>("create_segment_point_cloud: invalid segment or missing cluster");
        }
        
        auto& cluster = *segment->cluster();
        
        // Create point-plane pairs
        std::vector<std::pair<geo_point_t, WirePlaneId>> point_plane_pairs;
        for (const auto& point : path_points) {
            WirePlaneId wpid = dv->contained_by(point);
            point_plane_pairs.emplace_back(point, wpid);
        }
        
        // Get wpid_params (from detector configuration)
        const auto& wpids = cluster.grouping()->wpids();
        std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> wpid_params;
        std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_U_dir;
        std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_V_dir;
        std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_W_dir;
        std::set<int> apas;
        Facade::compute_wireplane_params(wpids, dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);
        
        // Create DynamicPointCloud
        auto dpc = std::make_shared<Facade::DynamicPointCloud>(wpid_params);
        
        // Create DPCPoints using factory function
        auto dpc_points = Facade::make_points_direct(&cluster, dv, wpid_params, point_plane_pairs, true);
        
        // Add points to DynamicPointCloud
        dpc->add_points(dpc_points);
        
        // Remove existing point cloud if it exists
        if (segment->dpcloud(cloud_name)) {
            segment->dpcloud(cloud_name, nullptr);
        }
        
        // Associate with segment
        segment->dpcloud(cloud_name, dpc);
        
        // Store global indices if provided
        if (!global_indices.empty()) {
            segment->set_global_indices(cloud_name, global_indices);
        }
    }

    void create_segment_fit_point_cloud(SegmentPtr segment,
                                    const IDetectorVolumes::pointer& dv,
                                    const std::string& cloud_name){
        std::vector<geo_point_t> fit_points;
        
        if (!segment || !segment->cluster()) {
            raise<RuntimeError>("create_segment_fit_point_cloud: invalid segment or missing cluster");
        }
        
        // Include all fit points unconditionally.  The origin-rejection guard was overly
        // cautious: fits from clear_fit carry valid wcpt-derived positions, and any fit
        // with index>=0 (fit.valid()) should be included even if it happens to land at
        // (0,0,0).  Downstream create_segment_point_cloud handles empty inputs safely.
        const auto& fits = segment->fits();
        fit_points.reserve(fits.size());
        for (const auto& fit : fits) {
            fit_points.push_back(fit.point);
        }
        create_segment_point_cloud(segment, fit_points, dv, cloud_name);
  
    }


    std::pair<double, WireCell::Point> segment_get_closest_point(SegmentPtr seg, const WireCell::Point& point, const std::string& cloud_name, const std::string& base_cloud_name){
        double min_dist = 1e9;
        WireCell::Point closest_point(0,0,0);
        
        if (!seg) {
            raise<RuntimeError>("get_closest_point: invalid segment");
        }
        
        auto dpc = seg->dpcloud(cloud_name);
        if (!dpc || dpc->npoints() == 0) {
            // Fall back to base_cloud_name (typically "main") if the requested cloud is
            // absent or empty (e.g. fitting produced no valid-index fit points).
            dpc = seg->dpcloud(base_cloud_name);
        }

        // If both clouds are absent or empty (e.g. segments from lightly-processed
        // other_clusters that share the same graph), treat the segment as infinitely
        // far away rather than crashing.  The caller's minimum-distance search will
        // simply skip this segment.
        if (!dpc || dpc->npoints() == 0) {
            return {min_dist, closest_point};
        }

        if (dpc->npoints() == 0) {
            return {min_dist, closest_point};
        }
        
        // Use KD-tree to find the closest point
        auto& kd_tree = dpc->kd3d();
        auto knn_results = kd_tree.knn(1, point);
        
        if (!knn_results.empty()) {
            size_t closest_index = knn_results[0].first;
            min_dist = std::sqrt(knn_results[0].second); // knn returns squared distance
            
            // Get the actual point from the DynamicPointCloud
            closest_point = dpc->point3d(closest_index);
        }
        
        return {min_dist, closest_point};
    }

    // NOTE (multi-APA semantics): this function intentionally takes a single (apa,face)
    // and returns the 2D distances from the segment's points that share that face.
    // For cross-APA segments, points on the other face are excluded — this is correct
    // behaviour because the query point is a 2D wire measurement at a specific face and
    // the comparison must use the same wire-coordinate system.  Callers are responsible
    // for passing the face that matches the query point.
    std::tuple<double, double, double> segment_get_closest_2d_distances(SegmentPtr seg, const WireCell::Point& point, int apa, int face, const std::string& cloud_name) {
        if (!seg) {
            raise<RuntimeError>("segment_get_closest_2d_distances: invalid segment");
        }

        auto dpc = seg->dpcloud(cloud_name);
        if (!dpc || dpc->npoints() == 0) {
            // Fall back to "main" cloud if requested cloud is absent or empty
            dpc = seg->dpcloud("main");
        }

        // If both clouds are absent or empty, return infinite distances so the
        // caller's minimum-distance search simply skips this segment.
        if (!dpc || dpc->npoints() == 0) {
            return {1e9, 1e9, 1e9};
        }

        if (dpc->npoints() == 0) {
            return {1e9, 1e9, 1e9};
        }
        
        // Use DynamicPointCloud's optimized method to get 2D distances for each plane
        auto closest_2d_u = dpc->get_closest_2d_point_info(point, 0, face, apa);  // U plane
        auto closest_2d_v = dpc->get_closest_2d_point_info(point, 1, face, apa);  // V plane  
        auto closest_2d_w = dpc->get_closest_2d_point_info(point, 2, face, apa);  // W plane
        
        // Extract distances for each plane (U=0, V=1, W=2)
        double min_dist_u = std::get<0>(closest_2d_u);
        double min_dist_v = std::get<0>(closest_2d_v);
        double min_dist_w = std::get<0>(closest_2d_w);
        
        return std::make_tuple(min_dist_u, min_dist_v, min_dist_w);
    }

    double segment_get_closest_2d_distance(SegmentPtr seg, const WireCell::Point& point, int apa, int face, int plane, const std::string& cloud_name) {
        if (!seg) {
            raise<RuntimeError>("segment_get_closest_2d_distance: invalid segment");
        }
        auto dpc = seg->dpcloud(cloud_name);
        if (!dpc || dpc->npoints() == 0) {
            // Fall back to "main" cloud if requested cloud is absent or empty
            dpc = seg->dpcloud("main");
            if (!dpc) {
                raise<RuntimeError>("segment_get_closest_2d_distance: segment missing DynamicPointCloud with name " + cloud_name + " and fallback 'main'");
            }
        }
        if (dpc->npoints() == 0) {
            raise<RuntimeError>("segment_get_closest_2d_distance: DynamicPointCloud has no points");
        }
        return std::get<0>(dpc->get_closest_2d_point_info(point, plane, face, apa));
    }

    std::tuple<WireCell::Point, WireCell::Vector, WireCell::Vector, bool> segment_search_kink(SegmentPtr seg, WireCell::Point& start_p, const std::string& cloud_name, double dQ_dx_threshold, double cathode_x, double cathode_kink_xcut){
        auto tmp_results = segment_get_closest_point(seg, start_p, cloud_name);
        WireCell::Point test_p = tmp_results.second;

        // Drift direction used to compute para_angles[i] = |angle_to_drift - 90°|.
        // The formula |acos(v·d/|v|) - 90°| is symmetric under negation of d, so
        // (1,0,0) is correct for both +x and -x drift (opposing APA faces in SBND/MicroBooNE).
        // For detectors with a non-x drift axis this would need to be resolved from
        // IDetectorVolumes::contained_by(seg start/end) — a future extension.
        WireCell::Vector drift_dir_abs(1,0,0);
        
        const auto& fits = seg->fits();
        if (fits.empty()) {
            WireCell::Point p1 = WireCell::Point(0,0,0);
            WireCell::Vector dir(0,0,0);
            return std::make_tuple(p1, dir, dir, false);
        }

        std::vector<double> refl_angles(fits.size(), 0);
        std::vector<double> para_angles(fits.size(), 0);
        
        // Start the angle search
        for (size_t i = 0; i < fits.size(); i++) {
            double angle1 = 0;
            double angle2 = 0;
            
            for (int j = 0; j < 6; j++) {
                WireCell::Vector v10(0,0,0);
                WireCell::Vector v20(0,0,0);
                
                const size_t offset = static_cast<size_t>((j + 1) * 2);
                if (i >= offset) {
                    v10 = fits[i].point - fits[i-offset].point;
                } else {
                    v10 = fits[i].point - fits.front().point;
                }
                
                if (i + offset < fits.size()) {
                    v20 = fits[i+offset].point - fits[i].point;
                } else {
                    v20 = fits.back().point - fits[i].point;
                }
                
                if (j == 0) {
                    double dot_product = v10.dot(v20);
                    double mag_product = v10.magnitude() * v20.magnitude();
                    if (mag_product > 0) {
                        angle1 = std::acos(std::max(-1.0, std::min(1.0, dot_product / mag_product))) / M_PI * 180.0;
                    }
                    
                    double drift_dot1 = v10.dot(drift_dir_abs);
                    double drift_dot2 = v20.dot(drift_dir_abs);
                    double drift_mag1 = v10.magnitude();
                    double drift_mag2 = v20.magnitude();
                    
                    double drift_angle1 = 90.0, drift_angle2 = 90.0;
                    if (drift_mag1 > 0) drift_angle1 = std::acos(std::max(-1.0, std::min(1.0, drift_dot1 / drift_mag1))) / M_PI * 180.0;
                    if (drift_mag2 > 0) drift_angle2 = std::acos(std::max(-1.0, std::min(1.0, drift_dot2 / drift_mag2))) / M_PI * 180.0;
                    
                    angle2 = std::max(std::abs(drift_angle1 - 90.0), std::abs(drift_angle2 - 90.0));
                } else {
                    if (v10.magnitude() != 0 && v20.magnitude() != 0) {
                        double dot_product = v10.dot(v20);
                        double mag_product = v10.magnitude() * v20.magnitude();
                        double current_angle1 = std::acos(std::max(-1.0, std::min(1.0, dot_product / mag_product))) / M_PI * 180.0;
                        angle1 = std::min(current_angle1, angle1);
                        
                        double drift_dot1 = v10.dot(drift_dir_abs);
                        double drift_dot2 = v20.dot(drift_dir_abs);
                        double drift_mag1 = v10.magnitude();
                        double drift_mag2 = v20.magnitude();
                        
                        double drift_angle1 = 90.0, drift_angle2 = 90.0;
                        if (drift_mag1 > 0) drift_angle1 = std::acos(std::max(-1.0, std::min(1.0, drift_dot1 / drift_mag1))) / M_PI * 180.0;
                        if (drift_mag2 > 0) drift_angle2 = std::acos(std::max(-1.0, std::min(1.0, drift_dot2 / drift_mag2))) / M_PI * 180.0;
                        
                        double current_angle2 = std::max(std::abs(drift_angle1 - 90.0), std::abs(drift_angle2 - 90.0));
                        angle2 = std::min(current_angle2, angle2);
                    }
                }
            }
            
            refl_angles[i] = angle1;
            para_angles[i] = angle2;
        }

        bool flag_check = false;
        int save_i = -1;
        bool flag_switch = false;
        bool flag_search = false;

        for (size_t i = 0; i < fits.size(); i++) {
            // Check if close to test point
            double dist_to_test = (test_p - fits[i].point).magnitude();
            if (dist_to_test < 0.1 * units::cm) flag_check = true;

            // Check distance constraints
            double dist_to_front = (fits[i].point - fits.front().point).magnitude();
            double dist_to_back = (fits[i].point - fits.back().point).magnitude();
            double dist_to_start = (fits[i].point - start_p).magnitude();

            if (dist_to_front < 1*units::cm ||
                dist_to_back < 1*units::cm ||
                dist_to_start < 1*units::cm) continue;

            if (flag_check) {
                // Cathode veto (doc pr/20 Part II, B0).  Gate ONLY the four accept
                // tests below: refl_angles/para_angles and the windowed sum_angles
                // are untouched, and the loop breaks at the first qualifying index,
                // so skipping a cathode index simply lets the scan continue and a
                // genuine kink elsewhere on the segment sees identical arithmetic.
                // The `> 0` guard and the strict `<` are both load-bearing for
                // byte-identicality when the knob is off.
                if (cathode_kink_xcut > 0 &&
                    std::abs(fits[i].point.x() - cathode_x) < cathode_kink_xcut) continue;

                // Calculate average and max dQ/dx in local region
                double ave_dQ_dx = 0;
                int ave_count = 0;
                double max_dQ_dx = fits[i].dQ / (fits[i].dx + 1e-9);

                for (int j = -2; j <= 2; j++) {
                    int idx = i + j;
                    if (idx >= 0 && idx < static_cast<int>(fits.size())) {
                        double local_dQ_dx = fits[idx].dQ / (fits[idx].dx + 1e-9);
                        ave_dQ_dx += local_dQ_dx;
                        ave_count++;
                        if (local_dQ_dx > max_dQ_dx) max_dQ_dx = local_dQ_dx;
                    }
                }
                if (ave_count != 0) ave_dQ_dx /= ave_count;

                // Calculate angle sums
                double sum_angles = 0;
                double nsum = 0;
                double sum_angles1 = 0;
                double nsum1 = 0;

                for (int j = -2; j <= 2; j++) {
                    int idx = i + j;
                    if (idx >= 0 && idx < static_cast<int>(fits.size())) {
                        if (para_angles[idx] > 10) {
                            sum_angles += pow(refl_angles[idx], 2);
                            nsum++;
                        }
                        if (para_angles[idx] > 7.5) {
                            sum_angles1 += pow(refl_angles[idx], 2);
                            nsum1++;
                        }
                    }
                }
                if (nsum != 0) sum_angles = sqrt(sum_angles / nsum);
                if (nsum1 != 0) sum_angles1 = sqrt(sum_angles1 / nsum1);

                // Apply kink detection criteria
                if (para_angles[i] > 10 && refl_angles[i] > 30 && sum_angles > 15) {
                    save_i = i;
                    break;
                } else if (para_angles[i] > 7.5 && refl_angles[i] > 45 && sum_angles1 > 25) {
                    save_i = i;
                    break;
                } else if (para_angles[i] > 15 && refl_angles[i] > 27 && sum_angles > 12.5) {
                    save_i = i;
                    break;
                } else if (para_angles[i] > 15 && refl_angles[i] > 22 && sum_angles > 19 &&
                          max_dQ_dx > dQ_dx_threshold*1.5 && ave_dQ_dx > dQ_dx_threshold) {
                    save_i = i;
                    flag_search = true;
                    break;
                }
            }
        }
        
        // Return results
        if (save_i > 0 && save_i+1 < static_cast<int>(fits.size())) {
            WireCell::Point p = fits[save_i].point;
            
            WireCell::Point prev_p(0,0,0);
            WireCell::Point next_p(0,0,0);
            int num_p = 0;
            int num_p1 = 0;
            
            double length1 = 0;
            double length2 = 0;
            WireCell::Point last_p1, last_p2;
            
            // Calculate direction vectors by averaging nearby points
            for (int i = 1; i < 10; i++) {
                if (save_i >= i) {
                    length1 += (fits[save_i-i].point - fits[save_i-i+1].point).magnitude();
                    prev_p = prev_p + fits[save_i-i].point;
                    last_p1 = fits[save_i-i].point;
                    num_p++;
                }
                if (save_i+i < static_cast<int>(fits.size())) {
                    length2 += (fits[save_i+i].point - fits[save_i+i-1].point).magnitude();
                    next_p = next_p + fits[save_i+i].point;
                    last_p2 = fits[save_i+i].point;
                    num_p1++;
                }
            }
            
            double length1_1 = (last_p1 - fits[save_i].point).magnitude();
            double length2_1 = (last_p2 - fits[save_i].point).magnitude();
            
            // Check for direction switch.
            // Guard: require the full 9-point window AND an absolute chord > 3 cm on the
            // "straight" side before trusting the straightness ratio.  When the kink sits
            // near a segment endpoint only 3-4 post-kink (or pre-kink) points exist, making
            // the path-length ≈ chord trivially regardless of true geometry.  Allowing
            // flag_switch to fire on such a degenerate window dispatches proto_extend_point
            // in the wrong direction and produces spurious near-endpoint tail segments.
            const double min_straight_chord  = 3.0 * units::cm;
            const int    min_straight_points = 9;

            if (num_p1 >= min_straight_points &&
                length2_1 > min_straight_chord &&
                std::abs(length2 - length2_1) < 0.03 * length2_1 &&
                length1 * length2_1 > 1.06 * length2 * length1_1) {
                flag_switch = true;
                flag_search = true;
            } else if (num_p >= min_straight_points &&
                       length1_1 > min_straight_chord &&
                       std::abs(length1 - length1_1) < 0.03 * length1_1 &&
                       length2 * length1_1 > 1.06 * length1 * length2_1) {
                flag_search = true;
            }

            prev_p = prev_p * (1.0/num_p);
            next_p = next_p * (1.0/num_p1);
            
            WireCell::Vector dir = (p - prev_p).norm();
            WireCell::Vector dir1 = (p - next_p).norm();
            
            // Calculate local charge density
            double sum_dQ = 0, sum_dx = 0;
            for (int i = -2; i <= 2; i++) {
                int idx = save_i + i;
                if (idx >= 0 && idx < static_cast<int>(fits.size())) {
                    sum_dQ += fits[idx].dQ;
                    sum_dx += fits[idx].dx;
                }
            }
            
            double local_dQdx = sum_dQ / (sum_dx + 1e-9);

            if (flag_search) {
                if (flag_switch) {
                    return std::make_tuple(p, dir1, dir, true);
                } else {
                    return std::make_tuple(p, dir, dir1, true);
                }
            // prototype (ProtoSegment.cxx:937): sum_dQ/sum_dx > 2500 e/mm =
            // 25000 e/cm = 25/43 of the uBooNE MIP median.  Expressed on the
            // routed dQ_dx_threshold (= m_mip_dqdx_median) so it scales with
            // the detector; exactly 25000/units::cm at the uBooNE default
            // (4300*25 = 107500, /43 = 2500, both FP-exact).  docs/pr/2 sec 7.3.
            } else if (local_dQdx > dQ_dx_threshold * 25.0 / 43.0) { //not too low ...
                if (flag_switch) {
                    return std::make_tuple(p, dir1, dir, false);
                } else {
                    return std::make_tuple(p, dir, dir1, false);
                }
            } else {
                if (flag_switch) {
                    return std::make_tuple(p, dir1, dir, true);
                } else {
                    return std::make_tuple(p, dir, dir1, true);
                }
            }
        } else {
            WireCell::Point p1 = fits.back().point;
            WireCell::Vector dir(0,0,0);
            return std::make_tuple(p1, dir, dir, false);
        }
    }




    std::tuple<bool, std::pair<SegmentPtr, SegmentPtr>, VertexPtr> break_segment(Graph& graph, SegmentPtr seg, Point point, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const IDetectorVolumes::pointer& dv, double max_dist/*=1e9*/)
    {
        /// sanity checks
        if (! seg->descriptor_valid()) {
            raise<RuntimeError>("break_segment: segment has invalid descriptor\n");
        }
        auto ed = seg->get_descriptor();
        auto vd1 = boost::source(ed, graph);
        auto vd2 = boost::target(ed, graph);
        auto [_, ingraph] = boost::edge(vd1, vd2, graph);
        if (! ingraph) {
            raise<RuntimeError>("break_segment: segment not in graph\n");
        }

        const auto& fits = seg->fits();
        auto itfits = closest_point(fits, point, owp_to_point<Fit>);
        
        // std::cout << "Break point in system units: " << point << std::endl;
        // std::cout << "Break point in cm: (" << point.x()/units::cm << ", " << point.y()/units::cm << ", " << point.z()/units::cm << ")" << std::endl;
        
        // for (size_t i = 0; i < fits.size(); ++i) {
        //     const auto& fit = fits[i];
        //     double dist = (fit.point - point).magnitude();
        //     std::cout << "  Point " << i << ": position=(" 
        //                 << fit.point.x()/units::cm << ", " << fit.point.y()/units::cm << ", " << fit.point.z()/units::cm
        //                 << "), dQ=" << fit.dQ << ", dx=" << fit.dx/units::cm 
        //                 << " | dist=" << dist/units::cm << " cm" << std::endl;
        // }


        // reject if test point is at begin or end of fits: would create a
        // degenerate 1-point segment.  Mirrors prototype nbreak_fit check.
        if (itfits == fits.begin() || itfits+1 == fits.end()) {
            return std::make_tuple(false, std::pair<SegmentPtr, SegmentPtr>(), VertexPtr());
        }

        const auto& wcpts = seg->wcpts();        
        auto itwcpts = closest_point(wcpts, point, owp_to_point<WCPoint>);

        // clamp the wcpts to not be first/last
        if (itwcpts == wcpts.begin()) {
            ++itwcpts;
        }
        else if (itwcpts+1 == wcpts.end()) {
            --itwcpts;
        }

        SPDLOG_LOGGER_TRACE(s_log, "break_segment: Closest point found: {} / {} {} / {} points in fits", itfits - fits.begin(), fits.size(), itwcpts - wcpts.begin(), wcpts.size());

        
        // update graph
        remove_segment(graph, seg);

        auto vtx1 = graph[vd1].vertex;
        auto vtx2 = graph[vd2].vertex;
        auto vtx = make_vertex(graph);

        // WARNING: Boost graph edges have no inherent orientation — source(e)/target(e)
        // do NOT necessarily correspond to wcpts.front()/wcpts.back().  Callers that need
        // an oriented (start-vertex, end-vertex) pair should use find_vertices() in
        // PRGraph.cxx, which disambiguates by comparing the wcpts.front() distance to each
        // candidate vertex and returns (vertex nearest front, vertex nearest back).
        auto seg1 = make_segment(graph, vtx1, vtx);
        auto seg2 = make_segment(graph, vtx, vtx2);


        // fill in the new objects.  All three get the middle thing
        // Split wcpts - break point included in both
        seg1->wcpts(std::vector<WCPoint>(wcpts.begin(), itwcpts+1));
        seg2->wcpts(std::vector<WCPoint>(itwcpts, wcpts.end()));
        vtx->wcpt(*itwcpts);

        seg1->cluster(seg->cluster()); 
        seg2->cluster(seg->cluster());

        // Split fits - break point included in both
        if (fits.size()>0){
            seg1->fits(std::vector<Fit>(fits.begin(), itfits+1));
            seg2->fits(std::vector<Fit>(itfits, fits.end()));
            vtx->fit(*itfits);
            vtx->fit_range(1*units::cm);  // prototype sets fit_range=1cm on new vertex
        }

        // Copy segment properties from original to both new segments (matching WCPPID)
        seg1->dir_weak(seg->dir_weak());
        seg2->dir_weak(seg->dir_weak());

        seg1->dirsign(seg->dirsign());
        seg2->dirsign(seg->dirsign());

        // Copy all flags
        seg1->flags_set(seg->flags());
        seg2->flags_set(seg->flags());

        // Copy particle_score (separate field, not in flags)
        seg1->particle_score(seg->particle_score());
        seg2->particle_score(seg->particle_score());

            
        if (seg->has_particle_info()) {
            // Copy particle info with 4-momentum for seg1
            int pdg = seg->particle_info()->pdg();
            auto four_momentum1 = segment_cal_4mom(seg1, pdg, particle_data, recomb_model);
            auto pinfo1 = std::make_shared<Aux::ParticleInfo>(
                pdg,
                particle_data->get_particle_mass(pdg),
                particle_data->pdg_to_name(pdg),
                four_momentum1
            );
            seg1->particle_info(pinfo1);
            
            // Copy particle info with 4-momentum for seg2
            auto four_momentum2 = segment_cal_4mom(seg2, pdg, particle_data, recomb_model);
            auto pinfo2 = std::make_shared<Aux::ParticleInfo>(
                pdg,
                particle_data->get_particle_mass(pdg),
                particle_data->pdg_to_name(pdg),
                four_momentum2
            );
            seg2->particle_info(pinfo2);
        }

        // Copy dynamic point clouds if they exist
        // The dpcloud() method returns the DynamicPointCloud associated with a given name
        if (seg->dpcloud("fit")) {
            create_segment_fit_point_cloud(seg1, dv, "fit");
            create_segment_fit_point_cloud(seg2, dv, "fit");
        }
        
        if (seg->dpcloud("main")) {
            // Convert WCPoint to Point for create_segment_point_cloud
            const auto& wcpts1 = seg1->wcpts();
            const auto& wcpts2 = seg2->wcpts();
            std::vector<geo_point_t> points1, points2;
            points1.reserve(wcpts1.size());
            points2.reserve(wcpts2.size());
            for (const auto& wcp : wcpts1) {
                points1.push_back(wcp.point);
            }
            for (const auto& wcp : wcpts2) {
                points2.push_back(wcp.point);
            }
            create_segment_point_cloud(seg1, points1, dv, "main");
            create_segment_point_cloud(seg2, points2, dv, "main");
        }
        
        if (seg->dpcloud("associate_points")) {
            // Redistribute associated points based on closest distance to seg1 vs seg2
            // This matches WCPPID lines 123-140
            
            auto orig_dpc = seg->dpcloud("associate_points");
            const size_t orig_npts = orig_dpc->npoints();
            
            // Separate points based on which segment they're closer to
            // (row indices into orig_dpc; values copied column-wise below)
            std::vector<size_t> rows1, rows2;
            rows1.reserve(orig_npts / 2);  // estimate
            rows2.reserve(orig_npts / 2);
            
            // Determine which cloud to use for distance calculations
            // Prefer "fit" points, but fall back to "main" if "fit" is not available
            std::string ref_cloud_name = "fit";
            if (!seg1->dpcloud("fit") || !seg2->dpcloud("fit")) {
                ref_cloud_name = "main";
                // Ensure main clouds exist for both segments
                if (!seg1->dpcloud("main") || !seg2->dpcloud("main")) {
                    raise<RuntimeError>("break_segment: cannot redistribute associate_points - neither 'fit' nor 'main' clouds available");
                }
            }
            
            // Iterate through all associated points
            for (size_t ir = 0; ir < orig_npts; ++ir) {
                WireCell::Point point = orig_dpc->point3d(ir);
                
                // Compute closest distance to seg1 and seg2 using reference cloud
                auto [dist1, _1] = segment_get_closest_point(seg1, point, ref_cloud_name);
                auto [dist2, _2] = segment_get_closest_point(seg2, point, ref_cloud_name);
                
                // Add point to closer segment
                if (dist1 < dist2) {
                    rows1.push_back(ir);
                } else {
                    rows2.push_back(ir);
                }
            }
            
            // Get wpid_params from the original cloud
            auto& cluster = *seg->cluster();
            const auto& wpids = cluster.grouping()->wpids();
            std::map<WirePlaneId, std::tuple<geo_point_t, double, double, double>> wpid_params;
            std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_U_dir;
            std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_V_dir;
            std::map<WirePlaneId, std::pair<geo_point_t, double>> wpid_W_dir;
            std::set<int> apas;
            Facade::compute_wireplane_params(wpids, dv, wpid_params, wpid_U_dir, wpid_V_dir, wpid_W_dir, apas);

            // Create new DynamicPointClouds for associated points if we have any
            if (!rows1.empty()) {
                // Create and populate seg1's associate_points cloud
                auto dpc1 = std::make_shared<Facade::DynamicPointCloud>(wpid_params);
                dpc1->add_points(*orig_dpc, rows1);
                seg1->dpcloud("associate_points", dpc1);
            }
            
            if (!rows2.empty()) {
                // Create and populate seg2's associate_points cloud
                auto dpc2 = std::make_shared<Facade::DynamicPointCloud>(wpid_params);
                dpc2->add_points(*orig_dpc, rows2);
                seg2->dpcloud("associate_points", dpc2);
            }
            
            // Note: KD-tree indices are automatically rebuilt when add_points() is called
        }
        
       
            
        return std::make_tuple(true,std::make_pair(seg1, seg2), vtx);
    }


    double segment_track_length(SegmentPtr seg, int flag, int n1, int n2, WireCell::Vector dir_perp)
    {
        double length = 0;
        
        if (flag == 1) {
            // Sum dx values from fits (equivalent to original flag==1 case)
            auto& fits = seg->fits();
            if (n1>=0 && n2 >=0){
                n1 = std::max(0, n1);
                n2 = std::min(static_cast<int>(fits.size())-1, n2);
                for (int i = n1; i+1 <= n2; ++i) {
                    auto& fit = fits[i];
                    if (fit.valid() && fit.dx > 0) {
                        length += fit.dx;
                    }
                }
            }else{
                for (auto& fit : fits) {
                    if (fit.valid() && fit.dx > 0) {
                        length += fit.dx;
                    }
                }
            }
        } else {
            // Calculate geometric length from fits (equivalent to original flag==0 case)
            const auto& fits = seg->fits();
            if (fits.size() < 2) {
                return 0.0;
            }
            if (n1 >=0 && n2 >=0){
                n1 = std::max(0, n1);
                n2 = std::min(static_cast<int>(fits.size())-1, n2);
                for (int i = n1; i + 1 <= n2; i++) {
                    const Point& p1 = fits[i].point;
                    const Point& p2 = fits[i + 1].point;
                    WireCell::Vector segment_vec = p2 - p1;
                    if (dir_perp.magnitude() > 0) {
                        double mag_sq = segment_vec.magnitude() * segment_vec.magnitude();
                        double dot_sq = std::pow(segment_vec.dot(dir_perp.norm()), 2);
                        length += std::sqrt(mag_sq - dot_sq);
                    }else{
                        length += segment_vec.magnitude();
                    }
                }
            }else{
                for (size_t i = 0; i + 1 < fits.size(); i++) {
                    const Point& p1 = fits[i].point;
                    const Point& p2 = fits[i + 1].point;
                    WireCell::Vector segment_vec = p2 - p1;
                    if (dir_perp.magnitude() > 0) {
                        double mag_sq = segment_vec.magnitude() * segment_vec.magnitude();
                        double dot_sq = std::pow(segment_vec.dot(dir_perp.norm()), 2);
                        length += std::sqrt(mag_sq - dot_sq);
                    }else{
                        length += segment_vec.magnitude();
                    }
                }
            }
        }

        // if (std::isnan(length)) {
        //     const auto& dbg_fits = seg->fits();
        //     std::cout << "segment_track_length: NaN! nfits=" << dbg_fits.size() << std::endl;
        //     for (size_t i = 0; i < dbg_fits.size(); i++) {
        //         std::cout << "  fit[" << i << "] point=("
        //                   << dbg_fits[i].point.x() << "," << dbg_fits[i].point.y() << "," << dbg_fits[i].point.z()
        //                   << ") dx=" << dbg_fits[i].dx/units::cm << std::endl;
        //     }
        // }

        return length;
    }

    double segment_track_direct_length(SegmentPtr seg, int n1, int n2, WireCell::Vector dir_perp){
        double length = 0;
        
        const auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }

        if (n1<0 && n2 <0){
            n1 = 0;
            n2 = static_cast<int>(fits.size()) - 1;
        }
        
        // Clamp indices to valid range (following WCPPID logic)
        if (n1 < 0) n1 = 0;
        if (n1 >= static_cast<int>(fits.size())) n1 = static_cast<int>(fits.size()) - 1;
        if (n2 < 0) n2 = static_cast<int>(fits.size()) - 1;
        if (n2 >= static_cast<int>(fits.size())) n2 = static_cast<int>(fits.size()) - 1;

        const Point& p1 = fits[n1].point;
        const Point& p2 = fits[n2].point;
        WireCell::Vector temp_dir = p1 - p2;
        
        if (dir_perp.magnitude() > 0) {
            // Calculate length with perpendicular direction subtracted
            double mag_sq = temp_dir.magnitude() * temp_dir.magnitude();
            double dot_sq = std::pow(temp_dir.dot(dir_perp.norm()), 2);
            length = std::sqrt(mag_sq - dot_sq);
        } else {
            // Simple direct distance
            length = temp_dir.magnitude();
        }

        return length;
    }

    double segment_track_max_deviation(SegmentPtr seg, int n1, int n2){
        double max_deviation = 0.0;

        const auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        if (n1<0 && n2 <0){
            n1 = 0;
            n2 = static_cast<int>(fits.size()) - 1;
        }

        // Handle default values and clamp indices (following WCPPID logic)
        if (n1 < 0) n1 = 0;
        if (n1 >= static_cast<int>(fits.size())) n1 = static_cast<int>(fits.size()) - 1;
        if (n2 < 0) n2 = static_cast<int>(fits.size()) - 1;
        if (n2 >= static_cast<int>(fits.size())) n2 = static_cast<int>(fits.size()) - 1;
        
        // Ensure n1 <= n2
        if (n1 > n2) std::swap(n1, n2);
        
        if (n1 != n2) {
            const Point& p1 = fits[n1].point;
            const Point& p2 = fits[n2].point;
            WireCell::Vector line_vec = p2 - p1;
            double line_length_sq = line_vec.magnitude2();
            
            for (int i = n1; i <= n2; i++) {
                const Point& test_point = fits[i].point;
                
                if (line_length_sq > 0) {
                    // Calculate distance from point to line using projection
                    WireCell::Vector point_vec = test_point - p1;
                    double projection = point_vec.dot(line_vec) / line_length_sq;
                    
                    // Clamp projection to line segment bounds
                    projection = std::max(0.0, std::min(1.0, projection));
                    
                    // Find closest point on line segment
                    Point closest_on_line = p1 + line_vec * projection;
                    
                    // Calculate distance
                    double distance = (test_point - closest_on_line).magnitude();
                    
                    if (distance > max_deviation) {
                        max_deviation = distance;
                    }
                } else {
                    // Line has zero length, distance is just point-to-point distance
                    double distance = (test_point - p1).magnitude();
                    if (distance > max_deviation) {
                        max_deviation = distance;
                    }
                }
            }
        }
        
        return max_deviation;
    }
        
        


    double segment_median_dQ_dx(SegmentPtr seg, int n1, int n2)
    {
        auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        // Handle default parameters (equivalent to get_medium_dQ_dx())
        if (n1 < 0 && n2 < 0) {
            n1 = 0;
            n2 = static_cast<int>(fits.size());
        }
        
        // Clamp indices to valid range (equivalent to WCPPID bounds checking)
        if (n1 < 0) n1 = 0;
        if (n1 + 1 > static_cast<int>(fits.size())) n1 = static_cast<int>(fits.size()) - 1;
        if (n2 < 0) n2 = 0;
        if (n2 + 1 > static_cast<int>(fits.size())) n2 = static_cast<int>(fits.size()) - 1;
        
        std::vector<double> vec_dQ_dx;
        vec_dQ_dx.reserve(n2 - n1 + 1);
        
        // Loop over specified range [n1, n2] (inclusive, matching WCPPID)
        for (int i = n1; i <= n2 && i < static_cast<int>(fits.size()); i++) {
            auto& fit = fits[i];
            if (fit.valid() && fit.dx > 0 && fit.dQ >= 0) {
                // Add small epsilon to avoid division by zero (same as WCPPID: 1e-9)
                vec_dQ_dx.push_back(fit.dQ / (fit.dx + 1e-9));
            }
        }
        
        if (vec_dQ_dx.empty()) {
            return 0.0;
        }
        
        // Use nth_element to find median (exact WCPPID algorithm)
        size_t median_index = vec_dQ_dx.size() / 2;
        std::nth_element(vec_dQ_dx.begin(), 
                        vec_dQ_dx.begin() + median_index, 
                        vec_dQ_dx.end());
        
        return vec_dQ_dx[median_index];
    }

    // doc sbnd_xin/docs/pr/40 -- see the header comment.
    bool segment_dqdx_spares_electron_reclass(SegmentPtr seg, double MIP_dQdx) {
        if (MIP_dQdx <= 0) return false;
        const double median = segment_median_dQ_dx(seg);
        if (median <= 0) return false;  // no evidence, not MIP-like evidence
        const double ratio = median / MIP_dQdx;
        return (ratio > 1.75) || (ratio < 1.2);
    }

    // doc sbnd_xin/docs/pr/40 round 2 F5 -- see the header comment.
    bool segment_has_proton_daughter(Graph& graph, SegmentPtr seg, VertexPtr main_vertex, double MIP_dQdx) {
        static const bool dbg = std::getenv("WCT_PROTON_DAUGHTER_DEBUG") != nullptr;
        if (!main_vertex || !seg || MIP_dQdx <= 0) {
            if (dbg) std::fprintf(stderr, "PROTON_DAUGHTER_DEBUG seg=%p: early false (main_vertex=%p seg=%p MIP=%g)\n",
                                   (void*)seg.get(), (void*)main_vertex.get(), (void*)seg.get(), MIP_dQdx);
            return false;
        }

        auto [v1, v2] = find_vertices(graph, seg);
        if (dbg) std::fprintf(stderr, "PROTON_DAUGHTER_DEBUG seg=%p (clus=%d idx=%zu): v1=%p v2=%p main_vertex=%p\n",
                               (void*)seg.get(), seg->cluster() ? seg->cluster()->get_cluster_id() : -1,
                               seg->get_graph_index(), (void*)v1.get(), (void*)v2.get(), (void*)main_vertex.get());
        if (!v1 || !v2) return false;

        VertexPtr far_v = nullptr;
        if (v1 == main_vertex) far_v = v2;
        else if (v2 == main_vertex) far_v = v1;
        else {
            if (dbg) std::fprintf(stderr, "PROTON_DAUGHTER_DEBUG seg=%p: neither endpoint IS main_vertex\n", (void*)seg.get());
            return false;  // does not emanate from the neutrino vertex
        }

        if (!far_v->descriptor_valid()) {
            if (dbg) std::fprintf(stderr, "PROTON_DAUGHTER_DEBUG seg=%p: far_v descriptor invalid\n", (void*)seg.get());
            return false;
        }

        int nnbr = 0;
        for (auto edesc : sorted_out_edges(far_v->get_descriptor(), graph)) {
            SegmentPtr nbr = graph[edesc].segment;
            if (!nbr || nbr == seg) continue;
            ++nnbr;
            const bool has_pi = nbr->has_particle_info();
            const int pdg = has_pi ? nbr->particle_info()->pdg() : 0;
            const double median = segment_median_dQ_dx(nbr);
            if (dbg) std::fprintf(stderr, "PROTON_DAUGHTER_DEBUG seg=%p: nbr=%p (clus=%d idx=%zu) has_pi=%d pdg=%d median/MIP=%g\n",
                                   (void*)seg.get(), (void*)nbr.get(),
                                   nbr->cluster() ? nbr->cluster()->get_cluster_id() : -1, nbr->get_graph_index(),
                                   has_pi, pdg, median / MIP_dQdx);
            if (!has_pi || pdg != 2212) continue;
            if (median > 0 && median / MIP_dQdx > 1.75) return true;
        }
        if (dbg) std::fprintf(stderr, "PROTON_DAUGHTER_DEBUG seg=%p: %d neighbour(s), none qualified\n", (void*)seg.get(), nnbr);
        return false;
    }

    // doc sbnd_xin/docs/pr/43 F1 -- a muon candidate is disqualified if a
    // charge-confirmed proton terminates its own decay chain, even when the
    // proton sits beyond the candidate's immediate far vertex.  The existing
    // 1-hop `n_proton` check in the muon-candidate loop
    // (NeutrinoVertexFinder.cxx examine_direction) and
    // segment_has_proton_daughter above both look only at the segment's own
    // far vertex; run 18255 evt 54351 has muon(54cm) -> muon-stub(2.7cm) ->
    // proton(13cm), i.e. the proton is TWO hops from the candidate's far
    // vertex, through a short collinear continuation segment neither check
    // reaches.
    //
    // `near_vertex` is the vertex `seg` is evaluated FROM (the vertex the
    // single-muon-per-cluster selection runs at).  The walk starts at seg's
    // OTHER endpoint and follows through at most `max_hops` further
    // segments, each of which must itself be a non-shower track
    // continuation (pdg 13, pdg 0/undetermined, or already pdg 211 -- never
    // a segment already flagged kShowerTrajectory/kShowerTopology or pdg
    // 11) at a simple degree-2 vertex (one segment in, one out).  A vertex
    // with any other branching (a genuine hadronic multi-prong vertex, or a
    // shower attaching) stops the walk rather than being treated as "the
    // same chain" -- this keeps the walk from reaching an unrelated proton
    // sitting past a real vertex activity.  Returns true on the first
    // charge-confirmed proton found (median dQ/dx > 1.75x MIP_dQdx, the
    // same threshold segment_has_proton_daughter uses).
    bool segment_chain_has_proton(Graph& graph, SegmentPtr seg, VertexPtr near_vertex, double MIP_dQdx, int max_hops) {
        if (!near_vertex || !seg || MIP_dQdx <= 0 || max_hops <= 0) return false;

        auto [v1, v2] = find_vertices(graph, seg);
        if (!v1 || !v2) return false;
        VertexPtr far_v = nullptr;
        if (v1 == near_vertex) far_v = v2;
        else if (v2 == near_vertex) far_v = v1;
        else return false;  // seg does not emanate from near_vertex

        SegmentPtr prev_seg = seg;
        VertexPtr cur_v = far_v;
        for (int hop = 0; hop < max_hops; ++hop) {
            if (!cur_v || !cur_v->descriptor_valid()) return false;

            SegmentPtr next_seg = nullptr;
            int nnbr = 0;
            for (auto edesc : sorted_out_edges(cur_v->get_descriptor(), graph)) {
                SegmentPtr nbr = graph[edesc].segment;
                if (!nbr || nbr == prev_seg) continue;
                ++nnbr;
                const bool has_pi = nbr->has_particle_info();
                const int pdg = has_pi ? nbr->particle_info()->pdg() : 0;
                if (has_pi && pdg == 2212) {
                    const double median = segment_median_dQ_dx(nbr);
                    if (median > 0 && median / MIP_dQdx > 1.75) return true;
                }
                const bool is_shower = nbr->flags_any(SegmentFlags::kShowerTrajectory) ||
                                        nbr->flags_any(SegmentFlags::kShowerTopology) ||
                                        (has_pi && std::abs(pdg) == 11);
                if (!is_shower && (!has_pi || pdg == 13 || pdg == 211)) next_seg = nbr;
            }
            // Only a simple degree-2 continuation (exactly one other
            // segment at this vertex, and it looks like more chain) is
            // followed further; anything else -- a branching vertex, a
            // dead end, a shower -- stops the walk here.
            if (nnbr != 1 || !next_seg) return false;

            auto [w1, w2] = find_vertices(graph, next_seg);
            VertexPtr next_v = nullptr;
            if (w1 == cur_v) next_v = w2;
            else if (w2 == cur_v) next_v = w1;
            else return false;

            prev_seg = next_seg;
            cur_v = next_v;
        }
        return false;
    }

    // doc sbnd_xin/docs/pr/43 F1 -- collect the same short, non-shower,
    // degree-2 continuation chain segment_chain_has_proton walks, for
    // relabeling every segment in a disqualified muon candidate's own chain
    // to pion.  Without this, demoting only the head segment (e.g. 17007 in
    // evt 54351) leaves an inconsistent muon stub (17005) between it and
    // the proton that disqualified it -- the owner's report explicitly
    // wants the WHOLE chain read as pion, not just its head.  The proton
    // itself (and anything past it) is never included: `next_seg`
    // deliberately excludes pdg 2212, so the walk stops one hop short of
    // relabeling the proton.  Returns segments in walk order (nearest
    // first); empty if `seg` does not emanate from `near_vertex` or no
    // continuation exists.
    std::vector<SegmentPtr> segment_chain_continuation(Graph& graph, SegmentPtr seg, VertexPtr near_vertex, int max_hops) {
        std::vector<SegmentPtr> out;
        if (!near_vertex || !seg || max_hops <= 0) return out;

        auto [v1, v2] = find_vertices(graph, seg);
        if (!v1 || !v2) return out;
        VertexPtr far_v = nullptr;
        if (v1 == near_vertex) far_v = v2;
        else if (v2 == near_vertex) far_v = v1;
        else return out;

        SegmentPtr prev_seg = seg;
        VertexPtr cur_v = far_v;
        for (int hop = 0; hop < max_hops; ++hop) {
            if (!cur_v || !cur_v->descriptor_valid()) return out;

            SegmentPtr next_seg = nullptr;
            int nnbr = 0;
            for (auto edesc : sorted_out_edges(cur_v->get_descriptor(), graph)) {
                SegmentPtr nbr = graph[edesc].segment;
                if (!nbr || nbr == prev_seg) continue;
                ++nnbr;
                const bool has_pi = nbr->has_particle_info();
                const int pdg = has_pi ? nbr->particle_info()->pdg() : 0;
                const bool is_shower = nbr->flags_any(SegmentFlags::kShowerTrajectory) ||
                                        nbr->flags_any(SegmentFlags::kShowerTopology) ||
                                        (has_pi && std::abs(pdg) == 11);
                if (!is_shower && pdg != 2212 && (!has_pi || pdg == 13 || pdg == 211)) next_seg = nbr;
            }
            if (nnbr != 1 || !next_seg) return out;

            out.push_back(next_seg);

            auto [w1, w2] = find_vertices(graph, next_seg);
            VertexPtr next_v = nullptr;
            if (w1 == cur_v) next_v = w2;
            else if (w2 == cur_v) next_v = w1;
            else return out;

            prev_seg = next_seg;
            cur_v = next_v;
        }
        return out;
    }

    // doc sbnd_xin/docs/pr/40 round 4 F8 -- see the header comment.
    bool segment_at_multi_proton_vertex(Graph& graph, SegmentPtr seg, VertexPtr main_vertex, double MIP_dQdx, int min_protons) {
        if (!main_vertex || !seg || MIP_dQdx <= 0 || min_protons <= 0) return false;

        auto [v1, v2] = find_vertices(graph, seg);
        if (!v1 || !v2) return false;

        for (VertexPtr end_v : {v1, v2}) {
            if (end_v == main_vertex) continue;  // the neutrino vertex is exempt
            if (!end_v->descriptor_valid()) continue;

            int nprotons = 0;
            for (auto edesc : sorted_out_edges(end_v->get_descriptor(), graph)) {
                SegmentPtr nbr = graph[edesc].segment;
                if (!nbr || nbr == seg) continue;
                if (!nbr->has_particle_info() || nbr->particle_info()->pdg() != 2212) continue;
                const double median = segment_median_dQ_dx(nbr);
                if (median > 0 && median / MIP_dQdx > 1.75) ++nprotons;
            }
            if (nprotons >= min_protons) return true;
        }
        return false;
    }

    // doc sbnd_xin/docs/pr/40 round 5 F10/F11 -- see the header comment.
    bool segment_is_straight_long_track(SegmentPtr seg, double min_length, double min_direct, double straight_ratio) {
        if (!seg) return false;
        const double length = segment_track_length(seg);
        if (length <= min_length) return false;
        const double direct_length = segment_track_direct_length(seg);
        return direct_length >= min_direct || direct_length > straight_ratio * length;
    }

    double segment_rms_dQ_dx(SegmentPtr seg)
    {
        auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        std::vector<double> vec_dQ_dx;
        vec_dQ_dx.reserve(fits.size());
        
        for (auto& fit : fits) {
            if (fit.valid() && fit.dx > 0 && fit.dQ >= 0) {
                // Add small epsilon to avoid division by zero (same as original)
                vec_dQ_dx.push_back(fit.dQ / (fit.dx + 1e-9));
            }
        }
        
        if (vec_dQ_dx.size() <= 1) {
            return 0.0;
        }

        // Calculate mean
        double sum = std::accumulate(vec_dQ_dx.begin(), vec_dQ_dx.end(), 0.0);
        double mean = sum / vec_dQ_dx.size();

        // Calculate sample variance (Bessel-corrected, divide by N-1, matches prototype)
        double sq_sum = std::inner_product(vec_dQ_dx.begin(), vec_dQ_dx.end(), vec_dQ_dx.begin(), 0.0);
        double variance = (sq_sum - mean * mean * vec_dQ_dx.size()) / (vec_dQ_dx.size() - 1);

        return std::sqrt(variance);
    }


    double segment_track_length_threshold(SegmentPtr seg, double threshold)
    {
        auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        double length = 0;
        for (auto& fit : fits) {
            if (fit.valid() && fit.dx > 0 ) {
                double dQ_dx = fit.dQ / (fit.dx + 1e-9); // Add epsilon to avoid division by zero
                if (dQ_dx > threshold || threshold == 0) {
                    length += fit.dx;
                }
            }
        }
        
        return length;
    }




   

    double segment_geometric_length(SegmentPtr seg, int n1, int n2, WireCell::Vector dir_perp)
    {
        return segment_track_length(seg, 0, n1, n2, dir_perp); // Always use geometric calculation
    }


    bool eval_ks_ratio(double ks1, double ks2, double ratio1, double ratio2){
        //  std::cout << ks1 << " " << ks2 << " " << ratio1 << " " << ratio2 << " " << sqrt(pow(ks2/0.06,2)+pow((ratio2-1)/0.06,2)) << " " << ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 << " " <<  ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 << " " << std::endl;
        if (ks1-ks2 >= 0.0) return false;
        if (sqrt(pow(ks2/0.06,2)+pow((ratio2-1)/0.06,2))< 1.4 && ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 > -0.02) return false;

        if (ks1 - ks2 < -0.02 && ((ks2 > 0.09 && fabs(ratio2-1) >0.1) || ratio2 > 1.5 || ks2 > 0.2)) return true;
        if ( ks1-ks2 + (fabs(ratio1-1)-fabs(ratio2-1))/1.5*0.3 < 0) return true;

        return false;
    }

    /// doc sbnd_xin/docs/pr/32 §11 F2.  false = today's set-only behaviour.
    bool g_shower_traj_refresh_flag = false;

    bool segment_is_shower_trajectory(SegmentPtr seg, double step_size, double mip_dQ_dx, bool straight_guard){
        bool flag_shower_trajectory = false;
        // doc pr/32 §11 F2: the prototype opens is_shower_trajectory() with
        // `flag_shower_trajectory = false` (ProtoSegment.cxx:544), BEFORE its
        // own `length > 50 cm` early return, so even the early-out path clears
        // the label.  Mirror that placement exactly.
        // f2_flag_cleared counts NET demotions -- the bit was set on entry and
        // this call did not put it back -- not every clear, because the common
        // case is clear-then-set-again and that is not a behaviour change.
        // Cheap enough to evaluate unconditionally; used only under the knob.
        const bool traj_was_set = seg->flags_any(SegmentFlags::kShowerTrajectory);
        if (g_shower_traj_refresh_flag) {
            seg->unset_flags(SegmentFlags::kShowerTrajectory);
        }
        // Count a net demotion at every exit that returns false.
        struct DemotionCounter {
            bool armed;
            const bool* result;
            ~DemotionCounter() {
                if (armed && !*result) {
                    g_pr32_audit.f2_flag_cleared.fetch_add(1, std::memory_order_relaxed);
                }
            }
        } demotion_counter{g_shower_traj_refresh_flag && traj_was_set, &flag_shower_trajectory};
        double length = segment_track_length(seg, 0);

        // Too long
        if (length > 50 * units::cm) return flag_shower_trajectory;
        
        const auto& fits = seg->fits();
        if (fits.empty()) return flag_shower_trajectory;
        
        int ncount = std::round(length / step_size);
        if (ncount == 0) ncount = 1;
        
        std::vector<std::pair<int,int>> sections(ncount);
        for (int i = 0; i < ncount; i++) {
            sections[i] = std::make_pair(
                std::round(fits.size() / ncount * i),
                std::round(fits.size() / ncount * (i + 1))
            );
        }
        sections.back().second = fits.size() - 1;
        
        int n_shower_like = 0;
        WireCell::Vector drift_dir(1, 0, 0);
        
        for (int j = 0; j < ncount; j++) {
            int first_idx = sections[j].first;
            int second_idx = sections[j].second;
            
            if (first_idx >= static_cast<int>(fits.size())) first_idx = fits.size() - 1;
            if (second_idx >= static_cast<int>(fits.size())) second_idx = fits.size() - 1;
            
            WireCell::Vector dir_1 = fits[first_idx].point - fits[second_idx].point;
            if (dir_1.magnitude() > 0) {
                dir_1 = dir_1.norm();
            }
            
            double tmp_dQ_dx = segment_median_dQ_dx(seg, first_idx, second_idx) / (mip_dQ_dx);
            
            // Calculate angle difference
            double dot_product = drift_dir.dot(dir_1);
            double angle_rad = std::acos(std::max(-1.0, std::min(1.0, dot_product)));
            double angle_diff = std::abs(angle_rad / M_PI * 180.0 - 90.0);
            
            if (angle_diff > 10) { // Not parallel case
                double direct_length = segment_track_direct_length(seg, first_idx, second_idx, WireCell::Vector(0,0,0));
                double integrated_length = segment_track_length(seg, 0, first_idx, second_idx, WireCell::Vector(0,0,0));
                double max_dev = segment_track_max_deviation(seg, first_idx, second_idx);
                
                double length_ratio;
                if (direct_length == 0) length_ratio = 1;
                else length_ratio = direct_length / integrated_length;
                
                if (tmp_dQ_dx * 0.11 + 2 * length_ratio < 2.03 && 
                    tmp_dQ_dx < 2 && 
                    length_ratio < 0.95 && 
                    (angle_diff < 60 || integrated_length < 10 * units::cm || 
                     (integrated_length >= 10 * units::cm && max_dev > 0.75 * units::cm))) {
                    n_shower_like++;
                }
            } else { // Parallel case
                WireCell::Vector dir_2 = drift_dir.cross(dir_1);
                if (dir_2.magnitude() > 0) {
                    dir_2 = dir_2.norm();
                }
                
                double direct_length = segment_track_direct_length(seg, first_idx, second_idx, dir_2);
                double integrated_length = segment_track_length(seg, 0, first_idx, second_idx, dir_2);
                double max_dev = segment_track_max_deviation(seg, first_idx, second_idx);
                
                double length_ratio;
                if (direct_length == 0) length_ratio = 1;
                else length_ratio = direct_length / integrated_length;
                
                if (tmp_dQ_dx * 0.11 + 2 * length_ratio < 2.06 && 
                    tmp_dQ_dx < 2 && 
                    length_ratio < 0.97 && 
                    (integrated_length < 10 * units::cm || 
                     (integrated_length >= 10 * units::cm && max_dev > 0.75 * units::cm))) {
                    n_shower_like++;
                }
            }
        }
        
        if (n_shower_like >= 0.5 * sections.size()) {
            flag_shower_trajectory = true;
        }

        // doc sbnd_xin/docs/pr/40 round 5 F11: the wiggliness test above has
        // no length/straightness cross-check (unlike segment_is_shower_
        // topology, which F3 already guards with a dQ/dx check).  A long,
        // straight track's own dQ/dx can still sit under the section-level
        // thresholds used above; straightness catches what dQ/dx alone does
        // not (SBND evt 55715 seg 15007: 1.70x MIP, just under the 1.75x
        // proton-like cut segment_dqdx_spares_electron_reclass uses).
        // false = legacy = byte-identical.
        if (flag_shower_trajectory && straight_guard && segment_is_straight_long_track(seg)) {
            flag_shower_trajectory = false;
        }

        // Set the flag on the segment if it's identified as shower trajectory
        if (flag_shower_trajectory) {
            seg->set_flags(SegmentFlags::kShowerTrajectory);
        }

        return flag_shower_trajectory;
    }

    bool segment_is_dir_weak(SegmentPtr seg)
    {
        // Check particle-type-based score thresholds (matches prototype is_dir_weak())
        auto pinfo = seg->particle_info();
        if (pinfo) {
            int pdg = std::abs(pinfo->pdg());
            double score = seg->particle_score();
            double length = segment_track_length(seg, 0);  // geometric length from fit points
            if (pdg == 13) {  // muon/antimuon
                if (score > 0.07 && length >= 5*units::cm) return true;
                if (score > 0.15 && length <  5*units::cm) return true;
            }
            if (pdg == 2212) {  // proton
                if (score > 0.13 && length >= 5*units::cm) return true;
                if (score > 0.27 && length <  5*units::cm) return true;
            }
        }
        // Fall through to static flag (set for electrons, very short tracks, etc.)
        return seg->dir_weak();
    }

    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg){
        const auto& fits = seg->fits();
        int flag_dir = seg->dirsign();

        if (fits.size() < 2) {
            SPDLOG_LOGGER_TRACE(s_log,
                "segment_cal_dir_3vector: seg id={} nfits={} dirsign={} — too few fits, returning (0,0,0)",
                seg->id(), fits.size(), flag_dir);
            return WireCell::Vector(0, 0, 0);
        }

        WireCell::Point p(0, 0, 0);

        if (flag_dir == 1) {
            // Forward direction: from first point using next few points
            for (size_t i = 1; i < 5 && i < fits.size(); i++) {
                p = p + (fits[i].point - fits[0].point);
            }
        } else if (flag_dir == -1) {
            // Backward direction: from last point using previous few points
            for (size_t i = 1; i < 5 && (fits.size() - i - 1) < fits.size(); i++) {
                if (fits.size() - i - 1 < fits.size()) {
                    p = p + (fits[fits.size() - i - 1].point - fits.back().point);
                }
            }
        } else {
            // flag_dir == 0: direction undetermined, return zero vector (matches prototype)
            SPDLOG_LOGGER_TRACE(s_log,
                "segment_cal_dir_3vector: seg id={} nfits={} dirsign=0 — direction undetermined, returning (0,0,0)",
                seg->id(), fits.size());
            return WireCell::Vector(0, 0, 0);
        }
        
        WireCell::Vector v1(p.x(), p.y(), p.z());
        if (v1.magnitude() > 0) {
            v1 = v1.norm();
        }
        return v1;
    }
    
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg, WireCell::Point& p, double dis_cut){
        const auto& fits = seg->fits();
        if (fits.empty()) {
            return WireCell::Vector(0, 0, 0);
        }
        
        WireCell::Point p1(0, 0, 0);
        int ncount = 0;
        
        for (size_t i = 0; i < fits.size(); i++) {
            double dis = (fits[i].point - p).magnitude();
            if (dis < dis_cut) {
                p1 = p1 + fits[i].point;
                ncount++;
            }
        }
        
        if (ncount == 0) {
            return WireCell::Vector(0, 0, 0);
        }
        
        WireCell::Point avg_point = p1 * (1.0 / ncount);
        WireCell::Vector v1 = avg_point - p;
        if (v1.magnitude() > 0) {
            v1 = v1.norm();
        }
        return v1;
    }
    
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg, int direction, int num_points, int start){
        const auto& fits = seg->fits();
        if (fits.empty() || start >= static_cast<int>(fits.size()) || start <= 0) {
            SPDLOG_LOGGER_TRACE(s_log, "segment_cal_dir_3vector: bad start point in segment_cal_dir_3vector");
            return WireCell::Vector(0, 0, 0);
        }
        
        WireCell::Point p(0, 0, 0);
        
        if (direction == 1) {
            // Forward direction
            for (int i = start; i < start + num_points - 1 && i < static_cast<int>(fits.size()); i++) {
                p = p + (fits[i].point - fits[start - 1].point);
            }
        } else if (direction == -1) {
            // Backward direction
             for (int i = start; i < start + num_points - 1; i++) {
                // WCPPID's bounds check
                if (i + 1 > static_cast<int>(fits.size())) break;
                
                // Ensure backward indices are valid
                int back_idx = fits.size() - i - 1;
                int ref_idx = fits.size() - start;
                
                if (back_idx >= 0 && back_idx < static_cast<int>(fits.size()) && 
                    ref_idx >= 0 && ref_idx < static_cast<int>(fits.size())) {
                    p = p + (fits[back_idx].point - fits[ref_idx].point);
                }
            }
        }
        
        WireCell::Vector v1(p.x(), p.y(), p.z());
        if (v1.magnitude() > 0) {
            v1 = v1.norm();
        }
        return v1;
    }

    // NOTE (multi-APA): recomb_model is a single global instance shared across all fit points
    // regardless of their per-fit paf{apa,face}.  For uBooNE (single APA+face) this is correct.
    // For multi-APA detectors (SBND, DUNE) a per-face electron-lifetime correction and a
    // per-face recombination model would be needed; callers must pre-apply lifetime corrections
    // to dQ before this function is called.  This is a known limitation inherited from the
    // prototype and deliberately deferred — see pid_direction_kinematics_review.md §0 G3.
    double segment_cal_kine_dQdx(SegmentPtr seg, const IRecombinationModel::pointer& recomb_model){
        if (!seg || !recomb_model) {
            return 0.0;
        }
        
        auto& fits = seg->fits();
        if (fits.empty()) {
            return 0.0;
        }
        
        double kine_energy = 0.0;
        
        
        for (size_t i = 0; i < fits.size(); i++) {
            if (!fits[i].valid() || fits[i].dx <= 0) continue;
            
            double dX = fits[i].dx;  // path length used for energy accumulation (may be shortened)
            double dQ = fits[i].dQ;
            // For the first and last fit points the fitted dx can be anomalously large because
            // the fit window extends past the segment endpoint.  When dx > 1.5× the inter-point
            // distance, use the inter-point distance as the path length instead.
            // IMPORTANT: only shorten the *accumulation* path dX, not the path used for dQ/dx —
            // the Box model is non-linear so computing dE/dx from a shortened dx would inflate
            // the dQ/dx and thus bias the energy upward.  Match the prototype (ProtoSegment.cxx:1356):
            //   dEdx = f(dQ/fits[i].dx);  kine_energy += dEdx * dX
            double dx_for_dQdx = fits[i].dx;  // always the true fitted path for dQ/dx computation
            if (i == 0 && fits.size() > 1) {
                // First point: check against distance to next point
                double dis = (fits[1].point - fits[0].point).magnitude();
                if (dis > 0 && dX > dis * 1.5) {
                    dX = dis;
                }
            } else if (i + 1 == fits.size() && fits.size() > 1) {
                // Last point: check against distance to previous point
                double dis = (fits[i].point - fits[i-1].point).magnitude();
                if (dis > 0 && dX > dis * 1.5) {
                    dX = dis;
                }
            }
            if (dX <= 0) continue;
            // std::cout << i << " " << fits[i].dQ << " " << fits[i].dx/units::cm << " " << dX/units::cm << std::endl;
            // Filter out unreasonable values using the true fitted dx (matches prototype)
            if (dQ/dx_for_dQdx / (43e3/units::cm) > 1000) dQ = 0;

            // Compute dE using the true fitted dx; accumulate over the (possibly shortened) dX
            double dE_per_dx = recomb_model->dE(dQ, dx_for_dQdx) / dx_for_dQdx;
            double dE = dE_per_dx * dX;

            // std::cout << dQ << " " << dX << " " << dE << std::endl;

            // Clamp to [0, 50 MeV/cm * dX]
            if (dE < 0) dE = 0;
            if (dE > 50 * units::MeV / units::cm * dX) dE = 50 * units::MeV / units::cm * dX;

            kine_energy += dE;
        }

        return kine_energy;
    }

    double cal_kine_dQdx(std::vector<double>& vec_dQ, std::vector<double>& vec_dx, const IRecombinationModel::pointer& recomb_model){
        if (vec_dQ.size() != vec_dx.size() || vec_dQ.empty() || !recomb_model) {
            return 0.0;
        }
        
        double kine_energy = 0.0;
        
        for (size_t i = 0; i < vec_dQ.size(); i++) {
              // Calculate dQ/dx with units conversion (same as original)
            double dQ = vec_dQ[i];
            double dx = vec_dx[i];
            
            // Filter out unreasonable values (same threshold as original)
            if (dQ/dx / (43e3/units::cm) > 1000) dQ = 0;
            
            // Calculate dE/dx using Box model inverse formula from original code
            double dE = recomb_model->dE(dQ, dx);

            // double dQp = (*recomb_model)(dE, dx);

            // std::cout << dQ << " " << dx << " " << dE << " " << units::MeV << " " << dQp << std::endl;
            
            // Apply bounds (same as original)
            if (dE < 0) dE = 0;
            if (dE > 50 * units::MeV / units::cm * dx) dE = 50 * units::MeV / units::cm * dx;

            // Calculate path length with special handling for first and last points
            kine_energy += dE;
        }
        
        return kine_energy;
    }

    std::vector<double> do_track_comp(std::vector<double>& L , std::vector<double>& dQ_dx, double compare_range, double offset_length, const Clus::ParticleDataSet::pointer& particle_data,  double MIP_dQdx, int skip_stop_samples, bool empty_abstain){

        double end_L = L.back() + 0.15*units::cm - offset_length;

        int ncount = 0;
        std::vector<double> vec_x;
        std::vector<double> vec_y;

        for (size_t i = 0; i != L.size(); i++) {
            if (end_L - L.at(i) < compare_range && end_L - L.at(i) > 0) { // check up to compared range
                vec_x.push_back(end_L - L.at(i));
                vec_y.push_back(dQ_dx.at(i));
                ncount++;
            }
        }

        // Endpoint-trim retry (pr/9 F1): samples were pushed in increasing-L
        // order, so the last entries are the ones nearest the hypothesized
        // stop.  Drop them without moving end_L (the template anchor).  Only
        // trim when at least one sample would survive.
        if (skip_stop_samples > 0 && ncount > skip_stop_samples) {
            vec_x.resize(vec_x.size() - skip_stop_samples);
            vec_y.resize(vec_y.size() - skip_stop_samples);
            ncount -= skip_stop_samples;
        }

        // If no points fall inside the comparison window, return "no direction signal" defaults.
        // doc sbnd_xin/docs/pr/31 §12 (F6, was P7): element 0 is read as
        // flag_forward = round(result.at(0)), so the legacy 1.0 means "this
        // orientation PASSED the direction gate" -- a segment about which
        // nothing is known becomes a directed muon on the both-pass branch.
        // empty_abstain returns 0.0 instead, the prototype's degenerate
        // answer (executed, not inferred: zero-bin TH1F KolmogorovTest gives
        // ks1 == ks2 == 0, so eval_ks_ratio's first gate returns false).
        // Elements 1-3 keep the 1e9 filler either way -- the finding is the
        // direction flag, not the particle scores.  The same literal is
        // applied at the missing-dEdx-function return below: it is the same
        // copy-pasted filler with the same consumer.
        if (ncount == 0) {
            return {empty_abstain ? 0.0 : 1.0, 1e9, 1e9, 1e9};
        }

        auto muon_fn     = particle_data->get_dEdx_function("muon");
        auto proton_fn   = particle_data->get_dEdx_function("proton");
        auto electron_fn = particle_data->get_dEdx_function("electron");
        if (!muon_fn || !proton_fn || !electron_fn) {
            return {empty_abstain ? 0.0 : 1.0, 1e9, 1e9, 1e9};
        }

        // Create reference vectors for different particles
        const size_t count = static_cast<size_t>(ncount);
        std::vector<double> muon_ref(count);
        std::vector<double> const_ref(count, MIP_dQdx);  // MIP-like constant
        std::vector<double> proton_ref(count);
        std::vector<double> electron_ref(count);

        for (size_t i = 0; i < count; i++) {
            muon_ref[i]     = muon_fn->scalar_function((vec_x[i])/units::cm) / units::cm;
            proton_ref[i]   = proton_fn->scalar_function((vec_x[i])/units::cm) / units::cm;
            electron_ref[i] = electron_fn->scalar_function((vec_x[i])/units::cm) / units::cm;
        }
        
        // Perform KS-like tests using kslike_compare
        double ks1 = WireCell::kslike_compare(vec_y, muon_ref);
        double ratio1 = std::accumulate(muon_ref.begin(), muon_ref.end(), 0.0) / 
                        (std::accumulate(vec_y.begin(), vec_y.end(), 0.0) + 1e-9);
        
        double ks2 = WireCell::kslike_compare(vec_y, const_ref);
        double ratio2 = std::accumulate(const_ref.begin(), const_ref.end(), 0.0) / 
                        (std::accumulate(vec_y.begin(), vec_y.end(), 0.0) + 1e-9);
        
        double ks3 = WireCell::kslike_compare(vec_y, proton_ref);
        double ratio3 = std::accumulate(proton_ref.begin(), proton_ref.end(), 0.0) / 
                        (std::accumulate(vec_y.begin(), vec_y.end(), 0.0) + 1e-9);
        
        double ks4 = WireCell::kslike_compare(vec_y, electron_ref);
        double ratio4 = std::accumulate(electron_ref.begin(), electron_ref.end(), 0.0) / 
                        (std::accumulate(vec_y.begin(), vec_y.end(), 0.0) + 1e-9);
        
        // std::cout << ks1 << " " << ratio1 << " " << ks2 << " " << ratio2 << " " << ks3 << " " << ratio3 << " " << ks4 << " " << ratio4 << std::endl;

        std::vector<double> results;
        // Convert bool result to double (1.0 for true, 0.0 for false)
        results.push_back(eval_ks_ratio(ks1, ks2, ratio1, ratio2) ? 1.0 : 0.0); // direction metric
        results.push_back(sqrt(pow(ks1, 2) + pow(ratio1-1, 2))); // muon information
        results.push_back(sqrt(pow(ks3, 2) + pow(ratio3-1, 2))); // proton information  
        results.push_back(sqrt(pow(ks4, 2) + pow(ratio4-1, 2))); // electron information
        
        return results;
    }
   
    double cal_kine_range(double L, int pdg_code, const Clus::ParticleDataSet::pointer& particle_data){

        IScalarFunction::pointer range_function = nullptr;

        if (abs(pdg_code) == 11) {        // electron
            range_function = particle_data->get_range_function("electron");
        }
        else if (abs(pdg_code) == 13) {   // muon
            range_function = particle_data->get_range_function("muon");
        }
        else if (abs(pdg_code) == 211) {  // pion
            range_function = particle_data->get_range_function("pion");
        }
        else if (abs(pdg_code) == 321) {  // kaon
            range_function = particle_data->get_range_function("kaon");
        }
        else if (abs(pdg_code) == 2212) { // proton
            range_function = particle_data->get_range_function("proton");
        }
        
        if (!range_function) {
            range_function = particle_data->get_range_function("muon");
        }
        if (!range_function) {
            return 0.0;
        }
        double kine_energy = range_function->scalar_function(L/units::cm) * units::MeV;
        return kine_energy;
    }

    // success, flag_dir, particle_type, particle_score
    std::tuple<bool, int, int, double> segment_do_track_pid(SegmentPtr segment, std::vector<double>& L , std::vector<double>& dQ_dx, const Clus::ParticleDataSet::pointer& particle_data, double compare_range , double offset_length, bool flag_force, const TrackPidOptions& pid_opts){
        const double MIP_dQdx = pid_opts.mip_dqdx;
        
        if (L.size() != dQ_dx.size() || L.empty() || !segment) {
            return std::make_tuple(false, 0, 0, 0.0);
        }
        
        std::vector<double> rL(L.size(), 0);
        std::vector<double> rdQ_dx(L.size(), 0);
        
        // Get reverse vectors
        for (size_t i = 0; i != L.size(); i++) {
            rL.at(i) = L.back() - L.at(L.size() - 1 - i);
            rdQ_dx.at(i) = dQ_dx.at(L.size() - 1 - i);
        }
        
        std::vector<double> result_forward = do_track_comp(L, dQ_dx, compare_range, offset_length, particle_data, MIP_dQdx, 0, pid_opts.track_comp_empty_abstain);
        std::vector<double> result_backward = do_track_comp(rL, rdQ_dx, compare_range, offset_length, particle_data, MIP_dQdx, 0, pid_opts.track_comp_empty_abstain);
        
        // Direction determination
        bool flag_forward = static_cast<bool>(std::round(result_forward.at(0)));
        bool flag_backward = static_cast<bool>(std::round(result_backward.at(0)));
        
        // Calculate length from path (total walk length over fits or wcpts)
        double length = segment_track_length(segment, 0);


        // // Calculate straight-line distance between endpoints (length1 equivalent)
        // double length1 = 0.0;
        // auto& fits = segment->fits();
        // length1 = (fits.front().point - fits.back().point).magnitude();
        
        // Forward particle type determination
        int forward_particle_type = 13; // default muon
        double min_forward_val = result_forward.at(1);
        if (result_forward.at(2) < min_forward_val) {
            min_forward_val = result_forward.at(2);
            forward_particle_type = 2212; // proton
        }
        if (result_forward.at(3) < min_forward_val && length < 20*units::cm) {
            min_forward_val = result_forward.at(3);
            forward_particle_type = 11; // electron
        }
        
        // Backward particle type determination  
        int backward_particle_type = 13; // default muon
        double min_backward_val = result_backward.at(1);
        if (result_backward.at(2) < min_backward_val) {
            min_backward_val = result_backward.at(2);
            backward_particle_type = 2212; // proton
        }
        if (result_backward.at(3) < min_backward_val && length < 20*units::cm) {
            min_backward_val = result_backward.at(3);
            backward_particle_type = 11; // electron
        }
        
        // Decision logic
        int flag_dir = 0;
        int particle_type = 0;
        double particle_score = 100.0;
        
        if (flag_forward == 1 && flag_backward == 0) {
            flag_dir = 1;
            particle_type = forward_particle_type;
            particle_score = min_forward_val;
            return std::make_tuple(true, flag_dir, particle_type, particle_score);
        }
        else if (flag_forward == 0 && flag_backward == 1) {
            flag_dir = -1;
            particle_type = backward_particle_type;
            particle_score = min_backward_val;
            return std::make_tuple(true, flag_dir, particle_type, particle_score);
        }
        else if (flag_forward == 1 && flag_backward == 1) {
            if (min_forward_val < min_backward_val) {
                flag_dir = 1;
                particle_type = forward_particle_type;
                particle_score = min_forward_val;
            }
            else {
                flag_dir = -1;
                particle_type = backward_particle_type;
                particle_score = min_backward_val;
            }
            return std::make_tuple(true, flag_dir, particle_type, particle_score);
        }
        else if (flag_forward == 0 && flag_backward == 0 && flag_force) {
            if (min_forward_val < min_backward_val) {
                particle_score = min_forward_val;
                particle_type = forward_particle_type;
                flag_dir = 1;
            }
            else {
                particle_score = min_backward_val;
                particle_type = backward_particle_type;
                flag_dir = -1;
            }
            return std::make_tuple(true, flag_dir, particle_type, particle_score);
        }
        
        // Improvement beyond the prototype (doc sbnd_xin/docs/pr/9 sec. 6 F1):
        // endpoint-trim retry.  Entered ONLY where the legacy logic abstains in
        // both orientations (no flag_force), so every decision the prototype
        // path makes is untouched.  The stopping tip's dQ/dx is unreliable when
        // the endpoint is ill-defined (deposit ends mid-bin => diluted, or
        // out-of-track charge piles into the last bin => inflated), yet it is
        // compared against the template's Bragg maximum — one bad tip sample
        // can veto an otherwise clean decision (SBND evt 172230: monotone
        // Bragg rise 98k->193k e/cm, tip collapses to 121k, both orientations
        // abstain; with the tip excluded the forward proton passes at score
        // 0.077).  Dynamic by construction: trim happens only on abstention,
        // exactly 1 sample, per-orientation at that orientation's hypothesized
        // stop end, template anchor unchanged, value-agnostic.
        if (pid_opts.endpoint_trim_retry && !flag_force) {
            std::vector<double> retry_forward  = do_track_comp(L,  dQ_dx,  compare_range, offset_length, particle_data, MIP_dQdx, 1, pid_opts.track_comp_empty_abstain);
            std::vector<double> retry_backward = do_track_comp(rL, rdQ_dx, compare_range, offset_length, particle_data, MIP_dQdx, 1, pid_opts.track_comp_empty_abstain);

            const bool rf = static_cast<bool>(std::round(retry_forward.at(0)));
            const bool rb = static_cast<bool>(std::round(retry_backward.at(0)));

            if (rf || rb) {
                // Same mu/p/e competition as the primary attempt (incl. the
                // <20 cm electron rule), evaluated on the trimmed results.
                int r_fwd_type = 13;
                double r_fwd_val = retry_forward.at(1);
                if (retry_forward.at(2) < r_fwd_val) { r_fwd_val = retry_forward.at(2); r_fwd_type = 2212; }
                if (retry_forward.at(3) < r_fwd_val && length < 20*units::cm) { r_fwd_val = retry_forward.at(3); r_fwd_type = 11; }

                int r_bwd_type = 13;
                double r_bwd_val = retry_backward.at(1);
                if (retry_backward.at(2) < r_bwd_val) { r_bwd_val = retry_backward.at(2); r_bwd_type = 2212; }
                if (retry_backward.at(3) < r_bwd_val && length < 20*units::cm) { r_bwd_val = retry_backward.at(3); r_bwd_type = 11; }

                if (rf && !rb) {
                    return std::make_tuple(true, 1, r_fwd_type, r_fwd_val);
                }
                if (rb && !rf) {
                    return std::make_tuple(true, -1, r_bwd_type, r_bwd_val);
                }
                // Both orientations pass on the trimmed data: pick the better
                // score, mirroring the primary attempt's both-pass branch.
                if (r_fwd_val < r_bwd_val) {
                    return std::make_tuple(true, 1, r_fwd_type, r_fwd_val);
                }
                return std::make_tuple(true, -1, r_bwd_type, r_bwd_val);
            }
            // Retry also abstained: fall through to the proton_dir_vote (which
            // votes on the UNTRIMMED results, unchanged pr/8 semantics).
        }

        // Improvement beyond the prototype (doc sbnd_xin/docs/pr/8): proton-template
        // direction vote.  Entered ONLY where the legacy logic abstains (both
        // orientations failed the muon-vs-flat gate, no flag_force), so every
        // decision the prototype path makes is untouched.  A guarded extension of
        // the flag_force mechanism above: proton-only, because the Bragg rise is
        // the one real directional dQ/dx signature (an abstaining flat muon's
        // fwd/bwd scores differ only by noise).  Guards:
        //   G1 good absolute proton fit — the score mixes shape and magnitude
        //      (sqrt(ks^2+(ratio-1)^2)), so MIP-like or shower profiles score >~1;
        //   G2 proton must win that orientation's existing mu/p/e competition;
        //   G3 significant fwd/bwd score asymmetry (a real direction claim);
        //   G4 the stop end must be topologically free — physical sanity, and it
        //      guarantees the result survives the caller's persistence condition.
        if (pid_opts.proton_dir_vote) {
            const double sp_f = result_forward.at(2);
            const double sp_b = result_backward.at(2);
            const bool fwd_wins = sp_f <= sp_b;
            const double sp_best  = fwd_wins ? sp_f : sp_b;
            const double sp_other = fwd_wins ? sp_b : sp_f;
            const int    win_type = fwd_wins ? forward_particle_type : backward_particle_type;
            const bool   free_end = fwd_wins ? (pid_opts.end_n == 1) : (pid_opts.start_n == 1);
            if (win_type == 2212
                && sp_best < pid_opts.proton_dir_score_max
                && sp_other > pid_opts.proton_dir_asym_min * sp_best
                && free_end) {
                return std::make_tuple(true, fwd_wins ? 1 : -1, 2212, sp_best);
            }
        }

        // Failure path: neither forward nor backward PID succeeded (and flag_force==false).
        // Do NOT write to segment->dirsign() here — this function has no segment side-effects;
        // callers (segment_determine_dir_track) set dirsign from the returned flag_dir.
        return std::make_tuple(false, 0, 0, 100.0);
    }

    // 4-momentum: E, px, py, pz,
    WireCell::D4Vector<double> segment_cal_4mom(SegmentPtr segment, int pdg_code, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx){
        double length = segment_track_length(segment, 0);
        double kine_energy = 0;

        // std::cout << "segment_cal_4mom: pdg=" << pdg_code
        //           << " length=" << length/units::cm << " cm"
        //           << " nfits=" << segment->fits().size() << std::endl;

        WireCell::D4Vector<double> results(0.0, 0.0, 0.0, 0.0); // 4-momentum: E, px, py, pz

        if (length < 4*units::cm){
            kine_energy = segment_cal_kine_dQdx(segment, recomb_model); // short track
        }else if (segment->flags_any(PR::SegmentFlags::kShowerTrajectory)){
            kine_energy = segment_cal_kine_dQdx(segment, recomb_model);
        }else{
            kine_energy = cal_kine_range(length, pdg_code, particle_data);
        }
        // std::cout << "segment_cal_4mom: kine_energy=" << kine_energy/units::MeV << " MeV" << std::endl;
        // results[4] = kine_energy;

        double particle_mass = particle_data->get_particle_mass(pdg_code);

        results[0]= kine_energy + particle_mass;
        double mom = sqrt(pow(results[0],2) - pow(particle_mass,2));
        auto v1 = segment_cal_dir_3vector(segment);

        results[1] = mom * v1.x();
        results[2] = mom * v1.y();
        results[3] = mom * v1.z();

        return results;
    }

    void segment_determine_dir_track(SegmentPtr segment, int start_n, int end_n, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx, bool flag_print, const TrackPidOptions& pid_opts) {
        if (!segment || !particle_data) {
            return;
        }

        // Per-segment copy: the vote's free-end guard (G4) needs the vertex
        // connectivity of THIS segment.
        TrackPidOptions opts = pid_opts;
        opts.start_n = start_n;
        opts.end_n = end_n;
        
        // Reset direction flag
        segment->dirsign(0);
        
        const auto& fits = segment->fits();
        int npoints = fits.size();
        int start_n1 = 0, end_n1 = npoints - 1;
        
        // If more than one point, exclude the vertex
        if (end_n != 1) {
            end_n1 = npoints - 2;
            npoints -= 1;
        }
        if (start_n != 1) {
            npoints -= 1;
            start_n1 = 1;
        }
        
        if (npoints == 0 || end_n1 < start_n1) return;
        
        std::vector<double> L(npoints, 0);
        std::vector<double> dQ_dx(npoints, 0);
        
        double dis = 0;
        // Unit convention (docs/pr/2 sec 7.1c): fits[i].dx is INTERNAL length
        // (mm = 1, cm = 10 -- WCP and WCT share this system), so dQ_dx[i] =
        // dQ / dx is [charge / internal unit].  do_track_comp's reference
        // templates are converted to the same scale (scalar_function(cm) /
        // units::cm, MIP_dQdx = m_mip_dqdx), so data and templates are
        // self-consistent; the KS and sum-ratio comparisons only need the two
        // sides to share a scale.  Do NOT divide dx by units::cm here: that
        // would make the data e/cm, 10x the templates -- exactly the SSM
        // get_scores bug fixed in the sec 7.8 round.  (The prototype's
        // do_track_pid path was also self-consistent, just in e/cm on both
        // sides: ProtoSegment.cxx data dQ_vec/(dx_vec/units::cm + 1e-9)
        // against e/cm TGraph tables and the 50e3 constant.)
        for (int i = start_n1; i <= end_n1; i++) {
            L.at(i - start_n1) = dis;
            if (fits[i].dx > 0) {
                dQ_dx.at(i - start_n1) = fits[i].dQ / (fits[i].dx + 1e-9);
            }
            if (i + 1 < static_cast<int>(fits.size())) {
                dis += (fits[i+1].point - fits[i].point).magnitude();
            }
        }
        
        int pdg_code = 0;
        double particle_score = 100.0;
        
        if (npoints >= 2) { // reasonably long
            bool tmp_flag_pid = false;
            
            if (start_n == 1 && end_n == 1 && npoints >= 15) {
                // Can use the dQ/dx to do PID and direction
                auto result = segment_do_track_pid(segment, L, dQ_dx, particle_data, 35*units::cm, 1*units::cm, true, opts);
                tmp_flag_pid = std::get<0>(result);
                if (tmp_flag_pid) {
                    segment->dirsign(std::get<1>(result));
                    pdg_code = std::get<2>(result);
                    particle_score = std::get<3>(result);
                }
                
                if (!tmp_flag_pid) {
                    result = segment_do_track_pid(segment, L, dQ_dx, particle_data, 15*units::cm, 1*units::cm, true, opts);
                    tmp_flag_pid = std::get<0>(result);
                    if (tmp_flag_pid) {
                        segment->dirsign(std::get<1>(result));
                        pdg_code = std::get<2>(result);
                        particle_score = std::get<3>(result);
                    }
                }
            } else {
                // Can use the dQ/dx to do PID and direction
                auto result = segment_do_track_pid(segment, L, dQ_dx, particle_data, 35*units::cm, 0*units::cm, false, opts);
                tmp_flag_pid = std::get<0>(result);
                if (tmp_flag_pid) {
                    segment->dirsign(std::get<1>(result));
                    pdg_code = std::get<2>(result);
                    particle_score = std::get<3>(result);
                }
                
                if (!tmp_flag_pid) {
                    result = segment_do_track_pid(segment, L, dQ_dx, particle_data, 15*units::cm, 0*units::cm, false, opts);
                    tmp_flag_pid = std::get<0>(result);
                    if (tmp_flag_pid) {
                        segment->dirsign(std::get<1>(result));
                        pdg_code = std::get<2>(result);
                        particle_score = std::get<3>(result);
                    }
                }
            }
        }
        
        double length = segment_track_length(segment, 0);

        // Compute median dQ/dx once over the trimmed range — reused in three branches below.
        // doc sbnd_xin/docs/pr/31 §12 (F4, was P8).  The prototype takes this
        // median over the SAME dQ_dx vector it hands to do_track_pid, at all
        // three of its sites (ProtoSegment.cxx:1574-1576 pdg==0 branch,
        // :1602-1611 both vertex-activity branches: nth_element over a copy of
        // dQ_dx, element size()/2).  segment_median_dQ_dx instead rebuilds a
        // FILTERED sample from fits() (drops !valid()/dx<=0/dQ<0), so one bad
        // fit is "a zero" to the PID and "not there" to the median, and the
        // filter changes the vector length so nth_element selects a different
        // order statistic.  dir_track_median_local restores the prototype's
        // self-consistency; the filtered helper itself is kept (its exclusion
        // of the dQ/1e-9 sentinel is judged an improvement) for every other
        // caller.  Same internal-unit scale on both paths.  false =>
        // filtered helper => byte-identical.
        double medium_dQ_dx = 0;
        if (opts.dir_track_median_local) {
            std::vector<double> tmp = dQ_dx;
            std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
            medium_dQ_dx = *std::next(tmp.begin(), tmp.size()/2);
        }
        else {
            medium_dQ_dx = segment_median_dQ_dx(segment, start_n1, end_n1);
        }

        // Short track what to do???
        if (pdg_code == 0) {
            if (medium_dQ_dx > MIP_dQdx * 1.75) {
                pdg_code = 2212; // proton
            } else if (medium_dQ_dx < MIP_dQdx * 1.2) {
                pdg_code = 13; // muon
            } else if (medium_dQ_dx < MIP_dQdx * 1.5 && length < 4*units::cm) {
                pdg_code = 13;
            }
        }

        // Electron and both end contain stuff
        if (abs(pdg_code) == 11 && (start_n > 1 && end_n > 1)) {
            segment->dir_weak(true);
            segment->dirsign(0);
            if (particle_score < 0.15) pdg_code = 13;
        } else if (abs(pdg_code) == 11 && ((start_n > 1 && segment->dirsign() == -1) || (end_n > 1 && segment->dirsign() == 1))) {
            segment->dir_weak(true);
            segment->dirsign(0);
            if (particle_score < 0.15) pdg_code = 13;
        } else if (length < 1.5*units::cm) {
            segment->dir_weak(true);
        }

        // Vertex activities
        if (length < 1.5*units::cm && (start_n == 1 || end_n == 1)) {
            if (start_n == 1 && end_n > 2) {
                segment->dirsign(-1);
                if (medium_dQ_dx > MIP_dQdx * 1.75) {
                    pdg_code = 2212;
                } else if (medium_dQ_dx < MIP_dQdx * 1.2) {
                    pdg_code = 211;
                }
            } else if (end_n == 1 && start_n > 2) {
                segment->dirsign(1);
                if (medium_dQ_dx > MIP_dQdx * 1.75) {
                    pdg_code = 2212;
                } else if (medium_dQ_dx < MIP_dQdx * 1.2) {
                    pdg_code = 211;
                }
            }
        }
        
        // If the particle score is really bad, make it a shower
        if (length > 10*units::cm && particle_score > 1.0 && particle_score < 100) {
            pdg_code = 11;
            particle_score = 200;
            segment->dirsign(0);
        }
        
        // Set particle mass and calculate 4-momentum
        // Only calculate if direction points toward a free end (matching WCPPID logic)
        const bool free_end_dir = (segment->dirsign() == 1 && end_n == 1) || (segment->dirsign() == -1 && start_n == 1);
        // doc sbnd_xin/docs/pr/40 F1 (= doc pr/7 sec 5 / pr/31 P14-F8):
        // track_pid_persist_dqdx persists type+mass whenever pdg_code != 0,
        // matching the prototype (ProtoSegment.cxx:1637-1639), and gates only
        // the 4-momentum/energy recompute on free_end_dir -- the existing
        // toolkit test, standing in for the prototype's
        // get_particle_4mom(3)>0 ("an energy was already computed") guard;
        // not independently re-verified against prototype source this round
        // (prototype_base symlink unavailable), see doc pr/40 for the
        // measured consequence either way.  false = legacy = byte-identical.
        //
        // doc sbnd_xin/docs/pr/40 round 5 F9: track_pid_persist_dqdx_electron_
        // guard narrows the rescue above -- it does not widen it.  An
        // UNDIRECTED (free_end_dir false) pdg_code==11 conclusion is never
        // F1's target population (proton/muon rescues from the median-dQ/dx
        // fallback or a confident-but-not-free-ended template match); letting
        // it persist anyway lets a single weak electron guess poison a
        // NEIGHBORING segment's flag_shower_in test (NeutrinoVertexFinder.cxx
        // :1320) into demoting an independently-confident muon/proton call.
        // false = today's shipped F1 behaviour = byte-identical.
        const bool f1_rescue = pid_opts.track_pid_persist_dqdx &&
            !(pid_opts.track_pid_persist_dqdx_electron_guard && pdg_code == 11 && !free_end_dir);
        if (pdg_code != 0 && (free_end_dir || f1_rescue)) {
            // Calculate 4-momentum using the identified particle type.  When
            // the knob rescues a non-free-end store, the momentum direction
            // is only as good as segment_cal_dir_3vector's no-argument
            // overload gives for a dirsign that may be 0 -- acceptable here
            // because downstream consumers gate on pdg/flag_shower, not on
            // this vector's orientation, for exactly the segments this knob
            // targets (see doc pr/40 Part 0).
            //
            // doc sbnd_xin/docs/pr/40 round 2 F4: the non-free-end branch below used
            // to store a rest-mass-only 4-vector (E = m, zero momentum) as a
            // literal stand-in for "the prototype's zero-momentum stub".
            // That makes Aux::ParticleInfo's kinetic_energy() (= E - mass)
            // exactly ZERO for every segment F1 rescues this way -- measured
            // as a real defect (SBND evt 174637 seg 9050, Bee PF node
            // "mu- 0 MeV") once track_pid_persist_dqdx went SBND-default-ON.
            // segment_cal_4mom itself has NO free-end dependence (its only
            // direction coupling is segment_cal_dir_3vector, which already
            // degrades gracefully to a zero 3-vector when dirsign()==0) --
            // the free-end gate here was purely external and stricter than
            // it needed to be.  track_pid_persist_4mom calls it
            // unconditionally instead, giving a correct KE with (at worst) a
            // zero-direction momentum.  Independent of
            // track_pid_persist_dqdx: false = legacy stub, byte-identical.
            WireCell::D4Vector<double> four_momentum(0.0, 0.0, 0.0, 0.0);
            if (free_end_dir || pid_opts.track_pid_persist_4mom) {
                four_momentum = segment_cal_4mom(segment, pdg_code, particle_data, recomb_model, MIP_dQdx);
            } else {
                const double mass = particle_data->get_particle_mass(pdg_code);
                four_momentum[0] = mass;  // at rest: matches the prototype's zero-momentum stub
            }

            // Create ParticleInfo with the identified particle
            auto pinfo = std::make_shared<Aux::ParticleInfo>(
                pdg_code,                                    // PDG code
                particle_data->get_particle_mass(pdg_code), // mass
                particle_data->pdg_to_name(pdg_code),       // name
                four_momentum                                // 4-momentum (E, px, py, pz)
            );

            // Store particle info in segment
            segment->particle_info(pinfo);
            segment->particle_score(particle_score);
        }
                
        if (flag_print) {
            // Match WCPPID output format: id, length, "Track", flag_dir, is_dir_weak, particle_type, mass, KE, particle_score
            double particle_mass = pdg_code != 0 ? particle_data->get_particle_mass(pdg_code) : 0.0;
            double kinetic_energy = 0.0;

            if (segment->has_particle_info()) {
                kinetic_energy = segment->particle_info()->kinetic_energy();
            }

            SPDLOG_LOGGER_TRACE(s_log, "segment_determine_dir_track: Seg {} cm Track {} {} {} {} MeV {} MeV {}",
                length/units::cm, segment->dirsign(), (segment->dir_weak() ? 1 : 0),
                pdg_code, particle_mass/units::MeV, kinetic_energy/units::MeV, particle_score);
            // doc sbnd_xin/docs/pr/40: the line above has no segment id, so it
            // cannot be joined to a Bee/calib segment.  Add the encoded id
            // (cluster_id*1000+graph_index, same scheme as PrDisplayDump.cxx)
            // on a separate WCT_PID_TRACE_DEBUG-gated tag rather than change
            // the format above (its comment claims WCPPID-format parity).
            SPDLOG_LOGGER_TRACE(s_log, "PID_TRACE_DEBUG id={} clus={} gidx={} L={:.1f}cm dir={} weak={} pdg={} score={:.4f}",
                (segment->cluster() ? segment->cluster()->get_cluster_id() : -1) * 1000 + static_cast<int>(segment->get_graph_index()),
                segment->cluster() ? segment->cluster()->get_cluster_id() : -1, segment->get_graph_index(),
                length/units::cm, segment->dirsign(), (segment->dir_weak() ? 1 : 0), pdg_code, particle_score);
        }
    }

     void segment_determine_shower_direction_trajectory(SegmentPtr segment, int start_n, int end_n, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx, bool flag_print, const TrackPidOptions& pid_opts){
        segment->dirsign(0);
        double length = segment_track_length(segment, 0);

        // For shower-trajectory segments PDG is always forced to electron (11).
        // particle_score = 100.0 is a sentinel meaning "PID not performed; score not applicable."
        // Callers must not interpret this score as a PID quality metric.
        // The prototype (ProtoSegment.cxx:1647) also did not compute a meaningful particle_score
        // for this code path; the sentinel is intentional, not a placeholder.
        int pdg_code = 11; // electron
        double particle_score = 100.0;
        
        if (start_n == 1 && end_n > 1){
            segment->dirsign(-1);
        } else if (start_n > 1 && end_n == 1){
            segment->dirsign(1);
        } else {
            // Try track PID first
            segment_determine_dir_track(segment, start_n, end_n, particle_data, recomb_model, MIP_dQdx, false, pid_opts);
            
            // Check if particle info was set and if it's not an electron
            if (segment->has_particle_info() && segment->particle_info()->pdg() != 11) {
                // Reset to electron and no direction
                pdg_code = 11;
                segment->dirsign(0);
            } else if (segment->has_particle_info()) {
                // Keep the electron identification from track PID
                pdg_code = segment->particle_info()->pdg();
            } else {
                // No particle info set, default to electron with no direction
                pdg_code = 11;
                segment->dirsign(0);
            }
        }
        
        // Always calculate 4-momentum for shower trajectories (matching WCPPID)
        auto four_momentum = segment_cal_4mom(segment, pdg_code, particle_data, recomb_model, MIP_dQdx);

        // Create ParticleInfo with the identified particle
        auto pinfo = std::make_shared<Aux::ParticleInfo>(
            pdg_code,                                    // PDG code (electron)
            particle_data->get_particle_mass(pdg_code), // mass
            particle_data->pdg_to_name(pdg_code),       // name
            four_momentum                                // 4-momentum (E, px, py, pz)
        );
                
        // Store particle info in segment
        segment->particle_info(pinfo);
        segment->particle_score(particle_score);

        if (flag_print) {
            // Match WCPPID output format: id, length, "S_traj", flag_dir, is_dir_weak, particle_type, mass, KE, particle_score
            double particle_mass = particle_data->get_particle_mass(pdg_code);
            double kinetic_energy = pinfo->kinetic_energy();
            
            SPDLOG_LOGGER_TRACE(s_log, "segment_determine_shower_direction_trajectory: Seg {} cm S_traj {} {} {} {} MeV {} MeV {}",
                length/units::cm, segment->dirsign(), (segment->dir_weak() ? 1 : 0),
                pdg_code, particle_mass/units::MeV, kinetic_energy/units::MeV, particle_score);
        }
     }

    void clustering_points_segments(std::vector<SegmentPtr> segments, const IDetectorVolumes::pointer& dv, const std::string& cloud_name, double search_range, double scaling_2d){
        // using Clock = std::chrono::steady_clock;
        // using MS = std::chrono::duration<double, std::milli>;
        // auto t_total = Clock::now();

        // Use cluster_id-based comparator so map_cluster_segs is iterated in a
        // deterministic, run-to-run stable order (not pointer-address order).
        struct ClusterIdCmp {
            bool operator()(Facade::Cluster* a, Facade::Cluster* b) const {
                return a->get_cluster_id() < b->get_cluster_id();
            }
        };
        std::map<Facade::Cluster*, std::vector<SegmentPtr>, ClusterIdCmp> map_cluster_segs;
        for (auto seg : segments){
            if (seg->cluster()){
                map_cluster_segs[seg->cluster()].push_back(seg);
            }
        }

        // Pre-build fit-point caches for ALL input segments upfront.
        // This enables cross-cluster ghost removal in the loop below, matching the
        // prototype's behaviour where all segments (across all clusters) are compared
        // during ghost-removal — not just those belonging to the cluster being processed.
        using ApFaceKey = std::tuple<int, int, int>;  // (plane, apa, face)
        struct Pt2D { double x, y; };
        // Use SegmentIndexCmp so iteration order is deterministic (graph index), not pointer-address.
        std::map<SegmentPtr, std::shared_ptr<Facade::DynamicPointCloud>, SegmentIndexCmp> seg_dpc_cache;
        std::map<SegmentPtr, std::vector<std::array<double, 3>>,         SegmentIndexCmp> seg_pts3d;
        std::map<SegmentPtr, std::map<ApFaceKey, std::vector<Pt2D>>,     SegmentIndexCmp> seg_pts2d;
        for (auto seg : segments) {
            auto dpc = seg->dpcloud("fit");
            if (!dpc) continue;
            seg_dpc_cache[seg] = dpc;
            auto& pts3d     = seg_pts3d[seg];
            auto& pts2d_map = seg_pts2d[seg];
            const auto& cols = dpc->points();
            for (size_t ip = 0; ip < cols.size(); ++ip) {
                pts3d.push_back({cols.x[ip], cols.y[ip], cols.z[ip]});
                for (int pind = 0; pind < 3; ++pind) {
                    const auto [jb, je] = cols.proj_range(ip, pind);
                    for (size_t j = jb; j < je; ++j) {
                        WirePlaneId wpid_2d(cols.p2d_wpid[j]);
                        pts2d_map[{pind, wpid_2d.apa(), wpid_2d.face()}]
                            .push_back({cols.p2d_x[j], cols.p2d_y[j]});
                    }
                }
            }
        }

        // F17: Build global per-(plane,apa,face) 2D KD-trees from all segment fit-point projections.
        // Querying these once per cluster point replaces the O(S) inner loop over segments,
        // reducing ghost-removal from O(P×S×F) to O(P×log(total_fits)) overall.
        using nfkd2d_t = NFKDVec::Tree<double, NFKDVec::IndexDynamic>;
        std::map<ApFaceKey, nfkd2d_t> global_kd2d;
        for (const auto& [seg, pts2d_map] : seg_pts2d) {
            for (const auto& [apfkey, pts] : pts2d_map) {
                auto [it, inserted] = global_kd2d.try_emplace(apfkey, 2);
                (void)inserted;
                std::vector<double> xs, ys;
                xs.reserve(pts.size());
                ys.reserve(pts.size());
                for (const auto& p : pts) { xs.push_back(p.x); ys.push_back(p.y); }
                it->second.append(nfkd2d_t::points_type{xs, ys});
            }
        }

        // Pre-populate the angle cache for every (apa, face) pair that actually appears in
        // seg_pts2d.  The old lazy approach iterated seg_dpc_cache at query time and could
        // silently store (0,0,0) when the first segment in cache order lacked that face —
        // giving wrong 2D projections for all subsequent points on that face.
        std::map<std::pair<int,int>, std::array<double, 3>> ang_cache;
        for (const auto& [seg, pts2d_map] : seg_pts2d) {
            auto dpc_it = seg_dpc_cache.find(seg);
            if (dpc_it == seg_dpc_cache.end()) continue;
            const auto& dpc = dpc_it->second;
            for (const auto& [apfkey, pts] : pts2d_map) {
                const auto& [pind, apa, face] = apfkey;
                auto key = std::make_pair(apa, face);
                if (ang_cache.count(key)) continue;  // already populated for this (apa,face)
                auto a = dpc->get_angles(face, apa);
                if (a[0] != 0.0 || a[1] != 0.0 || a[2] != 0.0) {
                    ang_cache[key] = a;
                }
            }
        }
        // O(1) lookup; unknown (apa,face) inserts (0,0,0) by default.
        auto get_angles_cached = [&](int apa, int face) -> const std::array<double, 3>& {
            auto [it, inserted] = ang_cache.emplace(
                std::make_pair(apa, face), std::array<double, 3>{0.0, 0.0, 0.0});
            return it->second;
        };

        // Fast 2D closest squared-distance: linear scan over precomputed fit-point projections.
        // Returns squared distances {d_u^2, d_v^2, d_w^2}.
        // Sentinel 1e18 for planes with no matching (apa, face) data (e.g. cross-APA).
        // Returns squared distances (not sqrt) so that the equality comparison with the
        // global KD-tree result (also squared) is exact — no sqrt round-trip ULP mismatch.
        auto get_2d_dist2_fast = [&](SegmentPtr seg,
                                     double qx, const std::array<double, 3>& qy,
                                     int apa, int face) -> std::tuple<double, double, double> {
            auto it = seg_pts2d.find(seg);
            if (it == seg_pts2d.end()) return {1e18, 1e18, 1e18};
            double min_d2[3] = {1e18, 1e18, 1e18};
            for (int pind = 0; pind < 3; ++pind) {
                auto it2 = it->second.find({pind, apa, face});
                if (it2 == it->second.end()) continue;
                for (const auto& [px, py] : it2->second) {
                    double dx = qx - px, dy = qy[pind] - py;
                    double d2 = dx*dx + dy*dy;
                    if (d2 < min_d2[pind]) min_d2[pind] = d2;
                }
            }
            return {min_d2[0], min_d2[1], min_d2[2]};  // squared distances
        };

        // Fast 3D closest distance: linear scan over precomputed fit-point coordinates.
        auto get_3d_dist_fast = [&](SegmentPtr seg, const geo_point_t& gp) -> double {
            auto it = seg_pts3d.find(seg);
            if (it == seg_pts3d.end()) return 1e9;
            double min_d2 = 1e18;
            for (const auto& [px, py, pz] : it->second) {
                double dx = gp.x()-px, dy = gp.y()-py, dz = gp.z()-pz;
                double d2 = dx*dx + dy*dy + dz*dz;
                if (d2 < min_d2) min_d2 = d2;
            }
            return std::sqrt(min_d2);
        };

        for (auto it : map_cluster_segs){
            auto clus = it.first;
            auto& segs = it.second;

            // get the default point cloud from cluster
            const auto& points = clus->points();

             // Get the graph directly
            const auto& graph = clus->find_graph("basic_pid");

            std::map<SegmentPtr, std::vector<geo_point_t>> map_segment_points;
            std::map<SegmentPtr, std::vector<size_t>> map_segment_global_indices;
            std::map<int, std::pair<SegmentPtr, double>> map_pindex_segment;

            // core algorithms

            // define steiner terminal for segments ...
            for (auto seg: segs){
                auto& fits = seg->fits();
                if (fits.size() > 2){
                    for (size_t i = 1; i+1 < fits.size(); i++){
                        geo_point_t gp = {fits[i].point.x(), fits[i].point.y(), fits[i].point.z()};
                        // use cluster to get the indices of the closest 5 points
                        auto closest_results = clus->kd_knn(5, gp);
                        for (const auto& [point_index, distance] : closest_results) {
                            if (map_pindex_segment.find(point_index) == map_pindex_segment.end()) {
                                map_pindex_segment[point_index] = std::make_pair(seg, distance);
                                break;
                            }
                        }
                    }
                }else if (fits.size() == 2){
                    geo_point_t gp = {(fits[0].point.x()+fits[1].point.x())/2., (fits[0].point.y()+fits[1].point.y())/2., (fits[0].point.z()+fits[1].point.z())/2.};
                    auto closest_results = clus->kd_knn(5, gp);
                    for (const auto& [point_index, distance] : closest_results) {
                        if (map_pindex_segment.find(point_index) == map_pindex_segment.end()) {
                            map_pindex_segment[point_index] = std::make_pair(seg, distance);
                            break;
                        }
                    }
                }
            }


            // these are terminals ...
            // auto t_ghost = Clock::now();
            if (map_pindex_segment.size()>0){

                // Convert terminals from int to vertex_type
                std::vector<WireCell::Clus::Graphs::Weighted::vertex_type> terminals;
                for (auto it = map_pindex_segment.begin(); it!=map_pindex_segment.end(); it++){
                    terminals.push_back(static_cast<WireCell::Clus::Graphs::Weighted::vertex_type>(it->first));
                }

                // auto t_voronoi = Clock::now();
                auto vor = WireCell::Clus::Graphs::Weighted::voronoi(graph, terminals);
                // std::cout << "[clustering_points_segments] voronoi took " << MS(Clock::now() - t_voronoi).count() << " ms" << std::endl;

                const int num_graph_vertices = boost::num_vertices(graph);
                // Use the flat vector directly for O(1) access by vertex index
                const auto& nearest_terminal = vor.terminal;

                // Squared 2D threshold — constant across all points in this cluster; hoisted out of loop.
                const double sq_2d_thr = (scaling_2d * search_range) * (scaling_2d * search_range);

                // now examine to remove ghost points ....
                for (int i = 0; i < num_graph_vertices; i++){
                    const int nt_i = static_cast<int>(nearest_terminal[i]);
                    if (map_pindex_segment.find(nt_i) == map_pindex_segment.end()) continue;
                    geo_point_t gp(points[0][i], points[1][i], points[2][i]);
                    auto main_sg = map_pindex_segment[nt_i].first;

                    auto point_wpid = clus->wire_plane_id(i);
                    auto apa = point_wpid.apa();
                    auto face = point_wpid.face();

                    // Compute projected query coordinates once per point (shared by all segments)
                    const auto& ang = get_angles_cached(apa, face);
                    const double qx = gp.x();
                    const std::array<double, 3> qy = {
                        std::cos(ang[0]) * gp.z() - std::sin(ang[0]) * gp.y(),
                        std::cos(ang[1]) * gp.z() - std::sin(ang[1]) * gp.y(),
                        std::cos(ang[2]) * gp.z() - std::sin(ang[2]) * gp.y()
                    };

                    // 3D closest distance to main segment (linear scan, used only for < search_range check)
                    std::pair<double, WireCell::Point> closest_dis_point = {get_3d_dist_fast(main_sg, gp), WireCell::Point(0,0,0)};

                    // All 2D comparisons use squared distances to avoid sqrt round-trip ULP mismatch
                    // between the KD-tree global query and the linear-scan per-segment result.
                    // "closest_2d_dis2" = min squared 2D distance from this point to main_sg's fit projections.
                    // "min_2d_dis2"     = min squared 2D distance from this point to ANY segment's projections.
                    // Both are computed without sqrt, so equality is exact when the same underlying point wins.
                    std::tuple<double, double, double> closest_2d_dis2 = get_2d_dist2_fast(main_sg, qx, qy, apa, face);

                    // F17: Replace O(S) inner loop with O(log N) global KD-tree queries.
                    // The global trees include main_sg's own points, so global_min2 <= closest_2d_dis2
                    // per plane.  Equality (exact, since both use the same squared-distance formula)
                    // means main_sg achieves the global minimum on that plane.
                    std::tuple<double, double, double> min_2d_dis2 = closest_2d_dis2;
                    {
                        const std::array<double, 2> query_u = {qx, qy[0]};
                        const std::array<double, 2> query_v = {qx, qy[1]};
                        const std::array<double, 2> query_w = {qx, qy[2]};
                        auto q2 = [&](int pind, const std::array<double,2>& q) {
                            auto kit = global_kd2d.find({pind, apa, face});
                            if (kit == global_kd2d.end()) return;
                            auto res = kit->second.knn(1, q);
                            if (!res.empty()) {
                                double d2 = res[0].second;  // already squared — no sqrt
                                if (pind == 0) std::get<0>(min_2d_dis2) = d2;
                                else if (pind == 1) std::get<1>(min_2d_dis2) = d2;
                                else               std::get<2>(min_2d_dis2) = d2;
                            }
                        };
                        q2(0, query_u);
                        q2(1, query_v);
                        q2(2, query_w);
                    }

                    // check against main_sg — all comparisons in squared-distance space
                    bool flag_change = true;

                    if (std::get<0>(min_2d_dis2) == std::get<0>(closest_2d_dis2) && std::get<1>(min_2d_dis2) == std::get<1>(closest_2d_dis2) && std::get<2>(min_2d_dis2) == std::get<2>(closest_2d_dis2)) // all closest
                    flag_change = false;
                    else if (std::get<0>(min_2d_dis2) == std::get<0>(closest_2d_dis2) && std::get<1>(min_2d_dis2) == std::get<1>(closest_2d_dis2) ) // 2 closest
                    flag_change = false;
                    else if (std::get<0>(min_2d_dis2) == std::get<0>(closest_2d_dis2) && std::get<2>(min_2d_dis2) == std::get<2>(closest_2d_dis2) )
                    flag_change = false;
                    else if (std::get<1>(min_2d_dis2) == std::get<1>(closest_2d_dis2) && std::get<2>(min_2d_dis2) == std::get<2>(closest_2d_dis2) )
                    flag_change = false;
                    else if (std::get<0>(min_2d_dis2) == std::get<0>(closest_2d_dis2) && (closest_dis_point.first < search_range || (std::get<1>(closest_2d_dis2) < sq_2d_thr && std::get<2>(closest_2d_dis2) < sq_2d_thr)) )
                    flag_change = false;
                    else if (std::get<1>(min_2d_dis2) == std::get<1>(closest_2d_dis2) && (closest_dis_point.first < search_range || (std::get<0>(closest_2d_dis2) < sq_2d_thr && std::get<2>(closest_2d_dis2) < sq_2d_thr) ))
                    flag_change = false;
                    else if (std::get<2>(min_2d_dis2) == std::get<2>(closest_2d_dis2) && (closest_dis_point.first < search_range || (std::get<1>(closest_2d_dis2) < sq_2d_thr && std::get<0>(closest_2d_dis2) < sq_2d_thr) ))
                    flag_change = false;

                    // deal with dead channels ...
                    if (!flag_change){
                        auto grouping = clus->grouping();
                        int ch_range = 0; // Default channel range for dead channel checking

                        // Check U plane (pind=0) for dead channels
                        if (grouping->get_closest_dead_chs(gp, ch_range, apa, face, 0) && std::get<0>(closest_2d_dis2) > sq_2d_thr){
                            if (std::get<1>(closest_2d_dis2) < sq_2d_thr || std::get<2>(closest_2d_dis2) < sq_2d_thr)
                                flag_change = true;
                        // Check V plane (pind=1) for dead channels
                        }else if (grouping->get_closest_dead_chs(gp, ch_range, apa, face, 1) && std::get<1>(closest_2d_dis2) > sq_2d_thr){
                            if (std::get<0>(closest_2d_dis2) < sq_2d_thr || std::get<2>(closest_2d_dis2) < sq_2d_thr)
                                flag_change = true;
                        // Check W plane (pind=2) for dead channels
                        }else if (grouping->get_closest_dead_chs(gp, ch_range, apa, face, 2) && std::get<2>(closest_2d_dis2) > sq_2d_thr){
                            if (std::get<1>(closest_2d_dis2) < sq_2d_thr || std::get<0>(closest_2d_dis2) < sq_2d_thr)
                                flag_change = true;
                        }
                   }
                
                    // change the point's clustering ...
                    if (!flag_change){
                        map_segment_global_indices[main_sg].push_back(i);
                        map_segment_points[main_sg].push_back(gp);
                    }
                }
            }


            // std::cout << "[clustering_points_segments] ghost removal took " << MS(Clock::now() - t_ghost).count() << " ms" << std::endl;

            // convert points to geo_point_t format
            // add points to segments ...
            // auto t_build = Clock::now();
            for (const auto& [seg, geo_points] : map_segment_points) {
                const auto& global_indices = map_segment_global_indices[seg];
                create_segment_point_cloud(seg, geo_points, dv, cloud_name, global_indices);
                // create_segment_point_cloud(seg, geo_points, dv, cloud_name);
            }

            // std::cout << "[clustering_points_segments] build point clouds took " << MS(Clock::now() - t_build).count() << " ms" << std::endl;

            // debug: print number of points assigned to each segment
            // for (const auto& [seg, geo_points] : map_segment_points) {
            //     std::cout << "[clustering_points_segments] cloud: " << cloud_name
            //               << "  seg id: " << seg->id()
            //               << "  npoints: " << geo_points.size() << std::endl;
            // }
        }
        // std::cout << "[clustering_points_segments] TOTAL took " << MS(Clock::now() - t_total).count() << " ms" << std::endl;
    }

    bool segment_determine_shower_direction(SegmentPtr segment, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const std::string& cloud_name, double MIP_dQdx, double rms_cut, double mip_dqdx, bool median_local){
        segment->dirsign(0);
        const auto& fits = segment->fits();
        
        if (fits.empty()) return false;
        
        // Get the fit point cloud for KD-tree queries
        auto dpcloud_fit = segment->dpcloud("fit");
        if (!dpcloud_fit) return false;
        
        // Get the associated point cloud 
        auto dpcloud_assoc = segment->dpcloud(cloud_name);
        if (!dpcloud_assoc) return false;
        
        const size_t assoc_npts = dpcloud_assoc->npoints();
        if (assoc_npts == 0) return false;
        
        // Initialize vectors to store analysis results for each fit point
        std::vector<std::vector<WireCell::Point>> local_points_vec(fits.size());
        std::vector<std::tuple<double, double, double>> vec_rms_vals(fits.size(), std::make_tuple(0,0,0));
        std::vector<double> vec_dQ_dx(fits.size(), 0);
        std::vector<WireCell::Vector> vec_dir(fits.size());
        
        // Build KD-tree index for fit points and associate points with nearest fit point
        auto& kd_tree_fit = dpcloud_fit->kd3d();
        
        for (size_t ia = 0; ia < assoc_npts; ++ia) {
            WireCell::Point test_p = dpcloud_assoc->point3d(ia);
            auto results = kd_tree_fit.knn(1, test_p);
            if (!results.empty()) {
                size_t closest_fit_idx = results.front().first;
                local_points_vec.at(closest_fit_idx).push_back(test_p);
            }
        }
        
        WireCell::Vector drift_dir_abs(1, 0, 0);  // Drift direction
        
        // Calculate local directions and RMS spreads for each fit point
        for (size_t i = 0; i < local_points_vec.size(); i++) {
            // Calculate local direction from neighboring fit points
            WireCell::Vector v1(0, 0, 0);
            for (size_t j = 1; j < 3; j++) {
                if (i + j < fits.size()) {
                    v1 += WireCell::Vector(
                        fits[i+j].point.x() - fits[i].point.x(),
                        fits[i+j].point.y() - fits[i].point.y(),
                        fits[i+j].point.z() - fits[i].point.z()
                    );
                }
                if (i >= j) {
                    v1 += WireCell::Vector(
                        fits[i].point.x() - fits[i-j].point.x(),
                        fits[i].point.y() - fits[i-j].point.y(),
                        fits[i].point.z() - fits[i-j].point.z()
                    );
                }
            }
            
            WireCell::Vector dir_1 = v1.magnitude() > 0 ? v1.norm() : WireCell::Vector(1, 0, 0);
            vec_dir.at(i) = dir_1;

            // Set up orthogonal coordinate system
            WireCell::Vector dir_2, dir_3;
            double angle_deg = std::acos(std::clamp(dir_1.dot(drift_dir_abs), -1.0, 1.0)) * 180.0 / M_PI;
            
            if (angle_deg < 7.5) {
                dir_1 = WireCell::Vector(1, 0, 0);
                dir_2 = WireCell::Vector(0, 1, 0);
                dir_3 = WireCell::Vector(0, 0, 1);
            } else {
                dir_2 = drift_dir_abs.cross(dir_1).norm();
                dir_3 = dir_1.cross(dir_2);
            }
            
            // Project associated points onto the local coordinate system
            std::vector<std::tuple<double, double, double>> vec_projs;
            for (const auto& pt : local_points_vec.at(i)) {
                double proj_1 = dir_1.dot(pt);
                double proj_2 = dir_2.dot(pt);
                double proj_3 = dir_3.dot(pt);
                vec_projs.push_back(std::make_tuple(proj_1, proj_2, proj_3));
            }
            
            // Calculate RMS spread in each direction
            int ncount = local_points_vec.at(i).size();
            if (ncount > 1) {
                WireCell::Point fit_pt(fits[i].point.x(), fits[i].point.y(), fits[i].point.z());
                std::tuple<double, double, double> means = std::make_tuple(
                    dir_1.dot(fit_pt),
                    dir_2.dot(fit_pt),
                    dir_3.dot(fit_pt)
                );
                
                for (const auto& proj : vec_projs) {
                    std::get<0>(vec_rms_vals.at(i)) += std::pow(std::get<0>(proj) - std::get<0>(means), 2);
                    std::get<1>(vec_rms_vals.at(i)) += std::pow(std::get<1>(proj) - std::get<1>(means), 2);
                    std::get<2>(vec_rms_vals.at(i)) += std::pow(std::get<2>(proj) - std::get<2>(means), 2);
                }
                
                std::get<0>(vec_rms_vals.at(i)) = std::sqrt(std::get<0>(vec_rms_vals.at(i)) / ncount);
                std::get<1>(vec_rms_vals.at(i)) = std::sqrt(std::get<1>(vec_rms_vals.at(i)) / ncount);
                std::get<2>(vec_rms_vals.at(i)) = std::sqrt(std::get<2>(vec_rms_vals.at(i)) / ncount);
            }
            
            // Calculate dQ/dx
            vec_dQ_dx.at(i) = fits[i].dQ / (fits[i].dx + 1e-9) / MIP_dQdx;
        }
        
        // Analyze spread characteristics
        double max_spread = 0;
        double large_spread_length = 0;
        double total_effective_length = 0;
        double total_length = 0;
        
        // bool flag_prev = false;
        for (size_t i = 0; i + 1 < local_points_vec.size(); i++) {
            double length = std::sqrt(
                std::pow(fits[i+1].point.x() - fits[i].point.x(), 2) +
                std::pow(fits[i+1].point.y() - fits[i].point.y(), 2) +
                std::pow(fits[i+1].point.z() - fits[i].point.z(), 2)
            );
            total_length += length;
            
            if (std::get<2>(vec_rms_vals.at(i)) != 0) {
                total_effective_length += length;
                if (std::get<2>(vec_rms_vals.at(i)) > rms_cut) {
                    large_spread_length += length;
                    // flag_prev = true;
                }
                if (std::get<2>(vec_rms_vals.at(i)) > max_spread) {
                    max_spread = std::get<2>(vec_rms_vals.at(i));
                }
            }
        }
        
        // Determine direction based on spread analysis
        int flag_dir = 0;
        

        // std::cout << "Shower Topology Direction: " << max_spread/units::cm << " " << large_spread_length/units::cm << " " << total_effective_length/units::cm << std::endl;

        // Check if this looks like a shower based on spread
        bool is_shower_like = (
            (max_spread > 0.7*units::cm && large_spread_length > 0.2 * total_effective_length && 
             total_effective_length > 3*units::cm && total_effective_length < 15*units::cm && 
             (large_spread_length > 2.7*units::cm || large_spread_length > 0.35 * total_effective_length)) ||
            (max_spread > 0.8*units::cm && large_spread_length > 0.3 * total_effective_length && 
             total_effective_length >= 15*units::cm) ||
            (max_spread > 1.0*units::cm && large_spread_length > 0.4 * total_effective_length)
        );
        
        if (is_shower_like) {
            WireCell::Vector main_dir1, main_dir2;
            bool flag_skip_angle1 = false;
            bool flag_skip_angle2 = false;
            
            // Create copies of points since segment_cal_dir_3vector expects non-const reference
            WireCell::Point front_pt = fits.front().point;
            WireCell::Point back_pt = fits.back().point;
            
            if (fits.front().point.z() < fits.back().point.z()) {
                main_dir1 = segment_cal_dir_3vector(segment, front_pt, 15*units::cm);
                main_dir2 = segment_cal_dir_3vector(segment, back_pt, 6*units::cm);
                double angle1 = std::acos(std::clamp(main_dir1.dot(drift_dir_abs), -1.0, 1.0)) * 180.0 / M_PI;
                if (std::fabs(angle1 - 90) < 10) flag_skip_angle1 = true;
            } else {
                main_dir1 = segment_cal_dir_3vector(segment, front_pt, 6*units::cm);
                main_dir2 = segment_cal_dir_3vector(segment, back_pt, 15*units::cm);
                double angle2 = std::acos(std::clamp(main_dir2.dot(drift_dir_abs), -1.0, 1.0)) * 180.0 / M_PI;
                if (std::fabs(angle2 - 90) < 10) flag_skip_angle2 = true;
            }
            
            // Group consecutive segments with large spread in forward direction
            // Each tuple: (start_index, end_index, max_rms_in_range)
            std::vector<std::tuple<int, int, double>> threshold_segs;
            const int sample_count = static_cast<int>(vec_dQ_dx.size());
            
            for (int i = 0; i < sample_count; i++) {
                double angle = std::acos(std::clamp(main_dir1.dot(vec_dir.at(i)), -1.0, 1.0)) * 180.0 / M_PI;
                if ((angle < 30 || (flag_skip_angle1 && angle < 60)) &&
                    (std::get<2>(vec_rms_vals.at(i))/units::cm > 0.4 || vec_dQ_dx.at(i) > 1.6)) {

                    if (threshold_segs.empty()) {
                        threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                    } else {
                        // Matches prototype: extend only when consecutive AND RMS is strictly
                        // larger than the current group max.  When consecutive but RMS dips, the
                        // point is silently skipped (neither extends nor starts a new group).
                        if (i == std::get<1>(threshold_segs.back()) + 1 &&
                            std::get<2>(threshold_segs.back()) < std::get<2>(vec_rms_vals.at(i))) {
                            std::get<1>(threshold_segs.back()) = i;
                            std::get<2>(threshold_segs.back()) = std::get<2>(vec_rms_vals.at(i));
                        } else {
                            if (i != std::get<1>(threshold_segs.back()) + 1)
                                threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                        }
                    }
                }
            }
            
            // Calculate total and max continuous length for forward direction
            double total_length1 = 0, max_length1 = 0;
            for (const auto& seg : threshold_segs) {
                int start_n = std::get<0>(seg);
                if (start_n > 0) start_n--;  // Include one segment before
                int end_n = std::get<1>(seg);
                
                double tmp_length = 0;
                for (int i = start_n; i < end_n && i + 1 < (int)fits.size(); i++) {
                    tmp_length += std::sqrt(
                        std::pow(fits[i+1].point.x() - fits[i].point.x(), 2) +
                        std::pow(fits[i+1].point.y() - fits[i].point.y(), 2) +
                        std::pow(fits[i+1].point.z() - fits[i].point.z(), 2)
                    );
                }
                total_length1 += tmp_length;
                if (tmp_length > max_length1) max_length1 = tmp_length;
            }
            
            // Group consecutive segments with large spread in backward direction
            threshold_segs.clear();
            
            for (int i = vec_dQ_dx.size() - 1; i >= 0; i--) {
                double angle = 180 - std::acos(std::clamp(main_dir2.dot(vec_dir.at(i)), -1.0, 1.0)) * 180.0 / M_PI;
                if ((angle < 30 || (flag_skip_angle2 && angle < 60)) &&
                    (std::get<2>(vec_rms_vals.at(i))/units::cm > 0.4 || vec_dQ_dx.at(i) > 1.6)) {

                    if (threshold_segs.empty()) {
                        threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                    } else {
                        if (i == std::get<1>(threshold_segs.back()) - 1 &&
                            std::get<2>(threshold_segs.back()) < std::get<2>(vec_rms_vals.at(i))) {
                            std::get<1>(threshold_segs.back()) = i;
                            std::get<2>(threshold_segs.back()) = std::get<2>(vec_rms_vals.at(i));
                        } else {
                            if (i != std::get<1>(threshold_segs.back()) - 1)
                                threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                        }
                    }
                }
            }
            
            // Calculate total and max continuous length for backward direction
            double total_length2 = 0, max_length2 = 0;
            for (const auto& seg : threshold_segs) {
                int start_n = std::get<0>(seg);
                if (start_n < (int)fits.size() - 1) start_n++;  // Include one segment after
                int end_n = std::get<1>(seg);
                
                double tmp_length = 0;
                for (int i = start_n; i > end_n && i > 0; i--) {
                    tmp_length += std::sqrt(
                        std::pow(fits[i-1].point.x() - fits[i].point.x(), 2) +
                        std::pow(fits[i-1].point.y() - fits[i].point.y(), 2) +
                        std::pow(fits[i-1].point.z() - fits[i].point.z(), 2)
                    );
                }
                total_length2 += tmp_length;
                if (tmp_length > max_length2) max_length2 = tmp_length;
            }
            
            // Compare using both total and max continuous lengths
            if (total_length1 + max_length1 > 1.1 * (total_length2 + max_length2)) {
                flag_dir = 1;
            } else if (1.1 * (total_length1 + max_length1) < total_length2 + max_length2) {
                flag_dir = -1;
            }
        } else {
            // Not shower-like, use simpler direction determination
            if (total_length < 5*units::cm) {
                // docs/pr/2 sec 8.1a: forward the configured MIP scales instead of the
                // uBooNE header defaults (50000/43000 e/cm).  uBooNE byte-identical:
                // the forwarded values equal the old defaults there.
                // doc pr/31 §12 (F4): this interior call passes a
                // default-constructed TrackPidOptions (NOT the caller's full
                // option set -- forwarding that here would unconditionally
                // switch proton_dir_vote et al. at this site), carrying ONLY
                // the median-source choice.  median_local=false => defaults
                // all the way => byte-identical.
                if (!segment_is_shower_trajectory(segment, 10*units::cm, mip_dqdx)) {
                    TrackPidOptions med_opts{};
                    med_opts.dir_track_median_local = median_local;
                    segment_determine_dir_track(segment, 0, fits.size(), particle_data, recomb_model, MIP_dQdx, false, med_opts);
                }
                // For short segments, could call determine_dir_track here if needed
            } else {
                // Count consistent directions at each end
                WireCell::Point front_pt = fits.front().point;
                WireCell::Point back_pt = fits.back().point;
                WireCell::Vector main_dir_front = segment_cal_dir_3vector(segment, front_pt, 6*units::cm); 
                int ncount_front = 0;
                for (size_t i = 0; i < vec_dQ_dx.size(); i++) {
                    double angle = std::acos(std::clamp(main_dir_front.dot(vec_dir.at(i)), -1.0, 1.0)) * 180.0 / M_PI;
                    if (angle < 30) ncount_front++;
                }

                WireCell::Vector main_dir_back = segment_cal_dir_3vector(segment, back_pt, 6*units::cm);
                int ncount_back = 0;
                for (int i = vec_dQ_dx.size() - 1; i >= 0; i--) {
                    double angle = 180 - std::acos(std::clamp(main_dir_back.dot(vec_dir.at(i)), -1.0, 1.0)) * 180.0 / M_PI;
                    if (angle < 30) ncount_back++;
                }
                
                if (1.2 * ncount_front < ncount_back) {
                    flag_dir = -1;
                } else if (ncount_front > 1.2 * ncount_back) {
                    flag_dir = 1;
                }
            }
        }
        
        segment->dirsign(flag_dir);
        return (flag_dir != 0);
    }

    bool segment_is_shower_topology(SegmentPtr segment, bool tmp_val, double MIP_dQ_dx, double demote_len, bool reset, bool dqdx_guard){
        int flag_dir = 0;
        bool flag_shower_topology = tmp_val;
        // doc sbnd_xin/docs/pr/31 §12 (F3, was P13).  The prototype's first two
        // statements are segment-state assignments made BEFORE any early return
        // (ProtoSegment.cxx:319-321): flag_shower_topology = tmp_val (every
        // caller passes false, so the flag is cleared on entry and re-set only
        // if the test still passes) and flag_dir = 0.  The toolkit port assigns
        // a LOCAL and the tail only ever set_flags, so a segment that qualified
        // on an earlier pass keeps a stale kShowerTopology flag -- which is what
        // routes determine_direction into the topology branch -- and a segment
        // whose clouds hit the four early returns below keeps a stale direction.
        // Re-entry is the normal path: this function runs in stage 3
        // (separate_track_shower) and at three stage-4 sites.  unset_flags
        // clears exactly the named bit (Flagged.h: other flags survive);
        // never clear_flags() here.  reset=false => set-only legacy path =>
        // byte-identical.
        if (reset) {
            if (!tmp_val) segment->unset_flags(SegmentFlags::kShowerTopology);
            segment->dirsign(0);
        }
        const auto& fits = segment->fits();

        if (fits.empty()) return false;
        
        // Get the fit point cloud for KD-tree queries
        auto dpcloud_fit = segment->dpcloud("fit");
        if (!dpcloud_fit) return false;
        
        // Get the associated point cloud 
        auto dpcloud_assoc = segment->dpcloud("associate_points");
        if (!dpcloud_assoc) return false;
        
        const size_t assoc_npts = dpcloud_assoc->npoints();
        if (assoc_npts == 0) return false;

        // sbnd_xin doc pr/25 sec 3: WCT_SHOWER_TOPO_DEBUG=1 dumps the decision
        // terms segment_is_shower_topology never otherwise logs (this function
        // had zero active log lines before this change), so the mechanism
        // behind a topology-vs-trajectory misclassification can be measured
        // directly instead of inferred from a 3-D geometric proxy (the trap
        // recorded after doc pr/25 sec 2 -- a proxy predicted "should not
        // fire" on SBND evt 321107/cluster 13, which did fire).  Read once,
        // env-gated, default off => zero behavior change either way.
        // Round 2 hoisted the flag above the per-fit-point loop so the loop can
        // record |dir_3.xhat|: dir_3 = dir_1 x (xhat x dir_1), so
        // dir_3.xhat = sin(angle-to-drift) and the sole axis this function
        // measures collapses onto the drift axis for an isochronous segment.
        static const bool s_shower_topo_dbg = (std::getenv("WCT_SHOWER_TOPO_DEBUG") != nullptr);
        std::vector<double> dbg_dir3x;  // |dir_3 . xhat| per populated bucket

        // Initialize vectors to store analysis results for each fit point
        std::vector<std::vector<WireCell::Point>> local_points_vec(fits.size());
        std::vector<std::tuple<double, double, double>> vec_rms_vals(fits.size(), std::make_tuple(0,0,0));
        std::vector<double> vec_dQ_dx(fits.size(), 0);
        
        // Build KD-tree index for fit points and associate points with nearest fit point
        auto& kd_tree_fit = dpcloud_fit->kd3d();
        
        for (size_t ia = 0; ia < assoc_npts; ++ia) {
            WireCell::Point test_p = dpcloud_assoc->point3d(ia);
            auto results = kd_tree_fit.knn(1, test_p);
            if (!results.empty()) {
                size_t closest_fit_idx = results.front().first;
                local_points_vec.at(closest_fit_idx).push_back(test_p);
            }
        }
        
        WireCell::Vector drift_dir_abs(1, 0, 0);  // Drift direction
        
        // Calculate local directions and RMS spreads for each fit point
        for (size_t i = 0; i < local_points_vec.size(); i++) {
            // Calculate local direction from neighboring fit points
            WireCell::Vector v1(0, 0, 0);
            for (size_t j = 1; j < 3; j++) {
                if (i + j < fits.size()) {
                    v1 += WireCell::Vector(
                        fits[i+j].point.x() - fits[i].point.x(),
                        fits[i+j].point.y() - fits[i].point.y(),
                        fits[i+j].point.z() - fits[i].point.z()
                    );
                }
                if (i >= j) {
                    v1 += WireCell::Vector(
                        fits[i].point.x() - fits[i-j].point.x(),
                        fits[i].point.y() - fits[i-j].point.y(),
                        fits[i].point.z() - fits[i-j].point.z()
                    );
                }
            }
            
            WireCell::Vector dir_1 = v1.magnitude() > 0 ? v1.norm() : WireCell::Vector(1, 0, 0);

            // Set up orthogonal coordinate system
            WireCell::Vector dir_2, dir_3;
            double angle_deg = std::acos(std::clamp(dir_1.dot(drift_dir_abs), -1.0, 1.0)) * 180.0 / M_PI;
            
            if (angle_deg < 7.5) {
                dir_1 = WireCell::Vector(1, 0, 0);
                dir_2 = WireCell::Vector(0, 1, 0);
                dir_3 = WireCell::Vector(0, 0, 1);
            } else {
                dir_2 = drift_dir_abs.cross(dir_1).norm();
                dir_3 = dir_1.cross(dir_2);
            }
            
            // Project associated points onto the local coordinate system
            std::vector<std::tuple<double, double, double>> vec_projs;
            for (const auto& pt : local_points_vec.at(i)) {
                double proj_1 = dir_1.dot(pt);
                double proj_2 = dir_2.dot(pt);
                double proj_3 = dir_3.dot(pt);
                vec_projs.push_back(std::make_tuple(proj_1, proj_2, proj_3));
            }
            
            // Calculate RMS spread in each direction
            int ncount = local_points_vec.at(i).size();
            if (ncount > 1) {
                WireCell::Point fit_pt(fits[i].point.x(), fits[i].point.y(), fits[i].point.z());
                std::tuple<double, double, double> means = std::make_tuple(
                    dir_1.dot(fit_pt),
                    dir_2.dot(fit_pt),
                    dir_3.dot(fit_pt)
                );
                
                for (const auto& proj : vec_projs) {
                    std::get<0>(vec_rms_vals.at(i)) += std::pow(std::get<0>(proj) - std::get<0>(means), 2);
                    std::get<1>(vec_rms_vals.at(i)) += std::pow(std::get<1>(proj) - std::get<1>(means), 2);
                    std::get<2>(vec_rms_vals.at(i)) += std::pow(std::get<2>(proj) - std::get<2>(means), 2);
                }
                
                std::get<0>(vec_rms_vals.at(i)) = std::sqrt(std::get<0>(vec_rms_vals.at(i)) / ncount);
                std::get<1>(vec_rms_vals.at(i)) = std::sqrt(std::get<1>(vec_rms_vals.at(i)) / ncount);
                std::get<2>(vec_rms_vals.at(i)) = std::sqrt(std::get<2>(vec_rms_vals.at(i)) / ncount);

                // doc pr/25 sec 3 round 2 diagnostic only -- see the hoisted
                // s_shower_topo_dbg comment above.
                if (s_shower_topo_dbg) dbg_dir3x.push_back(std::abs(dir_3.dot(drift_dir_abs)));
            }

            // Calculate dQ/dx
            vec_dQ_dx.at(i) = fits[i].dQ / (fits[i].dx + 1e-9) / MIP_dQ_dx;
        }
        
        // Analyze spread characteristics
        double max_spread = 0;
        double large_spread_length = 0;
        double total_effective_length = 0;
        
        double max_cont_length = 0;
        double max_cont_weighted_length = 0;
        double cont_length = 0;
        double cont_weighted_length = 0;
        bool flag_prev = false;
        
        for (size_t i = 0; i + 1 < local_points_vec.size(); i++) {
            double length = std::sqrt(
                std::pow(fits[i+1].point.x() - fits[i].point.x(), 2) +
                std::pow(fits[i+1].point.y() - fits[i].point.y(), 2) +
                std::pow(fits[i+1].point.z() - fits[i].point.z(), 2)
            );
            
            if (std::get<2>(vec_rms_vals.at(i)) != 0) {
                total_effective_length += length;
                if (std::get<2>(vec_rms_vals.at(i)) > 0.4 * units::cm) {
                    large_spread_length += length;
                    cont_length += length;
                    cont_weighted_length += length * std::get<2>(vec_rms_vals.at(i));
                    flag_prev = true;
                } else {
                    if (flag_prev && cont_length > max_cont_length) {
                        max_cont_length = cont_length;
                        max_cont_weighted_length = cont_weighted_length;
                    }
                    cont_length = 0;
                    cont_weighted_length = 0;
                    flag_prev = false;
                }
                if (std::get<2>(vec_rms_vals.at(i)) > max_spread) {
                    max_spread = std::get<2>(vec_rms_vals.at(i));
                }
            }
        }
        (void)max_cont_weighted_length; // Currently unused

        // std::cout << "Shower Topology Check: " << max_spread/units::cm << " " << large_spread_length/units::cm << " " << total_effective_length/units::cm << std::endl;

        int dbg_branch = 0;
        double dbg_frac_lsl = (total_effective_length > 0) ? large_spread_length / total_effective_length : 0;
        if (s_shower_topo_dbg) {
            if (max_spread > 0.7*units::cm && dbg_frac_lsl > 0.2 &&
                total_effective_length > 3*units::cm && total_effective_length < 15*units::cm &&
                (large_spread_length > 2.7*units::cm || dbg_frac_lsl > 0.35)) dbg_branch = 1;
            else if (max_spread > 0.8*units::cm && dbg_frac_lsl > 0.3 &&
                      total_effective_length >= 15*units::cm) dbg_branch = 2;
            else if (max_spread > 0.8*units::cm && large_spread_length > 8*units::cm &&
                      dbg_frac_lsl > 0.18) dbg_branch = 3;
            else if (max_spread > 1.0*units::cm && dbg_frac_lsl > 0.4) dbg_branch = 4;
            else if (max_spread > 1.0*units::cm && large_spread_length > 5*units::cm &&
                      dbg_frac_lsl > 0.23) dbg_branch = 5;

            // Per-bucket dir_3 RMS distribution over populated fit-point
            // buckets -- needed to tell "one bucket carries max_spread" from
            // "fifty do", which max_spread alone cannot.
            std::vector<double> dbg_rms;
            for (size_t i = 0; i + 1 < local_points_vec.size(); i++) {
                double r = std::get<2>(vec_rms_vals.at(i));
                if (r != 0) dbg_rms.push_back(r);
            }
            std::sort(dbg_rms.begin(), dbg_rms.end());
            auto n_over_cut = [&](double c) -> size_t {
                return std::count_if(dbg_rms.begin(), dbg_rms.end(),
                                     [c](double r) { return r > c*units::cm; });
            };
            auto pct = [&](double q) -> double {
                if (dbg_rms.empty()) return 0;
                size_t idx = std::min(dbg_rms.size() - 1, static_cast<size_t>(q * dbg_rms.size()));
                return dbg_rms[idx];
            };

            // Round 2: the longest CONTIGUOUS large-spread run.  The function
            // already computes max_cont_length above and discards it -- and
            // that copy only records a run when the run *ends*, so a run
            // reaching the last bucket is never counted.  Recomputed here
            // without that omission, for measurement only.
            double dbg_max_cont = 0, dbg_cur_cont = 0;
            for (size_t i = 0; i + 1 < local_points_vec.size(); i++) {
                if (std::get<2>(vec_rms_vals.at(i)) == 0) continue;
                double length = std::sqrt(
                    std::pow(fits[i+1].point.x() - fits[i].point.x(), 2) +
                    std::pow(fits[i+1].point.y() - fits[i].point.y(), 2) +
                    std::pow(fits[i+1].point.z() - fits[i].point.z(), 2));
                if (std::get<2>(vec_rms_vals.at(i)) > 0.4*units::cm) {
                    dbg_cur_cont += length;
                    if (dbg_cur_cont > dbg_max_cont) dbg_max_cont = dbg_cur_cont;
                } else {
                    dbg_cur_cont = 0;
                }
            }

            // Median |dir_3 . xhat| = median sin(angle-to-drift): 1.0 means the
            // measured axis IS the drift axis (isochronous segment).
            double dbg_dir3x_med = 0;
            if (!dbg_dir3x.empty()) {
                std::vector<double> tmp = dbg_dir3x;
                std::sort(tmp.begin(), tmp.end());
                dbg_dir3x_med = tmp[tmp.size()/2];
            }

            double dbg_geom_len = segment_track_length(segment, 0);
            SPDLOG_LOGGER_INFO(s_log,
                "shower_topo dbg: seg {} L {:.1f}cm assoc_npts {} nbuckets {} n_over0.4cm {} "
                "n_over0.7cm {} n_over0.8cm {} n_over1.0cm {} "
                "rms_p50 {:.2f}cm rms_p75 {:.2f}cm rms_p90 {:.2f}cm rms_p95 {:.2f}cm "
                "max_spread {:.2f}cm maxcont {:.1f}cm lsl {:.1f}cm tel {:.1f}cm "
                "lsl/tel {:.3f} tel/L {:.3f} dir3x {:.4f} branch {}",
                segment->id(), dbg_geom_len/units::cm, assoc_npts, dbg_rms.size(),
                n_over_cut(0.4), n_over_cut(0.7), n_over_cut(0.8), n_over_cut(1.0),
                pct(0.5)/units::cm, pct(0.75)/units::cm, pct(0.9)/units::cm, pct(0.95)/units::cm,
                max_spread/units::cm, dbg_max_cont/units::cm,
                large_spread_length/units::cm, total_effective_length/units::cm,
                dbg_frac_lsl, (dbg_geom_len > 0 ? total_effective_length/dbg_geom_len : 0),
                dbg_dir3x_med, dbg_branch);
        }

        // Determine if this is shower topology based on spread patterns
        if ((max_spread > 0.7*units::cm && large_spread_length > 0.2 * total_effective_length && 
             total_effective_length > 3*units::cm && total_effective_length < 15*units::cm && 
             (large_spread_length > 2.7*units::cm || large_spread_length > 0.35 * total_effective_length)) ||
            (max_spread > 0.8*units::cm && large_spread_length > 0.3 * total_effective_length && 
             total_effective_length >= 15*units::cm) ||
            (max_spread > 0.8*units::cm && large_spread_length > 8*units::cm && 
             large_spread_length > 0.18 * total_effective_length) ||
            (max_spread > 1.0*units::cm && large_spread_length > 0.4 * total_effective_length) ||
            (max_spread > 1.0*units::cm && large_spread_length > 5*units::cm && 
             large_spread_length > 0.23 * total_effective_length)) {
            
            flag_shower_topology = true;
        }
        
        // If identified as shower topology, determine direction
        if (flag_shower_topology) {
            // Group consecutive segments with large spread in forward direction
            std::vector<std::tuple<int, int, double>> threshold_segs;
            const int sample_count = static_cast<int>(vec_dQ_dx.size());
            
            for (int i = 0; i < sample_count; i++) {
                if (std::get<2>(vec_rms_vals.at(i))/units::cm > 0.4) {
                    if (threshold_segs.empty()) {
                        threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                    } else {
                        if (i == std::get<1>(threshold_segs.back()) + 1 &&
                            std::get<2>(threshold_segs.back()) < std::get<2>(vec_rms_vals.at(i))) {
                            std::get<1>(threshold_segs.back()) = i;
                            std::get<2>(threshold_segs.back()) = std::get<2>(vec_rms_vals.at(i));
                        } else {
                            if (i != std::get<1>(threshold_segs.back()) + 1)
                                threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                        }
                    }
                }
            }

            // Calculate total and max continuous length for forward direction
            double total_length1 = 0, max_length1 = 0;
            for (const auto& seg : threshold_segs) {
                int start_n = std::get<0>(seg);
                if (start_n > 0) start_n--;
                int end_n = std::get<1>(seg);
                
                double tmp_length = 0;
                for (int i = start_n; i < end_n && i + 1 < (int)fits.size(); i++) {
                    tmp_length += std::sqrt(
                        std::pow(fits[i+1].point.x() - fits[i].point.x(), 2) +
                        std::pow(fits[i+1].point.y() - fits[i].point.y(), 2) +
                        std::pow(fits[i+1].point.z() - fits[i].point.z(), 2)
                    );
                }
                total_length1 += tmp_length;
                if (tmp_length > max_length1) max_length1 = tmp_length;
            }
            
            // Group consecutive segments with large spread in backward direction
            threshold_segs.clear();
            
            for (int i = vec_dQ_dx.size() - 1; i >= 0; i--) {
                if (std::get<2>(vec_rms_vals.at(i))/units::cm > 0.4) {
                    if (threshold_segs.empty()) {
                        threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                    } else {
                        if (i == std::get<1>(threshold_segs.back()) - 1 &&
                            std::get<2>(threshold_segs.back()) < std::get<2>(vec_rms_vals.at(i))) {
                            std::get<1>(threshold_segs.back()) = i;
                            std::get<2>(threshold_segs.back()) = std::get<2>(vec_rms_vals.at(i));
                        } else {
                            if (i != std::get<1>(threshold_segs.back()) - 1)
                                threshold_segs.push_back(std::make_tuple(i, i, std::get<2>(vec_rms_vals.at(i))));
                        }
                    }
                }
            }
            
            // Calculate total and max continuous length for backward direction
            double total_length2 = 0, max_length2 = 0;
            for (const auto& seg : threshold_segs) {
                int start_n = std::get<0>(seg);
                if (start_n < (int)fits.size() - 1) start_n++;
                int end_n = std::get<1>(seg);
                
                double tmp_length = 0;
                for (int i = start_n; i > end_n && i > 0; i--) {
                    tmp_length += std::sqrt(
                        std::pow(fits[i-1].point.x() - fits[i].point.x(), 2) +
                        std::pow(fits[i-1].point.y() - fits[i].point.y(), 2) +
                        std::pow(fits[i-1].point.z() - fits[i].point.z(), 2)
                    );
                }
                total_length2 += tmp_length;
                if (tmp_length > max_length2) max_length2 = tmp_length;
            }
            
            // Determine direction based on spread comparison
            if (total_length1 + max_length1 > 1.1 * (total_length2 + max_length2)) {
                flag_dir = 1;
            } else if (1.1 * (total_length1 + max_length1) < total_length2 + max_length2) {
                flag_dir = -1;
            }
            
            // Override shower topology for very long segments with little spread
            double tmp_total_length = segment_track_length(segment, 0);
            bool dbg_guard_demoted = false;
            if (tmp_total_length > 50*units::cm &&
                total_length1 < 0.25 * tmp_total_length &&
                total_length2 < 0.25 * tmp_total_length) {
                flag_dir = 0;
                flag_shower_topology = false;
                dbg_guard_demoted = true;
            }

            // sbnd_xin doc pr/25 sec 3: make that guard unconditional above a
            // configurable length.  The legacy guard only demotes when the
            // large-spread runs are BOTH under 0.25*L, a condition that on
            // SBND evt 321107 missed by 0.008 (0.258/0.256) -- and the
            // discriminant behind it is drift-direction quantization noise
            // (sec 3.2-3.4: dir_3.xhat = sin(angle-to-drift), so for an
            // isochronous segment the only measured axis IS the drift axis,
            // where points are quantized on a 0.313cm time-slice lattice
            // against a 0.4cm cut).  An owner hand-scan of all 10 such
            // nu-candidate-main segments found 10/10 tracks, 0 showers.
            // Default 0 => never fires => byte-identical legacy path.
            if (demote_len > 0 && tmp_total_length > demote_len) {
                flag_dir = 0;
                flag_shower_topology = false;
                dbg_guard_demoted = true;
            }
            if (s_shower_topo_dbg) {
                SPDLOG_LOGGER_INFO(s_log,
                    "shower_topo dbg: seg {} guard branch {} L {:.1f}cm total_length1 {:.1f}cm({:.3f}) "
                    "total_length2 {:.1f}cm({:.3f}) demoted {} final_shower {}",
                    segment->id(), dbg_branch, tmp_total_length/units::cm,
                    total_length1/units::cm, (tmp_total_length > 0 ? total_length1/tmp_total_length : 0),
                    total_length2/units::cm, (tmp_total_length > 0 ? total_length2/tmp_total_length : 0),
                    dbg_guard_demoted, flag_shower_topology);
            }
        }

        // doc sbnd_xin/docs/pr/40 F3: the dQ/dx cross-check the 5-branch
        // geometric spread test never performs (vec_dQ_dx above is otherwise
        // dead -- pr/31 GOTCHA 5).  Runs regardless of which length-guard
        // branch fired above (both existing guards require length > their own
        // cut; this one does not).  MIP_dQ_dx is this function's own scale
        // parameter (mip_dqdx_median), consistent with every other threshold
        // here.  Does not touch flag_dir -- only whether the flag is set.
        // false = legacy = byte-identical.
        if (flag_shower_topology && dqdx_guard && segment_dqdx_spares_electron_reclass(segment, MIP_dQ_dx)) {
            flag_shower_topology = false;
        }

        if (flag_shower_topology) segment->set_flags(SegmentFlags::kShowerTopology);
        segment->dirsign(flag_dir);
        return flag_shower_topology;
    }

}
