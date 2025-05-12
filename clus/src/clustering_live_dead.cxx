#include "WireCellClus/IClusteringMethod.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"

class ClusteringLiveDead;
WIRECELL_FACTORY(ClusteringLiveDead, ClusteringLiveDead,
                 WireCell::IConfigurable, WireCell::Clus::IClusteringMethod)


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

    
class ClusteringLiveDead :  public IConfigurable, public Clus::IClusteringMethod, private NeedDV, private NeedScope {
    int dead_live_overlap_offset_{2};
public:
    ClusteringLiveDead() {}
    virtual ~ClusteringLiveDead() {}

    virtual void configure(const WireCell::Configuration& config) {
        NeedDV::configure(config);
        NeedScope::configure(config);
        
        dead_live_overlap_offset_ = get(config, "dead_live_overlap_offset", 2);
    }
    virtual Configuration default_configuration() const {
        Configuration cfg;
        return cfg;
    }
    

    virtual void clustering(Grouping& live_grouping, Grouping& dead_grouping) const {
        using spdlog::debug;

        // {
        //     auto wpids = live_grouping.wpids();
        //     for (const auto& wpid : wpids) {
        //         int apa = wpid.apa();
        //         int face = wpid.face();

        //         std::cout << "Test Number of Points: " << apa << " " << face << " " << live_grouping.get_num_points(apa, face, 0) << " " << live_grouping.get_num_points(apa, face, 1) << " " << live_grouping.get_num_points(apa, face, 2) << std::endl;

        //     }
        // }

        // check if the grouping's wpid ... 
        //std::cout << "Live: " << live_grouping.wpids().size() << " " << dead_grouping.wpids().size() << std::endl;
    
        // Check that groupings has less than one wpid
        if (live_grouping.wpids().size() > 1 || dead_grouping.wpids().size() > 1) {
            for (const auto& wpid : live_grouping.wpids()) {
                std::cout << "Live grouping wpid: " << wpid.name() << std::endl;
            }
            for (const auto& wpid : dead_grouping.wpids()) {
                std::cout << "Dead grouping wpid: " << wpid.name() << std::endl;
            }
            raise<ValueError>("Live %d > 1, Dead %d > 1", live_grouping.wpids().size(), dead_grouping.wpids().size());
        }
        auto [drift_dir, angle_u, angle_v, angle_w] = extract_geometry_params(live_grouping, m_dv);
    


        // form map from dead to set of live clusters ...
        std::map<const Cluster*, std::vector<Cluster*>> dead_live_cluster_mapping;
        std::vector<const Cluster*> dead_cluster_order;
        std::map<const Cluster*, std::vector<std::vector<const Blob*>>> dead_live_mcells_mapping;

        std::vector<Cluster*> live_clusters = live_grouping.children(); // copy

        for (auto& cluster : live_clusters) {
            if (cluster->get_default_scope().hash() != m_scope.hash()) {
                cluster->set_default_scope(m_scope);
                // std::cout << "Test: Set default scope: " << m_scope.pcname << " " << m_scope.coords[0] << " " << m_scope.coords[1] << " " << m_scope.coords[2] << " " << cluster->get_default_scope().hash() << " " << m_scope.hash() << std::endl;
            }
        }


        std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *cluster1, const Cluster *cluster2) {
            return cluster1->get_length() > cluster2->get_length();
        });
        // sort_clusters(live_clusters);

        auto dead_clusters = dead_grouping.children(); // copy
        sort_clusters(dead_clusters);

        for (size_t ilive = 0; ilive < live_clusters.size(); ++ilive) {
            const auto& live = live_clusters.at(ilive);
            if (!live->get_scope_filter(m_scope)) continue;
            for (size_t idead = 0; idead < dead_clusters.size(); ++idead) {
                const auto& dead = dead_clusters.at(idead);
	
                auto blobs = live->is_connected(*dead, dead_live_overlap_offset_);
                if (blobs.size() > 0) {
                    if (dead_live_cluster_mapping.find(dead) == dead_live_cluster_mapping.end()) {
                        dead_cluster_order.push_back(dead);
                    }
                    dead_live_cluster_mapping[dead].push_back(live);
                    dead_live_mcells_mapping[dead].push_back(blobs);
                }
            }
        }


        if (dead_live_cluster_mapping.empty()) {
            std::cerr
                << "WARNING: clustering_live: empty dead live cluster mapping,"
                << " ndead=" << dead_clusters.size()
                << " nlive=" << live_clusters.size() << std::endl;
        }

        // prepare a graph ...
        typedef cluster_connectivity_graph_t Graph;

        Graph g;
        std::unordered_map<int, int> ilive2desc;  // added live index to graph descriptor
        std::map<const Cluster*, int> map_cluster_index;
        for (const Cluster* live : live_grouping.children()) {
            size_t ilive = map_cluster_index.size();
            map_cluster_index[live] = ilive;
            ilive2desc[ilive] = boost::add_vertex(ilive, g);
        }

        std::set<std::pair<const Cluster*, const Cluster* > > tested_pairs;

        // start to form edges ...
        for (const auto& the_dead_cluster : dead_cluster_order) {
            const auto& connected_live_clusters = dead_live_cluster_mapping[the_dead_cluster];
            const auto& connected_live_mcells = dead_live_mcells_mapping[the_dead_cluster];

            if (connected_live_clusters.size() > 1) {

                for (size_t i = 0; i != connected_live_clusters.size(); i++) {
                    auto cluster_1 = connected_live_clusters.at(i);
                    const auto& blobs_1 = connected_live_mcells.at(i);
        //            cluster_connected_dead.insert(cluster_1);
                    cluster_1->set_flag(Flags::live_dead);
                    for (size_t j = i + 1; j < connected_live_clusters.size(); j++) {
                        const auto& cluster_2 = connected_live_clusters.at(j);
        

                        if (tested_pairs.find(std::make_pair(cluster_1, cluster_2)) == tested_pairs.end()) {
                            tested_pairs.insert(std::make_pair(cluster_1, cluster_2));
                            tested_pairs.insert(std::make_pair(cluster_2, cluster_1));

                            bool flag_merge = false;
                            const Blob* prev_mcell1 = 0;
                            const Blob* prev_mcell2 = 0;
                            const Blob* mcell1 = blobs_1.at(0);
                            const Blob* mcell2 = 0;

                            geo_point_t p1 = mcell1->center_pos();
                            std::tie(p1, mcell1) = cluster_1->get_closest_point_blob(p1);
                            //	      p1 = temp_pair.first;
                            //	      mcell1 = temp_pair.second;
                            geo_point_t p2(0, 0, 0);
                            while (mcell1 != prev_mcell1 || mcell2 != prev_mcell2) {
                                prev_mcell1 = mcell1;
                                prev_mcell2 = mcell2;

                                std::tie(p2, mcell2) = cluster_2->get_closest_point_blob(p1);
                                std::tie(p1, mcell1) = cluster_1->get_closest_point_blob(p2);
                            }
                            geo_point_t diff = p1 - p2;
                            double dis = diff.magnitude();

                            // std::cout << p1 << " " << p2 << " " << dis/units::cm << std::endl;

                            if (dis < 60 * units::cm) {
                                const double length_1 = cluster_1->get_length();
                                const double length_2 = cluster_2->get_length();

                                geo_point_t mcell1_center = cluster_1->calc_ave_pos(p1, 5 * units::cm);
                                geo_point_t dir1 = cluster_1->vhough_transform(mcell1_center, 30 * units::cm);
                                // protection against angles ...
                                geo_point_t dir5 = cluster_1->vhough_transform(p1, 30 * units::cm);
                                if (dir1.angle(dir5) > 120 / 180. * 3.1415926) dir1 = dir1 * (-1);

                                geo_point_t mcell2_center = cluster_2->calc_ave_pos(p2, 5 * units::cm);
                                geo_point_t dir3 = cluster_2->vhough_transform(mcell2_center, 30 * units::cm);
                                // Protection against angles
                                geo_point_t dir6 = cluster_2->vhough_transform(p2, 30 * units::cm);
                                if (dir3.angle(dir6) > 120 / 180. * 3.1415926) dir3 = dir3 * (-1);

                                geo_point_t dir2 = mcell2_center - mcell1_center;
                                geo_point_t dir4 = mcell1_center - mcell2_center;

                                double angle_diff1 = (3.1415926 - dir1.angle(dir2)) / 3.1415926 * 180.;  // 1 to 2
                                double angle_diff2 = (3.1415926 - dir3.angle(dir4)) / 3.1415926 * 180.;  // 2 to 1
                                double angle_diff3 = (3.1415926 - dir1.angle(dir3)) / 3.1415926 * 180.;  // 1 to 2

                                bool flag_para = false;

                                double angle1, angle2, angle3;
                                if (!flag_merge) {
                               
                                    angle1 = dir1.angle(drift_dir);
                                    angle2 = dir2.angle(drift_dir);
                                    angle3 = dir3.angle(drift_dir);

                                    if (fabs(angle1 - 3.1415926 / 2.) < 5 / 180. * 3.1415926 &&
                                        fabs(angle2 - 3.1415926 / 2.) < 5 / 180. * 3.1415926 &&
                                        fabs(angle3 - 3.1415926 / 2.) < 5 / 180. * 3.1415926) {
                                        if (dis < 10 * units::cm)  // if very parallel and close, merge any way
                                            flag_merge = true;
                                    }

                                    if (fabs(angle2 - 3.1415926 / 2.) < 7.5 / 180. * 3.1415926 &&
                                        (fabs(angle1 - 3.1415926 / 2.) < 7.5 / 180. * 3.1415926 ||
                                         fabs(angle3 - 3.1415926 / 2.) < 7.5 / 180. * 3.1415926) &&
                                        fabs(angle1 - 3.1415926 / 2.) + fabs(angle2 - 3.1415926 / 2.) +
                                        fabs(angle3 - 3.1415926 / 2.) <
                                        25 / 180. * 3.1415926) {
                                        flag_para = true;

                                        if (WireCell::Clus::Facade::is_angle_consistent(
                                                dir1, dir2, false, 15, angle_u, angle_v, angle_w, 3) &&
                                            WireCell::Clus::Facade::is_angle_consistent(
                                                dir3, dir2, true, 15, angle_u, angle_v, angle_w, 3))
                                            flag_merge = true;
                                    }
                                    else {
                                        bool flag_const1 = WireCell::Clus::Facade::is_angle_consistent(
                                            dir1, dir2, false, 10, angle_u, angle_v, angle_w, 2);
                                        bool flag_const2 = WireCell::Clus::Facade::is_angle_consistent(
                                            dir3, dir2, true, 10, angle_u, angle_v, angle_w, 2);

                                        if (flag_const1 && flag_const2) {
                                            flag_merge = true;
                                        }
                                        else if (flag_const1 && length_2 < 6 * units::cm && length_1 > 15 * units::cm) {
                                            if (WireCell::Clus::Facade::is_angle_consistent(
                                                    dir1, dir2, false, 5, angle_u, angle_v, angle_w, 3))
                                                flag_merge = true;
                                        }
                                        else if (flag_const2 && length_1 < 6 * units::cm && length_2 > 15 * units::cm) {
                                            if (WireCell::Clus::Facade::is_angle_consistent(
                                                    dir3, dir2, true, 5, angle_u, angle_v, angle_w, 3))
                                                flag_merge = true;
                                        }
                                    }
                                }

// This block of code comes from the prototype and should have parentheses
// applied to make the logic explicit but nobody wants to do that so we tell the
// compiler to be quiet about it.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

                                if (!flag_merge) {
                                    if (length_1 <= 12 * units::cm && length_2 <= 12 * units::cm) {
                                        // both are short
                                        if ((dis <= 3 * units::cm) &&
                                            ((angle_diff1 <= 45 || angle_diff2 <= 45) && (angle_diff3 < 60) ||
                                             (flag_para && (angle_diff1 <= 90 || angle_diff2 <= 90) &&
                                              angle_diff3 < 120)) ||
                                            (dis <= 5 * units::cm) && (angle_diff1 <= 30 || angle_diff2 <= 30) &&
                                            angle_diff3 < 45 ||
                                            (dis <= 15 * units::cm) && (angle_diff1 <= 15 || angle_diff2 <= 15) &&
                                            angle_diff3 < 20 ||
                                            (dis <= 60 * units::cm) && (angle_diff1 < 5 || angle_diff2 < 5) &&
                                            angle_diff3 < 10) {
                                            flag_merge = true;
                                        }
                                    }
                                    else if (length_1 > 12 * units::cm && length_2 <= 12 * units::cm) {
                                        // one is short
                                        if ((dis <= 3 * units::cm) &&
                                            ((angle_diff1 <= 45 || angle_diff2 <= 45) && (angle_diff3 < 60) ||
                                             (flag_para && (angle_diff1 <= 90 || angle_diff2 <= 90) &&
                                              angle_diff3 < 120)) ||
                                            dis <= 5 * units::cm && angle_diff1 <= 30 && angle_diff3 < 60 ||
                                            dis <= 15 * units::cm && (angle_diff1 <= 20) && angle_diff3 < 40 ||
                                            (angle_diff1 < 10 && dis <= 60 * units::cm && angle_diff3 < 15))
                                            flag_merge = true;
                                    }
                                    else if (length_2 > 12 * units::cm && length_1 <= 12 * units::cm) {
                                        // one is short
                                        if ((dis <= 3 * units::cm) &&
                                            ((angle_diff1 <= 45 || angle_diff2 <= 45) && (angle_diff3 < 60) ||
                                             (flag_para && (angle_diff1 <= 90 || angle_diff2 <= 90) &&
                                              angle_diff3 < 120)) ||
                                            dis <= 5 * units::cm && angle_diff2 <= 30 && angle_diff3 < 60 ||
                                            dis <= 15 * units::cm && (angle_diff2 <= 20) && angle_diff3 < 40 ||
                                            (angle_diff2 < 10 && dis <= 60 * units::cm && angle_diff3 < 15))
                                            flag_merge = true;
                                    }
                                    else {
                                        // both are long
                                        if ((dis <= 3 * units::cm) &&
                                            ((angle_diff1 <= 45 || angle_diff2 <= 45) && (angle_diff3 < 60) ||
                                             (flag_para && (angle_diff1 <= 90 || angle_diff2 <= 90) &&
                                              angle_diff3 < 120)) ||
                                            dis <= 5 * units::cm && (angle_diff1 <= 30 || angle_diff2 <= 30) &&
                                            angle_diff3 < 45 ||
                                            (dis <= 15 * units::cm) && (angle_diff1 <= 20 || angle_diff2 <= 20) &&
                                            angle_diff3 < 30 ||
                                            (angle_diff1 < 10 || angle_diff2 < 10) && (dis <= 60 * units::cm) &&
                                            angle_diff3 < 15)
                                            flag_merge = true;
                                    }
                                }
#pragma GCC diagnostic pop                            

                            }


                            if (flag_merge) {
                                boost::add_edge(ilive2desc[map_cluster_index[cluster_1]],
                                                ilive2desc[map_cluster_index[cluster_2]], g);
                            }
                        } // if (tested_pairs....)
                    } // j
                } // i
            } //if(connected_live_clusters.size()>1)
        }

        // new function to merge clusters ...
        merge_clusters(g, live_grouping);


        // set cluster id ... 
        int cluster_id = 1;
        for (auto* cluster : live_grouping.children()) {
            cluster->set_cluster_id(cluster_id++);
        }

        // live_clusters = live_grouping.children(); // copy
        // // Process each cluster
        // for (size_t iclus = 0; iclus < live_clusters.size(); ++iclus) {
        //     Cluster* cluster = live_clusters.at(iclus);
        //     auto& scope = cluster->get_default_scope();
        //     std::cout << "Test: " << iclus << " " << cluster->nchildren() << " " << scope.pcname << " " << scope.coords[0] << " " << scope.coords[1] << " " << scope.coords[2] << " " << cluster->get_scope_filter(scope)<< " " << cluster->get_center() << std::endl;
        // }
    }
};
#pragma GCC diagnostic pop
