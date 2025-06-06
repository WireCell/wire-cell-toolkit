#include <WireCellClus/ClusteringFuncs.h>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

// The original developers do not care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;


static
std::map<int, Cluster*> Separate_overclustering(Cluster *cluster)
{
    // can follow ToyClustering_separate to add clusters ...
    auto* grouping = cluster->grouping();
    const auto &tp = grouping->get_params();

    // copy the create_graph from the PR3D Cluster ...

    // PR3DClusterSelection new_clusters;

    // cluster->Create_point_cloud();
    const int N = cluster->npoints();
    std::shared_ptr<MCUGraph> graph = std::make_shared<MCUGraph>(N);

    // ToyPointCloud *point_cloud = cluster->get_point_cloud();
    std::vector<Blob*> mcells = cluster->children();

    // plane -> point -> wire index
    const auto& winds = cluster->wire_indices();
    const Cluster::time_blob_map_t &time_cells_set_map = cluster->time_blob_map();

    // WCP::WCPointCloud<double> &cloud = point_cloud->get_cloud();
    // WCP::WC2DPointCloud<double> &cloud_u = point_cloud->get_cloud_u();
    // WCP::WC2DPointCloud<double> &cloud_v = point_cloud->get_cloud_v();
    // WCP::WC2DPointCloud<double> &cloud_w = point_cloud->get_cloud_w();

    // blob -> wind -> points
    // std::map<Blob *, std::map<int, std::set<int>>> map_mcell_uindex_wcps;
    // std::map<Blob *, std::map<int, std::set<int>>> map_mcell_vindex_wcps;
    // std::map<Blob *, std::map<int, std::set<int>>> map_mcell_windex_wcps;
    std::map<const Blob *, std::map<int, std::set<int>>> map_mcell_wind_wcps[3];

    for (auto it = mcells.begin(); it != mcells.end(); it++) {
        Blob *mcell = (*it);
        // std::map<int, std::set<int>> map_uindex_wcps;
        // std::map<int, std::set<int>> map_vindex_wcps;
        // std::map<int, std::set<int>> map_windex_wcps;
        std::map<int, std::set<int>> map_wind_wcps[3];
        const std::vector<int> &wcps = cluster->get_blob_indices(mcell);
        for (const int point_index : wcps) {
            auto v = vertex(point_index, *graph);  // retrieve vertex descriptor
            (*graph)[v].index = point_index;
            // if (map_windex_wcps.find(wcp.index_w) == map_windex_wcps.end()) {
            //     std::set<int> wcps;
            //     wcps.insert(wcp.index);
            //     map_windex_wcps[wcp.index_w] = wcps;
            // }
            // else {
            //     map_windex_wcps[wcp.index_w].insert(wcp.index);
            // }
            for (size_t plane_ind=0; plane_ind!=3; ++plane_ind) {
                const int wind = winds[plane_ind][point_index];
                if (map_wind_wcps[plane_ind].find(wind) == map_wind_wcps[plane_ind].end()) {
                    std::set<int> wcps;
                    wcps.insert(point_index);
                    map_wind_wcps[plane_ind][wind] = wcps;
                }
                else {
                    map_wind_wcps[plane_ind][wind].insert(point_index);
                }
            }
        }
        // map_mcell_uindex_wcps[mcell] = map_uindex_wcps;
        // map_mcell_vindex_wcps[mcell] = map_vindex_wcps;
        // map_mcell_windex_wcps[mcell] = map_windex_wcps;
        for (size_t plane_ind=0; plane_ind!=3; ++plane_ind) {
            map_mcell_wind_wcps[plane_ind][mcell] = map_wind_wcps[plane_ind];
        }
    }

    int num_edges = 0;

    // create graph for points inside the same mcell
    for (auto it = mcells.begin(); it != mcells.end(); it++) {
        Blob *mcell = (*it);
        // std::vector<int> &wcps = point_cloud->get_mcell_indices(mcell);
        const std::vector<int> &wcps = cluster->get_blob_indices(mcell);
        int max_wire_interval = mcell->get_max_wire_interval();
        int min_wire_interval = mcell->get_min_wire_interval();
        std::map<int, std::set<int>> *map_max_index_wcps;
        std::map<int, std::set<int>> *map_min_index_wcps;
        // if (mcell->get_max_wire_type() == WirePlaneType_t(0)) {
        //     map_max_index_wcps = &map_mcell_uindex_wcps[mcell];
        // }
        // else if (mcell->get_max_wire_type() == WirePlaneType_t(1)) {
        //     map_max_index_wcps = &map_mcell_wind_wcps[plane_ind][mcell];
        // }
        // else {
        //     map_max_index_wcps = &map_mcell_windex_wcps[mcell];
        // }
        // if (mcell->get_min_wire_type() == WirePlaneType_t(0)) {
        //     map_min_index_wcps = &map_mcell_uindex_wcps[mcell];
        // }
        // else if (mcell->get_min_wire_type() == WirePlaneType_t(1)) {
        //     map_min_index_wcps = &map_mcell_vindex_wcps[mcell];
        // }
        // else {
        //     map_min_index_wcps = &map_mcell_windex_wcps[mcell];
        // }
        const int max_wire_type = mcell->get_max_wire_type();
        const int min_wire_type = mcell->get_min_wire_type();
        map_max_index_wcps = &map_mcell_wind_wcps[max_wire_type][mcell];
        map_min_index_wcps = &map_mcell_wind_wcps[min_wire_type][mcell];

        for (const int index1 : wcps) {
            // WCPointCloud<double>::WCPoint &wcp1 = cloud.pts[*it1];
            // int index1 = wcp1.index;
            int index_max_wire = winds[max_wire_type][index1];
            int index_min_wire = winds[min_wire_type][index1];
            // if (mcell->get_max_wire_type() == WirePlaneType_t(0)) {
            //     index_max_wire = wcp1.index_u;
            // }
            // else if (mcell->get_max_wire_type() == WirePlaneType_t(1)) {
            //     index_max_wire = wcp1.index_v;
            // }
            // else {
            //     index_max_wire = wcp1.index_w;
            // }
            // if (mcell->get_min_wire_type() == WirePlaneType_t(0)) {
            //     index_min_wire = wcp1.index_u;
            // }
            // else if (mcell->get_min_wire_type() == WirePlaneType_t(1)) {
            //     index_min_wire = wcp1.index_v;
            // }
            // else {
            //     index_min_wire = wcp1.index_w;
            // }

            std::vector<std::set<int> *> max_wcps_set;
            std::vector<std::set<int> *> min_wcps_set;

            // go through the first map and find the ones satisfying the condition
            for (auto it2 = map_max_index_wcps->begin(); it2 != map_max_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_max_wire) <= max_wire_interval) {
                    max_wcps_set.push_back(&(it2->second));
                }
            }
            // go through the second map and find the ones satisfying the condition
            for (auto it2 = map_min_index_wcps->begin(); it2 != map_min_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_min_wire) <= min_wire_interval) {
                    min_wcps_set.push_back(&(it2->second));
                }
            }

            std::set<int> wcps_set1;
            std::set<int> wcps_set2;

            for (auto it2 = max_wcps_set.begin(); it2 != max_wcps_set.end(); it2++) {
                wcps_set1.insert((*it2)->begin(), (*it2)->end());
            }
            for (auto it3 = min_wcps_set.begin(); it3 != min_wcps_set.end(); it3++) {
                wcps_set2.insert((*it3)->begin(), (*it3)->end());
            }

            // std::cout << max_wcps_set.size() << " " << min_wcps_set.size() << std::endl;
            // for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
            //	for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
            {
                std::set<int> common_set;
                set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),
                                 std::inserter(common_set, common_set.begin()));

                //	std::cout << "S0: " << common_set.size() << std::endl;

                for (const int index2 : common_set) {
                    // WCPointCloud<double>::WCPoint &wcp2 = cloud.pts[*it4];
                    if (index2 != index1) {
                        // int index2 = wcp2.index;
                        // std::cout << index1 << " " << index2 << std::endl;
                        // add edge ...
                        //auto edge = add_edge(index1, index2, *graph);
                        const geo_point_t wcp1 = cluster->point3d(index1);
                        const geo_point_t wcp2 = cluster->point3d(index2);
                        double dis = sqrt(pow(wcp1.x() - wcp2.x(), 2) + pow(wcp1.y() - wcp2.y(), 2) + pow(wcp1.z() - wcp2.z(), 2));
                        // if (edge.second) {
                            
                        //     (*graph)[edge.first].dist = dis;
                        //     num_edges++;
                        //     // std::cout << wcp1.x << " " << wcp1.y << " " << wcp1.z << " " << wcp1.index_u << " " <<
                        //     // wcp1.index_v << " " << wcp1.index_w << " " << wcp2.index_u << " " << wcp2.index_v << " "
                        //     // << wcp2.index_w << std::endl;
                        // }
                        auto edge = add_edge(index1, index2, WireCell::PointCloud::Facade::EdgeProp(dis),*graph);
                        if (edge.second) {
                            num_edges++;
                        }
                    }
                }
                //}
            }
        }
    }

    //  std::cout << "Xin: " << num_edges << " " << N << std::endl;

    std::vector<int> time_slices;
    for (auto it1 = time_cells_set_map.begin(); it1 != time_cells_set_map.end(); it1++) {
        time_slices.push_back((*it1).first);
    }

    std::vector<std::pair<const Blob *, const Blob *>> connected_mcells;

    for (size_t i = 0; i != time_slices.size(); i++) {
        const std::set<const Blob*, blob_less_functor> &mcells_set = time_cells_set_map.at(time_slices.at(i));

        // create graph for points in mcell inside the same time slice
        if (mcells_set.size() >= 2) {
            for (auto it2 = mcells_set.begin(); it2 != mcells_set.end(); it2++) {
                const Blob *mcell1 = *it2;
                auto it2p = it2;
                if (it2p != mcells_set.end()) {
                    it2p++;
                    for (auto it3 = it2p; it3 != mcells_set.end(); it3++) {
                        const Blob *mcell2 = *(it3);
                        // std::cout << mcell1 << " " << mcell2 << " " << mcell1->Overlap_fast(mcell2,2) << std::endl;
                        if (mcell1->overlap_fast(*mcell2, 2)) connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                    }
                }
            }
        }
        // create graph for points between connected mcells in adjacent time slices + 1, if not, + 2
        std::vector<std::set<const Blob*, blob_less_functor>> vec_mcells_set;
        if (i + 1 < time_slices.size()) {
            if (time_slices.at(i + 1) - time_slices.at(i) == 1*tp.nticks_live_slice) {
                vec_mcells_set.push_back(time_cells_set_map.at(time_slices.at(i + 1)));
                if (i + 2 < time_slices.size())
                    if (time_slices.at(i + 2) - time_slices.at(i) == 2*tp.nticks_live_slice)
                        vec_mcells_set.push_back(time_cells_set_map.at(time_slices.at(i + 2)));
            }
            else if (time_slices.at(i + 1) - time_slices.at(i) == 2*tp.nticks_live_slice) {
                vec_mcells_set.push_back(time_cells_set_map.at(time_slices.at(i + 1)));
            }
        }
        //    bool flag = false;
        for (size_t j = 0; j != vec_mcells_set.size(); j++) {
            //      if (flag) break;
            std::set<const Blob*, blob_less_functor> &next_mcells_set = vec_mcells_set.at(j);
            for (auto it1 = mcells_set.begin(); it1 != mcells_set.end(); it1++) {
                const Blob *mcell1 = (*it1);
                for (auto it2 = next_mcells_set.begin(); it2 != next_mcells_set.end(); it2++) {
                    const Blob *mcell2 = (*it2);
                    if (mcell1->overlap_fast(*mcell2, 2)) {
                        //	    flag = true; // correct???
                        connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                    }
                }
            }
        }
    }

    // establish edge ...
    const size_t max_num_nodes = 5;
    std::map<std::pair<int, int>, std::set<std::pair<double, int>>> closest_index;

    // std::cout << connected_mcells.size() << std::endl;
    for (auto it = connected_mcells.begin(); it != connected_mcells.end(); it++) {
        const Blob *mcell1 = (*it).first;
        const Blob *mcell2 = (*it).second;

        const std::vector<int>& wcps1 = cluster->get_blob_indices(mcell1);
        const std::vector<int>& wcps2 = cluster->get_blob_indices(mcell2);

        // test 2 against 1 ...
        int max_wire_interval = mcell1->get_max_wire_interval();
        int min_wire_interval = mcell1->get_min_wire_interval();
        std::map<int, std::set<int>> *map_max_index_wcps;
        std::map<int, std::set<int>> *map_min_index_wcps;

        // if (mcell1->get_max_wire_type() == WirePlaneType_t(0)) {
        //     map_max_index_wcps = &map_mcell_uindex_wcps[mcell2];
        // }
        // else if (mcell1->get_max_wire_type() == WirePlaneType_t(1)) {
        //     map_max_index_wcps = &map_mcell_vindex_wcps[mcell2];
        // }
        // else {
        //     map_max_index_wcps = &map_mcell_windex_wcps[mcell2];
        // }
        // if (mcell1->get_min_wire_type() == WirePlaneType_t(0)) {
        //     map_min_index_wcps = &map_mcell_uindex_wcps[mcell2];
        // }
        // else if (mcell1->get_min_wire_type() == WirePlaneType_t(1)) {
        //     map_min_index_wcps = &map_mcell_vindex_wcps[mcell2];
        // }
        // else {
        //     map_min_index_wcps = &map_mcell_windex_wcps[mcell2];
        // }
        map_max_index_wcps = &map_mcell_wind_wcps[mcell1->get_max_wire_type()][mcell2];
        map_min_index_wcps = &map_mcell_wind_wcps[mcell1->get_min_wire_type()][mcell2];

        for (const int index1 : wcps1) {
            // WCPointCloud<double>::WCPoint &wcp1 = cloud.pts[*it1];
            // int index1 = wcp1.index;
            int index_max_wire = winds[mcell1->get_max_wire_type()][index1];
            int index_min_wire = winds[mcell1->get_min_wire_type()][index1];
            // if (mcell1->get_max_wire_type() == WirePlaneType_t(0)) {
            //     index_max_wire = wcp1.index_u;
            // }
            // else if (mcell1->get_max_wire_type() == WirePlaneType_t(1)) {
            //     index_max_wire = wcp1.index_v;
            // }
            // else {
            //     index_max_wire = wcp1.index_w;
            // }
            // if (mcell1->get_min_wire_type() == WirePlaneType_t(0)) {
            //     index_min_wire = wcp1.index_u;
            // }
            // else if (mcell1->get_min_wire_type() == WirePlaneType_t(1)) {
            //     index_min_wire = wcp1.index_v;
            // }
            // else {
            //     index_min_wire = wcp1.index_w;
            // }
            std::vector<std::set<int> *> max_wcps_set;
            std::vector<std::set<int> *> min_wcps_set;
            // go through the first map and find the ones satisfying the condition
            for (auto it2 = map_max_index_wcps->begin(); it2 != map_max_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_max_wire) <= max_wire_interval) {
                    max_wcps_set.push_back(&(it2->second));
                }
            }
            // go through the second map and find the ones satisfying the condition
            for (auto it2 = map_min_index_wcps->begin(); it2 != map_min_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_min_wire) <= min_wire_interval) {
                    min_wcps_set.push_back(&(it2->second));
                }
            }

            std::set<int> wcps_set1;
            std::set<int> wcps_set2;

            for (auto it2 = max_wcps_set.begin(); it2 != max_wcps_set.end(); it2++) {
                wcps_set1.insert((*it2)->begin(), (*it2)->end());
            }
            for (auto it3 = min_wcps_set.begin(); it3 != min_wcps_set.end(); it3++) {
                wcps_set2.insert((*it3)->begin(), (*it3)->end());
            }

            //   for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
            //	for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
            {
                std::set<int> common_set;
                set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),
                                 std::inserter(common_set, common_set.begin()));

                //	std::cout << "S1: " << common_set.size() << std::endl;
                //	  std::cout << common_set.size() << std::endl;

                //	std::map<int,std::pair<int,double> > closest_index;

                // for (auto it4 = common_set.begin(); it4 != common_set.end(); it4++) {
                for (const int index2 : common_set) {
                    // WCPointCloud<double>::WCPoint &wcp2 = cloud.pts[*it4];
                    if (index2 != index1) {
                        // int index2 = wcp2.index;
                        const geo_point_t wcp1 = cluster->point3d(index1);
                        const geo_point_t wcp2 = cluster->point3d(index2);
                        double dis = sqrt(pow(wcp1.x() - wcp2.x(), 2) + pow(wcp1.y() - wcp2.y(), 2) + pow(wcp1.z() - wcp2.z(), 2));
                        const int time2 = cluster->blob_with_point(index2)->slice_index_min();
                        auto key = std::make_pair(index1, time2);

                        if (closest_index.find(key) == closest_index.end()) {
                            std::set<std::pair<double, int> > temp_sets;
                            temp_sets.insert(std::make_pair(dis,index2));
                            closest_index[key] = temp_sets;
                        }
                        else {
                            closest_index[key].insert(std::make_pair(dis,index2));
                            if (closest_index[key].size()>max_num_nodes){
                                auto it5 = closest_index[key].begin();
                                for (int qx = 0; qx!=max_num_nodes;qx++){
                                    it5++;
                                }
                                closest_index[key].erase(it5,closest_index[key].end());
                            }
                            //if (dis < closest_index[key].second || (std::abs(dis - closest_index[key].second) < 1e-10 && pind2 < closest_index[key].first)) closest_index[key] = std::make_pair(pind2, dis);
                        }

                        // if (closest_index.find(std::make_pair(index1, time2)) ==
                        //     closest_index.end()) {
                        //     closest_index[std::make_pair(index1, time2)] =
                        //         std::make_pair(index2, dis);
                        // }
                        // else {
                        //     if (dis < closest_index[std::make_pair(index1, time2)].second)
                        //         closest_index[std::make_pair(index1, time2)] =
                        //             std::make_pair(index2, dis);
                        // }
                    }
                }

                //	std::cout << closest_index.size() << std::endl;
                // for (auto it4 = closest_index.begin(); it4!=closest_index.end(); it4++){
                //   int index2 = it4->second.first;
                //   double dis = it4->second.second;
                //   auto edge = add_edge(index1,index2,*graph);
                //   if (edge.second){
                //     (*graph)[edge.first].dist = dis;
                //     num_edges ++;
                //   }
                // }
            }
            //}
        }

        // test 1 against 2 ...
        max_wire_interval = mcell2->get_max_wire_interval();
        min_wire_interval = mcell2->get_min_wire_interval();
        // if (mcell2->get_max_wire_type() == WirePlaneType_t(0)) {
        //     map_max_index_wcps = &map_mcell_uindex_wcps[mcell1];
        // }
        // else if (mcell2->get_max_wire_type() == WirePlaneType_t(1)) {
        //     map_max_index_wcps = &map_mcell_vindex_wcps[mcell1];
        // }
        // else {
        //     map_max_index_wcps = &map_mcell_windex_wcps[mcell1];
        // }
        // if (mcell2->get_min_wire_type() == WirePlaneType_t(0)) {
        //     map_min_index_wcps = &map_mcell_uindex_wcps[mcell1];
        // }
        // else if (mcell2->get_min_wire_type() == WirePlaneType_t(1)) {
        //     map_min_index_wcps = &map_mcell_vindex_wcps[mcell1];
        // }
        // else {
        //     map_min_index_wcps = &map_mcell_windex_wcps[mcell1];
        // }
        map_max_index_wcps = &map_mcell_wind_wcps[mcell2->get_max_wire_type()][mcell1];
        map_min_index_wcps = &map_mcell_wind_wcps[mcell2->get_min_wire_type()][mcell1];

        // for (auto it1 = wcps2.begin(); it1 != wcps2.end(); it1++) {
        for (const int index1 : wcps2) {
            // WCPointCloud<double>::WCPoint &wcp1 = cloud.pts[*it1];
            // int index1 = wcp1.index;
            int index_max_wire = winds[mcell2->get_max_wire_type()][index1];
            int index_min_wire = winds[mcell2->get_min_wire_type()][index1];
            // if (mcell2->get_max_wire_type() == WirePlaneType_t(0)) {
            //     index_max_wire = wcp1.index_u;
            // }
            // else if (mcell2->get_max_wire_type() == WirePlaneType_t(1)) {
            //     index_max_wire = wcp1.index_v;
            // }
            // else {
            //     index_max_wire = wcp1.index_w;
            // }
            // if (mcell2->get_min_wire_type() == WirePlaneType_t(0)) {
            //     index_min_wire = wcp1.index_u;
            // }
            // else if (mcell2->get_min_wire_type() == WirePlaneType_t(1)) {
            //     index_min_wire = wcp1.index_v;
            // }
            // else {
            //     index_min_wire = wcp1.index_w;
            // }
            std::vector<std::set<int> *> max_wcps_set;
            std::vector<std::set<int> *> min_wcps_set;
            // go through the first map and find the ones satisfying the condition
            for (auto it2 = map_max_index_wcps->begin(); it2 != map_max_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_max_wire) <= max_wire_interval) {
                    max_wcps_set.push_back(&(it2->second));
                }
            }
            // go through the second map and find the ones satisfying the condition
            for (auto it2 = map_min_index_wcps->begin(); it2 != map_min_index_wcps->end(); it2++) {
                if (fabs(it2->first - index_min_wire) <= min_wire_interval) {
                    min_wcps_set.push_back(&(it2->second));
                }
            }

            std::set<int> wcps_set1;
            std::set<int> wcps_set2;

            for (auto it2 = max_wcps_set.begin(); it2 != max_wcps_set.end(); it2++) {
                wcps_set1.insert((*it2)->begin(), (*it2)->end());
            }
            for (auto it3 = min_wcps_set.begin(); it3 != min_wcps_set.end(); it3++) {
                wcps_set2.insert((*it3)->begin(), (*it3)->end());
            }

            // for (auto it2 = max_wcps_set.begin(); it2!=max_wcps_set.end(); it2++){
            // 	for (auto it3 = min_wcps_set.begin(); it3!=min_wcps_set.end(); it3++){
            {
                std::set<int> common_set;
                set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),
                                 std::inserter(common_set, common_set.begin()));

                //	std::cout << "S2: " << common_set.size() << std::endl;

                //	std::map<int,std::pair<int,double> > closest_index;

                // for (auto it4 = common_set.begin(); it4 != common_set.end(); it4++) {
                for (const int index2 : common_set) {
                    // WCPointCloud<double>::WCPoint &wcp2 = cloud.pts[*it4];
                    if (index2 != index1) {
                        // int index2 = wcp2.index;
                        const geo_point_t wcp1 = cluster->point3d(index1);
                        const geo_point_t wcp2 = cluster->point3d(index2);
                        double dis = sqrt(pow(wcp1.x() - wcp2.x(), 2) + pow(wcp1.y() - wcp2.y(), 2) + pow(wcp1.z() - wcp2.z(), 2));

                        const int time2 = cluster->blob_with_point(index2)->slice_index_min();
                        auto key = std::make_pair(index1, time2);

                        if (closest_index.find(key) == closest_index.end()) {
                            std::set<std::pair<double, int> > temp_sets;
                            temp_sets.insert(std::make_pair(dis,index2));
                            closest_index[key] = temp_sets;
                        }
                        else {
                            closest_index[key].insert(std::make_pair(dis,index2));
                            if (closest_index[key].size()>max_num_nodes){
                                auto it5 = closest_index[key].begin();
                                for (int qx = 0; qx!=max_num_nodes;qx++){
                                    it5++;
                                }
                                closest_index[key].erase(it5,closest_index[key].end());
                            }
                            //if (dis < closest_index[key].second || (std::abs(dis - closest_index[key].second) < 1e-10 && pind2 < closest_index[key].first)) closest_index[key] = std::make_pair(pind2, dis);
                        }

                        // if (closest_index.find(std::make_pair(index1, time2)) ==
                        //     closest_index.end()) {
                        //     closest_index[std::make_pair(index1, time2)] =
                        //         std::make_pair(index2, dis);
                        // }
                        // else {
                        //     if (dis < closest_index[std::make_pair(index1, time2)].second)
                        //         closest_index[std::make_pair(index1, time2)] =
                        //             std::make_pair(index2, dis);
                        // }
                    }
                }

                // std::cout << closest_index.size() << std::endl;
                //  for (auto it4 = closest_index.begin(); it4!=closest_index.end(); it4++){
                //    int index2 = it4->second.first;
                //    double dis = it4->second.second;
                //    auto edge = add_edge(index1,index2,*graph);
                //    if (edge.second){
                //      (*graph)[edge.first].dist = dis;
                //      num_edges ++;
                //    }
                //  }

                // for (auto it4 = common_set.begin(); it4!=common_set.end(); it4++){
                //   WCPointCloud<double>::WCPoint& wcp2 = cloud.pts[*it4];
                //   if (wcp2.index != wcp1.index){
                //     int index2 = wcp2.index;
                //     auto edge = add_edge(index1,index2,*graph);
                //     if (edge.second){
                //       (*graph)[edge.first].dist =
                //       sqrt(pow(wcp1.x-wcp2.x,2)+pow(wcp1.y-wcp2.y,2)+pow(wcp1.z-wcp2.z,2)); num_edges ++;
                //     }
                //   }
                // }
            }
            //      }
        }
    }

    for (auto it4 = closest_index.begin(); it4 != closest_index.end(); it4++) {
        int index1 = it4->first.first;
        for (auto it5 = it4->second.begin(); it5!=it4->second.end(); it5++){
            int index2 = (*it5).second;
            double dis = (*it5).first;
            auto edge = add_edge(index1,index2,WireCell::PointCloud::Facade::EdgeProp(dis),*graph);
            if (edge.second){
                //      (*graph)[edge.first].dist = dis;
                num_edges ++;
            }
            // protect against dead cells ...
            //std::cout << dis/units::cm << std::endl;
            if (it5 == it4->second.begin() && dis > 0.25*units::cm)
                break;
        }

        // int index2 = it4->second.first;
        // double dis = it4->second.second;
        // // auto edge = add_edge(index1, index2, *graph);
        // // if (edge.second) {
        // //     (*graph)[edge.first].dist = dis;
        // //     num_edges++;
        // // }
        // auto edge = add_edge(index1, index2, WireCell::PointCloud::Facade::EdgeProp(dis),*graph);
        // if (edge.second) {
        //     num_edges++;
        // }
    }
    // end of copying ...

    // now form the connected components, point -> component
    std::vector<int> component(num_vertices(*graph));
    const int num = connected_components(*graph, &component[0]);

    // Create ordered components
    std::vector<ComponentInfo> ordered_components;
    ordered_components.reserve(component.size());
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components.emplace_back(i);
    }

    // Assign vertices to components
    for (size_t i = 0; i < component.size(); ++i) {
        ordered_components[component[i]].add_vertex(i);
    }

    // Sort components by minimum vertex index
    std::sort(ordered_components.begin(), ordered_components.end(), 
        [](const ComponentInfo& a, const ComponentInfo& b) {
            return a.min_vertex < b.min_vertex;
        });

    if (num <= 1) return {};



    // if (num > 1) {
        // For each component, create a point cloud
    std::vector<std::shared_ptr<Simple3DPointCloud>> pt_clouds;
    std::vector<std::vector<size_t>> pt_clouds_global_indices;
    for (const auto& comp : ordered_components) {
        auto pt_cloud = std::make_shared<Simple3DPointCloud>();
        std::vector<size_t> global_indices;
        
        for (size_t vertex_idx : comp.vertex_indices) {
            geo_point_t pt = cluster->point3d(vertex_idx);
            pt_cloud->add({pt.x(), pt.y(), pt.z()});
            global_indices.push_back(vertex_idx);
        }
        
        pt_clouds.push_back(pt_cloud);
        pt_clouds_global_indices.push_back(global_indices);
    }

        // for (int j = 0; j != num; j++) {
        //     pt_clouds.push_back(std::make_shared<Simple3DPointCloud>());
        // }

        // // std::vector<int>::size_type i;
        // // for (i = 0; i != component.size(); ++i) {
        // //     pt_clouds.at(component[i])->AddPoint(cloud.pts[i], cloud_u.pts[i], cloud_v.pts[i], cloud_w.pts[i]);
        // //     //   std::cout << "Vertex " << i << " " << cloud.pts[i].x << " " << cloud.pts[i].y << " " << cloud.pts[i].z
        // //     //   << " " << cloud.pts[i].index_u << " " << cloud.pts[i].index_v << " " << cloud.pts[i].index_w << " " <<
        // //     //   cloud.pts[i].mcell << " " << cloud.pts[i].mcell->GetTimeSlice()  << " is in component " << component[i]
        // //     //   << std::endl;
        // // }
        // // for (int j = 0; j != num; j++) {
        // //     pt_clouds.at(j)->build_kdtree_index();
        // // }
        // std::vector<int>::size_type i;
        // for (i = 0; i != component.size(); ++i) {
        //     geo_point_t pt = cluster->point3d(i);
        //     pt_clouds.at(component[i])->add({pt.x(), pt.y(), pt.z()});
        //     pt_clouds_global_indices.at(component[i]).push_back(i);
        // }

        std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis(
            num, std::vector<std::tuple<int, int, double>>(num));
        std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_mst(
            num, std::vector<std::tuple<int, int, double>>(num));

        std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir1(
            num, std::vector<std::tuple<int, int, double>>(num));
        std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir2(
            num, std::vector<std::tuple<int, int, double>>(num));
        std::vector<std::vector<std::tuple<int, int, double>>> index_index_dis_dir_mst(
            num, std::vector<std::tuple<int, int, double>>(num));

        for (int j = 0; j != num; j++) {
            for (int k = 0; k != num; k++) {
                index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
                index_index_dis_mst[j][k] = std::make_tuple(-1, -1, 1e9);
                index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                index_index_dis_dir_mst[j][k] = std::make_tuple(-1, -1, 1e9);
            }
        }

        // check against the closest distance ...
        // no need to have MST ...
        for (int j = 0; j != num; j++) {
            for (int k = j + 1; k != num; k++) {
                index_index_dis[j][k] = pt_clouds.at(j)->get_closest_points(*pt_clouds.at(k));

                if (num < 100 && pt_clouds.at(j)->get_num_points() > 100 && pt_clouds.at(k)->get_num_points() > 100 &&
                        (pt_clouds.at(j)->get_num_points() + pt_clouds.at(k)->get_num_points()) > 400 ||
                    pt_clouds.at(j)->get_num_points() > 500 && pt_clouds.at(k)->get_num_points() > 500) {
                    // WCPointCloud<double>::WCPoint wp1 = cloud.pts.at(std::get<0>(index_index_dis[j][k]));
                    // WCPointCloud<double>::WCPoint wp2 = cloud.pts.at(std::get<1>(index_index_dis[j][k]));
                    // Point p1(wp1.x, wp1.y, wp1.z);
                    // Point p2(wp2.x, wp2.y, wp2.z);
                    geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                    geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));

                    // TVector3 dir1 = cluster->VHoughTrans(p1, 30 * units::cm, pt_clouds.at(j));
                    // TVector3 dir2 = cluster->VHoughTrans(p2, 30 * units::cm, pt_clouds.at(k));
                    // dir1 *= -1;
                    // dir2 *= -1;

                    geo_point_t dir1 = cluster->vhough_transform(p1, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(j), pt_clouds_global_indices.at(j));
                    geo_point_t dir2 = cluster->vhough_transform(p2, 30 * units::cm, Cluster::HoughParamSpace::theta_phi, pt_clouds.at(k), pt_clouds_global_indices.at(k));
                    dir1 = dir1 * -1;
                    dir2 = dir2 * -1;

                    std::pair<int, double> result1 = pt_clouds.at(k)->get_closest_point_along_vec(
                        p1, dir1, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                    if (result1.first >= 0) {
                        index_index_dis_dir1[j][k] =
                            std::make_tuple(std::get<0>(index_index_dis[j][k]), result1.first, result1.second);
                    }

                    std::pair<int, double> result2 = pt_clouds.at(j)->get_closest_point_along_vec(
                        p2, dir2, 80 * units::cm, 5 * units::cm, 7.5, 3 * units::cm);

                    if (result2.first >= 0) {
                        index_index_dis_dir2[j][k] =
                            std::make_tuple(result2.first, std::get<1>(index_index_dis[j][k]), result2.second);
                    }
                }

                // Now check the path ...
                {
                    geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis[j][k]));
                    geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis[j][k]));
                    double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));

                    double step_dis = 1.0 * units::cm;
                    int num_steps = dis / step_dis + 1;
                    int num_bad = 0;
                    geo_point_t test_p;
                    for (int ii = 0; ii != num_steps; ii++) {
                        test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                                   p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                                   p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                        if (true) {
                            /// FIXME: how to add face information?
                            const bool good_point = cluster->grouping()->is_good_point(test_p, tp.face);
                            if (!good_point) num_bad++;
                        }
                    }

                    if (num_bad > 7 || num_bad > 2 && num_bad >= 0.75 * num_steps) {
                        index_index_dis[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }

                // Now check the path ...
                if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                    geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir1[j][k]));
                    geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir1[j][k]));

                    double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                    double step_dis = 1.0 * units::cm;
                    int num_steps = dis / step_dis + 1;
                    int num_bad = 0;
                    geo_point_t test_p;
                    for (int ii = 0; ii != num_steps; ii++) {
                        test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                                   p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                                   p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                        // if (!ct_point_cloud.is_good_point(test_p)) num_bad++;
                        if (true) {
                            /// FIXME: how to add face information?
                            const bool good_point = cluster->grouping()->is_good_point(test_p, tp.face);
                            if (!good_point) num_bad++;
                        }
                    }

                    if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                        index_index_dis_dir1[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }

                // Now check the path ...
                if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                    geo_point_t p1 = pt_clouds.at(j)->point(std::get<0>(index_index_dis_dir2[j][k]));
                    geo_point_t p2 = pt_clouds.at(k)->point(std::get<1>(index_index_dis_dir2[j][k]));

                    double dis = sqrt(pow(p1.x() - p2.x(), 2) + pow(p1.y() - p2.y(), 2) + pow(p1.z() - p2.z(), 2));
                    double step_dis = 1.0 * units::cm;
                    int num_steps = dis / step_dis + 1;
                    int num_bad = 0;
                    geo_point_t test_p;
                    for (int ii = 0; ii != num_steps; ii++) {
                        test_p.set(p1.x() + (p2.x() - p1.x()) / num_steps * (ii + 1),
                                   p1.y() + (p2.y() - p1.y()) / num_steps * (ii + 1),
                                   p1.z() + (p2.z() - p1.z()) / num_steps * (ii + 1));
                        // if (!ct_point_cloud.is_good_point(test_p)) num_bad++;
                        if (true) {
                            /// FIXME: how to add face information?
                            const bool good_point = cluster->grouping()->is_good_point(test_p, tp.face);
                            if (!good_point) num_bad++;
                        }
                    }

                    if (num_bad > 7 || (num_bad > 2 && num_bad >= 0.75 * num_steps)) {
                        index_index_dis_dir2[j][k] = std::make_tuple(-1, -1, 1e9);
                    }
                }
            }
        }

        // deal with MST
        {
            const int N = num;
            boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property,
                                  boost::property<boost::edge_weight_t, double>>
                temp_graph(N);

            for (int j = 0; j != num; j++) {
                for (int k = j + 1; k != num; k++) {
                    int index1 = j;
                    int index2 = k;
                    if (std::get<0>(index_index_dis[j][k]) >= 0)
                        /*auto edge =*/ add_edge(index1, index2, std::get<2>(index_index_dis[j][k]), temp_graph);
                }
            }

             // Process MST
            process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_mst);

            // {
            //     std::vector<int> possible_root_vertex;
            //     std::vector<int> component(num_vertices(temp_graph));
            //     const int num1 = connected_components(temp_graph, &component[0]);
            //     possible_root_vertex.resize(num1);
            //     std::vector<int>::size_type i;
            //     for (i = 0; i != component.size(); ++i) {
            //         possible_root_vertex.at(component[i]) = i;
            //     }

            //     for (size_t i = 0; i != possible_root_vertex.size(); i++) {
            //         std::vector<boost::graph_traits<MCUGraph>::vertex_descriptor> predecessors(
            //             num_vertices(temp_graph));

            //         prim_minimum_spanning_tree(temp_graph, &predecessors[0],
            //                                    boost::root_vertex(possible_root_vertex.at(i)));

            //         for (size_t j = 0; j != predecessors.size(); ++j) {
            //             if (predecessors[j] != j) {
            //                 if (j < predecessors[j]) {
            //                     index_index_dis_mst[j][predecessors[j]] = index_index_dis[j][predecessors[j]];
            //                 }
            //                 else {
            //                     index_index_dis_mst[predecessors[j]][j] = index_index_dis[predecessors[j]][j];
            //                 }
            //                 // std::cout << j << " " << predecessors[j] << " " << std::endl;
            //             }
            //             else {
            //                 // std::cout << j << " " << std::endl;
            //             }
            //         }
            //     }
            // }
        }

        // deal with MST for directionality
        {
            const int N = num;
            boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, boost::no_property,
                                  boost::property<boost::edge_weight_t, double>>
                temp_graph(N);

            for (int j = 0; j != num; j++) {
                for (int k = j + 1; k != num; k++) {
                    int index1 = j;
                    int index2 = k;
                    if (std::get<0>(index_index_dis_dir1[j][k]) >= 0 || std::get<0>(index_index_dis_dir2[j][k]) >= 0)
                        /*auto edge =*/ add_edge(
                            index1, index2,
                            std::min(std::get<2>(index_index_dis_dir1[j][k]), std::get<2>(index_index_dis_dir2[j][k])),
                            temp_graph);
                }
            }

            process_mst_deterministically(temp_graph, index_index_dis, index_index_dis_dir_mst);
            // {
            //     std::vector<int> possible_root_vertex;
            //     std::vector<int> component(num_vertices(temp_graph));
            //     const int num1 = connected_components(temp_graph, &component[0]);
            //     possible_root_vertex.resize(num1);
            //     std::vector<int>::size_type i;
            //     for (i = 0; i != component.size(); ++i) {
            //         possible_root_vertex.at(component[i]) = i;
            //     }

            //     for (size_t i = 0; i != possible_root_vertex.size(); i++) {
            //         std::vector<boost::graph_traits<MCUGraph>::vertex_descriptor> predecessors(
            //             num_vertices(temp_graph));
            //         prim_minimum_spanning_tree(temp_graph, &predecessors[0],
            //                                    boost::root_vertex(possible_root_vertex.at(i)));
            //         for (size_t j = 0; j != predecessors.size(); ++j) {
            //             if (predecessors[j] != j) {
            //                 if (j < predecessors[j]) {
            //                     index_index_dis_dir_mst[j][predecessors[j]] = index_index_dis[j][predecessors[j]];
            //                 }
            //                 else {
            //                     index_index_dis_dir_mst[predecessors[j]][j] = index_index_dis[predecessors[j]][j];
            //                 }
            //                 // std::cout << j << " " << predecessors[j] << " " << std::endl;
            //             }
            //             else {
            //                 // std::cout << j << " " << std::endl;
            //             }
            //         }
            //     }
            // }
        }

        for (int j = 0; j != num; j++) {
            for (int k = j + 1; k != num; k++) {
                if (std::get<2>(index_index_dis[j][k]) < 3 * units::cm) {
                    index_index_dis_mst[j][k] = index_index_dis[j][k];
                }

                // establish the path ...
                if (std::get<0>(index_index_dis_mst[j][k]) >= 0) {
                    // auto edge = add_edge(std::get<0>(index_index_dis_mst[j][k]), std::get<1>(index_index_dis_mst[j][k]),
                    //                      *graph);
                    // if (edge.second) {
                    const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_mst[j][k]));
                    const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_mst[j][k]));
                    float dis;
                    if (std::get<2>(index_index_dis_mst[j][k]) > 5 * units::cm) {
                        dis = std::get<2>(index_index_dis_mst[j][k]);
                    }
                    else {
                        dis = std::get<2>(index_index_dis_mst[j][k]);
                    }
                    // }

                    /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis),*graph);
                }

                if (std::get<0>(index_index_dis_dir_mst[j][k]) >= 0) {
                    if (std::get<0>(index_index_dis_dir1[j][k]) >= 0) {
                        // auto edge = add_edge(std::get<0>(index_index_dis_dir1[j][k]),
                        //                      std::get<1>(index_index_dis_dir1[j][k]), *graph);
                        // if (edge.second) {
                        const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir1[j][k]));
                        const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir1[j][k]));
                        float dis;
                        if (std::get<2>(index_index_dis_dir1[j][k]) > 5 * units::cm) {
                            dis = std::get<2>(index_index_dis_dir1[j][k]);
                        }
                        else {
                            dis = std::get<2>(index_index_dis_dir1[j][k]);
                        }
                        // }
                        /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis),*graph);
                    }
                    if (std::get<0>(index_index_dis_dir2[j][k]) >= 0) {
                        // auto edge = add_edge(std::get<0>(index_index_dis_dir2[j][k]),
                        //                      std::get<1>(index_index_dis_dir2[j][k]), *graph);
                        // if (edge.second) {
                        const int gind1 = pt_clouds_global_indices.at(j).at(std::get<0>(index_index_dis_dir2[j][k]));
                        const int gind2 = pt_clouds_global_indices.at(k).at(std::get<1>(index_index_dis_dir2[j][k]));
                        float dis;
                        if (std::get<2>(index_index_dis_dir2[j][k]) > 5 * units::cm) {
                            dis = std::get<2>(index_index_dis_dir2[j][k]);
                        }
                        else {
                            dis = std::get<2>(index_index_dis_dir2[j][k]);
                        }
                        // }
                        /*auto edge =*/ add_edge(gind1, gind2, WireCell::PointCloud::Facade::EdgeProp(dis), *graph);
                    }
                }
                // end check ...
            }
        }

        // std::cout << "Check: " << cluster->nchildren() << " " << num << std::endl;
        // study the independent component again ...
        {
            // point -> component
            std::vector<int> component1(num_vertices(*graph));
            const int num1 = connected_components(*graph, &component1[0]);

            if (num1 > 1) {

                // std::cout << "Check1: " << cluster->nchildren() << " " << num << " " << num1 << std::endl;

                // form new clusters ...
                // WCP logic ...
                // std::vector<std::set<const Blob*>> vec_mcells_set;
                // for (int ii = 0; ii != num1; ii++) {
                //     std::set<const Blob*> mcells_set;
                //     vec_mcells_set.push_back(mcells_set);
                // }
                // std::vector<int>::size_type i;
                // for (i = 0; i != component1.size(); ++i) {
                //     vec_mcells_set.at(component1[i]).insert(cloud.pts[i].mcell);
                // }
                // for (int ii = 0; ii != num1; ii++) {
                //     Cluster *cluster1 = new Cluster(1);
                //     for (auto it = vec_mcells_set.at(ii).begin(); it != vec_mcells_set.at(ii).end(); it++) {
                //         Blob *mcell = (*it);
                //         cluster1->AddCell(mcell, mcell->GetTimeSlice());
                //     }
                //     new_clusters.push_back(cluster1);
                // }
                std::vector<int> b2groupid(cluster->nchildren(), -1);
                std::vector<int>::size_type i;
                for (i = 0; i != component1.size(); ++i) {
                    const int bind = cluster->kd3d().major_index(i);
                    b2groupid.at(bind) = component1[i];
                }
                return grouping->separate(cluster, b2groupid, true); 
            }
        }

        // for (int i = 0; i != num; i++) {
        //     delete pt_clouds.at(i);
        // }

    // }  // if (num > 1)

    // delete graph;
    // return new_clusters;
    return {};
}

void WireCell::PointCloud::Facade::clustering_protect_overclustering(Grouping& live_grouping)
{
    std::vector<Cluster *> live_clusters = live_grouping.children();  // copy
    // sort the clusters by length using a lambda function
    // std::sort(live_clusters.begin(), live_clusters.end(), [](const Cluster *cluster1, const Cluster *cluster2) {
    //     return cluster1->get_length() > cluster2->get_length();
    // });
    for (size_t i = 0; i != live_clusters.size(); i++) {
        Cluster *cluster = live_clusters.at(i);
        // std::cout << "Cluster: " << i << " " << cluster->npoints() << std::endl;
        Separate_overclustering(cluster);
    }
}
