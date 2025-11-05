#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Grouping.h"

#include "connect_graphs.h"

using namespace WireCell;
using namespace WireCell::Clus;


void Graphs::connect_graph_closely(const Facade::Cluster& cluster, Weighted::Graph& graph, int num_neighbors)
{
    // What follows used to be in Cluster::Establish_close_connected_graph().
    // It is/was called from examine_graph() and Create_graph().

    using mcell_wire_wcps_map_t = std::map<const Facade::Blob*, std::map<int, std::set<int>>, Facade::BlobLess>;
    mcell_wire_wcps_map_t map_mcell_uindex_wcps, map_mcell_vindex_wcps, map_mcell_windex_wcps;

    std::map<Facade::Blob*, std::set<int>, Facade::BlobLess> map_mcell_indices;

    const auto& points = cluster.points();
    const auto& winds = cluster.wire_indices();

    for (Facade::Blob* mcell : cluster.children()) {
        std::map<int, std::set<int>> map_uindex_wcps;
        std::map<int, std::set<int>> map_vindex_wcps;
        std::map<int, std::set<int>> map_windex_wcps;

        std::vector<int> pinds = cluster.get_blob_indices(mcell);
        for (const int pind : pinds) {
            // auto v = vertex(pind, graph);  // retrieve vertex descriptor
            // (graph)[v].ident = pind;
            if (map_uindex_wcps.find(winds[0][pind]) == map_uindex_wcps.end()) {
                std::set<int> wcps;
                wcps.insert(pind);
                map_uindex_wcps[winds[0][pind]] = wcps;
            }
            else {
                map_uindex_wcps[winds[0][pind]].insert(pind);
            }

            if (map_vindex_wcps.find(winds[1][pind]) == map_vindex_wcps.end()) {
                std::set<int> wcps;
                wcps.insert(pind);
                map_vindex_wcps[winds[1][pind]] = wcps;
            }
            else {
                map_vindex_wcps[winds[1][pind]].insert(pind);
            }

            if (map_windex_wcps.find(winds[2][pind]) == map_windex_wcps.end()) {
                std::set<int> wcps;
                wcps.insert(pind);
                map_windex_wcps[winds[2][pind]] = wcps;
            }
            else {
                map_windex_wcps[winds[2][pind]].insert(pind);
            }
        }
        map_mcell_uindex_wcps[mcell] = map_uindex_wcps;
        map_mcell_vindex_wcps[mcell] = map_vindex_wcps;
        map_mcell_windex_wcps[mcell] = map_windex_wcps;
    }

    int num_edges = 0;

    // create graph for points inside the same mcell
    for (Facade::Blob* mcell : cluster.children()) {
        std::vector<int> pinds = cluster.get_blob_indices(mcell);
        int max_wire_interval = mcell->get_max_wire_interval();
        int min_wire_interval = mcell->get_min_wire_interval();
        std::map<int, std::set<int>>* map_max_index_wcps;
        std::map<int, std::set<int>>* map_min_index_wcps;
        if (mcell->get_max_wire_type() == 0) {
            map_max_index_wcps = &map_mcell_uindex_wcps[mcell];
        }
        else if (mcell->get_max_wire_type() == 1) {
            map_max_index_wcps = &map_mcell_vindex_wcps[mcell];
        }
        else {
            map_max_index_wcps = &map_mcell_windex_wcps[mcell];
        }
        if (mcell->get_min_wire_type() == 0) {
            map_min_index_wcps = &map_mcell_uindex_wcps[mcell];
        }
        else if (mcell->get_min_wire_type() == 1) {
            map_min_index_wcps = &map_mcell_vindex_wcps[mcell];
        }
        else {
            map_min_index_wcps = &map_mcell_windex_wcps[mcell];
        }

        for (const int pind1 : pinds) {
            int index_max_wire;
            int index_min_wire;
            if (mcell->get_max_wire_type() == 0) {
                index_max_wire = winds[0][pind1];
            }
            else if (mcell->get_max_wire_type() == 1) {
                index_max_wire = winds[1][pind1];
            }
            else {
                index_max_wire = winds[2][pind1];
            }
            if (mcell->get_min_wire_type() == 0) {
                index_min_wire = winds[0][pind1];
            }
            else if (mcell->get_min_wire_type() == 1) {
                index_min_wire = winds[1][pind1];
            }
            else {
                index_min_wire = winds[2][pind1];
            }
            std::vector<std::set<int>*> max_wcps_set;
            std::vector<std::set<int>*> min_wcps_set;
            // go through the first map and find the ones satisfying the condition
            for (auto it2 = map_max_index_wcps->begin(); it2 != map_max_index_wcps->end(); it2++) {
                if (std::abs(it2->first - index_max_wire) <= max_wire_interval) {
                    max_wcps_set.push_back(&(it2->second));
                }
            }
            // go through the second map and find the ones satisfying the condition
            for (auto it2 = map_min_index_wcps->begin(); it2 != map_min_index_wcps->end(); it2++) {
                if (std::abs(it2->first - index_min_wire) <= min_wire_interval) {
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

            {
                std::set<int> common_set;
                set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),
                                 std::inserter(common_set, common_set.begin()));

                for (auto it4 = common_set.begin(); it4 != common_set.end(); it4++) {
                    const int pind2 = *it4;
                    if (pind1 != pind2) {
                        // avoid duplicated edge addition
                        if (!boost::edge(pind1, pind2, graph).second) {

                            auto edge = add_edge(pind1,pind2,(sqrt(pow(points[0][pind1] - points[0][pind2], 2) +
                                                                                        pow(points[1][pind1] - points[1][pind2], 2) +
                                                                                        pow(points[2][pind1] - points[2][pind2], 2))),graph);
                            if (edge.second){
                                num_edges ++;
                            }
                        }
                    }
                }
            }
        }
    }

    // create graph for points between connected mcells, need to separate apa, face, and then ...
    std::map<int, std::map<int, std::vector<int> > > af_time_slices; // apa,face --> time slices 
    for (auto it = cluster.time_blob_map().begin(); it != cluster.time_blob_map().end(); it++) {
        int apa = it->first;
        for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
            int face = it1->first;
            std::vector<int> time_slices_vec;
            for (auto it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
                time_slices_vec.push_back(it2->first);
            }
            af_time_slices[apa][face] = time_slices_vec;
        }
    }
    
    std::vector<std::pair<const Facade::Blob*, const Facade::Blob*>> connected_mcells;

    for (auto it = af_time_slices.begin(); it != af_time_slices.end(); it++) {
        int apa = it->first;
        for (auto it1 = it->second.begin(); it1 != it->second.end(); it1++) {
            int face = it1->first;
            std::vector<int>& time_slices = it1->second;
            for (size_t i = 0; i != time_slices.size(); i++) {
                const auto& mcells_set = cluster.time_blob_map().at(apa).at(face).at(time_slices.at(i));

                // create graph for points in mcell inside the same time slice
                if (mcells_set.size() >= 2) {
                    for (auto it2 = mcells_set.begin(); it2 != mcells_set.end(); it2++) {
                        auto mcell1 = *it2;
                        auto it2p = it2;
                        if (it2p != mcells_set.end()) {
                            it2p++;
                            for (auto it3 = it2p; it3 != mcells_set.end(); it3++) {
                                auto mcell2 = *(it3);
                                if (mcell1->overlap_fast(*mcell2, 2))
                                    connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                            }
                        }
                    }
                }
                // create graph for points between connected mcells in adjacent time slices + 1, if not, + 2
                std::vector<Facade::BlobSet> vec_mcells_set;
                if (i + 1 < time_slices.size()) {
                    if (time_slices.at(i + 1) - time_slices.at(i) == 1*cluster.grouping()->get_nticks_per_slice().at(apa).at(face)) {
                        vec_mcells_set.push_back(cluster.time_blob_map().at(apa).at(face).at(time_slices.at(i + 1)));
                        if (i + 2 < time_slices.size())
                            if (time_slices.at(i + 2) - time_slices.at(i) == 2*cluster.grouping()->get_nticks_per_slice().at(apa).at(face))
                                vec_mcells_set.push_back(cluster.time_blob_map().at(apa).at(face).at(time_slices.at(i + 2)));
                    }
                    else if (time_slices.at(i + 1) - time_slices.at(i) == 2*cluster.grouping()->get_nticks_per_slice().at(apa).at(face)) {
                        vec_mcells_set.push_back(cluster.time_blob_map().at(apa).at(face).at(time_slices.at(i + 1)));
                    }
                }
                //    bool flag = false;
                for (size_t j = 0; j != vec_mcells_set.size(); j++) {
                    //      if (flag) break;
                    auto& next_mcells_set = vec_mcells_set.at(j);
                    for (auto it1 = mcells_set.begin(); it1 != mcells_set.end(); it1++) {
                        auto mcell1 = (*it1);
                        for (auto it2 = next_mcells_set.begin(); it2 != next_mcells_set.end(); it2++) {
                            auto mcell2 = (*it2);
                            if (mcell1->overlap_fast(*mcell2, 2)) {
                                //	    flag = true; // correct???
                                connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                            }
                        }
                    }
                }
            }
        }
    }

    // establish edge ...
    const int max_num_nodes = num_neighbors;
    std::map<std::pair<int, int>, std::set<std::pair<double, int>>> closest_index;

    for (auto it = connected_mcells.begin(); it != connected_mcells.end(); it++) {
        auto mcell1 = (*it).first;
        auto mcell2 = (*it).second;

        std::vector<int> pinds1 = cluster.get_blob_indices(mcell1);
        std::vector<int> pinds2 = cluster.get_blob_indices(mcell2);

        // test 2 against 1 ...
        int max_wire_interval = mcell1->get_max_wire_interval();
        int min_wire_interval = mcell1->get_min_wire_interval();
        std::map<int, std::set<int>>* map_max_index_wcps;
        std::map<int, std::set<int>>* map_min_index_wcps;

        if (mcell1->get_max_wire_type() == 0) {
            map_max_index_wcps = &map_mcell_uindex_wcps.at(mcell2);
        }
        else if (mcell1->get_max_wire_type() == 1) {
            map_max_index_wcps = &map_mcell_vindex_wcps.at(mcell2);
        }
        else {
            map_max_index_wcps = &map_mcell_windex_wcps.at(mcell2);
        }
        if (mcell1->get_min_wire_type() == 0) {
            map_min_index_wcps = &map_mcell_uindex_wcps.at(mcell2);
        }
        else if (mcell1->get_min_wire_type() == 1) {
            map_min_index_wcps = &map_mcell_vindex_wcps.at(mcell2);
        }
        else {
            map_min_index_wcps = &map_mcell_windex_wcps.at(mcell2);
        }

        for (const int pind1 : pinds1) {
            int index_max_wire;
            int index_min_wire;
            if (mcell1->get_max_wire_type() == 0) {
                index_max_wire = winds[0][pind1];
            }
            else if (mcell1->get_max_wire_type() == 1) {
                index_max_wire = winds[1][pind1];
            }
            else {
                index_max_wire = winds[2][pind1];
            }
            if (mcell1->get_min_wire_type() == 0) {
                index_min_wire = winds[0][pind1];
            }
            else if (mcell1->get_min_wire_type() == 1) {
                index_min_wire = winds[1][pind1];
            }
            else {
                index_min_wire = winds[2][pind1];
            }
            std::vector<std::set<int>*> max_wcps_set;
            std::vector<std::set<int>*> min_wcps_set;
            // go through the first map and find the ones satisfying the condition
            for (auto it2 = map_max_index_wcps->begin(); it2 != map_max_index_wcps->end(); it2++) {
                if (std::abs(it2->first - index_max_wire) <= max_wire_interval) {
                    max_wcps_set.push_back(&(it2->second));
                }
            }
            // go through the second map and find the ones satisfying the condition
            for (auto it2 = map_min_index_wcps->begin(); it2 != map_min_index_wcps->end(); it2++) {
                if (std::abs(it2->first - index_min_wire) <= min_wire_interval) {
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

            {
                std::set<int> common_set;
                set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),
                                 std::inserter(common_set, common_set.begin()));

                for (auto it4 = common_set.begin(); it4 != common_set.end(); it4++) {
                    const int pind2 = *it4;
                    if (pind1 != pind2) {
                        double dis = sqrt(pow(points[0][pind1] - points[0][pind2], 2) +
                                          pow(points[1][pind1] - points[1][pind2], 2) +
                                          pow(points[2][pind1] - points[2][pind2], 2));
                        auto b2 = cluster.blob_with_point(pind2);
                        auto key = std::make_pair(pind1, b2->slice_index_min());

                        if (closest_index.find(key) == closest_index.end()) {
                            std::set<std::pair<double, int> > temp_sets;
                            temp_sets.insert(std::make_pair(dis,pind2));
                            closest_index[key] = temp_sets;
                        }
                        else {
                            closest_index[key].insert(std::make_pair(dis,pind2));
                            if (closest_index[key].size() > static_cast<size_t>(max_num_nodes)) {
                                auto it5 = closest_index[key].begin();
                                for (int qx = 0; qx!=max_num_nodes;qx++){
                                    it5++;
                                }
                                closest_index[key].erase(it5,closest_index[key].end());
                            }
                        }
                    }
                }
            }
        }

        // test 1 against 2 ...
        max_wire_interval = mcell2->get_max_wire_interval();
        min_wire_interval = mcell2->get_min_wire_interval();
        if (mcell2->get_max_wire_type() == 0) {
            map_max_index_wcps = &map_mcell_uindex_wcps[mcell1];
        }
        else if (mcell2->get_max_wire_type() == 1) {
            map_max_index_wcps = &map_mcell_vindex_wcps[mcell1];
        }
        else {
            map_max_index_wcps = &map_mcell_windex_wcps[mcell1];
        }
        if (mcell2->get_min_wire_type() == 0) {
            map_min_index_wcps = &map_mcell_uindex_wcps[mcell1];
        }
        else if (mcell2->get_min_wire_type() == 1) {
            map_min_index_wcps = &map_mcell_vindex_wcps[mcell1];
        }
        else {
            map_min_index_wcps = &map_mcell_windex_wcps[mcell1];
        }
        for (const int pind1 : pinds2) {
            int index_max_wire;
            int index_min_wire;
            if (mcell2->get_max_wire_type() == 0) {
                index_max_wire = winds[0][pind1];
            }
            else if (mcell2->get_max_wire_type() == 1) {
                index_max_wire = winds[1][pind1];
            }
            else {
                index_max_wire = winds[2][pind1];
            }
            if (mcell2->get_min_wire_type() == 0) {
                index_min_wire = winds[0][pind1];
            }
            else if (mcell2->get_min_wire_type() == 1) {
                index_min_wire = winds[1][pind1];
            }
            else {
                index_min_wire = winds[2][pind1];
            }
            std::vector<std::set<int>*> max_wcps_set;
            std::vector<std::set<int>*> min_wcps_set;
            // go through the first map and find the ones satisfying the condition
            for (auto it2 = map_max_index_wcps->begin(); it2 != map_max_index_wcps->end(); it2++) {
                if (std::abs(it2->first - index_max_wire) <= max_wire_interval) {
                    max_wcps_set.push_back(&(it2->second));
                }
            }
            // go through the second map and find the ones satisfying the condition
            for (auto it2 = map_min_index_wcps->begin(); it2 != map_min_index_wcps->end(); it2++) {
                if (std::abs(it2->first - index_min_wire) <= min_wire_interval) {
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

            {
                std::set<int> common_set;
                set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),
                                 std::inserter(common_set, common_set.begin()));

                for (auto it4 = common_set.begin(); it4 != common_set.end(); it4++) {
                    const int pind2 = *it4;
                    if (pind1 != pind2) {
                        double dis = sqrt(pow(points[0][pind1] - points[0][pind2], 2) +
                                          pow(points[1][pind1] - points[1][pind2], 2) +
                                          pow(points[2][pind1] - points[2][pind2], 2));
                        auto b2 = cluster.blob_with_point(pind2);
                        auto key = std::make_pair(pind1, b2->slice_index_min());

                        if (closest_index.find(key) == closest_index.end()) {
                            std::set<std::pair<double, int> > temp_sets;
                            temp_sets.insert(std::make_pair(dis,pind2));
                            closest_index[key] = temp_sets;
                        }
                        else {
                            closest_index[key].insert(std::make_pair(dis,pind2));
                            if (closest_index[key].size() > static_cast<size_t>(max_num_nodes)) {
                                auto it5 = closest_index[key].begin();
                                for (int qx = 0; qx!=max_num_nodes;qx++){
                                    it5++;
                                }
                                closest_index[key].erase(it5,closest_index[key].end());
                            }
                        }
                    }
                }
            }
        }
    }

    for (auto it4 = closest_index.begin(); it4 != closest_index.end(); it4++) {
        int index1 = it4->first.first;
        for (auto it5 = it4->second.begin(); it5!=it4->second.end(); it5++){
            int index2 = (*it5).second;
            double dis = (*it5).first;

            if (!boost::edge(index1, index2, graph).second) {

                auto edge = add_edge(index1,index2,dis,graph);
                if (edge.second){
                    num_edges ++;
                }
            }
            // protect against dead cells ...
            if (it5 == it4->second.begin() && dis > 0.25*units::cm)
                break;
        }

    }

    (void)num_edges; // suppress unused variable warning
}

using namespace WireCell::Clus::Facade;


void Graphs::connect_graph_closely_pid(const Facade::Cluster& cluster, Weighted::Graph& graph)
{
    // PID-specific parameters (from prototype)
    const int max_num_nodes = 5;
    const double protection_distance = 0.25 * units::cm;
    
    // Build wire index maps for each blob (equivalent to mcell in prototype)
    using mcell_wire_wcps_map_t = std::map<const Facade::Blob*, std::map<int, std::set<int>>, Facade::BlobLess>;
    mcell_wire_wcps_map_t map_mcell_uindex_wcps, map_mcell_vindex_wcps, map_mcell_windex_wcps;

    const auto& points = cluster.points();
    const auto& winds = cluster.wire_indices();

    // Build wire index maps for each blob
    for (const Facade::Blob* mcell : cluster.children()) {
        std::map<int, std::set<int>> map_uindex_wcps;
        std::map<int, std::set<int>> map_vindex_wcps;
        std::map<int, std::set<int>> map_windex_wcps;

        std::vector<int> pinds = cluster.get_blob_indices(mcell);
        for (const int pind : pinds) {
            // Build U wire index map
            auto u_it = map_uindex_wcps.find(winds[0][pind]);
            if (u_it == map_uindex_wcps.end()) {
                std::set<int> wcps;
                wcps.insert(pind);
                map_uindex_wcps[winds[0][pind]] = wcps;
            } else {
                u_it->second.insert(pind);
            }

            // Build V wire index map
            auto v_it = map_vindex_wcps.find(winds[1][pind]);
            if (v_it == map_vindex_wcps.end()) {
                std::set<int> wcps;
                wcps.insert(pind);
                map_vindex_wcps[winds[1][pind]] = wcps;
            } else {
                v_it->second.insert(pind);
            }

            // Build W wire index map
            auto w_it = map_windex_wcps.find(winds[2][pind]);
            if (w_it == map_windex_wcps.end()) {
                std::set<int> wcps;
                wcps.insert(pind);
                map_windex_wcps[winds[2][pind]] = wcps;
            } else {
                w_it->second.insert(pind);
            }
        }
        
        map_mcell_uindex_wcps[mcell] = map_uindex_wcps;
        map_mcell_vindex_wcps[mcell] = map_vindex_wcps;
        map_mcell_windex_wcps[mcell] = map_windex_wcps;
    }

    int num_edges = 0;

    // Phase 1: Create graph for points inside the same mcell (blob)
    for (const Facade::Blob* mcell : cluster.children()) {
        std::vector<int> pinds = cluster.get_blob_indices(mcell);
        int max_wire_interval = mcell->get_max_wire_interval();
        int min_wire_interval = mcell->get_min_wire_interval();

        // Get appropriate wire index maps based on wire types
        std::map<int, std::set<int>>* map_max_index_wcps;
        std::map<int, std::set<int>>* map_min_index_wcps;
        
        if (mcell->get_max_wire_type() == 0) {
            map_max_index_wcps = &map_mcell_uindex_wcps[mcell];
        } else if (mcell->get_max_wire_type() == 1) {
            map_max_index_wcps = &map_mcell_vindex_wcps[mcell];
        } else {
            map_max_index_wcps = &map_mcell_windex_wcps[mcell];
        }

        if (mcell->get_min_wire_type() == 0) {
            map_min_index_wcps = &map_mcell_uindex_wcps[mcell];
        } else if (mcell->get_min_wire_type() == 1) {
            map_min_index_wcps = &map_mcell_vindex_wcps[mcell];
        } else {
            map_min_index_wcps = &map_mcell_windex_wcps[mcell];
        }

        for (const int pind1 : pinds) {
            int index_max_wire, index_min_wire;
            
            // Get wire indices for current point
            if (mcell->get_max_wire_type() == 0) {
                index_max_wire = winds[0][pind1];
            } else if (mcell->get_max_wire_type() == 1) {
                index_max_wire = winds[1][pind1];
            } else {
                index_max_wire = winds[2][pind1];
            }

            if (mcell->get_min_wire_type() == 0) {
                index_min_wire = winds[0][pind1];
            } else if (mcell->get_min_wire_type() == 1) {
                index_min_wire = winds[1][pind1];
            } else {
                index_min_wire = winds[2][pind1];
            }

            // Find candidate points within wire intervals
            std::vector<std::set<int>*> max_wcps_set;
            std::vector<std::set<int>*> min_wcps_set;

            // Find points within max wire interval
            for (auto& wire_pair : *map_max_index_wcps) {
                if (std::abs(wire_pair.first - index_max_wire) <= max_wire_interval) {
                    max_wcps_set.push_back(&wire_pair.second);
                }
            }

            // Find points within min wire interval
            for (auto& wire_pair : *map_min_index_wcps) {
                if (std::abs(wire_pair.first - index_min_wire) <= min_wire_interval) {
                    min_wcps_set.push_back(&wire_pair.second);
                }
            }

            // Create candidate sets
            std::set<int> wcps_set1, wcps_set2;
            for (const auto* wcp_set : max_wcps_set) {
                wcps_set1.insert(wcp_set->begin(), wcp_set->end());
            }
            for (const auto* wcp_set : min_wcps_set) {
                wcps_set2.insert(wcp_set->begin(), wcp_set->end());
            }

            // Find intersection of candidate sets
            std::set<int> common_set;
            std::set_intersection(wcps_set1.begin(), wcps_set1.end(), 
                                wcps_set2.begin(), wcps_set2.end(),
                                std::inserter(common_set, common_set.begin()));

            // Connect to all valid points in the same mcell
            for (const int pind2 : common_set) {
                if (pind2 != pind1) {
                    double distance = std::sqrt(
                        std::pow(points[0][pind1] - points[0][pind2], 2) +
                        std::pow(points[1][pind1] - points[1][pind2], 2) +
                        std::pow(points[2][pind1] - points[2][pind2], 2)
                    );
                    
                    if (!boost::edge(pind1, pind2, graph).second){
                        auto edge_result = add_edge(pind1, pind2, distance, graph);

                        // std::cout << mcell->slice_index_min() << " " << mcell->u_wire_index_min() << " " << mcell->v_wire_index_min() << " " 
                        //         << mcell->w_wire_index_min() << " " << pind1 << " " << pind2 
                        //         << " " << edge_result.second << std::endl;

                        if (edge_result.second) {
                            num_edges++;
                        }
                    }
                }
            }
        }
    }

    // Phase 2: Create graph for points between connected mcells across time slices
    const auto& time_cells_set_map = cluster.time_blob_map();
    
    // Build time slice structure
    std::map<int, std::map<int, std::vector<int>>> af_time_slices; // apa,face --> time slices
    for (const auto& apa_pair : time_cells_set_map) {
        int apa = apa_pair.first;
        for (const auto& face_pair : apa_pair.second) {
            int face = face_pair.first;
            std::vector<int> time_slices_vec;
            for (const auto& time_pair : face_pair.second) {
                time_slices_vec.push_back(time_pair.first);
            }
            af_time_slices[apa][face] = time_slices_vec;
        }
    }

    // Find connected mcells across time slices
    std::vector<std::pair<const Facade::Blob*, const Facade::Blob*>> connected_mcells;
    
    for (const auto& apa_pair : af_time_slices) {
        int apa = apa_pair.first;
        for (const auto& face_pair : apa_pair.second) {
            int face = face_pair.first;
            const std::vector<int>& time_slices = face_pair.second;
            
            for (size_t i = 0; i < time_slices.size(); i++) {
                const auto& mcells_set = time_cells_set_map.at(apa).at(face).at(time_slices[i]);
                
                // Connect mcells within the same time slice
                if (mcells_set.size() >= 2) {
                    for (auto it1 = mcells_set.begin(); it1 != mcells_set.end(); it1++) {
                        auto mcell1 = *it1;
                        auto it2 = it1;
                        it2++;
                        for (; it2 != mcells_set.end(); it2++) {
                            auto mcell2 = *it2;
                            if (mcell1->overlap_fast(*mcell2, 2)) {
                                connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                            }
                        }
                    }
                }
                
                // Connect mcells across adjacent time slices
                 std::vector<Facade::BlobSet> vec_mcells_set;
                if (i + 1 < time_slices.size()) {
                    if (time_slices.at(i + 1) - time_slices.at(i) == 1*cluster.grouping()->get_nticks_per_slice().at(apa).at(face)) {
                        vec_mcells_set.push_back(cluster.time_blob_map().at(apa).at(face).at(time_slices.at(i + 1)));
                        if (i + 2 < time_slices.size())
                            if (time_slices.at(i + 2) - time_slices.at(i) == 2*cluster.grouping()->get_nticks_per_slice().at(apa).at(face))
                                vec_mcells_set.push_back(cluster.time_blob_map().at(apa).at(face).at(time_slices.at(i + 2)));
                    }
                    else if (time_slices.at(i + 1) - time_slices.at(i) == 2*cluster.grouping()->get_nticks_per_slice().at(apa).at(face)) {
                        vec_mcells_set.push_back(cluster.time_blob_map().at(apa).at(face).at(time_slices.at(i + 1)));
                    }
                }
                
                
                // Find overlapping mcells between current and next time slices
                for (const auto& next_mcells_set : vec_mcells_set) {
                    for (const auto* mcell1 : mcells_set) {
                        for (const auto* mcell2 : next_mcells_set) {
                            if (mcell1->overlap_fast(*mcell2, 2)) {
                                connected_mcells.push_back(std::make_pair(mcell1, mcell2));
                            }
                        }
                    }
                }
            }
        }
    }

    // std::cout << "Connected mcells across time slices: " << connected_mcells.size() << std::endl;

    // Phase 3: Establish cross-mcell connections using PID-specific logic
    std::map<std::pair<int, int>, std::set<std::pair<double, int>>> closest_index;
    
    for (const auto& mcell_pair : connected_mcells) {
        const Facade::Blob* mcell1 = mcell_pair.first;
        const Facade::Blob* mcell2 = mcell_pair.second;
        
        // Process connections in both directions
        for (int direction = 0; direction < 2; direction++) {
            const Facade::Blob* source_mcell = (direction == 0) ? mcell1 : mcell2;
            const Facade::Blob* target_mcell = (direction == 0) ? mcell2 : mcell1;
            
            std::vector<int> source_pinds = cluster.get_blob_indices(source_mcell);
            
            int max_wire_interval = source_mcell->get_max_wire_interval();
            int min_wire_interval = source_mcell->get_min_wire_interval();
            
            std::map<int, std::set<int>>* map_max_index_wcps;
            std::map<int, std::set<int>>* map_min_index_wcps;
            
            // Select appropriate wire index maps based on source mcell's wire types
            if (source_mcell->get_max_wire_type() == 0) {
                map_max_index_wcps = &map_mcell_uindex_wcps[target_mcell];
            } else if (source_mcell->get_max_wire_type() == 1) {
                map_max_index_wcps = &map_mcell_vindex_wcps[target_mcell];
            } else {
                map_max_index_wcps = &map_mcell_windex_wcps[target_mcell];
            }
            
            if (source_mcell->get_min_wire_type() == 0) {
                map_min_index_wcps = &map_mcell_uindex_wcps[target_mcell];
            } else if (source_mcell->get_min_wire_type() == 1) {
                map_min_index_wcps = &map_mcell_vindex_wcps[target_mcell];
            } else {
                map_min_index_wcps = &map_mcell_windex_wcps[target_mcell];
            }
            
            for (const int pind1 : source_pinds) {
                int index_max_wire, index_min_wire;
                
                // Get wire indices for current point
                if (source_mcell->get_max_wire_type() == 0) {
                    index_max_wire = winds[0][pind1];
                } else if (source_mcell->get_max_wire_type() == 1) {
                    index_max_wire = winds[1][pind1];
                } else {
                    index_max_wire = winds[2][pind1];
                }
                
                if (source_mcell->get_min_wire_type() == 0) {
                    index_min_wire = winds[0][pind1];
                } else if (source_mcell->get_min_wire_type() == 1) {
                    index_min_wire = winds[1][pind1];
                } else {
                    index_min_wire = winds[2][pind1];
                }
                
                // Find candidate points within wire intervals
                std::vector<std::set<int>*> max_wcps_set;
                std::vector<std::set<int>*> min_wcps_set;
                
                for (auto& wire_pair : *map_max_index_wcps) {
                    if (std::abs(wire_pair.first - index_max_wire) <= max_wire_interval) {
                        max_wcps_set.push_back(&wire_pair.second);
                    }
                }
                
                for (auto& wire_pair : *map_min_index_wcps) {
                    if (std::abs(wire_pair.first - index_min_wire) <= min_wire_interval) {
                        min_wcps_set.push_back(&wire_pair.second);
                    }
                }
                
                // Create candidate sets
                std::set<int> wcps_set1, wcps_set2;
                for (const auto* wcp_set : max_wcps_set) {
                    wcps_set1.insert(wcp_set->begin(), wcp_set->end());
                }
                for (const auto* wcp_set : min_wcps_set) {
                    wcps_set2.insert(wcp_set->begin(), wcp_set->end());
                }
                
                // Find intersection
                std::set<int> common_set;
                std::set_intersection(wcps_set1.begin(), wcps_set1.end(), 
                                    wcps_set2.begin(), wcps_set2.end(),
                                    std::inserter(common_set, common_set.begin()));
                
                // Build closest index map for PID-specific connection limiting
                for (const int pind2 : common_set) {
                    if (pind2 != pind1) {
                        double distance = std::sqrt(
                            std::pow(points[0][pind1] - points[0][pind2], 2) +
                            std::pow(points[1][pind1] - points[1][pind2], 2) +
                            std::pow(points[2][pind1] - points[2][pind2], 2)
                        );
                        
                        // Use target mcell's time slice as key
                        int target_time_slice = target_mcell->slice_index_min();
                        auto key = std::make_pair(pind1, target_time_slice);
                        
                        auto it = closest_index.find(key);
                        if (it == closest_index.end()) {
                            std::set<std::pair<double, int>> temp_set;
                            temp_set.insert(std::make_pair(distance, pind2));
                            closest_index[key] = temp_set;
                        } else {
                            it->second.insert(std::make_pair(distance, pind2));
                            // Keep only the closest max_num_nodes connections
                            if (it->second.size() > static_cast<size_t>(max_num_nodes)) {
                                auto erase_it = it->second.begin();
                                std::advance(erase_it, max_num_nodes);
                                it->second.erase(erase_it, it->second.end());
                            }
                        }
                    }
                }
            }
        }
    }
    // std::cout << closest_index.size() << " closest index entries created" << std::endl;

    // Phase 4: Add the selected edges from closest_index map
    for (const auto& closest_pair : closest_index) {
        int pind1 = closest_pair.first.first;
        
        for (const auto& distance_pair : closest_pair.second) {
            int pind2 = distance_pair.second;
            double distance = distance_pair.first;
            
            if (!boost::edge(pind1, pind2, graph).second) {
                auto edge_result = add_edge(pind1, pind2, distance, graph);

                // std::cout << "Adding edge " << pind1 << " " << pind2 << " " << distance << " " << edge_result.second << std::endl;


                if (edge_result.second) {
                    num_edges++;
                }
            }
            
            // PID-specific protection: break if first connection is too far
            if (distance_pair == *closest_pair.second.begin() && distance > protection_distance) {
                break;
            }
        }
    }
    
    // Debug output (similar to prototype)
    // std::cout << "PID Graph: " << num_edges << " edges added" << std::endl;
    (void)num_edges; // suppress unused variable warning
}