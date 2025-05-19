#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Grouping.h"

using namespace WireCell;
using namespace WireCell::Clus;



namespace WireCell::Clus {
    Graph::Ident::graph_ptr close_connected_graph(const Facade::Cluster& cluster);
}


Graph::Ident::graph_ptr WireCell::Clus::close_connected_graph(const Facade::Cluster& cluster)
{
    auto the_graph = std::make_unique<Graph::Ident::graph_type>(cluster.npoints());
    
    // What follows used to be in Cluster::Establish_close_connected_graph().
    // It is/was called from deprecated examine_graph() and Create_graph().

    using mcell_wire_wcps_map_t = std::map<const Facade::Blob*, std::map<int, std::set<int>>, Facade::blob_less_functor>;
    mcell_wire_wcps_map_t map_mcell_uindex_wcps, map_mcell_vindex_wcps, map_mcell_windex_wcps;

    std::map<Facade::Blob*, std::set<int>, Facade::blob_less_functor> map_mcell_indices;

    const auto& points = cluster.points();
    const auto& winds = cluster.wire_indices();

    for (Facade::Blob* mcell : cluster.children()) {
        std::map<int, std::set<int>> map_uindex_wcps;
        std::map<int, std::set<int>> map_vindex_wcps;
        std::map<int, std::set<int>> map_windex_wcps;

        std::vector<int> pinds = cluster.get_blob_indices(mcell);
        for (const int pind : pinds) {
            auto v = vertex(pind, *the_graph);  // retrieve vertex descriptor
            (*the_graph)[v].ident = pind;
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

            {
                std::set<int> common_set;
                set_intersection(wcps_set1.begin(), wcps_set1.end(), wcps_set2.begin(), wcps_set2.end(),
                                 std::inserter(common_set, common_set.begin()));

                for (auto it4 = common_set.begin(); it4 != common_set.end(); it4++) {
                    const int pind2 = *it4;
                    if (pind1 != pind2) {
                        auto edge = add_edge(pind1,pind2,Graph::Ident::EdgeProp(sqrt(pow(points[0][pind1] - points[0][pind2], 2) +
                                                                                     pow(points[1][pind1] - points[1][pind2], 2) +
                                                                                     pow(points[2][pind1] - points[2][pind2], 2))),*the_graph);
                        if (edge.second){
                            num_edges ++;
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
                std::vector<Facade::Cluster::BlobSet> vec_mcells_set;
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
    const int max_num_nodes = 5;
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
                            if (closest_index[key].size()>max_num_nodes){
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
                            if (closest_index[key].size()>max_num_nodes){
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
            auto edge = add_edge(index1,index2,Graph::Ident::EdgeProp(dis),*the_graph);
            if (edge.second){
                num_edges ++;
            }
            // protect against dead cells ...
            if (it5 == it4->second.begin() && dis > 0.25*units::cm)
                break;
        }

    }

    return the_graph;
}
