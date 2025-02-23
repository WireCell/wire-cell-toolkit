#ifndef WIRECELLAUX_GROUPINGHELPER
#define WIRECELLAUX_GROUPINGHELPER

// #include "WireCellClus/MultiAlgBlobClustering.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"

#include <sstream>
#include <functional>
#include <unordered_set>
#include <unordered_map>

#include <string>

namespace WireCell::PointCloud::Facade {
  // create a new function to take in the original_grouping and the newly created shadow_grouping
    // the return of this new function should be a std::map<Cluster*, std::pair<Cluster*, Cluster*> >. 
    // The first Cluster should be the original cluster, the first of the pair Clusters should be the corresponding shadow_cluster, the second of the pair Cluster should be the main cluster due to separation. 
    // Inside this function, the algorithm should run as the following
    // 1. Loop over all clusters in the original Grouping, and find the corresponnding cluster in the shadow grouping. form a map between the two clusters.
    // 2. loop through these pairs, and doing separation using splits and using get_pcarray("isolated","perblob")
    // 3. form the output pairs using the key of the splits, and the corresponding the main cluster.
    // Use the mapping as needed
    // Process original and shadow groupings to create a mapping between clusters
    std::map<Cluster*, std::tuple<Cluster*, int, Cluster*>> process_groupings_helper(
        Grouping& original, 
        Grouping& shadow, 
        const std::string& aname = "isolated", 
        const std::string& pname = "perblob");  // Removed const here

}

#endif