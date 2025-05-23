#ifndef WIRECELLCLUS_FACADE_CLUSTERCACHE
#define WIRECELLCLUS_FACADE_CLUSTERCACHE

#include "WireCellClus/Facade_Blob.h"
#include "WireCellClus/Facade_Util.h"
#include "WireCellClus/Graphs.h"

// DO NOT #include Facade_Cluster.h itself. It depends on us, not vice versa.

#include <map>
#include <vector>

namespace WireCell::Clus::Facade {


    // ALL cached items for Cluster internal use go in this struct.
    //
    // DO NOT PLACE THEM BARE DIRECTLY IN Cluster.

    struct ClusterCache {
        // order is synchronized with children()

        using time_blob_map_t = std::map<int, std::map<int, std::map<int, BlobSet> > >; // apa, face, time, blobset

        time_blob_map_t time_blob_map;  // lazy, do not access directly.
        std::map<const Blob*, std::vector<int>> map_mcell_indices; // lazy, do not access directly.

        // Add to private members in Facade_Cluster.h:
        std::vector<geo_point_t> hull_points;
        // mutable bool m_hull_calculated{false};  // use hull_poits.size()

        // Cached and lazily calculated in get_length().
        // Getting a new node invalidates by setting to 0.
        double length{0};
        // Cached and lazily calculated in npoints()
        int npoints{0};


        // mutable bool m_pca_calculated{false}; // use pca_axis.size()
        struct PCA {
            geo_point_t center;
            std::vector<geo_vector_t> axis;
            std::vector<double> values;
            // if vectors are empty, PCA is invalid.
        };
        std::unique_ptr<PCA> pca;

        // A cache of named graphs held by unique_ptr.  For now, this cache is
        // internal to Cluster.
        using graph_type = Graphs::Weighted::Graph;        
        using graph_ptr = std::unique_ptr<graph_type>;

        std::map<std::string, graph_ptr> graphs;

        // Cluster makes its own graphs for the purpose of calculating shortest
        // paths.  This is the lazy cache.  For now, the key string is either
        // "basic" or "ctpc" to distinguish the type of graph.
        std::map<std::string, Graphs::Weighted::GraphAlgorithms> spgraphs;

        // wire plane IDs by point (3d scoped view) index
        std::vector<int> point_wpids;

        // wire plane IDs by blob (child node) index
        std::vector<WireCell::WirePlaneId> blob_wpids;
    };

}

#endif
