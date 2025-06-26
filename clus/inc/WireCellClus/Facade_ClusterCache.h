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

        // Maps apa/face/slice to child blob facade.  Depends on PC tree structure.
        time_blob_map_t time_blob_map;

        // Maps child blob node facade to the set of point indices for points in
        // that blob.  Depends on DEFAULT SCOPE.
        std::map<const Blob*, std::vector<int>> map_mcell_indices;

        // The subset of points() that make up a convex hull.  Depends on DEFAULT SCOPE.
        std::vector<geo_point_t> hull_points;

        // The "length" of a cluster estimated by time and wire extents.
        // A zero length is invalid.
        double length{0};
        // The number of points in the DEFAULT SCOPE.
        int npoints{0};


        // mutable bool m_pca_calculated{false}; // use pca_axis.size()
        struct PCA {
            geo_point_t center;
            std::vector<geo_vector_t> axis;
            std::vector<double> values;
            // if vectors are empty, PCA is invalid.
        };
        // Depends on DEFAULT SCOPE
        std::unique_ptr<PCA> pca;


        // Wire plane IDs by point (3d scoped view) index. Depends on the RAW SCOPE.
        std::vector<int> point_wpids;

        // Wire plane IDs by blob (child node) index. Depends on the DEFAULT SCOPE.
        std::vector<WireCell::WirePlaneId> blob_wpids;

        // Wire indices by point index (3 vectors for u,v,w)
        std::vector<int> point_u_wire_indices;
        std::vector<int> point_v_wire_indices; 
        std::vector<int> point_w_wire_indices;
        
        // Charge values by point index (3 vectors for u,v,w)
        std::vector<double> point_u_charges;
        std::vector<double> point_v_charges;
        std::vector<double> point_w_charges;
        
        // Charge uncertainties by point index (3 vectors for u,v,w) 
        std::vector<double> point_u_charge_uncs;
        std::vector<double> point_v_charge_uncs;
        std::vector<double> point_w_charge_uncs;
    };

}

#endif
