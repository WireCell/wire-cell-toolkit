/** This provides a "workspace" class for creating the steiner graph.

    It is meant to be equivalent to the Steiner-related SUBSET of WCP's
    PR3DCluster methods and data.

    See also SteinerFunctions.h for any free functions.
*/

#ifndef WIRECELLCLUS_STEINER
#define WIRECELLCLUS_STEINER

#include "WireCellClus/Graphs.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/Facade_Blob.h"

#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IDetectorVolumes.h"

#include <set>
#include <map>

namespace WireCell::Clus::Steiner {


    class Grapher {
    public:

        // This holds various "global" and const info sources.  See
        // CreateSteinerGraph for an example of how it is provided.
        struct Config {
            IBlobSampler::pointer sampler;
            IDetectorVolumes::pointer dv;
        };

        /// Construct with an existing cluster facade.  Caller must assure the
        /// underlying cluster node is kept live.
        Grapher(Facade::Cluster& cluster, const Config& cfg);
        Grapher() = delete;

        Facade::Cluster& cluster() { return m_cluster; }
        const Facade::Cluster& cluster() const { return m_cluster; }
        Config config() const { return m_config; }

        ///
        ///  Types
        ///
        
        /// This is top-level interface which is called in the visit() loop.  It
        /// is expected to create and add to our cluster using the given
        /// graph_name a "steiner graph" of C++ type:
        using graph_type = WireCell::Clus::Graphs::Weighted::Graph;

        /// A type that holds a set of graph vertices.  A vertex (descriptor) IS an index.
        using graph_vertex_set = std::set<size_t>;

        /// A type that maps blobs to graph vertices
        using blob_vertex_map = std::map<const Facade::Blob*, graph_vertex_set>;


        ///
        ///  The main entry method
        ///

        /// Create and return a steiner graph for the cluster.
        graph_type create_steiner_graph();


        ///
        ///  Intermediate methods
        ///

        /// Populate blob PCs with sampled points.  This was a free function in
        /// WCP but we make it a method to make use of the Config.  NOTE: this
        /// may be better replaced with a method that clones a cluster from an
        /// existing, already sampled cluster.
        void calc_sampling_points(/*, ...*/);

        graph_vertex_set find_peak_point_indices(bool disable_dead_mix_cell);
        blob_vertex_map form_cell_points_map();
        graph_vertex_set find_steiner_terminals(bool disable_dead_mix_cell=true);
        void establish_same_blob_steiner_edges(graph_type& graph, 
                                               bool disable_dead_mix_cell=true, int flag=1);
    
        graph_type create_steiner_tree(/*what type for point_cloud_steiner?*/);

    private:
        // The Grapher "wraps" a Cluster.  As the Cluster is a *facade* of an
        // underlying PC tree node, we do not own the Cluster and we rely on
        // whoever owns us to keep the underlying cluster node alive. as long as
        // we are alive.
        Facade::Cluster& m_cluster;     

        // This holds various "global" info sources
        const Config& m_config;


        // XIN: add any more data and methods you need here.  

    };


}

#endif
