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
#include "WireCellClus/IPCTransform.h"

#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IDetectorVolumes.h"

#include "WireCellUtil/Logging.h"

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
            WireCell::Clus::IPCTransformSet::pointer pcts;
        };
        Log::logptr_t log;

        /// Construct with an existing cluster facade.  Caller must assure the
        /// underlying cluster node is kept live.
        Grapher(Facade::Cluster& cluster, const Config& cfg, Log::logptr_t log);
        Grapher() = delete;

        ///
        ///  Types
        ///
        
        /// Forward some types from Graphs.h
        using graph_type = WireCell::Clus::Graphs::Weighted::graph_type;
        using vertex_type = WireCell::Clus::Graphs::Weighted::vertex_type;
        using edge_type = WireCell::Clus::Graphs::Weighted::edge_type;
        using vertex_set = WireCell::Clus::Graphs::Weighted::vertex_set;
        using edge_set = WireCell::Clus::Graphs::Weighted::edge_set;
        using edge_weight_type = WireCell::Clus::Graphs::Weighted::edge_weight_type;

        /// A type that maps blobs to graph vertices
        using blob_vertex_map = std::map<const Facade::Blob*, vertex_set>;


        ///
        /// Basic data accessors.
        ///
        
        Facade::Cluster& cluster() { return m_cluster; }
        const Facade::Cluster& cluster() const { return m_cluster; }
        Config config() const { return m_config; }


        ///
        /// Helper methods - these are general purpose, primitive.
        /// 

        ///  Get a graph, possibly making it on the fly if flavor is one of the
        ///  3 reserved names.
        graph_type& get_graph(const std::string& flavor = "basic");
        const graph_type& get_graph(const std::string& flavor = "basic") const ;
        
        /// Return a PC held by the cluster node of the given name.  If it does
        /// not exist, one is derived from the default scoped view, saved to
        /// that name, and returned.
        PointCloud::Dataset& get_point_cloud(const std::string& name = "steiner");

        /// Store a point cloud by std::move()
        void put_point_cloud(PointCloud::Dataset&& pc, const std::string& name = "steiner");
        /// Store a point cloud by copy
        void put_point_cloud(const PointCloud::Dataset& pc, const std::string& name = "steiner");


        ///
        /// Temporary algorithm to show how to do some things.
        ///
        graph_type fake_steiner_graph();        


        ///
        ///  The real main entry method
        ///

        /// Create and return a steiner graph for the cluster.
        graph_type create_steiner_graph();


        ///
        ///  Intermediate algorithm methods
        ///

        /// Populate blob PCs with sampled points.  This was a free function in
        /// WCP but we make it a method to make use of the Config.  NOTE: this
        /// may be better replaced with a method that clones a cluster from an
        /// existing, already sampled cluster.
        void calc_sampling_points(/*, ...*/);

        vertex_set find_peak_point_indices(bool disable_dead_mix_cell);
        blob_vertex_map form_cell_points_map();
        vertex_set find_steiner_terminals(bool disable_dead_mix_cell=true);
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
