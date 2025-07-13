// This provides RetileCluster aka "IPCTreeMutate".  
//
// Warning: this lives up to its name.  It may change the input cluster.
// 
// It requires the input cluster node to have a (grouping) node parent.
//
// The retiling of a cluster follows this general sequence:
//
// 1) constructs layers of "activity" from input grouping.
// 2) applies "hacks" to the activity.
// 3) runs WCT tiling to create blobs.
// 4) runs blobs sampling to make point clouds
// 5) produces clusters such that the new blobs formed from an old cluster form a new "shadow" cluster.
// 6) forms a PC-tree
// 7) outputs the new grouping
//


#ifndef WIRECELLCLUS_RETILE_CLUSTER_H
#define WIRECELLCLUS_RETILE_CLUSTER_H


#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/RayHelpers.h"

#include "WireCellIface/IBlob.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IPCTreeMutate.h"

#include "WireCellAux/PlaneTools.h"

#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/ClusteringFuncsMixins.h"


#include "WireCellClus/IEnsembleVisitor.h"
#include "WireCellClus/ClusteringFuncs.h"
#include "WireCellClus/ClusteringFuncsMixins.h"

#include "WireCellIface/IConfigurable.h"

#include "WireCellUtil/NamedFactory.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SamplingHelpers.h"

#include "WireCellUtil/PointTree.h"

#include "WireCellAux/SimpleSlice.h"
#include "WireCellClus/GroupingHelper.h"


#include <vector>
#include <map>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

namespace WireCell::Clus {



class RetileCluster : public IConfigurable, public IPCTreeMutate, protected Clus::NeedDV, protected Clus::NeedPCTS {

   

public:

    RetileCluster() {}
    virtual ~RetileCluster() {};

    // IConfigurable API
    void configure(const WireCell::Configuration& config);
    virtual Configuration default_configuration() const {
        Configuration cfg;
        return cfg;
    }

    // IPCTreeMutate API
    virtual std::unique_ptr<node_t> mutate(node_t& node) const;

protected:
     // Step 0. Collect grouping info
    Facade::Cluster* reinitialize(Points::node_type& node) const;

     // Cache
    mutable Grouping* m_grouping = nullptr;
    mutable std::map<WirePlaneId , std::vector<double> > m_wpid_angles;
    std::map<int, std::map<int, std::vector<Aux::WirePlaneInfo>>> m_plane_infos;
    /** Configuration: "sampler" (required)

        The type/name an IBlobSampler for producing the "3d" point cloud.

        If not given, the retailed blob tree nodes will not have point clouds.
    */

    bool m_verbose{false};

    
    std::map<int, std::map<int, WireCell::IBlobSampler::pointer>> m_samplers;

    // fixme: this restricts the retiling to single-anode-face clusters.
    // As such, it will likely freak out if fed clusters that have been
    // stitched across anode faces.  Since tiling is inherently a per-face
    // operation, this may be okay.
    /** Configuration "face" (optional, default is 0)

        The INDEX of the face in the anode's list of faces to use.
    */
    std::map<int, std::map<int, IAnodeFace::pointer>> m_face; // now apa/face --> m_face

    // Step 3. Form IBlobs from activities.
    std::vector<WireCell::IBlob::pointer> make_iblobs(std::map<std::pair<int, int>, std::vector<WireCell::RayGrid::measure_t> >& map_slices_measures, int apa, int face) const;

    // Step 1. Build activities from blobs in a cluster.
    void get_activity(const Cluster& cluster, std::map<std::pair<int, int>, std::vector<WireCell::RayGrid::measure_t> >& map_slices_measures, int apa, int face) const;


    // Step 2. Modify activity to suit.
    void hack_activity(const Cluster& cluster,
                       std::map<std::pair<int, int>, std::vector<WireCell::RayGrid::measure_t> >& map_slices_measures,
                       const std::vector<size_t>& path_wcps,
                       int apa, int face) const;

private:

   
    std::set<const Blob*> remove_bad_blobs(const Cluster& cluster, Cluster& shad_cluster, int tick_span, int apa, int face) const;


    // Remaining steps are done in the operator() directly.

   

    /** Configuration "cut_time_low" (optional, default is -1e9)
        Lower bound for time cut in nanoseconds
    */
    double m_cut_time_low;

    /** Configuration "cut_time_high" (optional, default is 1e9)
        Upper bound for time cut in nanoseconds
    */
    double m_cut_time_high;



    /** Configuration "anode" (required)

        The type/name of the anode.
    */


    // Wrap up getting the shortest path for the cluster high/low points.
    const std::vector<size_t>& cluster_path_wcps(const Cluster* cluster) const {
        // find the highest and lowest points
        std::pair<geo_point_t, geo_point_t> pair_points = cluster->get_highest_lowest_points();
        // std::cerr << "retile: hilo: " << pair_points.first << " " << pair_points.second << std::endl;
        int high_idx = cluster->get_closest_point_index(pair_points.first);
        int low_idx = cluster->get_closest_point_index(pair_points.second);
        return cluster->graph_algorithms().shortest_path(high_idx, low_idx);
    }
};                              // RetileCluster

}

#endif // WIRECELLCLUS_RETILE_CLUSTER_H