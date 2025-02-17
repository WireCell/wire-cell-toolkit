// This defines a class to be called by ClusteringFuncs.
//
// We do not provide a function like others but instead code to the
// ClusteringFuncs API because some info like RayGrid::Coordinates is needed for
// each call.
//
// When called it  "retiles the cluster" using this sequence.
//
// 1) constructs layers of "activity" from input grouping.
// 2) applies "hacks" to the activity.
// 3) runs WCT tiling to create blobs.
// 4) runs blobs sampling to make point clouds
// 5) produces clusters such that the new blobs formed from an old cluster form a new "shadow" cluster.
// 6) forms a PC-tree
// 7) outputs the new grouping
//
// Developer notes, delete when no longer relevant:
// 
// Unlike usual ..., the many FIXME's here and in the .cxx can NOT be ignored.
// Remove them as they are fixed.
//
// FIXME: (2) above is omitted, to start with, making this a very complex type
// of "copy".  Xin will add the "hacks".
//
// FIXME: (7) is currently impossible until issue #377 provides a solution.  We
// start by making and then dropping the "shadow" grouping.
//


#ifndef WIRECELLCLUS_CLUSTERINGRETILE
#define WIRECELLCLUS_CLUSTERINGRETILE

#include "WireCellUtil/RayTiling.h"
#include "WireCellIface/IBlob.h"
#include "WireCellIface/IBlobSampler.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellClus/Facade_Grouping.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellUtil/RayHelpers.h"
#include "WireCellAux/PlaneTools.h"


#include <vector>

// This is not an appropriate name but I propagate it because that is what is
// used in ClusteringFuncs.h.  Better that it be WireCell::Clus.
namespace WireCell::PointCloud::Facade {

    // There is a hard-wired factory method in ClusteringFuncs to which this
    // class is added.
    class ClusteringRetile {
    public:
        ClusteringRetile(const WireCell::Configuration& config);
        
        // fixme: we can not satisfy this type by including CluteringFuncs.h
        // because that header must include this one.  A refactoring would ease
        // this problem.
        using cluster_set_t = std::set<const Cluster*>;

        // FIXME: #377.
        void operator()(Grouping& live_clusters, Grouping& dead_clusters,
                        cluster_set_t& cluster_connected_dead) const;

    private:

        // Step 1. Build activities from blobs in a cluster.
        void get_activity(const Cluster& cluster, std::map<std::pair<int, int>, std::vector<WireCell::RayGrid::measure_t> >& map_slices_measures) const;


        // Step 2. Modify activity to suit.
        void hack_activity(const Cluster& cluster, std::map<std::pair<int, int>, std::vector<WireCell::RayGrid::measure_t> >& map_slices_measures) const;

        // Step 3. Form IBlobs from activities.
        std::vector<WireCell::IBlob::pointer> make_iblobs(std::map<std::pair<int, int>, std::vector<WireCell::RayGrid::measure_t> >& map_slices_measures) const;

        std::set<const WireCell::PointCloud::Facade::Blob*> remove_bad_blobs(const Cluster& cluster, Cluster& shad_cluster, int tick_span) const;

        // Remaining steps are done in the operator() directly.

        /** Configuration: "sampler" (required)

            The type/name an IBlobSampler for producing the "3d" point cloud.

            If not given, the retailed blob tree nodes will not have point clouds.
        */
        WireCell::IBlobSampler::pointer m_sampler;

        // fixme: this restricts ClusteringRetile to single-anode-face clusters.
        // As such, it will likely freak out if fed clusters that have been
        // stitched across anode faces.  Since tiling is inherently a per-face
        // operation, this may be okay.
        /** Configuration "face" (optional, default is 0)

            The INDEX of the face in the anode's list of faces to use.
        */
        IAnodeFace::pointer m_face;

        /** Configuration "cut_time_low" (optional, default is -1e9)
            Lower bound for time cut in nanoseconds
        */
        double m_cut_time_low;

        /** Configuration "cut_time_high" (optional, default is 1e9)
            Upper bound for time cut in nanoseconds
        */
        double m_cut_time_high;

        std::vector<Aux::WirePlaneInfo> m_plane_infos;

        /** Configuration "anode" (required)

            The type/name of the anode.
        */
        // nothing to store.

        // Process original and shadow groupings to create a mapping between clusters
        std::map<Cluster*, std::tuple<Cluster*, int, Cluster*>> process_groupings(Grouping& original, Grouping& shadow, const std::string& aname = "isolated", const std::string& pname = "perblob") const;
    };
}


#endif
