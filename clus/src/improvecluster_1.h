// ImproveCluster_1 - First level cluster improvement using Steiner tree methods
//
// This class inherits from RetileCluster and provides enhanced cluster
// improvement functionality by incorporating Steiner tree algorithms
// from the Wire-Cell Prototype.

#ifndef WIRECELLCLUS_IMPROVE_CLUSTER_1_H
#define WIRECELLCLUS_IMPROVE_CLUSTER_1_H

#include "retile_cluster.h"  // Include the RetileCluster header
#include "WireCellClus/BadBlobRuns.h"

#include "WireCellAux/Logger.h"
#include "WireCellUtil/NamedFactory.h"

#include <vector>

namespace WireCell::Clus {

    using namespace WireCell;
    using namespace WireCell::Clus;
    using namespace WireCell::Clus::Facade;
    using namespace WireCell::PointCloud::Tree;

    class ImproveCluster_1 : public RetileCluster, public Aux::Logger {

    public:

        ImproveCluster_1();
        virtual ~ImproveCluster_1();

        // IConfigurable API - extend the base configuration
        void configure(const WireCell::Configuration& config) override;
        virtual Configuration default_configuration() const override;

        // IPCTreeMutate API - override to add Steiner tree improvements
        virtual std::unique_ptr<node_t> mutate(node_t& node) const override;

    protected:
       void get_activity_improved(const Cluster& cluster, std::map<std::pair<int, int>,std::vector<WireCell::RayGrid::measure_t>>& map_slices_measures, int apa, int face) const;

       void hack_activity_improved(const Cluster& cluster, std::map<std::pair<int, int>, std::vector<WireCell::RayGrid::measure_t> >& map_slices_measures, const std::vector<size_t>& path_wcps, int apa, int face) const;

       std::vector<WireCell::IBlob::pointer> make_iblobs_improved(std::map<std::pair<int, int>, std::vector<WireCell::RayGrid::measure_t> >& map_slices_measures, int apa, int face) const;


       std::vector<const Blob*> remove_bad_blobs(const Cluster& cluster, Cluster& shad_cluster, int tick_span, int apa, int face) const;

       // doc pdvd/40 round 3: the knob-ON / census path of remove_bad_blobs
       // (BadBlobRuns.h holds the decision core).  all_new_blobs is the
       // deterministically ordered blob list remove_bad_blobs already builds.
       std::vector<const Blob*> remove_bad_blobs_runs(const Cluster& cluster,
                                                      const std::vector<const Blob*>& all_new_blobs,
                                                      const std::map<int, BlobSet>& orig_time_blob_map,
                                                      const std::map<int, BlobSet>& new_time_blob_map,
                                                      int tick_span, int apa, int face) const;

       // doc pdvd/40 round 3.  The retile fabricates blobs along a whole-cluster
       // shortest path and remove_bad_blobs is its only anti-ghost filter.  The
       // historical filter (a faithful port of the prototype, kept byte-for-byte
       // when both knobs are off) votes per connected component, by that
       // component's FIRST blob, and only when there is more than one component.
       //
       // bad_blob_max_run (length; C++ default 0 = OFF): when > 0 every new blob
       // is tested for support (overlap with an ORIGINAL blob within one slice),
       // a component is kept iff ANY blob is supported, and a connected run of
       // UNSUPPORTED blobs inside a kept component whose bounding-box diagonal
       // exceeds this length is removed whole.  Adjacency gains same-slice
       // overlap so a column at fixed drift time is one run, not N singletons.
       // bad_blob_report (C++ default false): log-only census, one DEBUG line
       // "BADBLOB ..." per (cluster, apa, face), independent of the bound.
       double m_bad_blob_max_run{0.0};
       bool m_bad_blob_report{false};

    private:
 
       
    };

}
#endif // WIRECELLCLUS_IMPROVE_CLUSTER_1_H