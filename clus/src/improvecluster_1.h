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

#include <set>
#include <tuple>
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

       // doc pdhd/08 stage 1.  bridge_cells is the set of (time_slice, plane,
       // wire) cells a bridge has painted.  It is owned by the caller so the two
       // hack_activity_improved calls of ImproveCluster_2::mutate share one set
       // (no mutable member, no data race), which is what lets the census tell a
       // "cell already occupied" skip caused by an EARLIER BRIDGE from one caused
       // by real charge or by an admitted dead channel.  nullptr = no accounting.
       // `which` only labels the BRIDGE report line.  Both are report-only.
       using bridge_cells_t = std::set<std::tuple<int, int, int>>;

       void hack_activity_improved(const Cluster& cluster, std::map<std::pair<int, int>, std::vector<WireCell::RayGrid::measure_t> >& map_slices_measures, const std::vector<size_t>& path_wcps, int apa, int face, const char* which = "", bridge_cells_t* bridge_cells = nullptr) const;

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

       // doc pdhd/08.  Two default-OFF knobs against the retiler's fabrication,
       // measured on PDHD 029107 events 991 and 1079 (doc pdhd/08 sec 2, 3).
       //
       // hack_max_bridge (length; C++ default 0 = uncapped = the prototype,
       // ImprovePR3DCluster.cxx:973-990): hack_activity_improved fills each gap
       // between consecutive shortest-path points with int(gap/0.3 cm)+1 invented
       // points, each painting a 7x7 (dt,dw) disc in all three planes, with no
       // length limit -- the longest single gap measured is 480 cm.  Above this
       // length the gap is left unfilled (both real endpoints are kept).  The
       // census showed 0.3-0.7 % of bridges, those over 20 cm, create 49-62 % of
       // every cell the retiler invents, and that only 0.3-0.8 % of the cells a
       // bridge touches were already covered by an admitted dead channel -- so a
       // length cap is well targeted and does not fight dead-region bridging.
       //
       // bad_blob_run_merge (length; C++ default 0 = OFF): before applying
       // bad_blob_max_run, transitively merge the runs of one component whose
       // bounding boxes lie within this distance, and bound the merged group
       // instead (BadBlobRuns::analyze step 4).  This closes the fragmentation
       // escape of doc pdhd/07: a 127.9 cm drift column cut into 100.8 / 18.5 /
       // 3.2 / 0.0 cm runs by 2.2 cm gaps, of which only the first exceeded the
       // 20 cm bound.  Requires bad_blob_max_run > 0 (it feeds that bound).
       // Sensitive: the extra removal runs 12 % -> 74 % of survivors as this
       // goes 1 -> 8 cm, and at 5 cm it demonstrably fuses two objects.
       double m_hack_max_bridge{0.0};
       double m_bad_blob_run_merge{0.0};

    private:
 
       
    };

}
#endif // WIRECELLCLUS_IMPROVE_CLUSTER_1_H