// ImproveCluster_1 - First level cluster improvement using Steiner tree methods
//
// This class inherits from RetileCluster and provides enhanced cluster
// improvement functionality by incorporating Steiner tree algorithms
// from the Wire-Cell Prototype.

#ifndef WIRECELLCLUS_IMPROVE_CLUSTER_1_H
#define WIRECELLCLUS_IMPROVE_CLUSTER_1_H

#include "retile_cluster.h"  // Include the RetileCluster header

#include "WireCellUtil/NamedFactory.h"

#include <vector>

namespace WireCell::Clus {

    using namespace WireCell;
    using namespace WireCell::Clus;
    using namespace WireCell::Clus::Facade;
    using namespace WireCell::PointCloud::Tree;

    class ImproveCluster_1 : public RetileCluster {

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


       std::set<const Blob*> remove_bad_blobs(const Cluster& cluster, Cluster& shad_cluster, int tick_span, int apa, int face) const;

    private:
 
       
    };

}
#endif // WIRECELLCLUS_IMPROVE_CLUSTER_1_H