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

    private:

       
    };

}
#endif // WIRECELLCLUS_IMPROVE_CLUSTER_1_H