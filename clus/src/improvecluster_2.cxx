// ImproveCluster_2 - Second level cluster improvement 
//
// This class inherits from ImproveCluster_1 and provides additional
// cluster improvement functionality, building upon the Steiner tree
// enhancements from the first level.

#include "improvecluster_1.h"  // Include the ImproveCluster_1 header

#include "WireCellUtil/NamedFactory.h"

#include <vector>

namespace WireCell::Clus {

    class ImproveCluster_2 : public ImproveCluster_1 {

    public:

        ImproveCluster_2();
        virtual ~ImproveCluster_2();

        // IConfigurable API - extend the base configuration
        void configure(const WireCell::Configuration& config) override;
        virtual Configuration default_configuration() const override;

        // IPCTreeMutate API - override to add second level improvements
        virtual std::unique_ptr<node_t> mutate(node_t& node) const override;

    private:

       

    };

} // namespace WireCell::Clus

WIRECELL_FACTORY(ImproveCluster_2, WireCell::Clus::ImproveCluster_2,
                 WireCell::IConfigurable, WireCell::IPCTreeMutate)

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using namespace WireCell::PointCloud::Tree;

namespace WireCell::Clus {

    ImproveCluster_2::ImproveCluster_2() 
    {
    }

    ImproveCluster_2::~ImproveCluster_2() 
    {
    }

    void ImproveCluster_2::configure(const WireCell::Configuration& cfg)
    {
        // Configure base class first
        ImproveCluster_1::configure(cfg);
        
  
    }

    Configuration ImproveCluster_2::default_configuration() const
    {
        Configuration cfg = ImproveCluster_1::default_configuration();
        
    
        return cfg;
    }

    std::unique_ptr<ImproveCluster_2::node_t> ImproveCluster_2::mutate(node_t& node) const
    {
        // First, apply the base class improvements (ImproveCluster_1) 
        auto improved_node = ImproveCluster_1::mutate(node);

        if (!improved_node) {
            return nullptr;
        }

        // Get the cluster from the improved node
        auto cluster = improved_node->value.facade<Cluster>();

        std::cout << m_grouping->get_name() << " " << m_grouping->children().size() << std::endl;

        // Move the improved cluster to m_grouping (reparent it)
        auto& new_cluster = m_grouping->make_child();
        new_cluster.take_children(*cluster);  // Move all blobs from improved cluster
        new_cluster.from(*cluster);       
        std::cout << "Test: " << m_grouping->children().size() << " " << new_cluster.grouping()->get_name() << std::endl;

        // Remove this cluster from the grouping
        auto* cluster_ptr = &new_cluster;
        m_grouping->destroy_child(cluster_ptr, true);
        std::cout << "Test: " << m_grouping->children().size() << " " << std::endl;

        std::cout << "Improving cluster with advanced Steiner tree methods..." << std::endl;

        // TODO: Add additional cluster improvements here
        // Since we destroyed the original cluster, we might need to:
        // 1. Create a new improved cluster, or
        // 2. Return nullptr if the cluster should be removed
        
        // For now, return nullptr since the cluster has been destroyed
        return nullptr;
    }

} // namespace WireCell::Clus
