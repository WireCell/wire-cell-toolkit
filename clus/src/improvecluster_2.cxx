// ImproveCluster_2 - Second level cluster improvement 
//
// This class inherits from ImproveCluster_1 and provides additional
// cluster improvement functionality, building upon the Steiner tree
// enhancements from the first level.

#include "improvecluster_1.h"  // Include the ImproveCluster_1 header
#include "SteinerGrapher.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Logging.h"

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

// Segregate this weird choice for namespace.
namespace WCF = WireCell::Clus::Facade;

// Nick name for less typing.
namespace WRG = WireCell::RayGrid;
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
        // get the original cluster
        auto* orig_cluster = reinitialize(node);

        // First: get the shortest path from the original cluster
        // Create a SteinerGrapher instance with the cluster
        // You'll need to provide appropriate configuration
        Steiner::Grapher::Config grapher_config;
        // Configure as needed - you may need to access member variables
        // that provide detector volumes and point cloud transform sets
        grapher_config.dv = m_dv;     // From NeedDV mixin
        grapher_config.pcts = m_pcts; // From NeedPCTS mixin
        
        // Create the Steiner::Grapher instance
        auto log = Log::logger("improve_cluster_2 mutate");
        Steiner::Grapher orig_steiner_grapher(*orig_cluster, grapher_config, log);
        auto& orig_graph = orig_steiner_grapher.get_graph("basic_pid");

        if (m_verbose) std::cout << "ImproveCluster_2 " << " Orig Graph vertices: " << boost::num_vertices(orig_graph) << ", edges: " << boost::num_edges(orig_graph) << std::endl;

        // Establish same blob steiner edges
        orig_steiner_grapher.establish_same_blob_steiner_edges("basic_pid", true);
        std::vector<size_t> orig_path_point_indices;
        {  
            auto pair_points = orig_cluster->get_two_boundary_wcps();
            auto first_index  =   orig_cluster->get_closest_point_index(pair_points.first);
            auto second_index =   orig_cluster->get_closest_point_index(pair_points.second);
            orig_path_point_indices = orig_cluster->graph_algorithms("basic_pid").shortest_path(first_index, second_index);
        }
        if (m_verbose) std::cout << "ImproveCluster_2 " << " Origi Shortest path indices: " << orig_path_point_indices.size() << " ; Graph vertices: " << boost::num_vertices(orig_graph) << ", edges: " << boost::num_edges(orig_graph)<< std::endl;
        
        orig_steiner_grapher.remove_same_blob_steiner_edges("basic_pid");
        if (m_verbose) std::cout << "ImproveCluster_2 " << " Orig Graph vertices: " << boost::num_vertices(orig_graph) << ", edges: " << boost::num_edges(orig_graph) << std::endl;



        // Second, make a temp_cluster based on the original cluster via ImproveCluster_1
        if (m_verbose) std::cout << "ImproveCluster_2: Grouping" << m_grouping->get_name() << " " << m_grouping->children().size() << std::endl;

        auto temp_node = ImproveCluster_1::mutate(node);
        auto temp_cluster_1 = temp_node->value.facade<Cluster>();
        auto& temp_cluster = m_grouping->make_child();
        temp_cluster.take_children(*temp_cluster_1);  // Move all blobs from improved cluster
        temp_cluster.from(*orig_cluster);
        if (m_verbose) std::cout << "ImproveCluster_2: Grouping" << m_grouping->get_name() << " " << m_grouping->children().size() << std::endl;

        Steiner::Grapher temp_steiner_grapher(temp_cluster, grapher_config, log);
        auto& temp_graph = temp_steiner_grapher.get_graph("basic_pid");
        if (m_verbose) std::cout << "ImproveCluster_2 " << " Temp Graph vertices: " << boost::num_vertices(temp_graph) << ", edges: " << boost::num_edges(temp_graph) << std::endl;
        temp_steiner_grapher.establish_same_blob_steiner_edges("basic_pid", false);
        std::vector<size_t> temp_path_point_indices;
        {  
            auto pair_points = temp_cluster.get_two_boundary_wcps();
            auto first_index  =   temp_cluster.get_closest_point_index(pair_points.first);
            auto second_index =   temp_cluster.get_closest_point_index(pair_points.second);
            temp_path_point_indices = temp_cluster.graph_algorithms("basic_pid").shortest_path(first_index, second_index);
        }
        if (m_verbose) std::cout << "ImproveCluster_2 " << " Temp Shortest path indices: " << temp_path_point_indices.size() << " ; Graph vertices: " << boost::num_vertices(temp_graph) << ", edges: " << boost::num_edges(temp_graph)<< std::endl;
        temp_steiner_grapher.remove_same_blob_steiner_edges("basic_pid");
        if (m_verbose) std::cout << "ImproveCluster_2 " << " Temp Graph vertices: " << boost::num_vertices(temp_graph) << ", edges: " << boost::num_edges(temp_graph) << std::endl;


        // star to construct a new cluster
        auto wpids = orig_cluster->wpids_blob();
        std::set<WirePlaneId> wpid_set(wpids.begin(), wpids.end());

        // make a new node from the existing grouping
        auto& new_cluster = m_grouping->make_child(); // make a new cluster inside 

        for (auto it = wpid_set.begin(); it != wpid_set.end(); ++it) {
            int apa = it->apa();
            int face = it->face();
            const auto& angles = m_wpid_angles.at(*it);

            std::map<std::pair<int, int>, std::vector<WRG::measure_t> > map_slices_measures;
            
            // get original activities ...
            get_activity_improved(*orig_cluster, map_slices_measures, apa, face);

            // hack activity according to original cluster
            hack_activity_improved(*orig_cluster, map_slices_measures, orig_path_point_indices, apa, face); // may need more args

            // hack activities according to the new cluster
            hack_activity_improved(temp_cluster, map_slices_measures, temp_path_point_indices, apa, face); // may need more args

            // Step 3.
            auto iblobs = make_iblobs(map_slices_measures, apa, face);

            if (m_verbose) std::cout << "ImproveCluster_2: new cluster " << iblobs.size() << " iblobs for apa " << apa << " face " << face << std::endl;

            auto niblobs = iblobs.size();
            // start to sampling points 
            int npoints = 0;
            for (size_t bind=0; bind<niblobs; ++bind) {
          
                const IBlob::pointer iblob = iblobs[bind];
                auto sampler = m_samplers.at(apa).at(face);
                const double tick = m_grouping->get_tick().at(apa).at(face);

                auto pcs = Aux::sample_live(sampler, iblob, angles, tick, bind);
                // DO NOT EXTEND FURTHER! see #426, #430

                if (pcs["3d"].size()==0) continue; // no points ...
                // Access 3D coordinates
                auto pc3d = pcs["3d"];  // Get the 3D point cloud dataset
                auto x_coords = pc3d.get("x")->elements<double>();  // Get X coordinates
                // auto y_coords = pc3d.get("y")->elements<double>();  // Get Y coordinates  
                // auto z_coords = pc3d.get("z")->elements<double>();  // Get Z coordinates
                // auto ucharge_val = pc3d.get("ucharge_val")->elements<double>();  // Get U charge
                // auto vcharge_val = pc3d.get("vcharge_val")->elements<double>();  // Get V charge
                // auto wcharge_val = pc3d.get("wcharge_val")->elements<double>();  // Get W charge
                // auto ucharge_err = pc3d.get("ucharge_unc")->elements<double>();  // Get U charge error
                // auto vcharge_err = pc3d.get("vcharge_unc")->elements<double>();  // Get V charge error
                // auto wcharge_err = pc3d.get("wcharge_unc")->elements<double>();  // Get W charge error

                // std::cout << "ImproveCluster_1 PCS: " << pcs.size() << " " 
                //           << pcs["3d"].size() << " " 
                //           << x_coords.size() << std::endl;     

                npoints +=x_coords.size();
                if (pcs.empty()) {
                    SPDLOG_DEBUG("ImproveCluster_1: skipping blob {} with no points", iblob->ident());
                    continue;
                }
                new_cluster.node()->insert(Tree::Points(std::move(pcs)));
            }
            if (m_verbose) std::cout << "ImproveCluster_2: " << npoints << " points sampled for apa " << apa << " face " << face << " Blobs " << niblobs << std::endl;

            // remove bad blobs
            int tick_span = map_slices_measures.begin()->first.second -  map_slices_measures.begin()->first.first;
            auto blobs_to_remove = remove_bad_blobs(*orig_cluster, new_cluster, tick_span, apa, face);
            for (const Blob* blob : blobs_to_remove) {
                Blob& b = const_cast<Blob&>(*blob);
                new_cluster.remove_child(b);
            }
            if (m_verbose) std::cout << "ImproveCluster_2: " << blobs_to_remove.size() << " blobs removed for apa " << apa << " face " << face << " " << new_cluster.children().size() << std::endl;
        }


        // Remove this cluster from the grouping
        auto* temp_cluster_ptr = &temp_cluster;
        m_grouping->destroy_child(temp_cluster_ptr, true);
        if (m_verbose) std::cout << "ImproveCluster_2: Grouping" << m_grouping->get_name() << " " << m_grouping->children().size() << std::endl;

        auto& default_scope = orig_cluster->get_default_scope();
        auto& raw_scope = orig_cluster->get_raw_scope();

        if (m_verbose) std::cout << "ImproveCluster_1: Scope: " << default_scope.hash() << " " << raw_scope.hash() << std::endl;
        if (default_scope.hash()!=raw_scope.hash()){
            auto correction_name = orig_cluster->get_scope_transform(default_scope);
            // std::vector<int> filter_results = c
            new_cluster.add_corrected_points(m_pcts, correction_name);
            // Get the new scope with corrected points
            const auto& correction_scope = new_cluster.get_scope(correction_name);
            // Set this as the default scope for viewing
            new_cluster.from(*orig_cluster); // copy state from original cluster
            // std::cout << "Test: Same:" << default_scope.hash() << " " << raw_scope.hash() << std::endl; 
        }

        // auto retiled_node = new_cluster.node();


        return m_grouping->remove_child(new_cluster);

    }

} // namespace WireCell::Clus
