// ImproveCluster_2 - Second level cluster improvement 
//
// This class inherits from ImproveCluster_1 and provides additional
// cluster improvement functionality, building upon the Steiner tree
// enhancements from the first level.

#include "improvecluster_1.h"  // Include the ImproveCluster_1 header
#include "SteinerGrapher.h"
#include "WireCellUtil/NamedFactory.h"
#include <chrono>

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
                 WireCell::INamed, WireCell::IConfigurable, WireCell::IPCTreeMutate)

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
        using Clock = std::chrono::steady_clock;
        using MS = std::chrono::duration<double, std::milli>;
        auto t_mutate_start = Clock::now();
        auto t0 = Clock::now();

        // get the original cluster
        auto* orig_cluster = reinitialize(node);
        SPDLOG_LOGGER_DEBUG(log, "timing: reinitialize took {} ms", MS(Clock::now()-t0).count());

        // First: get the shortest path from the original cluster
        // Create a SteinerGrapher instance with the cluster
        // You'll need to provide appropriate configuration
        Steiner::Grapher::Config grapher_config;
        // Configure as needed - you may need to access member variables
        // that provide detector volumes and point cloud transform sets
        grapher_config.dv = m_dv;     // From NeedDV mixin
        grapher_config.pcts = m_pcts; // From NeedPCTS mixin
        grapher_config.perf = m_verbose; // toggle timing printouts with verbose=true
        
        // Create the Steiner::Grapher instance
        t0 = Clock::now();
        Steiner::Grapher orig_steiner_grapher(*orig_cluster, grapher_config, log);
        auto& orig_graph = orig_steiner_grapher.get_graph("basic_pid"); // this is good for the original cluster
        SPDLOG_LOGGER_DEBUG(log, "timing: get_graph(basic_pid) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_DEBUG(log, "Orig Graph vertices: {}, edges: {}", boost::num_vertices(orig_graph), boost::num_edges(orig_graph));

        // Establish same blob steiner edges
        t0 = Clock::now();
        orig_steiner_grapher.establish_same_blob_steiner_edges("basic_pid", true);
        std::vector<size_t> orig_path_point_indices;
        {  
            auto pair_points = orig_cluster->get_two_boundary_wcps();
            auto first_index  =   orig_cluster->get_closest_point_index(pair_points.first);
            auto second_index =   orig_cluster->get_closest_point_index(pair_points.second);
            orig_path_point_indices = orig_cluster->graph_algorithms("basic_pid").shortest_path(first_index, second_index);
        }
        SPDLOG_LOGGER_DEBUG(log, "timing: establish+shortest_path(basic_pid) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_DEBUG(log, "Orig shortest path indices: {} ; Graph vertices: {}, edges: {}", orig_path_point_indices.size(), boost::num_vertices(orig_graph), boost::num_edges(orig_graph));
        
        t0 = Clock::now();
        orig_steiner_grapher.remove_same_blob_steiner_edges("basic_pid");
        SPDLOG_LOGGER_DEBUG(log, "timing: remove_same_blob_steiner_edges(basic_pid) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_DEBUG(log, "Orig Graph vertices: {}, edges: {}", boost::num_vertices(orig_graph), boost::num_edges(orig_graph));



        // Second, make a temp_cluster based on the original cluster via ImproveCluster_1
        SPDLOG_LOGGER_DEBUG(log, "Grouping {} {}", m_grouping->get_name(), m_grouping->children().size());

        t0 = Clock::now();
        auto temp_node = ImproveCluster_1::mutate(node);
        auto temp_cluster_1 = temp_node->value.facade<Cluster>();
        auto& temp_cluster = m_grouping->make_child();
        temp_cluster.take_children(*temp_cluster_1);  // Move all blobs from improved cluster
        temp_cluster.from(*orig_cluster);
        SPDLOG_LOGGER_DEBUG(log, "timing: ImproveCluster_1::mutate took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_DEBUG(log, "Grouping {} {}", m_grouping->get_name(), m_grouping->children().size());

        t0 = Clock::now();
        Steiner::Grapher temp_steiner_grapher(temp_cluster, grapher_config, log);
        
        // this requires CTPC and ref_point cloud of original cluster
        auto& temp_graph = temp_cluster.find_graph("ctpc_ref_pid", *orig_cluster, m_dv, m_pcts);
        SPDLOG_LOGGER_DEBUG(log, "timing: Grapher+find_graph(ctpc_ref_pid) took {} ms", MS(Clock::now()-t0).count());
        //temp_steiner_grapher.get_graph("basic_pid"); 
        

        SPDLOG_LOGGER_DEBUG(log, "Temp Graph vertices: {}, edges: {}", boost::num_vertices(temp_graph), boost::num_edges(temp_graph));
        t0 = Clock::now();
        temp_steiner_grapher.establish_same_blob_steiner_edges("ctpc_ref_pid", false);
        std::vector<size_t> temp_path_point_indices;
        {  
            auto pair_points = temp_cluster.get_two_boundary_wcps();
            auto first_index  =   temp_cluster.get_closest_point_index(pair_points.first);
            auto second_index =   temp_cluster.get_closest_point_index(pair_points.second);
            temp_path_point_indices = temp_cluster.graph_algorithms("ctpc_ref_pid").shortest_path(first_index, second_index);
        }
        SPDLOG_LOGGER_DEBUG(log, "timing: establish+shortest_path(ctpc_ref_pid) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_DEBUG(log, "Temp shortest path indices: {} ; Graph vertices: {}, edges: {}", temp_path_point_indices.size(), boost::num_vertices(temp_graph), boost::num_edges(temp_graph));
        t0 = Clock::now();
        temp_steiner_grapher.remove_same_blob_steiner_edges("ctpc_ref_pid");
        SPDLOG_LOGGER_DEBUG(log, "timing: remove_same_blob_steiner_edges(ctpc_ref_pid) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_DEBUG(log, "Temp Graph vertices: {}, edges: {}", boost::num_vertices(temp_graph), boost::num_edges(temp_graph));


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
            t0 = Clock::now();
            get_activity_improved(*orig_cluster, map_slices_measures, apa, face);
            SPDLOG_LOGGER_DEBUG(log, "timing: get_activity_improved (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());

            // hack activity according to original cluster
            t0 = Clock::now();
            hack_activity_improved(*orig_cluster, map_slices_measures, orig_path_point_indices, apa, face); // may need more args
            SPDLOG_LOGGER_DEBUG(log, "timing: hack_activity(orig) (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());

            // hack activities according to the new cluster
            t0 = Clock::now();
            hack_activity_improved(temp_cluster, map_slices_measures, temp_path_point_indices, apa, face); // may need more args
            SPDLOG_LOGGER_DEBUG(log, "timing: hack_activity(temp) (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());

            // Step 3.
            t0 = Clock::now();
            auto iblobs = make_iblobs_improved(map_slices_measures, apa, face);
            SPDLOG_LOGGER_DEBUG(log, "timing: make_iblobs_improved (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
            SPDLOG_LOGGER_DEBUG(log, "new cluster {} iblobs for apa {} face {}", iblobs.size(), apa, face);

            auto niblobs = iblobs.size();
            // start to sampling points 
            int npoints = 0;
            t0 = Clock::now();
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
                    SPDLOG_LOGGER_DEBUG(log, "skipping blob {} with no points", iblob->ident());
                    continue;
                }
                new_cluster.node()->insert(Tree::Points(std::move(pcs)));
            }
            SPDLOG_LOGGER_DEBUG(log, "timing: sample_live loop (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
            SPDLOG_LOGGER_DEBUG(log, "{} points sampled for apa {} face {} Blobs {}", npoints, apa, face, niblobs);

            // remove bad blobs
            t0 = Clock::now();
            int tick_span = map_slices_measures.begin()->first.second -  map_slices_measures.begin()->first.first;
            auto blobs_to_remove = remove_bad_blobs(*orig_cluster, new_cluster, tick_span, apa, face);
            for (const Blob* blob : blobs_to_remove) {
                Blob& b = const_cast<Blob&>(*blob);
                new_cluster.remove_child(b);
            }
            SPDLOG_LOGGER_DEBUG(log, "timing: remove_bad_blobs (apa={},face={}) took {} ms", apa, face, MS(Clock::now()-t0).count());
            SPDLOG_LOGGER_DEBUG(log, "{} blobs removed for apa {} face {} remaining {}", blobs_to_remove.size(), apa, face, new_cluster.children().size());
        }


        // Remove this cluster from the grouping
        t0 = Clock::now();
        auto* temp_cluster_ptr = &temp_cluster;
        m_grouping->destroy_child(temp_cluster_ptr, true);
        SPDLOG_LOGGER_DEBUG(log, "timing: destroy_child(temp_cluster) took {} ms", MS(Clock::now()-t0).count());
        SPDLOG_LOGGER_DEBUG(log, "Grouping {} {}", m_grouping->get_name(), m_grouping->children().size());

        auto& default_scope = orig_cluster->get_default_scope();
        auto& raw_scope = orig_cluster->get_raw_scope();

        SPDLOG_LOGGER_DEBUG(log, "Scope: {} {}", default_scope.hash(), raw_scope.hash());
        if (default_scope.hash()!=raw_scope.hash()){
            t0 = Clock::now();
            auto correction_name = orig_cluster->get_scope_transform(default_scope);
            // std::vector<int> filter_results = c
            new_cluster.add_corrected_points(m_pcts, correction_name);
            // Get the new scope with corrected points
            // const auto& correction_scope = new_cluster.get_scope(correction_name);
            // Set this as the default scope for viewing
            new_cluster.from(*orig_cluster); // copy state from original cluster
            SPDLOG_LOGGER_DEBUG(log, "timing: add_corrected_points took {} ms", MS(Clock::now()-t0).count());
            // std::cout << "Test: Same:" << default_scope.hash() << " " << raw_scope.hash() << std::endl; 
        }

        // auto retiled_node = new_cluster.node();

        SPDLOG_LOGGER_DEBUG(log, "timing: mutate() TOTAL took {} ms", MS(Clock::now()-t_mutate_start).count());
        return m_grouping->remove_child(new_cluster);

    }

} // namespace WireCell::Clus
