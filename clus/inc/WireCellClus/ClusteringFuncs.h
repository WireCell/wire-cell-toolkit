#include "WireCellClus/MultiAlgBlobClustering.h"
#include "WireCellClus/Facade.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/Persist.h"
#include "WireCellAux/TensorDMpointtree.h"
#include "WireCellAux/TensorDMdataset.h"
#include "WireCellAux/TensorDMcommon.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/Graph.h"


#include <fstream>

namespace WireCell::PointCloud::Facade {
    using namespace WireCell::PointCloud::Tree;

    using cluster_set_t = std::set<const Cluster*>;

    using cluster_connectivity_graph_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, int>;

    using cluster_vector_t = std::vector<Cluster*>;


    // clustering_util.cxx
    // merging clustering function
    void merge_clusters(cluster_connectivity_graph_t& g, // 
			Grouping& live_clusters,
			cluster_set_t& cluster_connected_dead); // in/out
    
    
    // clustering_live_dead.cxx
    // first function ...
    void clustering_live_dead(Grouping& live_clusters,
                              const Grouping& dead_clusters,
                              cluster_set_t& cluster_connected_dead,  // in/out
                              const int dead_live_overlap_offset      // specific params
    );
    class ClusteringLiveDead {
       public:
        ClusteringLiveDead(const WireCell::Configuration& config)
        {
            // FIXME: throw if not found?
            dead_live_overlap_offset_ = get(config, "dead_live_overlap_offset", 2);
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            clustering_live_dead(live_clusters, dead_clusters, cluster_connected_dead, dead_live_overlap_offset_);
            return {};
        }

       private:
        int dead_live_overlap_offset_{2};
    };

    // clustering_extend.cxx
    //helper function ..
    double Find_Closest_Points(const Cluster& cluster1,
			       const Cluster& cluster2,
			       double length_1,
			       double length_2,
			       double length_cut,
			       geo_point_t& p1_save, // output
			       geo_point_t& p2_save,  // output
                   bool flag_print = false
			       );
			       
    
    // clustering_extend.cxx
    // second function ...
    void clustering_extend(Grouping& live_clusters,
                           cluster_set_t& cluster_connected_dead,            // in/out
                           const int flag,                                                //
                           const double length_cut = 150*units::cm,                       //
                           const int num_try = 0,                                         //
                           const double length_2_cut = 3*units::cm,                       //
                           const int num_dead_try =3                                      //
			   );
    class ClusteringExtend {
       public:
        ClusteringExtend(const WireCell::Configuration& config)
        {
            // FIXME: throw if not found?
            flag_ = get(config, "flag", 0);
            length_cut_ = get(config, "length_cut", 150*units::cm);
            num_try_ = get(config, "num_try", 0);
            length_2_cut_ = get(config, "length_2_cut", 3*units::cm);
            num_dead_try_ = get(config, "num_dead_try", 3);
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            clustering_extend(live_clusters, cluster_connected_dead, flag_, length_cut_, num_try_, length_2_cut_, num_dead_try_);
            return {};
        }

       private:
        int flag_{0};
        double length_cut_{150*units::cm};
        int num_try_{0};
        double length_2_cut_{3*units::cm};
        int num_dead_try_{3};
    };
    class ClusteringExtendLoop {
       public:
        ClusteringExtendLoop(const WireCell::Configuration& config)
        {
            // FIXME: throw if not found?
            num_try_ = get(config, "num_try", 0);
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            // for very busy events do less ...
            int num_try = num_try_;
            if (live_clusters.nchildren() > 1100) num_try = 1;
            for (int i = 0; i != num_try; i++) {
                // deal with prolong case
                clustering_extend(live_clusters, cluster_connected_dead, 1, 150*units::cm, 0);
                // deal with parallel case
                clustering_extend(live_clusters, cluster_connected_dead, 2, 30*units::cm, 0);
                // extension regular case
                clustering_extend(live_clusters, cluster_connected_dead, 3, 15*units::cm, 0);
                // extension ones connected to dead region ...
                if (i == 0) {
                    clustering_extend(live_clusters, cluster_connected_dead, 4, 60 * units::cm, i);
                }
                else {
                    clustering_extend(live_clusters, cluster_connected_dead, 4, 35 * units::cm, i);
                }
            }
            return {};
        }

       private:
        int num_try_{0};
    };

    bool Clustering_4th_prol(const Cluster& cluster1,
			     const Cluster& cluster2,
			     double length_2,
			     geo_point_t& earliest_p,
			     geo_point_t& dir_earlp,
			     double length_cut);
    
    bool Clustering_4th_para(const Cluster& cluster1,
			     const Cluster& cluster2,
			     double length_1, double length_2,
			     geo_point_t& earliest_p,
			     geo_point_t& dir_earlp,
			     double length_cut);

    bool Clustering_4th_reg(const Cluster& cluster1,
			    const Cluster& cluster2,
			    double length_1, double length_2,
			    geo_point_t p1, double length_cut);

    bool Clustering_4th_dead(const Cluster& cluster1,
			     const Cluster& cluster2,
			     double length_1, double length_2, double length_cut, int num_dead_try=3);
      

    // clustering_regular.cxx
    // third function 
    void clustering_regular(Grouping& live_clusters,
                            cluster_set_t& cluster_connected_dead,            // in/out
                            const double length_cut = 45*units::cm,                                       //
                            bool flag_enable_extend = true                                       //
    );
    class ClusteringRegular {
       public:
        ClusteringRegular(const WireCell::Configuration& config)
        {
            // FIXME: throw if not found?
            length_cut_ = get(config, "length_cut", 45*units::cm);
            flag_enable_extend_ = get(config, "flag_enable_extend", true);
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            clustering_regular(live_clusters, cluster_connected_dead, length_cut_, flag_enable_extend_);
            return {};
        }

       private:
        double length_cut_{45*units::cm};
        bool flag_enable_extend_{true};
    };

    bool Clustering_1st_round(const Cluster& cluster1,
			      const Cluster& cluster2,
			      double length_1,
			      double length_2,
			      double length_cut = 45*units::cm,
			      bool flag_enable_extend = true);

    // clustering_parallel_prolong.cxx:
    void clustering_parallel_prolong(Grouping& live_clusters,
                                     cluster_set_t& cluster_connected_dead, // in/out
                                     const double length_cut = 35*units::cm
    );
    class ClusteringParallelProlong {
       public:
        ClusteringParallelProlong(const WireCell::Configuration& config)
        {
            // FIXME: throw if not found?
            length_cut_ = get(config, "length_cut", 35*units::cm);
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            clustering_parallel_prolong(live_clusters, cluster_connected_dead, length_cut_);
            return {};
        }

       private:
        double length_cut_{35*units::cm};
    };

    bool Clustering_2nd_round(const Cluster& cluster1,
			      const Cluster& cluster2,
			      double length_1,
			      double length_2,
			      double length_cut = 35*units::cm);
    
    // clustering_close.cxx
    void clustering_close(Grouping& live_clusters,           // 
                          cluster_set_t& cluster_connected_dead, // in/out
                          const double length_cut = 1*units::cm //
    );
    class ClusteringClose {
       public:
        ClusteringClose(const WireCell::Configuration& config)
        {
            // FIXME: throw if not found?
            length_cut_ = get(config, "length_cut", 1*units::cm);
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            clustering_close(live_clusters, cluster_connected_dead, length_cut_);
            return {};
        }

       private:
        double length_cut_{1*units::cm};
    };

    bool Clustering_3rd_round( const Cluster& cluster1,
			       const Cluster& cluster2,
			       double length_1,
			       double length_2,
			       double length_cut = 1*units::cm);

    // void clustering_separate(Grouping& live_grouping,
    //                          std::map<int, std::pair<double, double>>& dead_u_index,
    //                          std::map<int, std::pair<double, double>>& dead_v_index,
    //                          std::map<int, std::pair<double, double>>& dead_w_index);

    void clustering_separate(Grouping& live_grouping,
                             const bool use_ctpc);
    class ClusteringSeparate {
       public:
        ClusteringSeparate(const WireCell::Configuration& config)
        {
            // FIXME: throw if not found?
            use_ctpc_ = get(config, "use_ctpc", true);
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            clustering_separate(live_clusters, use_ctpc_);
            return {};
        }

       private:
        double use_ctpc_{true};
    };

    void clustering_connect1(Grouping& live_grouping);
    class ClusteringConnect1 {
       public:
        ClusteringConnect1(const WireCell::Configuration& config)
        {
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            clustering_connect1(live_clusters);
            return {};
        }

       private:
    };

    void clustering_deghost(Grouping& live_grouping,
                            const bool use_ctpc,
                            double length_cut = 0);
    class ClusteringDeGhost {
       public:
        ClusteringDeGhost(const WireCell::Configuration& config)
        {
            // FIXME: throw if not found?
            use_ctpc_ = get(config, "use_ctpc", true);
            length_cut_ = get(config, "length_cut", 0);
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            clustering_deghost(live_clusters, use_ctpc_, length_cut_);
            return {};
        }

       private:
        double use_ctpc_{true};
        double length_cut_{0};
    };

    void clustering_examine_x_boundary(Grouping& live_grouping);
    class ClusteringExamineXBoundary {
       public:
        ClusteringExamineXBoundary(const WireCell::Configuration& config)
        {
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            clustering_examine_x_boundary(live_clusters);
            return {};
        }

       private:
    };

    void clustering_protect_overclustering(Grouping& live_grouping);
    class ClusteringProtectOverClustering {
       public:
        ClusteringProtectOverClustering(const WireCell::Configuration& config)
        {
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            clustering_protect_overclustering(live_clusters);
            return {};
        }

       private:
    };

    void clustering_neutrino(Grouping &live_grouping, int num_try);
    class ClusteringNeutrino {
       public:
        ClusteringNeutrino(const WireCell::Configuration& config)
        {
            // FIXME: throw if not found?
            num_try_ = get(config, "num_try", 1);
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            for (int i = 0; i != num_try_; i++) {
                clustering_neutrino(live_clusters, i);
            }
            return {};
        }

       private:
        int num_try_{1};
    };

    map_cluster_cluster_vec clustering_isolated(Grouping& live_grouping);
    class ClusteringIsolated {
       public:
        ClusteringIsolated(const WireCell::Configuration& config)
        {
        }

        map_cluster_cluster_vec operator()(Grouping& live_clusters, Grouping& dead_clusters, cluster_set_t& cluster_connected_dead) const
        {
            return clustering_isolated(live_clusters);
        }

       private:
    };


    // time_slice_length is length span for a slice
    bool JudgeSeparateDec_1(const Cluster* cluster, const geo_point_t& drift_dir, const double length, const double time_slice_length);
    /// @attention contains hard-coded distance cuts
    /// @param boundary_points return the boundary points
    /// @param independent_points return the independent points
    bool JudgeSeparateDec_2(const Cluster* cluster, const geo_point_t& drift_dir,
                               std::vector<geo_point_t>& boundary_points, std::vector<geo_point_t>& independent_points,
                               const double cluster_length);
    




    std::vector<Cluster *> Separate_1(const bool use_ctpc, Cluster *cluster,
                                                         std::vector<geo_point_t> &boundary_points,
                                                         std::vector<geo_point_t> &independent_points,
                                                         std::map<int, std::pair<double, double>> &dead_u_index,
                                                         std::map<int, std::pair<double, double>> &dead_v_index,
                                                         std::map<int, std::pair<double, double>> &dead_w_index,
                                                         double length);

    std::vector<int> Separate_2(Cluster *cluster, const double dis_cut =  5*units::cm, const size_t ticks_per_slice = 4);


    inline std::function<map_cluster_cluster_vec(Grouping&, Grouping&, cluster_set_t&)> getClusteringFunction(const WireCell::Configuration& config) {
        std::string function_name = config["name"].asString();

        if (function_name == "clustering_live_dead") {
            return ClusteringLiveDead(config);
        }
        else if (function_name == "clustering_extend") {
            return ClusteringExtend(config);
        }
        else if (function_name == "clustering_regular") {
            return ClusteringRegular(config);
        }
        else if (function_name == "clustering_parallel_prolong") {
            return ClusteringParallelProlong(config);
        }
        else if (function_name == "clustering_close") {
            return ClusteringClose(config);
        }
        else if (function_name == "clustering_extend_loop") {
            return ClusteringExtendLoop(config);
        }
        else if (function_name == "clustering_separate") {
            return ClusteringSeparate(config);
        }
        else if (function_name == "clustering_connect1") {
            return ClusteringConnect1(config);
        }
        else if (function_name == "clustering_deghost") {
            return ClusteringDeGhost(config);
        }
        else if (function_name == "clustering_examine_x_boundary") {
            return ClusteringExamineXBoundary(config);
        }
        else if (function_name == "clustering_protect_overclustering") {
            return ClusteringProtectOverClustering(config);
        }
        else if (function_name == "clustering_neutrino") {
            return ClusteringNeutrino(config);
        }
        else if (function_name == "clustering_isolated") {
            return ClusteringIsolated(config);
        }
        else {
            throw std::invalid_argument("Unknown function name in configuration");
        }
    }
}  // namespace WireCell::PointCloud::Facade
