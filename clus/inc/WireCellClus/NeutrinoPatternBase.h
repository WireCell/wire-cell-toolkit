#include "WireCellClus/PRGraph.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/PRShower.h"
#include "WireCellClus/NeutrinoTaggerInfo.h"
#include "WireCellClus/IClusGeomHelper.h"

namespace WireCell::Clus::PR {

    /** BDT input features for the best pi0 candidate found by id_pi0_with_vertex (flag==1)
     *  or id_pi0_without_vertex (flag==2).  All angles in radians, energies and distances
     *  in WireCell internal units (multiply by 1/units::MeV, 1/units::cm for output).
     *  Default-initialised to zero (flag==0 means no pi0 found).
     */
    /// Type aliases for the pre-collected 2D charge maps used by cal_kine_charge.
    /// Collect once per shower_clustering_with_nv call via TrackFitting::collect_2D_charge()
    /// and reuse across all showers/merges to avoid O(N_hits) re-collection per shower.
    using ChargeMap = std::map<TrackFitting::CoordReadout, TrackFitting::ChargeMeasurement>;
    using WireMap   = std::map<std::pair<int,int>, std::vector<std::tuple<int,int,int>>>;

    struct Pi0KineFeatures {
        int    flag{0};      ///< 0=none, 1=with_vertex, 2=without_vertex
        double mass{0};      ///< reconstructed pi0 invariant mass
        double vtx_dis{0};   ///< distance from pi0 decay vertex to main vertex
        double energy_1{0};  ///< kine_charge of first shower
        double theta_1{0};   ///< polar angle of first shower direction (from z-axis)
        double phi_1{0};     ///< azimuthal angle of first shower direction
        double dis_1{0};     ///< distance from pi0 vertex to first shower start point
        double energy_2{0};  ///< kine_charge of second shower
        double theta_2{0};   ///< polar angle of second shower direction
        double phi_2{0};     ///< azimuthal angle of second shower direction
        double dis_2{0};     ///< distance from pi0 vertex to second shower start point
        double angle{0};     ///< opening angle between the two shower directions
    };

    class PatternAlgorithms{
        public:
        bool m_perf{false};  // if true, print per-step timing to stdout

        // 2D charge maps cached for the duration of shower_clustering_with_nv.
        // Populated once by collect_charge_maps(); reused by calculate_shower_kinematics
        // and all cal_kine_charge call sites to avoid O(N_hits) re-collection per shower.
        ChargeMap m_charge_2d_u, m_charge_2d_v, m_charge_2d_w;
        WireMap   m_map_apa_ch_plane_wires;

        // Populate the cached charge maps from track_fitter.
        // Call once at the start of shower_clustering_with_nv; the maps are then
        // valid for the entire call tree beneath it.
        void collect_charge_maps(TrackFitting& track_fitter);
        std::vector<VertexPtr> find_cluster_vertices(Graph& graph, const Facade::Cluster& cluster);
        std::vector<SegmentPtr> find_cluster_segments(Graph& graph, const Facade::Cluster& cluster);
        bool clean_up_graph(Graph& graph, const Facade::Cluster& cluster);

        SegmentPtr init_first_segment(Graph& graph, Facade::Cluster& cluster, Facade::Cluster* main_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_back_search = true);

        // find the shortest path using steiner graph
        std::vector<Facade::geo_point_t> do_rough_path(const Facade::Cluster& cluster,Facade::geo_point_t& first_point, Facade::geo_point_t& last_point);
        std::vector<Facade::geo_point_t> do_rough_path_reg_pc(const Facade::Cluster& cluster, Facade::geo_point_t& first_point, Facade::geo_point_t& last_point,  std::string graph_name = "relaxed_pid");
        // create a segment given a path
        SegmentPtr create_segment_for_cluster(WireCell::Clus::Facade::Cluster& cluster, IDetectorVolumes::pointer dv, const std::vector<Facade::geo_point_t>& path_points, int dir = 0);
        // create a segment given two vertices, null, if failed
        SegmentPtr create_segment_from_vertices(Graph& graph, Facade::Cluster& cluster, VertexPtr v1, VertexPtr v2, IDetectorVolumes::pointer dv);
        // replace a segment and vertex with another segment and vertex, assuming the original vertex only connect to this segment
        bool replace_segment_and_vertex(Graph& graph, SegmentPtr& seg, VertexPtr& vtx, std::list<Facade::geo_point_t>& path_point_list, Facade::geo_point_t& break_point, IDetectorVolumes::pointer dv);
        bool replace_segment_and_vertex(Graph& graph, SegmentPtr& seg, VertexPtr old_vertex, VertexPtr new_vertex, IDetectorVolumes::pointer dv);
        bool break_segment_into_two(Graph& graph, VertexPtr vtx1, SegmentPtr seg, VertexPtr vtx2, std::list<Facade::geo_point_t>& path_point_list1, Facade::geo_point_t& break_point, std::list<Facade::geo_point_t>& path_point_list2, IDetectorVolumes::pointer dv, SegmentPtr& out_seg2);


        // return the point and its index in the steiner tree as a pair
        std::pair<Facade::geo_point_t, size_t> proto_extend_point(const Facade::Cluster& cluster, Facade::geo_point_t& p, Facade::geo_vector_t& dir, Facade::geo_vector_t& dir_other, bool flag_continue);
        // return Steiner Graph path in wcps_list1 and wcps_list2
        bool proto_break_tracks(const Facade::Cluster& cluster, const Facade::geo_point_t& first_wcp, Facade::geo_point_t& curr_wcp, const Facade::geo_point_t& last_wcp, std::list<Facade::geo_point_t>& wcps_list1, std::list<Facade::geo_point_t>& wcps_list2, bool flag_pass_check);
        // breaking segments ...
        bool break_segments(Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, std::vector<SegmentPtr>& remaining_segments, float dis_cut = 0);
        // merge vertices within 0.1 cm after break_segments, then refit if changed
        bool merge_nearby_vertices(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        // merge two segments to one
        bool merge_two_segments_into_one(Graph& graph, SegmentPtr& seg1, VertexPtr& vtx, SegmentPtr& seg2, IDetectorVolumes::pointer dv);
        // merge vertex into another
        bool merge_vertex_into_another(Graph& graph, VertexPtr& vtx_from, VertexPtr& vtx_to, IDetectorVolumes::pointer dv);

        // get direction with  distance cut ... 
        Facade::geo_vector_t vertex_get_dir(VertexPtr& vertex, Graph& graph, double dis_cut = 5*units::cm);
        Facade::geo_vector_t vertex_segment_get_dir(VertexPtr& vertex, SegmentPtr& segment, Graph& graph, double dis_cut = 2*units::cm);


        // Structure examination
        void examine_structure(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv); // call examine_structure_1 and examine_structure_2      
        bool examine_structure_1(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_2(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        bool examine_structure_3(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_4(VertexPtr vertex, bool flag_final_vertex, Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // identify other segments giving the graph ...
        void find_other_segments(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv , bool flag_break_track =true, double search_range = 1.5*units::cm, double scaling_2d = 0.8);
        VertexPtr find_vertex_other_segment(Graph& graph, Facade::Cluster& cluster, SegmentPtr seg, bool flag_forwrard, Facade::geo_point_t& wcp, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );
        std::tuple<VertexPtr, SegmentPtr, Facade::geo_point_t> check_end_point(Graph& graph, Facade::Cluster& cluster, std::vector<Facade::geo_point_t>& tracking_path, bool flag_front = true, double vtx_cut1 = 0.9 * units::cm, double vtx_cut2 = 2.0 * units::cm, double sg_cut1  = 2.0 * units::cm, double sg_cut2  = 1.2 * units::cm);
        bool modify_vertex_isochronous(Graph& graph, Facade::Cluster& cluster, VertexPtr vtx, VertexPtr v1, SegmentPtr sg, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool modify_segment_isochronous(Graph& graph, Facade::Cluster& cluster, SegmentPtr sg1, VertexPtr v1, SegmentPtr sg, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double dis_cut = 6*units::cm, double angle_cut = 15, double extend_cut = 15*units::cm);

        // examine segment
        void examine_segment(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool crawl_segment(Graph& graph, Facade::Cluster& cluster, SegmentPtr seg, VertexPtr vertex, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );
        void examine_partial_identical_segments(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        //examine vertices
        void examine_vertices(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_1(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_1p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_vertices_2(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_4(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, VertexPtr main_vertex = nullptr);
        bool examine_vertices_4p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::geo_point_t get_local_extension(Facade::Cluster& cluster, const Facade::geo_point_t& wcp);
        void examine_vertices_3(Graph& graph, Facade::Cluster& main_cluster, std::pair<VertexPtr, VertexPtr> main_cluster_initial_pair_vertices, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // master pattern recognition function
        bool find_proto_vertex(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_break_track = true, int nrounds_find_other_tracks = 2, bool flag_back_search = true);
        
        void init_point_segment(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // examine the structure of the patterns ... 
        bool examine_structure_final_1(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final_1p(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final_2(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final_3(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_structure_final(Graph& graph, VertexPtr main_vertex, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // EM shower related
        void clustering_points(Graph& graph, Facade::Cluster& cluster, const IDetectorVolumes::pointer& dv, const std::string& cloud_name = "associate_points", double search_range = 1.2*units::cm, double scaling_2d = 0.7);
        void separate_track_shower(Graph&graph, Facade::Cluster& cluster);
        // Direction
        void determine_direction(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        std::pair<int, double> calculate_num_daughter_showers(Graph& graph, VertexPtr vertex, SegmentPtr segment, bool flag_count_shower = true);
        std::pair<int, double> calculate_num_daughter_tracks(Graph& graph, VertexPtr vtx, SegmentPtr sg, bool flag_count_shower = false, double length_cut = 0);
        // find_cont_muon_segment_nue: find adjacent segment continuing in same direction (MIP-like)
        std::pair<SegmentPtr, VertexPtr> find_cont_muon_segment_nue(Graph& graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx = false);
        void examine_good_tracks(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data);
        // about fix maps
        void fix_maps_multiple_tracks_in(Graph& graph, Facade::Cluster& cluster);
        void fix_maps_shower_in_track_out(Graph& graph, Facade::Cluster& cluster);
        void improve_maps_one_in(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_strong_check = true);
        void improve_maps_shower_in_track_out(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_strong_check = true);
        void improve_maps_no_dir_tracks(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void improve_maps_multiple_tracks_in(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void judge_no_dir_tracks_close_to_showers(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, IDetectorVolumes::pointer dv);
        bool examine_maps(Graph&graph, Facade::Cluster& cluster);
        void examine_all_showers(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data);
        void shower_determining_in_main_cluster(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, IDetectorVolumes::pointer dv);
        void set_default_shower_particle_info(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);

        // PCA calculation
        std::pair<Facade::geo_point_t, Facade::geo_vector_t> calc_PCA_main_axis(std::vector<Facade::geo_point_t>& points);

        // vertex related functions 
        bool search_for_vertex_activities(Graph& graph, VertexPtr vertex, std::vector<SegmentPtr>& segments_set, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double search_range = 1.5*units::cm);
        bool eliminate_short_vertex_activities(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex, std::set<SegmentPtr>& existing_segments, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        std::tuple<bool, int, int> examine_main_vertex_candidate(Graph& graph, VertexPtr vertex);
        VertexPtr compare_main_vertices_all_showers(Graph& graph, Facade::Cluster& cluster, std::vector<VertexPtr>& vertex_candidates, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        float calc_conflict_maps(Graph& graph, VertexPtr vertex);
        VertexPtr compare_main_vertices(Graph& graph, Facade::Cluster& cluster, std::vector<VertexPtr>& vertex_candidates);
        std::pair<SegmentPtr, VertexPtr> find_cont_muon_segment(Graph &graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx = false);
        bool examine_direction(Graph& graph, VertexPtr vertex, VertexPtr main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_final = false);

        bool fit_vertex(Facade::Cluster& cluster, VertexPtr vertex, VertexPtr main_vertex, std::vector<SegmentPtr>& sg_set, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void improve_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_search_vertex_activity = true , bool flag_final_vertex = false);
        void determine_main_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr& main_vertex, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void change_daughter_type(Graph& graph, VertexPtr vertex, SegmentPtr segment, int particle_type, double mass, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_main_vertices_local(Graph& graph, std::vector<VertexPtr>& vertices, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);

        // cluster functions ...
        Facade::geo_vector_t calc_dir_cluster(Graph& graph, Facade::Cluster& cluster, const Facade::geo_point_t& orig_p, double dis_cut);
        Facade::Cluster* swap_main_cluster(Facade::Cluster& new_main_cluster, Facade::Cluster& old_main_cluster, std::vector<Facade::Cluster*>& other_clusters);
        void examine_main_vertices(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters);

        VertexPtr compare_main_vertices_global(Graph& graph, std::vector<VertexPtr>& vertex_candidates, Facade::Cluster& main_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::Cluster* check_switch_main_cluster(Graph& graph, ClusterVertexMap map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::Cluster* check_switch_main_cluster_2(Graph& graph, VertexPtr temp_main_vertex, Facade::Cluster* max_length_cluster, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters);
        VertexPtr determine_overall_main_vertex(Graph& graph, ClusterVertexMap map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_dev_chain = true);

        // Deep-learning vertex refinement.  Returns true if the DL network changed
        // the selected vertex (in which case the traditional determine_overall_main_vertex
        // should NOT be called).  Returns false if the DL network was unavailable, failed,
        // or did not improve on the candidate vertices (fall back to traditional).
        bool determine_overall_main_vertex_DL(Graph& graph, ClusterVertexMap& map_cluster_main_vertices, Facade::Cluster*& main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const std::string& dl_weights, double dl_vtx_cut, double dQdx_scale = 0.1, double dQdx_offset = -1000.0);

        // global information transfer
        void transfer_info_from_segment_to_cluster(Graph& graph, Facade::Cluster& cluster, const std::string& cloud_name = "associated_points");

        // print information
        void print_segs_info(Graph& graph, Facade::Cluster& cluster, VertexPtr vertex= nullptr);

        // shower related functions
        void update_shower_maps(IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters);
        void shower_clustering_with_nv_in_main_cluster(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_connecting_to_main_vertex(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters);
        void shower_clustering_with_nv_from_main_cluster(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters);
        void shower_clustering_with_nv_from_vertices(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, IndexedShowerSet& showers, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void calculate_shower_kinematics(IndexedShowerSet& showers, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_merge_showers(IndexedShowerSet& showers, VertexPtr main_vertex, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_in_other_clusters(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_save = true);
        void examine_shower_1(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_showers(Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void id_pi0_with_vertex(int acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void id_pi0_without_vertex(int acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, IndexedSegmentSet& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_with_nv(int acc_segment_id, IndexedShowerSet& pi0_showers, ShowerIntMap& map_shower_pio_id, std::map<int, std::vector<ShowerPtr > >& map_pio_id_showers, std::map<int, std::pair<double, int> >& map_pio_id_mass, std::map<int, std::pair<int, int> >& map_pio_id_saved_pair, Pi0KineFeatures& pio_kine, IndexedVertexSet& vertices_in_long_muon, IndexedSegmentSet& segments_in_long_muon, Graph& graph, VertexPtr main_vertex, IndexedShowerSet& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, ClusterVertexMap map_cluster_main_vertices, ShowerVertexMap& map_vertex_in_shower, ShowerSegmentMap& map_segment_in_shower, VertexShowerSetMap& map_vertex_to_shower, ClusterPtrSet& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);


        // deghost related functions ...
        void order_clusters(Graph& graph, std::vector<Facade::Cluster*>& ordered_clusters, std::map<Facade::Cluster*, std::vector<SegmentPtr> >& map_cluster_to_segments, std::map<Facade::Cluster*, double>& map_cluster_total_length);
        void deghost_clusters(Graph& graph, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void order_segments(std::vector<SegmentPtr>& ordered_segments, std::vector<SegmentPtr>& segments);
        void deghost_segments(Graph& graph, ClusterVertexMap map_cluster_main_vertices, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void deghosting(Graph& graph, ClusterVertexMap map_cluster_main_vertices, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );

        // energy calculation ...
        double cal_corr_factor(WireCell::Point& pt, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        // kinematics output ...
        void init_tagger_info(TaggerInfo& ti);

        // cosmic tagger
        bool bad_reconstruction(Graph& graph, VertexPtr main_vertex, ShowerPtr shower,
                                bool flag_fill = false, TaggerInfo* ti = nullptr);
        bool cosmic_tagger(Graph& graph, VertexPtr main_vertex,
                           IndexedShowerSet& showers,
                           ShowerSegmentMap& map_segment_in_shower,
                           VertexShowerSetMap& map_vertex_to_shower,
                           IndexedSegmentSet& segments_in_long_muon,
                           Facade::Cluster* main_cluster,
                           std::vector<Facade::Cluster*>& all_clusters,
                           IDetectorVolumes::pointer dv,
                           TaggerInfo& ti);

        // numu tagger
        // count_daughters: counts track branches and total branches at the far end of a muon.
        // For a segment: BFS from the vertex closer to main_vertex through sg, count what is beyond.
        // For a shower:  finds the last segment of the long-muon chain, then counts from its near end.
        std::pair<int, int> count_daughters(Graph& graph, SegmentPtr sg, VertexPtr main_vertex);
        std::pair<int, int> count_daughters(Graph& graph, ShowerPtr shower, VertexPtr main_vertex,
                                            std::set<SegmentPtr>& segments_in_long_muon_plain);
        // numu_tagger: fills TaggerInfo numu_cc_* BDT features and returns
        //   {flag_long_muon, max_muon_length}.
        // Prototype: WCPPID::NeutrinoID::numu_tagger() in NeutrinoID_numu_tagger.h.
        std::pair<bool, double> numu_tagger(Graph& graph,
                                            VertexPtr main_vertex,
                                            IndexedShowerSet& showers,
                                            IndexedSegmentSet& segments_in_long_muon,
                                            Facade::Cluster* main_cluster,
                                            TaggerInfo& ti);

        // nue_tagger: fills TaggerInfo NuE BDT features and returns flag_nue.
        // Prototype: WCPPID::NeutrinoID::nue_tagger() in NeutrinoID_nue_tagger.h.
        bool nue_tagger(Graph& graph,
                        Facade::Cluster* main_cluster,
                        VertexPtr main_vertex,
                        int apa, int face,
                        IndexedShowerSet& showers,
                        VertexShowerSetMap& map_vertex_to_shower,
                        IndexedShowerSet& pi0_showers,
                        ShowerIntMap& map_shower_pio_id,
                        std::map<int, std::vector<ShowerPtr>>& map_pio_id_showers,
                        std::map<int, std::pair<double,int>>& map_pio_id_mass,
                        IDetectorVolumes::pointer dv,
                        ParticleDataSet::pointer particle_data,
                        double muon_length,
                        TaggerInfo& ti);

        // singlephoton_tagger: fills TaggerInfo shw_sp_* BDT features and returns flag_sp.
        // Prototype: WCPPID::NeutrinoID::singlephoton_tagger() in NeutrinoID_singlephoton_tagger.h.
        // apa/face are derived internally from dv->contained_by(main_vertex_pt) so callers
        // do not need to know detector geometry details.
        bool singlephoton_tagger(Graph& graph,
                                 Facade::Cluster* main_cluster,
                                 VertexPtr main_vertex,
                                 IndexedShowerSet& showers,
                                 VertexShowerSetMap& map_vertex_to_shower,
                                 ShowerIntMap& map_shower_pio_id,
                                 std::map<int, std::vector<ShowerPtr>>& map_pio_id_showers,
                                 std::map<int, std::pair<double,int>>& map_pio_id_mass,
                                 IDetectorVolumes::pointer dv,
                                 TaggerInfo& ti);

        // ssm_tagger: fills TaggerInfo ssm_* and ssmsp_* features and returns flag_ssm.
        // Prototype: WCPPID::NeutrinoID::ssm_tagger() in NeutrinoID_ssm_tagger.h.
        bool ssm_tagger(Graph& graph,
                        VertexPtr main_vertex,
                        IndexedShowerSet& showers,
                        ShowerVertexMap& map_vertex_in_shower,
                        ShowerSegmentMap& map_segment_in_shower,
                        const Pi0KineFeatures& pio_kine,
                        int flag_ssmsp,
                        int& acc_segment_id,
                        const ParticleDataSet::pointer& particle_data,
                        const IRecombinationModel::pointer& recomb_model,
                        TaggerInfo& ti);

        KineInfo fill_kine_tree(VertexPtr main_vertex, IndexedShowerSet& showers,
                                const Pi0KineFeatures& pio_kine,
                                Graph& graph, TrackFitting& track_fitter,
                                IDetectorVolumes::pointer dv,
                                WireCell::IClusGeomHelper::pointer geom_helper,
                                const Clus::ParticleDataSet::pointer& particle_data,
                                const IRecombinationModel::pointer& recomb_model);
        // Convenience overloads: collect 2D charge maps internally (safe for isolated calls).
        double cal_kine_charge(ShowerPtr Shower, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        double cal_kine_charge(SegmentPtr segment, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        // Fast overload: reuse pre-collected 2D charge maps (avoids O(N_hits) collection per call).
        // Collect maps once with track_fitter.collect_2D_charge() and pass here when calling in a loop.
        double cal_kine_charge(ShowerPtr shower,
            const ChargeMap& charge_2d_u, const ChargeMap& charge_2d_v, const ChargeMap& charge_2d_w,
            const WireMap& map_apa_ch_plane_wires,
            TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

    };
}
