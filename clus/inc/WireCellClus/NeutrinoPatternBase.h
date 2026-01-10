#include "WireCellClus/PRGraph.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/PRShower.h"

namespace WireCell::Clus::PR {
    class PatternAlgorithms{
        public:
        std::set<VertexPtr> find_cluster_vertices(Graph& graph, const Facade::Cluster& cluster);
        std::set<SegmentPtr> find_cluster_segments(Graph& graph, const Facade::Cluster& cluster);
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
        bool break_segment_into_two(Graph& graph, VertexPtr vtx1, SegmentPtr seg, VertexPtr vtx2, std::list<Facade::geo_point_t>& path_point_list1, Facade::geo_point_t& break_point, std::list<Facade::geo_point_t>& path_point_list2, IDetectorVolumes::pointer dv);


        // return the point and its index in the steiner tree as a pair
        std::pair<Facade::geo_point_t, size_t> proto_extend_point(const Facade::Cluster& cluster, Facade::geo_point_t& p, Facade::geo_vector_t& dir, Facade::geo_vector_t& dir_other, bool flag_continue);
        // return Steiner Graph path in wcps_list1 and wcps_list2
        bool proto_break_tracks(const Facade::Cluster& cluster, const Facade::geo_point_t& first_wcp, const Facade::geo_point_t& curr_wcp, const Facade::geo_point_t& last_wcp, std::list<Facade::geo_point_t>& wcps_list1, std::list<Facade::geo_point_t>& wcps_list2, bool flag_pass_check);
        // breaking segments ...
        bool break_segments(Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, std::vector<SegmentPtr>& remaining_segments, float dis_cut = 0);
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

        // examine segment
        void examine_segment(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool crawl_segment(Graph& graph, Facade::Cluster& cluster, SegmentPtr seg, VertexPtr vertex, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );
        void examine_partial_identical_segments(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

        //examine vertices
        void examine_vertices(Graph& graph, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_vertices_1(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_vertices_1p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_vertices_2(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_vertices_4(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_vertices_4p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::geo_point_t get_local_extension(Facade::Cluster& cluster, Facade::geo_point_t& wcp);
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
        void clustering_points(Graph& graph, Facade::Cluster& cluster, const IDetectorVolumes::pointer& dv, const std::string& cloud_name = "associated_points", double search_range = 1.2*units::cm, double scaling_2d = 0.7);
        void separate_track_shower(Graph&graph, Facade::Cluster& cluster);
        // Direction
        void determine_direction(Graph& graph, Facade::Cluster& cluster, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        std::pair<int, double> calculate_num_daughter_showers(Graph& graph, VertexPtr vertex, SegmentPtr segment, bool flag_count_shower = true); 
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

        // PCA calculation
        std::pair<Facade::geo_point_t, Facade::geo_vector_t> calc_PCA_main_axis(std::vector<Facade::geo_point_t>& points);

        // vertex related functions 
        bool search_for_vertex_activities(Graph& graph, VertexPtr vertex, std::set<SegmentPtr>& segments_set, Facade::Cluster& cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, double search_range = 1.5*units::cm);
        bool eliminate_short_vertex_activities(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex, std::set<SegmentPtr>& existing_segments, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        std::tuple<bool, int, int> examine_main_vertex_candidate(Graph& graph, VertexPtr vertex);
        VertexPtr compare_main_vertices_all_showers(Graph& graph, Facade::Cluster& cluster, std::vector<VertexPtr>& vertex_candidates, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        float calc_conflict_maps(Graph& graph, VertexPtr vertex);
        VertexPtr compare_main_vertices(Graph& graph, Facade::Cluster& cluster, std::vector<VertexPtr>& vertex_candidates);
        std::pair<SegmentPtr, VertexPtr> find_cont_muon_segment(Graph &graph, SegmentPtr sg, VertexPtr vtx, bool flag_ignore_dQ_dx = false);
        bool examine_direction(Graph& graph, VertexPtr vertex, VertexPtr main_vertex, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_final = false);

        bool fit_vertex(Facade::Cluster& cluster, VertexPtr vertex, VertexPtr main_vertex, std::set<SegmentPtr>& sg_set, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void improve_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_search_vertex_activity = true , bool flag_final_vertex = false);
        void determine_main_vertex(Graph& graph, Facade::Cluster& cluster, VertexPtr main_vertex, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool flag_print = true);
        void change_daughter_type(Graph& graph, VertexPtr vertex, SegmentPtr segment, int particle_type, double mass, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_main_vertices_local(Graph& graph, std::vector<VertexPtr>& vertices, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);

        // cluster functions ...
        Facade::geo_vector_t calc_dir_cluster(Graph& graph, Facade::Cluster& cluster, const Facade::geo_point_t& orig_p, double dis_cut);
        Facade::Cluster* swap_main_cluster(Facade::Cluster& new_main_cluster, Facade::Cluster& old_main_cluster, std::vector<Facade::Cluster*>& other_clusters);
        void examine_main_vertices(Graph& graph, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters);

        VertexPtr compare_main_vertices_global(Graph& graph, std::vector<VertexPtr>& vertex_candidates, Facade::Cluster& main_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::Cluster* check_switch_main_cluster(Graph& graph, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        Facade::Cluster* check_switch_main_cluster_2(Graph& graph, VertexPtr temp_main_vertex, Facade::Cluster* max_length_cluster, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters);
        VertexPtr determine_overall_main_vertex(Graph& graph, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model );

        // global information transfer
        void transfer_info_from_segment_to_cluster(Graph& graph, Facade::Cluster& cluster, const std::string& cloud_name = "associated_points");

        // print information
        void print_segs_info(Graph& graph, Facade::Cluster& cluster, VertexPtr vertex= nullptr);

        // shower related functions
        void update_shower_maps(std::set<ShowerPtr>& showers,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters);
        void shower_clustering_with_nv_in_main_cluster(Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon);
        void shower_clustering_connecting_to_main_vertex(Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters);
        void shower_clustering_with_nv_from_main_cluster(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, std::set<ShowerPtr>& showers,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters);
        void shower_clustering_with_nv_from_vertices(Graph& graph, VertexPtr main_vertex, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, std::set<ShowerPtr>& showers,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void calculate_shower_kinematics(std::set<ShowerPtr>& showers, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void examine_merge_showers(std::set<ShowerPtr>& showers, VertexPtr main_vertex, std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, std::set<VertexPtr>& vertices_in_long_muon, std::set<SegmentPtr>& segments_in_long_muon, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);
        void shower_clustering_in_other_clusters(Graph& graph, VertexPtr main_vertex, std::set<ShowerPtr>& showers, Facade::Cluster* main_cluster, std::vector<Facade::Cluster*>& other_clusters, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices,  std::map<VertexPtr, ShowerPtr>& map_vertex_in_shower,  std::map<SegmentPtr, ShowerPtr>& map_segment_in_shower, std::map<VertexPtr, std::set<ShowerPtr> >& map_vertex_to_shower, std::set<Facade::Cluster*>& used_shower_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model);


        // deghost related functions ...
        void order_clusters(Graph& graph, std::vector<Facade::Cluster*>& ordered_clusters, std::map<Facade::Cluster*, std::vector<SegmentPtr> >& map_cluster_to_segments, std::map<Facade::Cluster*, double>& map_cluster_total_length);
        void deghost_clusters(Graph& graph, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void order_segments(std::vector<SegmentPtr>& ordered_segments, std::vector<SegmentPtr>& segments);
        void deghost_segments(Graph& graph, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        void deghosting(Graph& graph, std::map<Facade::Cluster*, VertexPtr> map_cluster_main_vertices, std::vector<Facade::Cluster*>& all_clusters, TrackFitting& track_fitter, IDetectorVolumes::pointer dv );

        // energy calculation ...
        double cal_corr_factor(WireCell::Point& pt, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        double cal_kine_charge(ShowerPtr Shower, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        double cal_kine_charge(SegmentPtr segment, Graph& graph, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);

    };
}
