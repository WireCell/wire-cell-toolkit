#include "WireCellClus/PRGraph.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellClus/TrackFitting.h"

namespace WireCell::Clus::PR {
    class PatternAlgorithms{
        public:
        std::set<VertexPtr> find_cluster_vertices(Graph& graph, const Facade::Cluster& cluster);
        std::set<SegmentPtr> find_cluster_segments(Graph& graph, const Facade::Cluster& cluster);
        bool clean_up_graph(Graph& graph, const Facade::Cluster& cluster);

        SegmentPtr init_first_segment(Graph& graph, Facade::Cluster& cluster, Facade::Cluster* main_cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv, bool flag_back_search = true);

        // find the shortest path using steiner graph
        std::vector<Facade::geo_point_t> do_rough_path(const Facade::Cluster& cluster,Facade::geo_point_t& first_point, Facade::geo_point_t& last_point);
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

        //examine vertices
        bool examine_vertices_1(Graph&graph, Facade::Cluster&cluster, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
        bool examine_vertices_1p(Graph&graph, VertexPtr v1, VertexPtr v2, TrackFitting& track_fitter, IDetectorVolumes::pointer dv);
    };
}
