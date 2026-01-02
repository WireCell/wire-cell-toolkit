#ifndef WIRECELLCLUS_MYFCN_H
#define WIRECELLCLUS_MYFCN_H

#include "WireCellClus/PRGraph.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellUtil/Units.h"

#include <vector>
#include <tuple>

namespace WireCell::Clus::PR {

    /// MyFCN: Vertex fitting class for pattern recognition
    /// Performs vertex position fitting by minimizing distances from segments
    class MyFCN {
    public:
        /// Constructor
        /// @param vtx The vertex to fit
        /// @param flag_vtx_constraint Whether to constrain the vertex position
        /// @param vtx_constraint_range Range for vertex constraint
        /// @param vertex_protect_dis Protection distance for vertex
        /// @param vertex_protect_dis_short_track Protection distance for short tracks
        /// @param fit_dis Fitting distance
        MyFCN(VertexPtr vtx, 
              bool flag_vtx_constraint = false, 
              double vtx_constraint_range = 1*units::cm, 
              double vertex_protect_dis = 1.5*units::cm, 
              double vertex_protect_dis_short_track = 0.9*units::cm, 
              double fit_dis = 6*units::cm);
        
        ~MyFCN();

        /// Update the fitting range parameters
        void update_fit_range(double tmp_vertex_protect_dis = 1.5*units::cm, 
                             double tmp_vertex_protect_dis_short_track = 0.9*units::cm, 
                             double tmp_fit_dis = 6*units::cm);
        
        /// Add a segment to the fitting
        void AddSegment(SegmentPtr sg);
        
        /// Fit the vertex position
        /// @return pair of (success flag, fitted position)
        std::pair<bool, Facade::geo_point_t> FitVertex();
        
        /// Update information after fitting
        /// @param fit_pos The fitted position
        /// @param temp_cluster The cluster being processed
        /// @param default_dis_cut Default distance cut
        void UpdateInfo(Facade::geo_point_t fit_pos, 
                       Facade::Cluster* temp_cluster, 
                       double default_dis_cut = 4.0*units::cm);
        
        /// Get segment information at index i
        /// @return pair of (segment pointer, index)
        std::pair<SegmentPtr, int> get_seg_info(int i);
        
        /// Get number of fittable tracks
        int get_fittable_tracks();
        
        /// Get vertex constraint flag
        bool get_flag_vtx_constraint() { return flag_vtx_constraint; }
        
        /// Set vertex constraint flag
        void set_flag_vtx_constraint(bool val) { flag_vtx_constraint = val; }
        
        /// Set vertex constraint range
        void set_vtx_constraint_range(double val) { vtx_constraint_range = val; }
        
        /// Get the vector of segments
        std::vector<SegmentPtr>& get_segments() { return segments; }
        
        /// Get the vector of point vectors
        std::vector<std::vector<Facade::geo_point_t>>& get_vec_points() { return vec_points; }
        
        /// Print points for debugging
        void print_points();
        
        /// Set enforce two track fit flag
        void set_enforce_two_track_fit(bool val) { enforce_two_track_fit = val; }
        
        /// Get enforce two track fit flag
        bool get_enforce_two_track_fit() { return enforce_two_track_fit; }
        
    private:
        VertexPtr vtx;
        bool enforce_two_track_fit;
        bool flag_vtx_constraint;
        double vtx_constraint_range;
        
        double vertex_protect_dis;
        double vertex_protect_dis_short_track;
        double fit_dis;
        
        std::vector<SegmentPtr> segments;
        std::vector<std::vector<Facade::geo_point_t>> vec_points;
        
        // PCA directions: tuple of (dir1, dir2, dir3)
        std::vector<std::tuple<Facade::geo_point_t, Facade::geo_point_t, Facade::geo_point_t>> vec_PCA_dirs;
        
        // PCA eigenvalues: tuple of (val1, val2, val3)
        std::vector<std::tuple<double, double, double>> vec_PCA_vals;
        
        // Centers for each segment
        std::vector<Facade::geo_point_t> vec_centers;
    };

} // namespace WireCell::Clus::PR

#endif // WIRECELLCLUS_MYFCN_H
