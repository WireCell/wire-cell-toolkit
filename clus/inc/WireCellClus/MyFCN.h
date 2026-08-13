#ifndef WIRECELLCLUS_MYFCN_H
#define WIRECELLCLUS_MYFCN_H

#include "WireCellClus/PRGraph.h"
#include "WireCellClus/Facade_Cluster.h"
#include "WireCellUtil/Units.h"
#include "WireCellClus/TrackFitting.h"


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
        /// @return true if the vertex and its segments were actually updated.
        ///         The prototype counterpart (NeutrinoID_improve_vertex.h:816)
        ///         is void and cannot fail; the toolkit adds three guards for
        ///         conditions the prototype's single-APA geometry cannot reach
        ///         (fit_pos outside every detector volume, missing wpid
        ///         offset/slope entry, missing steiner_pc).  Each aborts before
        ///         any write, so callers must not treat that as a fit.
        ///         See sbnd_xin/docs/pr/28 sec. 3.3.
        bool UpdateInfo(Facade::geo_point_t fit_pos,
                       Facade::Cluster& temp_cluster,
                       TrackFitting& track_fitter, IDetectorVolumes::pointer dv,
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

        /// doc sbnd_xin/docs/pr/51 round 7 (robust vertex fit): parameters of
        /// the per-leg dynamic direction-window substitution.  With on=false
        /// (the default) AddSegment's epilogue never runs and every code path
        /// is byte-identical to the legacy fit.  When on, a non-shower leg
        /// whose fits-chord is >= min_len gets a second, re-seat-free
        /// direction estimate over the annulus (rin, rout] with
        ///   rin  = max(vertex_protect_dis, reseat + rin_margin)
        ///          (reseat mirrors UpdateInfo: default_dis_cut iff the leg's
        ///           far wcpt end is > 2*default_dis_cut from the vertex)
        ///   rout = clamp(rout_frac * chord, rout_min, rout_max)
        /// and its PCA/center are SUBSTITUTED (vec_PCA_dirs / vec_PCA_vals /
        /// vec_centers only -- vec_points, hence ntracks and the npoints
        /// prior weight, stay production) iff the folded inner-vs-outer axis
        /// angle exceeds angle_deg, the outer window holds >= min_pts points
        /// and its anisotropy sqrt(l0/l1) >= min_aniso.  Substituted axes are
        /// hemisphere-oriented toward their production counterparts (the
        /// FitVertex pair-angle census is sign-sensitive).  prior_range
        /// replaces vtx_constraint_range in FitVertex iff >= 1 leg was
        /// substituted AND exactly 2 legs are fittable (a distorted 2-track
        /// vertex needs more corrective authority; >= 3-track vertices are
        /// already well braced).  All lengths internal units, angle deg.
        struct RobustParams {
            bool on{false};
            double min_len{10 * units::cm};
            double rin_margin{2 * units::cm};
            double rout_frac{0.5};
            double rout_min{9 * units::cm};
            double rout_max{18 * units::cm};
            double angle{20};
            int min_pts{5};
            double min_aniso{3.0};
            double prior_range{1.0 * units::cm};
        };

        /// Enable the robust per-leg substitution (call BEFORE AddSegment).
        void set_robust(const RobustParams& p) { m_robust = p; }

        /// Number of legs whose direction estimate was substituted.
        int robust_substituted() const { return m_n_substituted; }

        /// Undo every substitution (production PCA/centers restored).  Called
        /// by the fit_vertex charge-veto path so a rejected robust fit
        /// re-seats exactly as the production fit would have.
        void restore_production_pca();

    private:
        /// AddSegment epilogue: compute the outer-window estimate for the leg
        /// just pushed and substitute if the disagreement gate fires.
        void robust_maybe_substitute_leg(SegmentPtr sg);

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

        // doc sbnd_xin/docs/pr/51 round 7: robust-substitution state.  Backups
        // are index-aligned tuples (leg index, dirs, vals, center) so
        // restore_production_pca() is exact; std::vector only (determinism:
        // no pointer-keyed containers).
        RobustParams m_robust;
        int m_n_substituted{0};
        std::vector<std::tuple<size_t,
                               std::tuple<Facade::geo_point_t, Facade::geo_point_t, Facade::geo_point_t>,
                               std::tuple<double, double, double>,
                               Facade::geo_point_t>> m_robust_backup;
    };

} // namespace WireCell::Clus::PR

#endif // WIRECELLCLUS_MYFCN_H
