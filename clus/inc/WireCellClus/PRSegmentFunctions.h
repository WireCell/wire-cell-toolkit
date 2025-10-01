#ifndef WIRECELL_CLUS_PRSEGMENTFUNCTIONS
#define WIRECELL_CLUS_PRSEGMENTFUNCTIONS

#include "WireCellClus/PRGraph.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/D4Vector.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellIface/IRecombinationModel.h"
#include "WireCellClus/ParticleDataSet.h"

namespace WireCell::Clus::PR {

    using geo_point_t = WireCell::Point;

    /// Replace the segment in the graph with two new segments that meet at a
    /// new vertex nearest to the point.
    ///
    /// The input segment is removed from the graph.
    ///
    /// The point must be withing max_dist of the segment.
    ///
    /// Returns true if the graph was modified.
    bool break_segment(Graph& graph, SegmentPtr seg, Point point,
                       double max_dist=1e9*units::cm);
    // patter recognition
    std::tuple<WireCell::Point, WireCell::Vector, WireCell::Vector, bool> segment_search_kink(SegmentPtr seg, WireCell::Point& start_p, const std::string& cloud_name = "fit", double dQ_dx_threshold = 43000/units::cm );

    /// Calculate track length from segment
    ///
    /// If flag == 1 and segment has fitted dx values, sum the dx values from fits.
    /// If flag == 0, calculate geometric length from wcpts.
    ///
    /// @param seg The segment to calculate length for
    /// @param flag Calculation method: 0=geometric from points, 1=from fit dx values
    /// @return Track length
    double segment_track_length(SegmentPtr seg, int flag = 0, int n1 = -1, int n2 = -1, WireCell::Vector dir_perp = WireCell::Vector(0,0,0));
    double segment_track_direct_length(SegmentPtr seg, int n1 = -1, int n2 = -1, WireCell::Vector dir_perp = WireCell::Vector(0,0,0));
    double segment_track_max_deviation(SegmentPtr seg, int n1 = -1, int n2 = -1);
    /// Calculate track length above dQ/dx threshold
    ///
    /// Extracts dQ and dx from segment's fits and calculates length above threshold.
    ///
    /// @param seg The segment containing fit data
    /// @param threshold dQ/dx threshold value
    /// @return Length of track segments above threshold
    double segment_track_length_threshold(SegmentPtr seg, double threshold = 75000./units::cm);
    /// Calculate track length from segment using geometric distance between points
    ///
    /// This is a convenience function that always uses geometric calculation
    /// regardless of available dx data.
    ///
    /// @param seg The segment to calculate length for
    /// @return Geometric track length
    double segment_geometric_length(SegmentPtr seg, int n1 = -1, int n2 = -1, WireCell::Vector dir_perp = WireCell::Vector(0,0,0));


    /// Calculate median dQ/dx for a segment
    ///
    /// Extracts dQ and dx from segment's fits and calculates median dQ/dx.
    ///
    /// @param seg The segment containing fit data
    /// @return Median dQ/dx value (0 if no valid fits)
    double segment_median_dQ_dx(SegmentPtr seg);
    double segment_rms_dQ_dx(SegmentPtr seg);
    
    
    /// Create and associate a DynamicPointCloud with a segment from path points
    ///
    /// @param segment The segment to associate the DynamicPointCloud with
    /// @param path_points Vector of 3D points to process
    /// @param dv Detector volume for wire plane ID determination
    /// @param cloud_name Name for the DynamicPointCloud (default: "main")
    void create_segment_point_cloud(SegmentPtr segment,
                                    const std::vector<geo_point_t>& path_points,
                                    const IDetectorVolumes::pointer& dv,
                                    const std::string& cloud_name = "main");

    void create_segment_fit_point_cloud(SegmentPtr segment,
                                    const IDetectorVolumes::pointer& dv,
                                    const std::string& cloud_name = "fit");

    std::pair<double, WireCell::Point> segment_get_closest_point(SegmentPtr seg, const WireCell::Point& point, const std::string& cloud_name = "fit");
    std::tuple<double, double, double> segment_get_closest_2d_distances(SegmentPtr seg, const WireCell::Point& point, int apa, int face, const std::string& cloud_name = "fit");


    // PID related 
    bool eval_ks_ratio(double ks1, double ks2, double ratio1, double ratio2);
    std::vector<double> do_track_comp(std::vector<double>& L , std::vector<double>& dQ_dx, double compare_range, double offset_length, const Clus::ParticleDataSet::pointer& particle_data, double MIP_dQdx = 50000/units::cm);
    // success, flag_dir, pdg_code, particle_score
    std::tuple<bool, int, int, double> segment_do_track_pid(SegmentPtr segment, std::vector<double>& L , std::vector<double>& dQ_dx, double compare_range , double offset_length, bool flag_force, const Clus::ParticleDataSet::pointer& particle_data, double MIP_dQdx = 50000/units::cm);

    // direction calculation ...
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg);
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg, WireCell::Point& p, double dis_cut);
    WireCell::Vector segment_cal_dir_3vector(SegmentPtr seg, int direction, int num_points, int start);
    void segment_determine_dir_track(SegmentPtr segment, int start_n, int end_n, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx = 43000/units::cm, bool flag_print = false);

    // kinemiatics calculations ...
    double segment_cal_kine_dQdx(SegmentPtr seg, const IRecombinationModel::pointer& recomb_model);
    double cal_kine_dQdx(std::vector<double>& vec_dQ, std::vector<double>& vec_dx, const IRecombinationModel::pointer& recomb_model);
    double cal_kine_range(double L, int pdg_code, const Clus::ParticleDataSet::pointer& particle_data);
    // 4-momentum: E, px, py, pz
    WireCell::D4Vector<double> segment_cal_4mom(SegmentPtr segment, int pdg_code, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx = 50000/units::cm);

    // EMshower PID
    void clustering_points_segments(std::set<SegmentPtr> segments, const IDetectorVolumes::pointer& dv, const std::string& cloud_name = "associate_points", double search_range = 1.2*units::cm, double scaling_2d = 0.7);

    bool segment_is_shower_trajectory(SegmentPtr seg, double step_size = 10*units::cm, double mip_dQ_dx = 50000 / units::cm);
    void segment_determine_shower_direction_trajectory(SegmentPtr segment, int start_n, int end_n, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double MIP_dQdx = 43000/units::cm, bool flag_print = false);
    
    bool segment_determine_shower_direction(SegmentPtr segment, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, const std::string& cloud_name = "associate_points", double MIP_dQdx = 43000/units::cm, double rms_cut= 0.4*units::cm);
    bool segment_is_shower_topology(SegmentPtr seg, bool tmp_val=false, double MIP_dQ_dx = 43000/units::cm);
}

#endif
