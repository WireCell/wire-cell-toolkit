#ifndef WIRECELL_CLUS_PRSEGMENTFUNCTIONS
#define WIRECELL_CLUS_PRSEGMENTFUNCTIONS

#include "WireCellClus/PRGraph.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"
#include "WireCellIface/IDetectorVolumes.h"

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


    /// Calculate median dQ/dx for a segment
    ///
    /// Extracts dQ and dx from segment's fits and calculates median dQ/dx.
    ///
    /// @param seg The segment containing fit data
    /// @return Median dQ/dx value (0 if no valid fits)
    double segment_median_dQ_dx(SegmentPtr seg);
    double segment_rms_dQ_dx(SegmentPtr seg);
    
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

    std::tuple<WireCell::Point, WireCell::Vector, WireCell::Vector, bool> segment_search_kink(SegmentPtr seg, WireCell::Point& start_p, const std::string& cloud_name = "fit", double dQ_dx_threshold = 43000/units::cm );

    bool eval_ks_ratio(double ks1, double ks2, double ratio1, double ratio2);

    bool segment_is_shower_trajectory(SegmentPtr seg, double step_size, double mip_dQ_dx = 50000 / units::cm);


}

#endif
