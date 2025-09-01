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
    double segment_track_length(SegmentPtr seg, int flag = 0);
    
    /// Calculate median dQ/dx for a segment
    ///
    /// Extracts dQ and dx from segment's fits and calculates median dQ/dx.
    ///
    /// @param seg The segment containing fit data
    /// @return Median dQ/dx value (0 if no valid fits)
    double segment_median_dQ_dx(SegmentPtr seg);
    
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
    double segment_geometric_length(SegmentPtr seg);

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

    // Helper functions for external data (compatibility with existing workflows)

    // /// Calculate median dQ/dx from external data vectors
    // ///
    // /// @param dQ Vector of charge values
    // /// @param dx Vector of distance values
    // /// @return Median dQ/dx value
    // double calculate_median_dQ_dx(const std::vector<double>& dQ, 
    //                               const std::vector<double>& dx);
    
    // /// Calculate track length above dQ/dx threshold from external data vectors
    // ///
    // /// @param dQ Vector of charge values
    // /// @param dx Vector of distance values  
    // /// @param threshold dQ/dx threshold value
    // /// @return Length of track segments above threshold
    // double calculate_track_length_threshold(const std::vector<double>& dQ,
    //                                        const std::vector<double>& dx,
    //                                        double threshold);

}

#endif
