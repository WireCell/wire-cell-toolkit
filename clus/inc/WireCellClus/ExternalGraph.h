/** pr_ml doc 79 (P3): a segment graph built OUTSIDE the toolkit (for example from the E-1 network's outputs), read from a JSON file
 * and turned into a PR::Graph for the trajectory + dQ/dx fit.  Used only by TaggerCheckNeutrino's pr_fit_only = "external" mode
 * (default off); nothing else calls it.
 *
 * File format (version 1), positions in mm in the reconstruction's point frame:
 *
 *   { "version": 1, "units": "mm",
 *     "vertices": [ {"id": 0, "p": [x, y, z]}, ... ],
 *     "segments": [ {"id": 0, "v": [v0, v1], "pts": [[x, y, z], ...], "type": "track"}, ... ],
 *     "clumps":   [ {"id": 0, "pts": [[x, y, z], ...], "type": "shower"}, ... ] }
 *
 * A segment's first and last points must equal its two vertices' positions (to 0.01 mm); a segment has at least 2 points;
 * vertex ids are unique.  read_external_graph() refuses a file that breaks any of these (ValueError).
 */
#ifndef WIRECELLCLUS_EXTERNALGRAPH
#define WIRECELLCLUS_EXTERNALGRAPH

#include "WireCellClus/PRGraph.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellIface/IDetectorVolumes.h"
#include "WireCellUtil/Point.h"

#include <string>
#include <vector>

namespace WireCell::Clus::PR {

    class PatternAlgorithms;   // NeutrinoPatternBase.h (no include guard; included by the .cxx)

    struct ExternalGraph {
        struct Vtx { int id{-1}; WireCell::Point p; };
        struct Seg { int id{-1}; int v0{-1}; int v1{-1}; std::vector<WireCell::Point> pts; std::string type; };
        struct Clump { int id{-1}; std::vector<WireCell::Point> pts; std::string type; };
        std::vector<Vtx> vertices;
        std::vector<Seg> segments;
        std::vector<Clump> clumps;
    };

    /// Parse and check a version-1 file (see the header comment).  Positions are returned in WCT units.
    ExternalGraph read_external_graph(const std::string& path);

    /// What build_external_graph() did, for the EXTGRAPH log line.
    struct ExternalGraphStats {
        int seg_read{0}, seg_used{0}, seg_skipped{0}, seg_other_nearer{0}, seg_refused{0};
        int clump_read{0}, clump_used{0}, clump_other_nearer{0};
        std::vector<double> seg_med_dist;   // per used segment, median distance of its points to the cluster (WCT units)
        std::vector<std::pair<int, SegmentPtr>> used;   // pr_ml doc 92 (O9): (file segment id, segment) per used segment, file order
    };

    /// pr_ml doc 92 (O9 / O13, default off): ONE assignment of every segment and clump of the file to its nearest cluster.
    /// For each item, the median over its points of the distance to a cluster's default-scope points is computed for every
    /// cluster of `clusters` (the first is always computed; an empty other cluster is passed over); the item goes to the
    /// cluster with the smallest median, ties to the earlier cluster in the list (so to the main when it is first), when that
    /// median is at most max_med_dist; else to none (-1).  With the main first and the other clusters after it, the main's
    /// items are exactly those build_external_graph() takes for the main with the same others (used or refused).
    struct ExternalAssignment {
        std::vector<int> seg_cluster, clump_cluster;   // index into `clusters`, -1 = none
        std::vector<int> seg_nearest, clump_nearest;   // index of the nearest cluster even when farther than max_med_dist (-1: none)
        std::vector<double> seg_md, clump_md;          // that nearest median (WCT units; +inf when no cluster)
    };
    ExternalAssignment assign_external_graph(const ExternalGraph& eg, const std::vector<const Facade::Cluster*>& clusters,
                                             double max_med_dist);

    /// pr_ml doc 92 (O13, default off): build_external_graph() for the items assigned to clusters[index] (duplicated, not
    /// shared: the main keeps build_external_graph).  Segments are created, added to the graph and the fitter, and clumps
    /// handed to the fitter, exactly as there; seg_skipped counts items whose nearest cluster is this one but too far,
    /// seg_other_nearer items whose nearest is another cluster (likewise for clumps).
    ExternalGraphStats build_external_graph_assigned(PatternAlgorithms& pa, Graph& graph, Facade::Cluster& cluster, TrackFitting& tf,
                                                     IDetectorVolumes::pointer dv, const ExternalGraph& eg,
                                                     const ExternalAssignment& as, int index);

    /// Build the graph for ONE cluster: a segment (or clump) is used when this cluster is its nearest cluster -- the median over
    /// its points of the distance to a cluster's default-scope points is smallest for this cluster among it and `others` -- and
    /// that median is at most max_med_dist; the others belong to another cluster (or a wrong frame) and are skipped.
    /// Vertices are created on first use and shared by every segment that names them.  Used segments are added to the fitter
    /// (TrackFitting::add_segment); used clumps are handed to it with TrackFitting::set_charge_clumps (they act only when the
    /// fitter's dqdx_clump_term parameter is on).
    ExternalGraphStats build_external_graph(PatternAlgorithms& pa, Graph& graph, Facade::Cluster& cluster, TrackFitting& tf,
                                            IDetectorVolumes::pointer dv, const ExternalGraph& eg, double max_med_dist,
                                            const std::vector<const Facade::Cluster*>& others = {});

}

#endif
