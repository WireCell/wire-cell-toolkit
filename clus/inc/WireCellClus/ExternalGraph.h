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
    };

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
