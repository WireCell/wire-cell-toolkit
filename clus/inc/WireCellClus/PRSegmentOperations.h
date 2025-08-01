#ifndef WIRECELL_CLUS_PRSEGMENTOPERATIONS
#define WIRECELL_CLUS_PRSEGMENTOPERATIONS

#include "WireCellClus/PRGraph.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

namespace WireCell::Clus::PR {

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


}

#endif
