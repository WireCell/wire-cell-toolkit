#ifndef WIRECELL_CLUS_PR_TRAJECTORY
#define WIRECELL_CLUS_PR_TRAJECTORY

#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRVertex.h"

namespace WireCell::Clus::PR {

    /** A PR::Trajectory is a base class for a model of a particle trajectory.
     *
     * A "trajectory" describes one sub-graph of PR::Vertex nodes and
     * PR::Segment edges that is in some "covering" the reconstructed activity
     * for a single identified object.
     *
     * This base class provides features and behavior common to all types of.
     * Subclasses may be define to provide additional features or variant
     * behavior.  In particular, see subclass PR::Shower or PR::Track.
     *
     * The PR::Trajectory knows its sub-graph 
     * Through a PR::Trajectory the local sub-graph may be queried or mutated.  
     */
    class Trajectory {
    public:

        // FIXME: open question: do we have Vertex and Segment hold their
        // descriptors, forcing graph monogamy or do we add maps from these pointers
        // to their descriptors as graph properties and write a bunch of boost::*
        // graph function wrappers to assure maps and graph are synced?

        
        graph_type& graph() { return m_graph; }
        const graph_type& graph() const { return m_graph; }
        

    private:
        graph_type& m_graph;
    };

}

#endif
