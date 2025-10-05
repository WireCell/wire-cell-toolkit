#ifndef WIRECELL_CLUS_PR_VERTEX
#define WIRECELL_CLUS_PR_VERTEX

#include "WireCellClus/PRCommon.h"
#include "WireCellUtil/Flagged.h"
#include "WireCellClus/PRGraphType.h"

namespace WireCell::Clus::Facade {
    class Cluster;
}

namespace WireCell::Clus::PR {

    /** The flags used to categorize a Vertex
     *
     * These are used in Vertex via the "Flagged" base class (from util).
     */ 
    enum class VertexFlags {
        /// The vertex has no particular category
        kUndefined = 0,
        /// The vertex is determined to location of neutrino interaction.
        kNeutrinoVertex = 1<<1,
        /// Same meaning as WCP's "flag_fit_fix".
        kFitFix = 1<<2,

    };

    /** A PR::Vertex instance represents a connection with one or more PR::Segment intances.
     
        A PR::Vertex has:
     
        - an associated 3D point as well as scalar "fit" information which
          includes a potentially different 3D point.
        
        - a set of possible FLAGS (see VertexFlags and Flagged base class)
     
        - a bare pointer to a Facade::Cluster.  This may be nullptr.  And it can
          be invalid if the user does something to destroy the cluster while
          this object still lives.

        - an optional graph node descriptor (can be "graph_type::null_vertex()")
      
        A PR::Vertex must be constructed through `make_vertex()`.
      
        Note, a PR::Vertex is analogous to the ProtoVertex of WCP.
     */
    class Vertex
        : public Flagged<VertexFlags> // can set flags
        , public Graphed<node_descriptor> // may live in a graph
        , public HasCluster<Segment>      // has an associated Cluster*.
    {
    public:        

        /// Getters

        /// The "initial" point information.
        const WCPoint& wcpt() const { return m_wcpt; }
        WCPoint& wcpt() { return m_wcpt; }
        
        /// Information about this vertex provided by some "fit".
        const Fit& fit() const { return m_fit; }
        Fit& fit() { return m_fit; }

        

        /// Chainable setters
        Vertex& wcpt(const WCPoint& wcpt) { m_wcpt = wcpt; return *this; }
        Vertex& fit(const Fit& fit) { m_fit = fit; return *this; }

        // index and range ...
        void set_fit_index(int idx) { m_fit.index = idx; }
        int get_fit_index() const { return m_fit.index; }
        void set_fit_range(double range) { m_fit.range = range; }
        double get_fit_range() const { return m_fit.range; }

        //
        void reset_fit_prop(){
            // Clear the kFitFix flag using keep_flags with inverted mask
            this->keep_flags(static_cast<VertexFlags>(~static_cast<int>(VertexFlags::kFitFix)));
            m_fit.reset();
        }

        // Distance from "initial point" to fit point.
        double fit_distance() {
            return m_fit.distance(m_wcpt.point);
        };
        
        double get_dis(WireCell::Point point){
            return m_fit.distance(point);
        };

    private:

        WCPoint m_wcpt;         // "initial" point.
        Fit m_fit;


    };

    /// See PRGraph.h for related functions.
}
#endif

