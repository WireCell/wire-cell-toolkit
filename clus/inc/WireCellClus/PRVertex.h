#ifndef WIRECELL_CLUS_PR_VERTEX
#define WIRECELL_CLUS_PR_VERTEX

#include "WireCellClus/PRCommon.h"
#include "WireCellUtil/Flagged.h"
#include "WireCellClus/PRGraphType.h"

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
     
        - an individual ID and the ID of an associated cluster (see Identities
          base class).  Note, the ID is not the same as the index for the
          associated graph node.
     
        - an optional graph node descriptor (can be "graph_type::null_vertex()")
      
        A PR::Vertex must be constructed through `make_vertex()`.
      
        Note, a PR::Vertex is analogous to the ProtoVertex of WCP.
     */
    class Vertex
        : public Flagged<VertexFlags> // can set flags
        , public Identities           // hold id and cluster_id
        , public Graphed<node_descriptor> // may live in a graph
    {
    public:        

        /// Getters

        /// The "initial" point information.
        const WCPoint& wcpt() const { return m_wcpt; }
        WCPoint& wcpt() { return m_wcpt; }
        
        /// Information about this vertex provided by some "fit".
        const Fit& fit() const { return m_fit; }
        Fit& fit() { return m_fit; }

        /// Base provides these getters:
        using Identities::ident;
        using Identities::cluster_id;

        /// Chainable setters

        Vertex& wcpt(const WCPoint& wcpt) { m_wcpt = wcpt; return *this; }
        Vertex& fit(const Fit& fit) { m_fit = fit; return *this; }
        Vertex& ident(int id) { Identities::ident(id); return *this; }
        Vertex& cluster_id(int cid) { Identities::cluster_id(cid); return *this; }

        // Distance from "initial point" to fit point.
        double fit_distance() {
            return m_fit.distance(m_wcpt.point);
        };

    private:

        WCPoint m_wcpt;         // "initial" point.
        Fit m_fit;

    };

    /// See PRGraph.h for related functions.
}
#endif

