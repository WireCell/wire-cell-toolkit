#ifndef WIRECELL_CLUS_PR_SEGMENT
#define WIRECELL_CLUS_PR_SEGMENT

#include "WireCellClus/PRCommon.h"
#include "WireCellClus/PRGraphType.h"
#include "WireCellUtil/Flagged.h"

namespace WireCell::Clus::PR {

    /** The flags used to categorize a segment.
     *
     * These are used by the Segment class via its Flagged base (from util/).
     *
     */
    enum class SegmentFlags {
        /// The segment has no particular category.
        kUndefined = 0,
        /// The segment is part of a shower trajectory
        kShowerTrajectory = 1<<1,
        /// The segment is part of a shower topology
        kShowerTopology = 1<<2,
        /// The segment should not have a muon check.
        kAvoidMuonCheck = 1<<3,
        /// The fits are provided.
        kFit = 1<<4,
    };


    /// The segment caries two types of PCs each of which span 3d and 3 2D PCs
    struct SegmentPointClouds {
    };

    /** A segment represents a connection between vertices in a larger trajectory.
     *
     * A segment has:
     *
     * - a vector of associated 3D point and corresponding "fit" information
     *   which includes a potentially different 3D point.
     *
     * - a set of possible FLAGS (see SegmentFlags and Flagged base class)
     *
     * - an ID and the ID of an associated cluster (see Identities base class)
     *
     * - a generic graph edge descriptor (see Graphed base class and PR::Vertex
     *   commentary for more information on the nature this descriptor).
     * 
     * Caution, although the points held by the segment should be ordered
     * "along" the segment, the graph edge representing the segment is NOT
     * ORDERED.  This means that the point in the vertex found at the `source()`
     * node for the segments edge is not necessarily closest to the segment's
     * first point.
     *
     * Note, a PR::Segment is essentially the ProtoSegment of WCP.
     */
    class Segment
    : public Flagged<SegmentFlags> // can set flags
    , public Identities            // hold id and cluster_id
    , public Graphed<edge_descriptor> // may live in a graph
    {
    public:

        // Getters

        /// Get the const vector of fits.
        const std::vector<Fit>& fits() const { return m_fits; }

        /// Get the mutable vector of fits.
        std::vector<Fit>& fits() { return m_fits; }
        
        /// Get the const original points.
        const std::vector<WCPoint>& wcpts() const { return m_wcpts; }

        /// Get the mutable original points.
        std::vector<WCPoint>& wcpts() { return m_wcpts; }

        /// Get the sign +1/0/-1 (was "flag_dir" in WCT).
        int dirsign() const { return m_dirsign; }

        // Base provides these getters:
        using Identities::ident;
        using Identities::cluster_id;

        // Chainable setters

        /// Chainable setter of fits vector.
        Segment& fits(const std::vector<Fit>& f) { m_fits = f; return *this; };
        /// Chainable setter of original points vector.
        Segment& wcpts(const std::vector<WCPoint>& pts) { m_wcpts = pts; return *this; }
        /// Chainable setter of dirsign.
        Segment& dirsign(int dirsign) {
            if (dirsign == 0) m_dirsign = 0;
            else m_dirsign = dirsign > 0 ? 1 : -1;
            return *this;
        }

        // Override base so we return self type.

        /// Set our ident. 
        Segment& ident(int id) { Identities::ident(id); return *this; }

        /// Set ID of related cluster.
        Segment& cluster_id(int cid) { Identities::cluster_id(cid); return *this; }
        

    private:

        std::vector<WCPoint> m_wcpts;
        std::vector<Fit> m_fits;

        int m_dirsign{0};

        // Still must consider adding:
        // + pcloud_fit
        // + pcloud_associated
        // - but NOT pcloud_associated_steiner as it is never used
    };

}
#endif
