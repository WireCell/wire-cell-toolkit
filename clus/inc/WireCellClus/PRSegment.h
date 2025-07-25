#ifndef WIRECELL_CLUS_PR_SEGMENT
#define WIRECELL_CLUS_PR_SEGMENT

#include "WireCellClus/PRCommon.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellUtil/Flagged.h"

namespace WireCell::Clus::PR {

    /** The flags used to categorize a segment.
     *
     * These are used by the Segment class via its Flagged base (from util/).
     *
     */
    enum SegmentFlags {
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


    /** A segment represents a connection between vertices in a larger trajectory.
     *
     * A segment has:
     *
     * - a vector of associated 3D point and corresponding "fit" information
     * which includes a potentially different 3D point.
     *
     * - a set of possible FLAGS (see SegmentFlags and Flagged base class)
     *
     * - an ID and the ID of an associated cluster (see Identities base class)
     *
     * - can live in a graph as an edge (see Graphed base class).
     * 
     * A PR::Segment is essentially the ProtoSegment of WCP.
     */
    class Segment
    : public Flagged<SegmentFlags>         // can set flags
    , public Identities                    // hold id and cluster_id
    , public Graphed<edge_descriptor_type> // live in a graph
    {
    public:


        /// Getters

        /** Return the fits for associated points.
         *
         * Note, check the kFit flag or the valid() on individual Fit value to
         * determine if they are provide.
         */
        const std::vector<Fit>& fits() const { return m_fits; }
        std::vector<Fit>& fits() { return m_fits; }

        /// The "initial" point information.
        const std::vector<WCPoint>& wcpts() const { return m_wcpts; }
        std::vector<WCPoint>& wcpts() { return m_wcpts; }

        /// The directional sign (was "flag_dir" in WCT).
        int dirsign() const { return m_dirsign; }

        /// Base provides these getters:
        using Identities::ident;
        using Identities::cluster_id;

        /// Chainable setters

        Segment& wcpts(const std::vector<WCPoint>& wcpts) { m_wcpts = wcpts; return *this; }
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
