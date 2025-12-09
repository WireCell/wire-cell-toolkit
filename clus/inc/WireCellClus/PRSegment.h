#ifndef WIRECELL_CLUS_PR_SEGMENT
#define WIRECELL_CLUS_PR_SEGMENT

#include "WireCellClus/PRCommon.h"
#include "WireCellClus/PRGraphType.h"
#include "WireCellUtil/Flagged.h"
#include "WireCellAux/ParticleInfo.h"
#include "WireCellIface/IDetectorVolumes.h"
#include <memory>


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


    /** A segment represents a connection between vertices in a larger trajectory.
     *
     * A segment has:
     *
     * - a vector of associated 3D point and corresponding "fit" information
     *   which includes a potentially different 3D point.
     *
     * - a set of possible FLAGS (see SegmentFlags and Flagged base class)
     *
     * - a bare pointer to a Facade::Cluster.  This may be nullptr.  And it can
     *   be invalid if the user does something to destroy the cluster while this
     *   object still lives.
     *
     * - a generic graph edge descriptor (see Graphed base class and PR::Vertex
     *   commentary for more information on the nature this descriptor).
     * 
     * Caution, although the points held by the segment should be ordered
     * "along" the segment, the graph edge representing the segment is NOT
     * ORDERED.  This means that the point in the vertex found at the `source()`
     * node for the segments edge is not necessarily closest to the segment's
     * first point.  See `find_endpoints()` for one way to resolve this
     * directional ambiguity.
     *
     * Note, a PR::Segment is essentially the ProtoSegment of WCP.
     */
    class Segment
    : public Flagged<SegmentFlags> // can set flags
    , public Graphed<edge_descriptor> // may live in a graph
    , public HasCluster<Segment>      // has an associated Cluster*.
    , public HasDPCs<Segment>      // has associated DynamicPointClouds.
    , public std::enable_shared_from_this<Segment>  // allows shared_from_this()
    {
    public:

        // Getters
        const std::shared_ptr<Aux::ParticleInfo>& particle_info() const { return m_particle_info; }
        std::shared_ptr<Aux::ParticleInfo>& particle_info() { return m_particle_info; }

        // Chainable setter
        Segment& particle_info(std::shared_ptr<Aux::ParticleInfo> pinfo) {
            m_particle_info = pinfo;
            return *this; 
        }
        
        // Convenience method to check if particle info is available
        bool has_particle_info() const { return m_particle_info != nullptr; }


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
        bool dir_weak() const { return m_dir_weak; }

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
        Segment& dir_weak(bool weak) { m_dir_weak = weak; return *this; }

        // reset fit properties ...
        void reset_fit_prop(); 
        int fit_index(int i);
        void fit_index(int i, int idx);
        bool fit_flag_skip(int i);
        void fit_flag_skip(int i, bool flag);

        void set_fit_associate_vec(std::vector<PR::Fit >& tmp_fit_vec, const IDetectorVolumes::pointer& dv,const std::string& cloud_name="fit");
        
        
    private:

        std::vector<WCPoint> m_wcpts;
        std::vector<Fit> m_fits;

        int m_dirsign{0};
        bool m_dir_weak{false};

        std::shared_ptr<Aux::ParticleInfo> m_particle_info{nullptr};



        // Still must consider adding:
        // + pcloud_fit
        // + pcloud_associated
        // - but NOT pcloud_associated_steiner as it is never used
    };

}
#endif
