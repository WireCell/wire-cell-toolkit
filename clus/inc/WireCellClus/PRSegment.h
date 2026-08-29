#ifndef WIRECELL_CLUS_PR_SEGMENT
#define WIRECELL_CLUS_PR_SEGMENT

#include "WireCellClus/PRCommon.h"
#include "WireCellClus/PRGraphType.h"
#include "WireCellUtil/Flagged.h"
#include "WireCellAux/ParticleInfo.h"
#include "WireCellIface/IDetectorVolumes.h"
#include <limits>
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
        /// The segment is an arm of a two-end dQ/dx break (doc
        /// sbnd_xin/docs/pr/48): its travel direction away from the break
        /// vertex was established by the break's own two-arm
        /// stopping-template accept.  determine_direction reconstructs that
        /// outward direction (from this flag plus which endpoint vertex
        /// carries VertexFlags::kProtectedBreak) and lets it stand over a
        /// WEAK direction recompute (a strong recompute still wins).  Set
        /// only by the default-OFF two_end_break pass => byte-identical
        /// when off.
        kTwoEndBreakArm = 1<<5,
        /// The segment was demoted from EM shower to stopping muon by the
        /// Michel-stem guard (doc sbnd_xin/docs/pr/74 round 4): it is the
        /// muon half of a muon+Michel pair at the neutrino vertex.  Read
        /// only by stem_backfill, which must not absorb such a segment back
        /// into the Michel shower it was just separated from.  Set only by
        /// the default-OFF shower_traj_michel_stem pass => byte-identical
        /// when off.  Nothing serialises the raw flags word (only named bits
        /// are ever tested), so a new bit is inert in every output.
        kMuonStemGuard = 1<<6,
        /// The segment was declined by the pass4_angle long-track absorb
        /// guard (doc sbnd_xin/docs/pr/123 round 2): a straight-long
        /// track-like segment (e.g. SBND 18255-171572's 125cm muon) that
        /// the cone would have absorbed into an EM shower.  Read only by
        /// the default-OFF pf_orphan_guard_freed (PF root node) and
        /// kine_count_guard_freed (kine push) passes, so such a track is
        /// not lost from the PF tree / kine_reco_Enu when nothing else
        /// claims it.  Set only when shower_pass4_track_guard_len > 0 =>
        /// byte-identical when that knob is off; nothing serialises the
        /// raw flags word, so the bit alone is inert in every output.
        kPass4GuardFreed = 1<<7,
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
        double particle_score() const { return m_particle_score; }
        void particle_score(double score) { m_particle_score = score; }

        // Chainable setter.  doc sbnd_xin/docs/pr/40: body moved to
        // PRSegment.cxx so it can see the complete Facade::Cluster type for
        // the WCT_PID_WRITE_DEBUG diagnostic (env-gated, same idiom as
        // WCT_SHOWER_TOPO_DEBUG / WCT_DET_DEBUG; proven physics-inert by
        // doc pr/40's G1 byte-identical gate).  __builtin_FILE/LINE capture
        // the caller with zero call-site churn across the ~40 writers.
        Segment& particle_info(std::shared_ptr<Aux::ParticleInfo> pinfo,
                                const char* _caller_file = __builtin_FILE(),
                                int _caller_line = __builtin_LINE());
        
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
        void clear_fit(const IDetectorVolumes::pointer& dv, const std::string& cloud_name="fit");

        int fit_index(int i);
        void fit_index(int i, int idx);
        bool fit_flag_skip(int i);
        void fit_flag_skip(int i, bool flag);

        void set_fit_associate_vec(std::vector<PR::Fit> tmp_fit_vec, const IDetectorVolumes::pointer& dv, const std::string& cloud_name="fit");
        
        // Global indices management for point clouds
        void set_global_indices(const std::string& cloud_name, std::vector<size_t> indices) {
            m_pc_global_indices[cloud_name] = std::move(indices);
        }
        
        const std::vector<size_t>& global_indices(const std::string& cloud_name) const {
            static const std::vector<size_t> empty_vec;
            auto it = m_pc_global_indices.find(cloud_name);
            return (it != m_pc_global_indices.end()) ? it->second : empty_vec;
        }
        
        bool has_global_indices(const std::string& cloud_name) const {
            return m_pc_global_indices.find(cloud_name) != m_pc_global_indices.end();
        }

        bool reset_global_indices();
        bool reset_global_indices(const std::string& cloud_name);

        // Add public accessor:
        int id() const { return m_id; }
        void set_id(int id) { m_id = id; }

        size_t get_graph_index() const { return m_graph_index; }
        void set_graph_index(size_t idx) { m_graph_index = idx; }

    private:
        // Add to PRSegment.h private section:
        int m_id{-1};
        size_t m_graph_index{std::numeric_limits<size_t>::max()};

        std::vector<WCPoint> m_wcpts;
        std::vector<Fit> m_fits;

        int m_dirsign{0};
        bool m_dir_weak{false};

        std::shared_ptr<Aux::ParticleInfo> m_particle_info{nullptr};
        double m_particle_score{100};
        
        // Mapping from local DPC index to global cluster point index
        std::map<std::string, std::vector<size_t>> m_pc_global_indices;



        // Still must consider adding:
        // + pcloud_fit
        // + pcloud_associated
        // - but NOT pcloud_associated_steiner as it is never used
    };

}
#endif
