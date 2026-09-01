#ifndef WIRECELL_CLUS_PR_SHOWER
#define WIRECELL_CLUS_PR_SHOWER

#include "WireCellClus/PRCommon.h"
#include "WireCellClus/PRTrajectoryView.h"
#include "WireCellClus/PRVertex.h"
#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/Facade_Cluster.h"

#include "WireCellUtil/Flagged.h"
#include "WireCellUtil/Point.h"

#include "WireCellIface/IRecombinationModel.h"
#include "WireCellClus/ParticleDataSet.h"

#include <map>
#include <set>

namespace WireCell::Clus::PR {

    /** The "flags" that may be set on a shower.
     */
    enum class ShowerFlags {
        /// The segment has no particular category.
        kUndefined = 0,
        /// The shower flag
        kShower = 1<<1,
        /// The kinematics flag
        kKinematics = 1<<2,
    };

    /** The data attributes of a PR::Shower
        
        The original WCShower has a huge number of attributes that are merely
        carried by the shower.  This struct factors out that data to facilitate
        writing stand-alone functions.

        Anything that is not part of the core shower-as-graph-view concept gets
        lumped into this ShowerData struct.

        Note, "flags" are set via Shower's Flagged base class.

    */
    struct ShowerData
    {
        int particle_type;
        double kenergy_range;
        double kenergy_dQdx;
        double kenergy_charge;
        double kenergy_best;
    
        WireCell::Point start_point;
        WireCell::Point end_point;
        WireCell::Vector init_dir;

        /// 1 for direct connection, 2 for indirection connection with a gap, 3
        /// for associations (not clear if this should be connected or not
        int start_connection_type;
    };


    /** Model a shower-like view of a trajectory.

        This is the WCT equivalent to a WCT WCShower.
     */
    class Shower
        : public TrajectoryView
        , public Flagged<ShowerFlags> // can set flags
        , public HasDPCs<Shower>      // has associated DynamicPointClouds.
    {
    public:

        Shower(Graph& graph);
        
        virtual ~Shower();

        // The bag of attributes is directly exposed to user.
        ShowerData data;


        // Getters

        /** Get the vertex that is considered the start of the shower.
         */
        VertexPtr start_vertex();

        /** Get the segment that is considered the start of the shower.
         */
        SegmentPtr start_segment();

        // Chainable setters

        /** Chainable setter of start vertex.

            The vertex must already be added to the underlying graph that this
            shower views.

            The vertex will be added to the view.

            The vertex will replace any existing start vertex and not remove the
            prior vertex from the shower's view.  Use `Shower::remove_vertex()`
            to explicitly remove from the view and `PR::remove_vertex()` to
            remove it from the underlying graph, if either are required.

            If the vertex lacks a valid descriptor, eg has yet to be added to
            the underlying graph, this function is a no-op and the stored
            start_vertex is nullified.
        */
        Shower& set_start_vertex(VertexPtr vtx, int type);
        std::pair<VertexPtr, int> get_start_vertex_and_type() {
            return std::make_pair(m_start_vertex, data.start_connection_type);
        }

        /** Chainable setter of start segment.

            This has the same semantics and caveats as the chainable setter:
            `start_vertex(VertexPtr)`.
        */
        Shower& set_start_segment(SegmentPtr seg, bool flag_include_vertices = false, const std::string& cloud_name_fit = "fit", const std::string& cloud_name_associate = "associate_points");

        Shower& set_start_point( WireCell::Point pt ) {
            data.start_point = pt;
            return *this;
        }
        WireCell::Point get_start_point() const {
            return data.start_point;
        }
        WireCell::Point get_end_point() const {
            return data.end_point;
        }

        void add_segment(SegmentPtr seg, bool flag_include_vertices = false, const std::string& cloud_name_fit = "fit", const std::string& cloud_name_associate = "associate_points");

        // stable shower ID (assigned at construction, unique per run)
        int get_shower_id() const { return m_shower_id; }

        // particle type
        int get_particle_type(){return data.particle_type;};
        void set_particle_type(int val){data.particle_type = val;};

        // access flags
        void set_flag_kinematics(bool val);
        bool get_flag_kinematics();
        bool get_flag_shower();

        // return kinematic energy estimation
        double get_kine_range(){return data.kenergy_range;};
        double get_kine_dQdx(){return data.kenergy_dQdx;};
        void set_kine_charge(double val){data.kenergy_charge = val;};
        double get_kine_charge(){return data.kenergy_charge;};
        // doc pr/101 (K3): hadronic showers publish sum dE/dx as best.
        void set_kine_best(double val){data.kenergy_best = val;};
        double get_kine_best(){  
            if (data.kenergy_best != 0) return data.kenergy_best; else return data.kenergy_charge;};

        // return initial direction ...
        WireCell::Vector& get_init_dir(){return data.init_dir;};

        // Get point cloud by name (convenience wrapper around dpcloud)
        std::shared_ptr<Facade::DynamicPointCloud> get_pcloud(const std::string& cloud_name = "fit") {
            return this->dpcloud(cloud_name);
        }
        std::shared_ptr<const Facade::DynamicPointCloud> get_pcloud(const std::string& cloud_name = "associate_points") const {
            return this->dpcloud(cloud_name);
        }

        // Add all segments and vertices from another shower to this one
        void add_shower(Shower& temp_shower, const std::string& cloud_name_fit = "fit", const std::string& cloud_name_associate = "associate_points");
        // absorb_track_guard (doc pr/40 round 6 F12): when true, the
        // flood-fill skips (and terminates the walk at) a confidently PID'd
        // non-electron segment that is long and straight
        // (segment_is_straight_long_track); long-muon pseudo-showers
        // (get_particle_type()==13) are exempt.  false = legacy flood-fill =
        // byte-identical.  See the comment in the implementation.
        // walk_visited_parity (doc pr/91 round 3): the flood-fill's frontier
        // test is `!has_node(v)` -- a vertex already present in the view is
        // never enqueued, even if it was added by set_start_vertex() and its
        // neighborhood was never scanned.  The prototype's
        // WCPPID::WCShower::set_start_vertex (WCShower.cxx:529-532) assigns
        // start_vertex and start_connection_type ONLY -- map_vtx_segs (the
        // toolkit's has_node() equivalent for this purpose) is untouched --
        // so a former start vertex is invisible to the prototype's own
        // frontier test (map_vtx_segs.find(vtx) != end(), WCShower.cxx:735)
        // and gets enqueued fresh on the very next complete_structure call.
        // Instrumented on SBND nueCC evt 168596: shower 2 (the 2039 MeV
        // electron) was seeded on vertex 14027 via set_start_vertex, then
        // re-seated onto the main vertex by examine_showers
        // (NeutrinoShowerClustering.cxx:3114-3119), which re-floods with a
        // FRESH used_segments set -- but 14027 is already a view node, so the
        // re-flood never crosses it, the 4.74 cm proton stub hanging off it
        // is never examined, and a 7.7 cm electron stub 96% inside this
        // shower's own point cloud is left as a separate shower that later
        // pairs into a spurious pi0.  When true, m_walked_nodes (populated
        // unconditionally below, see .cxx) replaces has_node() as the
        // frontier test: any node added to the view but never actually
        // scanned by a flood-fill worklist is re-enqueued.  The CURRENT
        // start vertex is still never walked (the `!= m_start_vertex` guards
        // are unchanged) -- only FORMER start vertices and other view-only
        // insertions become reachable.  Default false = legacy
        // has_node()-gated frontier, byte-identical.
        void complete_structure_with_start_segment(IndexedSegmentSet& used_segments, const std::string& cloud_name_fit = "fit", const std::string& cloud_name_associate = "associate_points", bool absorb_track_guard = false, bool walk_visited_parity = false);


        // get the information from the shower.
        // exclude_start_vertex (doc pr/38): omit m_start_vertex from
        // used_vertices.  The prototype's WCShower::fill_sets reads
        // map_vtx_segs, which NEVER holds the start vertex -- every write
        // path skips it (`if (vtx == start_vertex) continue;`,
        // WCShower.cxx:547 and complete_structure_with_start_segment
        // :708-716/:733-745) -- while the toolkit view holds it via
        // set_start_vertex -> add_vertex.  Default false = legacy view
        // semantics, byte-identical for all existing callers.
        void fill_sets(IndexedVertexSet& used_vertices, IndexedSegmentSet& used_segments, bool flag_exclude_start_segment = true,
                       bool exclude_start_vertex = false);
        void fill_point_vector(std::vector<WireCell::Point>& points, bool flag_main = true);
        TrajectoryView& fill_maps();

        int count_connected_segments(SegmentPtr seg);
        std::pair<IndexedVertexSet, IndexedSegmentSet> get_connected_pieces(SegmentPtr seg);

        std::pair<SegmentPtr, VertexPtr> get_last_segment_vertex_long_muon(IndexedSegmentSet& segments_in_muons);

        // some simple get functions
        int get_num_main_segments();
        int get_num_segments();

        double get_total_length();
        double get_total_length(Facade::Cluster* cluster);
        double get_total_track_length();

        std::vector<double> get_stem_dQ_dx(VertexPtr vertex, SegmentPtr segment, int limit = 20, double mip_dqdx_median = 43000/units::cm);

        // calculate the kinematics
        // doc sbnd_xin/docs/pr/40 round 3: main_vertex + protect_proton_daughter_pion
        // (C++ defaults nullptr/false = legacy = byte-identical) let this
        // function's own shower/track majority-vote reassignment (see .cxx)
        // decline to overwrite m_start_segment when set_default_shower_
        // particle_info (NeutrinoPatternBase.cxx) already relabelled it pion
        // via the same segment_has_proton_daughter test.  Without this guard,
        // this function silently reverted that relabel -- traced end-to-end
        // (SBND evt 256587 seg 11079) in round 2; see porting_dictionary.md.
        // proton_daughter_mip_dqdx is DELIBERATELY separate from mip_dqdx
        // (this function's own 4-momentum-calc scale, bound to the caller's
        // m_mip_dqdx=50000/units::cm): the guard must re-derive
        // segment_has_proton_daughter's verdict on the SAME scale F5 used
        // (m_mip_dqdx_median=43000/units::cm) or it can silently disagree
        // with F5's own decision at a different threshold.
        // vote_track_pid_counts: doc sbnd_xin/docs/pr/93 Cause C knob
        // (shower_vote_track_pid_counts).  C++ default false = legacy vote
        // (only confirmed protons count as track).
        // accept_pid_guard: doc sbnd_xin/docs/pr/93 Cause B knob
        // (shower_accept_pid_guard) -- the same knob that declines the forced
        // set_pdg(11) at the two acceptance sites also declines THIS vote's
        // final overwrite when the start segment carries a confident
        // non-electron template PID (segment_confident_nonelectron_pid);
        // without it the vote re-flips the segment one line after the
        // acceptance-site guard spared it.  C++ defaults false = legacy.
        void update_particle_type(const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double mip_dqdx = 50000/units::cm, VertexPtr main_vertex = nullptr, bool protect_proton_daughter_pion = false, double proton_daughter_mip_dqdx = 43000/units::cm, bool vote_track_pid_counts = false, bool accept_pid_guard = false, double accept_pid_min_len = 50*units::cm);
        // exclude_start_vertex_from_endpoint (doc pr/39): same prototype-parity
        // rule as fill_sets's exclude_start_vertex above, applied to the
        // farthest-vertex search that sets data.end_point.  The prototype's
        // WCShower::calculate_kinematics (WCShower.cxx:377-387) searches
        // map_vtx_segs, which never holds the start vertex; the toolkit search
        // here walks the full node set with no such exclusion, so a detached
        // (conn_type 2/3) shower's end_point can collapse onto its own start
        // vertex (e.g. the neutrino vertex) instead of growing away from it.
        // Default false = legacy search over every node, byte-identical.
        //
        // endpoint_skip_orphan_vertices (doc pr/91 round 1): the same search
        // must also skip a node that NO member segment of this shower touches.
        // Such "orphan" nodes exist because set_start_vertex() calls
        // add_vertex() (the toolkit-only divergence recorded at fill_sets
        // above), so a conn-2/3 shower's view holds a vertex belonging to
        // somebody else's cluster.  exclude_start_vertex_from_endpoint hides it
        // only while THIS shower still calls it its start vertex; add_shower()
        // then imports it wholesale into an absorber where that exclusion no
        // longer applies, and it wins the farthest-vertex search.  Measured on
        // SBND 169626/174752/347129/394532: 6 absorbs import a foreign vertex
        // and 5 become the reported end_point -- 394532's 30 MeV and 66 MeV
        // showers end on each other's charge.  end_point feeds the Bee PF node
        // and the cal_dir_3vector/angle consumers; kine_* sum over member
        // segments and are untouched either way.
        // Default false = legacy search, byte-identical.
        void calculate_kinematics(const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool exclude_start_vertex_from_endpoint = false, bool endpoint_skip_orphan_vertices = false);
        // doc pr/101 (K4): best_mode 0 = legacy dQdx, 1 = range over the
        // muon chain, 2 = range when the far muon vertex is a dead-end and
        // dQdx/range is within [1-ratio_lo, 1+ratio_hi], else dQdx.
        void calculate_kinematics_long_muon(IndexedSegmentSet& segments_in_muons, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool exclude_start_vertex_from_endpoint = false, int best_mode = 0, double ratio_lo = 0.3, double ratio_hi = 0.5, bool range_empty_chain_fallback = false, bool members_geometry = false);

        // doc sbnd_xin/docs/pr/93 round 4 (shower_detach_track_stem).
        // Remove `prefix` -- a connected chain of track segments beginning at
        // the CURRENT start segment -- from this shower's view, re-root the
        // shower at `new_start_vertex` with connection type 2 (the
        // pseudo-gamma / disconnected-shower rendering class), re-seat the
        // start segment on the remaining member closest to the new start
        // vertex (graph-index tie-break), and REBUILD the named point clouds
        // from the remaining members only.  The rebuild is required because
        // the shower clouds are add-only merges of member clouds and both
        // kine_charge (kine_charge_from_maps over shower pclouds) and the
        // conn-2 start_point derivation (closest "fit"-cloud point to the
        // start vertex) read them -- naive edge removal would leave the
        // detached track's charge inside the daughter's energy.
        // Prefix-only vertices (including the old root) are removed from the
        // view so they cannot win the farthest-vertex end_point search; any
        // vertex still touched by a remaining member is kept.
        // Refuses (returns 0) when the prefix is empty, would empty the
        // shower, or new_start_vertex is invalid.  Caller must re-run
        // update_particle_type / calculate_kinematics / set_kine_charge.
        // Returns the number of segments peeled.  m_shower_id is untouched.
        int detach_track_prefix(const std::vector<SegmentPtr>& prefix,
                                VertexPtr new_start_vertex,
                                const std::string& cloud_name_fit = "fit",
                                const std::string& cloud_name_associate = "associate_points");

        // doc sbnd_xin/docs/pr/99 round 2 (shower_ghost_member_drop).
        // Remove one interior/peripheral member -- a projective-ghost
        // segment (charge-starved 2D duplicate of a sibling member, the
        // pr/83 op1-proj class found INSIDE shower membership) -- from this
        // shower's view.  Leaf-only contract: refuses (returns 0) when the
        // ghost is the start segment, when removal would empty the shower,
        // or when any remaining member would be stranded (view-restricted
        // connectivity from the start segment must still reach every
        // remaining member -- the PRShower.cxx:1586 nsegments==nconnected
        // branch and the stranded-member kine sums are the hazards this
        // guard exists for).  On success: ghost-only vertices leave the
        // view (walked marks erased), the named point clouds are REBUILT
        // from the remaining members (same reason as detach_track_prefix:
        // kine_charge reads the clouds, and a ghost's negative charge must
        // leave the energy), caches invalidated, flag_kinematics cleared.
        // Caller owns any full-graph edit (the ghost segment itself) and
        // must re-run update_particle_type / calculate_kinematics /
        // set_kine_charge / update_shower_maps.  Bookkeeping forked BY
        // DUPLICATION from detach_track_prefix (production method stays
        // byte-untouched).  Returns 1 on removal, 0 on refusal.
        int drop_ghost_member(SegmentPtr ghost,
                              const std::string& cloud_name_fit = "fit",
                              const std::string& cloud_name_associate = "associate_points");

        // doc sbnd_xin/docs/pr/123 round 1 (shower_pass4_prune_detached).
        // Remove a SET of members -- a spatially detached component -- from
        // this shower's view WITHOUT re-rooting: start vertex, connection
        // type and start segment stay.  Refuses (returns 0) when the set is
        // empty, contains the start segment, includes a non-member, or would
        // empty the shower.  Vertices touched only by removed members leave
        // the view (walked marks erased); the named point clouds are REBUILT
        // from the remaining members (kine_charge reads the clouds).  No
        // connectivity requirement on the remainder -- unlike
        // drop_ghost_member's leaf-only contract, the caller removes a whole
        // spatial component and re-seeds it as its own shower.  Caller must
        // re-run calculate_kinematics / set_kine_charge / update_shower_maps.
        // Bookkeeping forked BY DUPLICATION from detach_track_prefix
        // (production method stays byte-untouched).  Returns the number of
        // segments removed.
        int detach_member_set(const std::vector<SegmentPtr>& members,
                              const std::string& cloud_name_fit = "fit",
                              const std::string& cloud_name_associate = "associate_points");

        // doc pr/99 round 3 (C1b, kine_charge_rebuild).  Build and return an
        // EPHEMERAL point cloud of the named kind ("fit" or
        // "associate_points") from the shower's CURRENT members only --
        // prototype parity with WCShower::rebuild_point_clouds(), which WCP
        // calls before every cal_kine_charge cloud read
        // (NeutrinoID_energy_reco.h:99).  The stored add-only clouds are
        // deliberately NOT touched: taggers and the pi0 start_point
        // derivation query them after the energy recompute, and a rebuild
        // changes row order hence kNN tie-breaks.  Loop forked BY
        // DUPLICATION from the detach_track_prefix rebuild (production
        // method stays byte-untouched).  Returns nullptr when no member
        // carries such a cloud.
        std::shared_ptr<Facade::DynamicPointCloud>
        rebuild_pcloud(const std::string& cloud_name);

    private:

        Graph& m_full_graph;
        VertexPtr m_start_vertex;
        SegmentPtr m_start_segment;
        int m_shower_id{-1};

        // doc pr/91 round 3: the prototype's map_vtx_segs key set, restated
        // as "has this vertex's neighborhood actually been scanned by a
        // flood-fill worklist" -- distinct from view MEMBERSHIP (has_node()),
        // which set_start_vertex()/set_start_segment()/add_segment()/
        // add_shower() can grow without ever visiting a vertex's edges.
        // Populated unconditionally wherever those methods already add a
        // vertex (pure bookkeeping, membership-tested only, never iterated --
        // no pointer-order dependence); read only by
        // complete_structure_with_start_segment's walk_visited_parity branch.
        std::set<node_descriptor> m_walked_nodes;

        // Lazy caches — invalidated whenever segments are added.
        mutable double m_total_length_cache{-1.0};
        mutable int    m_num_main_segs_cache{-1};

        void invalidate_segment_caches() {
            m_total_length_cache  = -1.0;
            m_num_main_segs_cache = -1;
        }

    };

    using ShowerPtr = std::shared_ptr<Shower>;

    /// Shower::m_shower_id before set_shower_id() has been called.  Same hazard
    /// as PR::kUnindexed but a DIFFERENT sentinel (-1, not SIZE_MAX): showers
    /// are numbered by the clustering pass, not by graph insertion.
    constexpr int kUnassignedShowerId = -1;

    struct ShowerIndexCmp {
        bool operator()(const ShowerPtr& a, const ShowerPtr& b) const {
            const int ia = a->get_shower_id();
            const int ib = b->get_shower_id();
            if (ia == kUnassignedShowerId || ib == kUnassignedShowerId) [[unlikely]] {
                static std::atomic<bool> warned{false};
                if (!warned.exchange(true)) warn_unindexed("PR::Shower");
            }
            return ia < ib;
        }
    };
    using IndexedShowerSet   = std::set<ShowerPtr, ShowerIndexCmp>;
    using ShowerVertexMap    = std::map<VertexPtr,  ShowerPtr, VertexIndexCmp>;
    using ShowerSegmentMap   = std::map<SegmentPtr, ShowerPtr, SegmentIndexCmp>;
    using VertexShowerSetMap = std::map<VertexPtr,  IndexedShowerSet, VertexIndexCmp>;
    using ShowerIntMap       = std::map<ShowerPtr,  int, ShowerIndexCmp>;

    struct ClusterPtrCmp {
        bool operator()(Facade::Cluster* a, Facade::Cluster* b) const {
            return a->get_cluster_id() < b->get_cluster_id();
        }
    };
    using ClusterPtrSet    = std::set<Facade::Cluster*, ClusterPtrCmp>;
    using ClusterVertexMap = std::map<Facade::Cluster*, VertexPtr, ClusterPtrCmp>;

}
#endif
