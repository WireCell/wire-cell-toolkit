#include "WireCellClus/PRShower.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellClus/PRShowerFunctions.h"
#include "WireCellClus/DynamicPointCloud.h"
#include "WireCellUtil/Logging.h"
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <unordered_set>

static auto s_log = WireCell::Log::logger("clus.NeutrinoPattern");
static std::atomic<int> s_shower_id_counter{0};

namespace WireCell::Clus::PR {

    namespace {
        /// Deep-copy a segment's dynamic point cloud so the shower owns its own.
        ///
        /// The seeding branches below used to hand the shower the segment's own
        /// shared_ptr (dpcloud(name, ptr) stores it verbatim, PRCommon.h:78).
        /// Every later merge into the shower then grew the *segment's* cloud:
        /// instrumented on SBND nueCC evt 172230, shower 0's start segment saw
        /// its own "fit" cloud go 44 -> 786 points via 67 merges, i.e. it ended
        /// up holding the entire shower.  Consumers that read a segment's cloud
        /// as a proximity mask (NeutrinoEnergyReco.cxx:266-267,
        /// PRSegmentFunctions.cxx, MultiAlgBlobClustering.cxx:861) would then be
        /// answering for the shower rather than the segment.  DynamicPointCloud
        /// holds a unique_ptr k-d tree so it is not copy-constructible; rebuild
        /// via the public API instead.
        std::shared_ptr<Facade::DynamicPointCloud>
        clone_dpc(const Facade::DynamicPointCloud& src)
        {
            auto out = std::make_shared<Facade::DynamicPointCloud>(src.get_wpid_params());
            out->add_points(src);
            return out;
        }

        // doc sbnd_xin/docs/pr/91 round 1 -- WCT_SHOWER_ENDPOINT_DEBUG.
        //
        // get_end_point() is not a point on the shower's charge: it is the
        // farthest *vertex* of the shower's TrajectoryView from start_point
        // (the two loops below).  The view's node set can therefore decide the
        // end point even when no member segment reaches that vertex, which is
        // how SBND 169626's `e- 107` (cluster 22, members spanning z 396-410)
        // ends up reporting an end at z=461 on a cluster-13 vertex owned by the
        // 567 MeV shower, and how 394532's 30 MeV and 66 MeV showers end up
        // reporting each other's vertices.  Both appeared only after
        // shower_dedup_start_seg (doc pr/84 round 3) started calling
        // add_shower(), which unions the absorbed shower's nodes wholesale.
        // This probe prints every candidate vertex and whether a member segment
        // of this shower actually touches it.  Log/stderr only: no effect on
        // emitted bytes.
        bool endpoint_dbg()
        {
            static const bool dbg = std::getenv("WCT_SHOWER_ENDPOINT_DEBUG") != nullptr;
            return dbg;
        }

        // doc pr/91 round 3 -- WCT_SHOWER_WALK_DEBUG.  Fires only when
        // walk_visited_parity re-enqueues a vertex the legacy !has_node()
        // test would have pruned (i.e. it was already present in the view but
        // never walked).  Census only -- no attempt to classify how the node
        // first entered the view (former start vertex vs add_shower import);
        // that can be read off by comparing against the WCT_SHOWER_CONTENT_
        // DEBUG / WCT_SHOWER_ENDPOINT_DEBUG dumps for the same event.
        // Log/stderr only, no effect on emitted bytes.
        bool walk_dbg()
        {
            static const bool dbg = std::getenv("WCT_SHOWER_WALK_DEBUG") != nullptr;
            return dbg;
        }
        // doc sbnd_xin/docs/pr/93 round 3 -- WCT_SHOWER_ABSORB_DEBUG.
        // Instruments every membership-growth path (flood-fill ADD/EXCLUDE
        // in complete_structure_with_start_segment, plus the call-site tags
        // and direct add_segment/add_shower probes in
        // NeutrinoShowerClustering.cxx) so the absorbing site of any one
        // segment can be read off a single grep.  stderr only, no effect on
        // emitted bytes.
        bool absorb_dbg()
        {
            static const bool dbg = std::getenv("WCT_SHOWER_ABSORB_DEBUG") != nullptr;
            return dbg;
        }
        int vtx_display_id(const VertexPtr& v)
        {
            if (!v) return -1;
            const auto* c = v->cluster();
            return (c ? c->get_cluster_id() : 0) * 1000 + static_cast<int>(v->get_graph_index());
        }
        int seg_display_id(const SegmentPtr& s)
        {
            if (!s) return -1;
            int sid = s->id();
            if (sid < 0) sid = static_cast<int>(s->get_graph_index());
            const auto* c = s->cluster();
            return c ? c->get_cluster_id() * 1000 + sid : sid;
        }
        void probe_endpoint(Shower& sh, Graph& g, const WireCell::Point& sp,
                            const VertexPtr& start_vtx, int conn, bool excl,
                            bool skip_orphans, const char* tag)
        {
            if (!endpoint_dbg()) return;
            std::set<node_descriptor> touched;
            std::set<int> member_clusters;
            for (auto edesc : ordered_edges(sh, g)) {
                touched.insert(boost::source(edesc, g));
                touched.insert(boost::target(edesc, g));
                SegmentPtr seg = g[edesc].segment;
                if (seg && seg->cluster()) member_clusters.insert(seg->cluster()->get_cluster_id());
            }
            std::fprintf(stderr,
                         "SHOWER_ENDPOINT tag=%s shower_id=%d start_seg=%d conn=%d nseg=%d "
                         "excl_start_vtx=%d skip_orphans=%d start=(%.3f,%.3f,%.3f) "
                         "member_clusters=%zu\n",
                         tag, sh.get_shower_id(), seg_display_id(sh.start_segment()),
                         conn, sh.get_num_segments(), excl ? 1 : 0, skip_orphans ? 1 : 0,
                         sp.x() / WireCell::units::cm, sp.y() / WireCell::units::cm,
                         sp.z() / WireCell::units::cm, member_clusters.size());
            double best = -1;
            int best_id = -1, best_touch = -1;
            for (auto vdesc : ordered_nodes(sh, g)) {
                VertexPtr v = g[vdesc].vertex;
                if (!v) continue;
                const int touch0 = touched.count(vdesc) ? 1 : 0;
                const bool skipped = (excl && v == start_vtx) ||
                                     (skip_orphans && !touch0);
                const double dis = (sp - v->fit().point).magnitude();
                const int touch = touch0;
                std::fprintf(stderr,
                             "SHOWER_ENDPOINT   shower_id=%d cand_vtx=%d cluster=%d "
                             "(%.3f,%.3f,%.3f) dis=%.3f touched_by_member=%d skipped=%d\n",
                             sh.get_shower_id(), vtx_display_id(v),
                             v->cluster() ? v->cluster()->get_cluster_id() : -1,
                             v->fit().point.x() / WireCell::units::cm,
                             v->fit().point.y() / WireCell::units::cm,
                             v->fit().point.z() / WireCell::units::cm,
                             dis / WireCell::units::cm, touch, skipped ? 1 : 0);
                if (!skipped && dis > best) { best = dis; best_id = vtx_display_id(v); best_touch = touch; }
            }
            std::fprintf(stderr,
                         "SHOWER_ENDPOINT   shower_id=%d WINNER vtx=%d dis=%.3f "
                         "touched_by_member=%d\n",
                         sh.get_shower_id(), best_id, best / WireCell::units::cm, best_touch);
        }
    }

 
    // Default initialization constructor following WCPPID WCShower logic
    Shower::Shower(Graph& graph)
        : TrajectoryView(graph)
        , m_full_graph(graph)
    {
        // Assign a stable, unique shower ID
        m_shower_id = s_shower_id_counter.fetch_add(1, std::memory_order_relaxed);

        // Initialize all ShowerData members to defaults
        data.particle_type = 0;
        data.kenergy_range = 0;
        data.kenergy_dQdx = 0;
        data.kenergy_charge = 0;
        data.kenergy_best = 0;
        
        data.start_point = Point(0, 0, 0);
        data.end_point = Point(0, 0, 0);
        data.init_dir = Vector(0, 0, 0);
        
        data.start_connection_type = 0;
        
        // Initialize start vertex and segment to nullptr
        m_start_vertex = nullptr;
        m_start_segment = nullptr;
        
        // Set shower flag (following WCPPID)
        set_flags(ShowerFlags::kShower);
        unset_flags(ShowerFlags::kKinematics);
    }

    Shower::~Shower()
    {
    }


    VertexPtr Shower::start_vertex()
    {
        return m_start_vertex;
    }

    SegmentPtr Shower::start_segment()
    {
        return m_start_segment;
    }


    // Chainable setters

    /// Chainable setter of start vertex
    Shower& Shower::set_start_vertex(VertexPtr vtx, int type)
    {
        if (!vtx || ! vtx->descriptor_valid()) {
            m_start_vertex = nullptr;
            return *this;
        }
        this->add_vertex(vtx);
        m_start_vertex = vtx;
        data.start_connection_type = type;


        return *this;
    }

    
    
    /// Chainable setter of start segment
    Shower& Shower::set_start_segment(SegmentPtr seg, bool flag_include_vertices, const std::string& cloud_name_fit, const std::string& cloud_name_associate)
    {
        if (! seg->descriptor_valid()) {
            m_start_segment = nullptr;
            return *this;
        }
        // Membership gate: if the segment is already part of this shower's
        // view its points are already in the shower clouds, and the merge
        // below would duplicate them.  Instrumented on SBND nueCC evt 163543:
        // set_start_segment() on an already-member segment re-added 9 "fit" /
        // 85 "associate_points" points (the examine_showers
        // add_segment-then-set_start_segment sequence), and the
        // fresh-used_segments complete_structure pass re-added another
        // 11 / 140 through add_segment().  Graph membership was always
        // idempotent (std::set); the cloud merge now is too.
        const bool was_member = this->has_edge(seg->get_descriptor());
        // doc pr/84 round 3: a RETARGET (start segment replaced by a different
        // one after construction) is one of the two ways two showers end up
        // sharing a start segment.  Env-gated stderr only, byte-neutral.
        if (std::getenv("WCT_SHOWER_CREATE_DEBUG") && m_start_segment && m_start_segment != seg) {
            std::fprintf(stderr, "SHOWER_CREATE_DEBUG retarget shower_id=%d seg %d -> %d\n",
                         m_shower_id, m_start_segment->id(), seg->id());
        }
        TrajectoryView::add_segment(seg);
        m_start_segment = seg;
        invalidate_segment_caches();

        // If flag_include_vertices is true, add all vertices connected to this segment
        if (flag_include_vertices) {
            // Get the two vertices connected to this segment from the full graph
            auto vertices = find_vertices(m_full_graph, seg);

            // Add each vertex to the view (skip the start_vertex)
            if (vertices.first && vertices.first != m_start_vertex) {
                this->add_vertex(vertices.first);
                m_walked_nodes.insert(vertices.first->get_descriptor());  // doc pr/91 r3
            }
            if (vertices.second && vertices.second != m_start_vertex) {
                this->add_vertex(vertices.second);
                m_walked_nodes.insert(vertices.second->get_descriptor());  // doc pr/91 r3
            }
        }
        
        // Merge dynamic point clouds from segment to shower with the provided
        // names.  Seeding an as-yet-cloudless shower is always safe (there is
        // nothing to duplicate); the merge branch is gated on !was_member.
        if (!cloud_name_fit.empty()) {
            auto seg_dpc_fit = seg->dpcloud(cloud_name_fit);
            if (seg_dpc_fit) {
                auto shower_dpc_fit = this->dpcloud(cloud_name_fit);
                if (!shower_dpc_fit) {
                    // First segment: seed the shower DPC with a COPY of this
                    // segment's DPC -- never the segment's own object (clone_dpc).
                    this->dpcloud(cloud_name_fit, clone_dpc(*seg_dpc_fit));
                } else if (!was_member && shower_dpc_fit != seg_dpc_fit) {
                    // F14: merge wpid_params so the shower DPC can answer queries for
                    // all (apa,face) pairs present in any constituent segment's DPC.
                    // The `!=` test mirrors add_segment(): a merge onto an aliased
                    // cloud would be a self-append (doc pr/11 sec 6.4).
                    shower_dpc_fit->merge_wpid_params(*seg_dpc_fit);
                    shower_dpc_fit->add_points(*seg_dpc_fit);
                }
            }
        }

        if (!cloud_name_associate.empty()) {
            auto seg_dpc_associate = seg->dpcloud(cloud_name_associate);
            if (seg_dpc_associate) {
                auto shower_dpc_associate = this->dpcloud(cloud_name_associate);
                if (!shower_dpc_associate) {
                    this->dpcloud(cloud_name_associate, clone_dpc(*seg_dpc_associate));
                } else if (!was_member && shower_dpc_associate != seg_dpc_associate) {
                    shower_dpc_associate->merge_wpid_params(*seg_dpc_associate);
                    shower_dpc_associate->add_points(*seg_dpc_associate);
                }
            }
        }

        return *this;
    }

    void Shower::add_segment(SegmentPtr seg, bool flag_include_vertices, const std::string& cloud_name_fit, const std::string& cloud_name_associate)
    {
        if (! seg->descriptor_valid()) {
            return;
        }
        // Membership gate -- same rationale as set_start_segment() above:
        // callers routinely re-add segments that are already in the shower
        // (fresh used_segments sets handed to
        // complete_structure_with_start_segment, the examine_showers merge
        // loop), and the DPC merge below was unconditional while graph
        // membership is a no-op-on-repeat std::set insert.  Instrumented on
        // SBND nueCC evt 163543: a re-added segment put 11 "fit" / 140
        // "associate_points" duplicate points into the shower cloud.
        const bool was_member = this->has_edge(seg->get_descriptor());
        TrajectoryView::add_segment(seg);
        invalidate_segment_caches();

        // If flag_include_vertices is true, add all vertices connected to this segment
        if (flag_include_vertices) {
            // Get the two vertices connected to this segment from the full graph
            auto vertices = find_vertices(m_full_graph, seg);

            // Unlike set_start_segment() above, this branch has no
            // `!= m_start_vertex` guard on add_vertex() itself -- untouched,
            // byte-neutral.  The m_walked_nodes bookkeeping DOES guard on it
            // (doc pr/91 r3): if vertices.first/second is the CURRENT start
            // vertex, recording it here would wrongly mark it "walked" while
            // it is still exempt from the frontier test by the unconditional
            // `!= m_start_vertex` checks in complete_structure_with_start_
            // segment -- and that mark would incorrectly survive a later
            // re-seat, when this vertex becomes a FORMER start vertex that
            // the fix is specifically meant to re-expand.
            if (vertices.first) {
                this->add_vertex(vertices.first);
                if (vertices.first != m_start_vertex) m_walked_nodes.insert(vertices.first->get_descriptor());
            }
            if (vertices.second) {
                this->add_vertex(vertices.second);
                if (vertices.second != m_start_vertex) m_walked_nodes.insert(vertices.second->get_descriptor());
            }
        }

        // Merge dynamic point clouds from segment to shower with the provided names.
        // The shower always owns its cloud (clone_dpc on the seeding branch), so a
        // merge can never append a DPCBatch to itself -- which is what used to
        // double the cloud on every pass and blow the heap to a 224 GB peak RSS
        // (SBND MCP2025C evt 399118).  The `!= seg_dpc` tests remain as a cheap
        // invariant guard for any future path that reintroduces sharing.
        if (!cloud_name_fit.empty()) {
            auto seg_dpc_fit = seg->dpcloud(cloud_name_fit);
            if (seg_dpc_fit) {
                auto shower_dpc_fit = this->dpcloud(cloud_name_fit);
                if (!shower_dpc_fit) {
                    this->dpcloud(cloud_name_fit, clone_dpc(*seg_dpc_fit));
                } else if (!was_member && shower_dpc_fit != seg_dpc_fit) {
                    shower_dpc_fit->merge_wpid_params(*seg_dpc_fit);  // F14
                    shower_dpc_fit->add_points(*seg_dpc_fit);
                }
            }
        }

        if (!cloud_name_associate.empty()) {
            auto seg_dpc_associate = seg->dpcloud(cloud_name_associate);
            if (seg_dpc_associate) {
                auto shower_dpc_associate = this->dpcloud(cloud_name_associate);
                if (!shower_dpc_associate) {
                    this->dpcloud(cloud_name_associate, clone_dpc(*seg_dpc_associate));
                } else if (!was_member && shower_dpc_associate != seg_dpc_associate) {
                    shower_dpc_associate->merge_wpid_params(*seg_dpc_associate);  // F14
                    shower_dpc_associate->add_points(*seg_dpc_associate);
                }
            }
        }
    }

    void Shower::set_flag_kinematics(bool val){
        if (val) {
            set_flags(ShowerFlags::kKinematics);
        } else {
            unset_flags(ShowerFlags::kKinematics);
        }
    }
    
    bool Shower::get_flag_kinematics(){
        return flags_any(ShowerFlags::kKinematics);
    }
    
    bool Shower::get_flag_shower(){
        return flags_any(ShowerFlags::kShower);
    }

    void Shower::add_shower(Shower& shower, const std::string& cloud_name_fit, const std::string& cloud_name_associate){

        invalidate_segment_caches();

        // ordered_nodes/ordered_edges rather than the raw view sets: those are
        // unordered_sets keyed on heap addresses (see get_total_length()).  The
        // node loop is set-insertion and so order-blind, but the edge loop
        // below appends points into batch_fit/batch_associate, and THAT order
        // becomes the shower point cloud's row order -- which decides kNN ties
        // in shower_get_closest_point().  Ordered for the same reason both
        // loops are here (doc pr/28 sec 15).
        // doc pr/91 round 1: the node loop is unconditional -- every vertex of
        // the absorbed shower joins this one's view, including vertices no
        // surviving member segment touches.  get_end_point() is a farthest-
        // VERTEX search over exactly this set, so this is candidate mechanism
        // (a) for the cross-shower end points on SBND 169626 / 394532.
        const size_t pr91_nodes_before = this->nodes().size();
        for (auto vdesc : ordered_nodes(shower, m_full_graph)) {
            VertexPtr vtx = m_full_graph[vdesc].vertex;
            if (vtx && vtx->descriptor_valid()) this->add_vertex(vtx);
        }
        // doc pr/91 round 3: mirror the prototype's add_shower (WCShower.cxx:
        // 681-692), which merges the absorbed shower's own map_vtx_segs keys
        // -- so a vertex the absorbed shower never walked (e.g. its own
        // former start vertex) stays unwalked here too, and remains eligible
        // for re-expansion under walk_visited_parity.  Same graph
        // (m_full_graph is shared across all showers in one clustering pass),
        // so descriptors carry over directly.  Membership-tested only.
        m_walked_nodes.insert(shower.m_walked_nodes.begin(), shower.m_walked_nodes.end());
        if (endpoint_dbg()) {
            std::fprintf(stderr,
                         "SHOWER_ENDPOINT tag=add_shower into_sid=%d into_seg=%d from_sid=%d "
                         "from_seg=%d nodes %zu -> %zu (from had %zu)\n",
                         this->get_shower_id(), seg_display_id(this->start_segment()),
                         shower.get_shower_id(), seg_display_id(shower.start_segment()),
                         pr91_nodes_before, this->nodes().size(), shower.nodes().size());
        }

        // Batch-collect all points before adding to avoid repeated vector reallocations.
        // For each cloud: the first segment sets the shower's cloud pointer; all subsequent
        // segments' points are gathered into a single vector and appended in one call.
        Facade::DPCBatch batch_fit;
        Facade::DPCBatch batch_associate;

        for (auto edesc : ordered_edges(shower, m_full_graph)) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (!seg || !seg->descriptor_valid()) continue;
            // Membership gate -- same rationale as add_segment(): a segment
            // already in this shower's view has its points in the shower
            // clouds already, and batching them again below would duplicate
            // them (both showers holding the same segment is exactly the
            // absorb-loop situation this function serves).
            if (this->has_edge(edesc)) continue;
            // Graph-only add: DPCs are handled by the batch logic below to avoid
            // double-adding points (Shower::add_segment would also merge DPCs).
            TrajectoryView::add_segment(seg);

            if (!cloud_name_fit.empty()) {
                auto seg_dpc = seg->dpcloud(cloud_name_fit);
                if (seg_dpc) {
                    if (!this->dpcloud(cloud_name_fit)) {
                        // Seed with a COPY (clone_dpc), never the segment's own
                        // shared_ptr: dpcloud(name, ptr) stores the pointer
                        // verbatim (PRCommon.h), so seeding by pointer ALIASES
                        // the shower's cloud to the segment's -- the defect-F
                        // shape (doc pr/11 sec 6.7) this file's other seeding
                        // branches were already cured of.  With an aliased
                        // cloud a later pass through this loop would append
                        // the accumulated cloud to itself.
                        this->dpcloud(cloud_name_fit, clone_dpc(*seg_dpc));
                    } else if (this->dpcloud(cloud_name_fit) != seg_dpc) {
                        // F14: merge wpid_params before batching so the shower DPC
                        // serves queries for all (apa,face) pairs across all segments.
                        this->dpcloud(cloud_name_fit)->merge_wpid_params(*seg_dpc);
                        batch_fit.append(seg_dpc->points());
                    }
                }
            }

            if (!cloud_name_associate.empty()) {
                auto seg_dpc = seg->dpcloud(cloud_name_associate);
                if (seg_dpc) {
                    if (!this->dpcloud(cloud_name_associate)) {
                        this->dpcloud(cloud_name_associate, clone_dpc(*seg_dpc));
                    } else if (this->dpcloud(cloud_name_associate) != seg_dpc) {
                        this->dpcloud(cloud_name_associate)->merge_wpid_params(*seg_dpc);
                        batch_associate.append(seg_dpc->points());
                    }
                }
            }
        }

        // Single bulk add_points call per cloud instead of one per segment
        if (!batch_fit.empty()) {
            if (auto dpc = this->dpcloud(cloud_name_fit)) dpc->add_points(std::move(batch_fit));
        }
        if (!batch_associate.empty()) {
            if (auto dpc = this->dpcloud(cloud_name_associate)) dpc->add_points(std::move(batch_associate));
        }

    }

    // doc sbnd_xin/docs/pr/93 round 4 (shower_detach_track_stem) -- see the
    // header comment for the contract.  In-place mutation: m_shower_id and
    // the Shower object identity are preserved (pr/92's
    // dropped_satellite_shower_ids and IndexedShowerSet keying depend on it).
    int Shower::detach_track_prefix(const std::vector<SegmentPtr>& prefix,
                                    VertexPtr new_start_vertex,
                                    const std::string& cloud_name_fit,
                                    const std::string& cloud_name_associate)
    {
        if (prefix.empty()) return 0;
        if (!new_start_vertex || !new_start_vertex->descriptor_valid()) return 0;
        // Refuse to empty the shower: at least one member must remain.
        if (this->edges().size() <= prefix.size()) return 0;
        for (const auto& sg : prefix) {
            if (!sg || !sg->descriptor_valid() || !this->has_edge(sg->get_descriptor())) return 0;
        }

        // Collect the prefix chain's vertices (candidates for view removal).
        // Vector + membership-set of descriptors; never iterated by pointer.
        std::vector<VertexPtr> prefix_vtxs;
        std::unordered_set<size_t> seen_vtx_idx;
        for (const auto& sg : prefix) {
            auto [va, vb] = find_vertices(m_full_graph, sg);
            for (VertexPtr v : {va, vb}) {
                if (!v || !v->descriptor_valid() || v == new_start_vertex) continue;
                const size_t idx = m_full_graph[v->get_descriptor()].index;
                if (seen_vtx_idx.insert(idx).second) prefix_vtxs.push_back(v);
            }
        }

        // Remove the prefix edges from the view.
        int peeled = 0;
        for (const auto& sg : prefix) {
            if (TrajectoryView::remove_segment(sg)) ++peeled;
        }

        // A prefix vertex leaves the view only if NO remaining member still
        // touches it -- otherwise (e.g. a branch point shared with an EM
        // member) it must stay.  Membership set of graph indices, order-free.
        std::unordered_set<size_t> keep_vtx_idx;
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (!seg || !seg->descriptor_valid()) continue;
            auto [ka, kb] = find_vertices(m_full_graph, seg);
            if (ka && ka->descriptor_valid()) keep_vtx_idx.insert(m_full_graph[ka->get_descriptor()].index);
            if (kb && kb->descriptor_valid()) keep_vtx_idx.insert(m_full_graph[kb->get_descriptor()].index);
        }
        for (const auto& v : prefix_vtxs) {
            const size_t idx = m_full_graph[v->get_descriptor()].index;
            if (keep_vtx_idx.count(idx)) continue;
            TrajectoryView::remove_vertex(v);
            // Without this, the old root (e.g. the MAIN vertex) stays a view
            // node and wins calculate_kinematics' farthest-vertex end_point
            // search.  Also forget its walked mark: if some later pass
            // re-adds it, it is genuinely un-walked again.
            m_walked_nodes.erase(v->get_descriptor());
        }

        // Re-root: conn type 2 on purpose -- the pseudo-gamma rendering
        // class (append_pseudo_shower), the pi0 disconnected_showers class,
        // and kine_energy_included stays 1 (only type 3 is penalized).
        set_start_vertex(new_start_vertex, 2);

        // Re-seat the start segment: remaining member closest to the new
        // start vertex; graph-index tie-break.  ordered_edges = stable order.
        const WireCell::Point nv_pt = new_start_vertex->fit().valid()
            ? new_start_vertex->fit().point : new_start_vertex->wcpt().point;
        SegmentPtr new_start = nullptr;
        double best_dis = -1;
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (!seg || !seg->descriptor_valid()) continue;
            const double dis = segment_get_closest_point(seg, nv_pt).first;
            if (!new_start || dis < best_dis) {
                new_start = seg;
                best_dis = dis;
            }
        }
        if (new_start) {
            // Already a member => the cloud-merge branch in
            // set_start_segment is a no-op (was_member gate): no duplication.
            set_start_segment(new_start);
        }

        // Rebuild the shower point clouds from the REMAINING members only.
        // The clouds are add-only merges, and both kine_charge and the
        // conn-2 start_point derivation read them -- without the rebuild the
        // detached track's charge stays inside the daughter's energy.
        // Same code shape as add_shower's batch merge; ordered_edges order.
        for (const std::string& cname : {cloud_name_fit, cloud_name_associate}) {
            if (cname.empty()) continue;
            this->dpcloud(cname, nullptr);
            Facade::DPCBatch batch;
            for (auto edesc : ordered_edges(*this, m_full_graph)) {
                SegmentPtr seg = m_full_graph[edesc].segment;
                if (!seg || !seg->descriptor_valid()) continue;
                auto seg_dpc = seg->dpcloud(cname);
                if (!seg_dpc) continue;
                if (!this->dpcloud(cname)) {
                    this->dpcloud(cname, clone_dpc(*seg_dpc));
                } else if (this->dpcloud(cname) != seg_dpc) {
                    this->dpcloud(cname)->merge_wpid_params(*seg_dpc);
                    batch.append(seg_dpc->points());
                }
            }
            if (!batch.empty()) {
                if (auto dpc = this->dpcloud(cname)) dpc->add_points(std::move(batch));
            }
        }

        invalidate_segment_caches();
        set_flag_kinematics(false);
        return peeled;
    }

    // doc sbnd_xin/docs/pr/99 round 2 (shower_ghost_member_drop).  Contract
    // and rationale in PRShower.h; bookkeeping forked BY DUPLICATION from
    // detach_track_prefix above -- that production method stays untouched.
    std::shared_ptr<Facade::DynamicPointCloud> Shower::rebuild_pcloud(const std::string& cloud_name)
    {
        // doc pr/99 round 3 (C1b).  Contract in PRShower.h.  Same code shape
        // as the detach_track_prefix rebuild loop (ordered_edges order), but
        // writing into a fresh cloud instead of this->dpcloud.
        std::shared_ptr<Facade::DynamicPointCloud> out = nullptr;
        Facade::DPCBatch batch;
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (!seg || !seg->descriptor_valid()) continue;
            auto seg_dpc = seg->dpcloud(cloud_name);
            if (!seg_dpc) continue;
            if (!out) {
                out = clone_dpc(*seg_dpc);
            } else {
                out->merge_wpid_params(*seg_dpc);
                batch.append(seg_dpc->points());
            }
        }
        if (out && !batch.empty()) out->add_points(std::move(batch));
        return out;
    }

    int Shower::drop_ghost_member(SegmentPtr ghost,
                                  const std::string& cloud_name_fit,
                                  const std::string& cloud_name_associate)
    {
        if (!ghost || !ghost->descriptor_valid() || !this->has_edge(ghost->get_descriptor())) return 0;
        if (ghost == m_start_segment) return 0;
        // Refuse to empty the shower: at least one member must remain.
        if (this->edges().size() <= 1) return 0;
        if (!m_start_segment || !m_start_segment->descriptor_valid()
            || !this->has_edge(m_start_segment->get_descriptor())) return 0;

        // Leaf-only guard: with the ghost edge filtered out, every remaining
        // member must stay reachable from the start segment.  Test-remove is
        // safe -- TrajectoryView::remove/add_segment are pure filter-set
        // edits (no cloud merge, unlike Shower::add_segment).
        TrajectoryView::remove_segment(ghost);
        const int n_reach = count_connected_segments(m_start_segment);
        if (n_reach != static_cast<int>(this->edges().size())) {
            TrajectoryView::add_segment(ghost);  // restore: removal would strand a member
            return 0;
        }

        // Ghost-only vertices leave the view (same rule as
        // detach_track_prefix: a vertex still touched by a remaining member
        // must stay; a departed one must not win the farthest-vertex
        // end_point search, and its walked mark must be forgotten).
        std::unordered_set<size_t> keep_vtx_idx;
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (!seg || !seg->descriptor_valid()) continue;
            auto [ka, kb] = find_vertices(m_full_graph, seg);
            if (ka && ka->descriptor_valid()) keep_vtx_idx.insert(m_full_graph[ka->get_descriptor()].index);
            if (kb && kb->descriptor_valid()) keep_vtx_idx.insert(m_full_graph[kb->get_descriptor()].index);
        }
        auto [gva, gvb] = find_vertices(m_full_graph, ghost);
        for (VertexPtr v : {gva, gvb}) {
            if (!v || !v->descriptor_valid()) continue;
            const size_t idx = m_full_graph[v->get_descriptor()].index;
            if (keep_vtx_idx.count(idx)) continue;
            TrajectoryView::remove_vertex(v);
            m_walked_nodes.erase(v->get_descriptor());
        }

        // Rebuild the shower point clouds from the REMAINING members only
        // (detach_track_prefix rationale: the clouds are add-only merges and
        // kine_charge reads them -- the ghost's charge must leave).
        for (const std::string& cname : {cloud_name_fit, cloud_name_associate}) {
            if (cname.empty()) continue;
            this->dpcloud(cname, nullptr);
            Facade::DPCBatch batch;
            for (auto edesc : ordered_edges(*this, m_full_graph)) {
                SegmentPtr seg = m_full_graph[edesc].segment;
                if (!seg || !seg->descriptor_valid()) continue;
                auto seg_dpc = seg->dpcloud(cname);
                if (!seg_dpc) continue;
                if (!this->dpcloud(cname)) {
                    this->dpcloud(cname, clone_dpc(*seg_dpc));
                } else if (this->dpcloud(cname) != seg_dpc) {
                    this->dpcloud(cname)->merge_wpid_params(*seg_dpc);
                    batch.append(seg_dpc->points());
                }
            }
            if (!batch.empty()) {
                if (auto dpc = this->dpcloud(cname)) dpc->add_points(std::move(batch));
            }
        }

        invalidate_segment_caches();
        set_flag_kinematics(false);
        return 1;
    }

    void Shower::complete_structure_with_start_segment(IndexedSegmentSet& used_segments, const std::string& cloud_name_fit, const std::string& cloud_name_associate, bool absorb_track_guard, bool walk_visited_parity) {
        if (!m_start_segment || !m_start_segment->descriptor_valid()) return;

        // doc sbnd_xin/docs/pr/40 round 6 F12: the flood-fill below has no
        // per-segment test at all, so one shower-flagged seed swallows every
        // connected segment regardless of its own PID (round-5 G2a/G2c
        // failure mechanism).  When absorb_track_guard is on, a confidently
        // PID'd non-electron that is long and straight is not absorbed and
        // the walk terminates there (its far vertex is never enqueued, so a
        // Michel candidate beyond a stopping muon stays unclaimed for the
        // later seeding passes).  Long-muon pseudo-showers are exempt: the
        // in_main_cluster seeder records particle_type 13 on the shower
        // BEFORE calling here (NeutrinoShowerClustering.cxx:132), and broken
        // muon tracks must keep reassembling through this flood-fill.  The
        // excluded segment is deliberately NOT inserted into used_segments
        // (the predicate is a pure segment property -- every other
        // flood-fill re-excludes it) and its shower flags are deliberately
        // NOT consulted (a stale flag must not defeat the exclusion).
        // false = legacy = byte-identical.
        const bool apply_guard = absorb_track_guard && this->get_particle_type() != 13;
        auto guard_excludes = [&](const SegmentPtr& seg) -> bool {
            if (!apply_guard) return false;
            if (!seg->has_particle_info() || !seg->particle_info()) return false;
            const int pdg = seg->particle_info()->pdg();
            if (pdg == 0 || std::abs(pdg) == 11) return false;
            return segment_is_straight_long_track(seg);
        };
        if (absorb_dbg()) {
            std::fprintf(stderr, "SHOWER_ABSORB walk_begin shower_start_seg=%d start_vtx=%d cached_type=%d guard_arg=%d apply_guard=%d\n",
                         seg_display_id(m_start_segment), vtx_display_id(m_start_vertex),
                         this->get_particle_type(), (int)absorb_track_guard, (int)apply_guard);
        }
        
        std::vector<SegmentPtr> new_segments;
        std::vector<VertexPtr> new_vertices;
        
        // Add start_segment to the view and mark as used
        used_segments.insert(m_start_segment);
        
      
        // Find vertices connected to start_segment (excluding start_vertex)
        auto vertices = find_vertices(m_full_graph, m_start_segment);
        if (vertices.first && vertices.first != m_start_vertex) {
            this->add_vertex(vertices.first);
            new_vertices.push_back(vertices.first);
            m_walked_nodes.insert(vertices.first->get_descriptor());  // doc pr/91 r3: about to be walked below
        }
        if (vertices.second && vertices.second != m_start_vertex) {
            this->add_vertex(vertices.second);
            new_vertices.push_back(vertices.second);
            m_walked_nodes.insert(vertices.second->get_descriptor());  // doc pr/91 r3
        }
        
        // Worklist algorithm: explore connected segments and vertices
        while (!new_vertices.empty() || !new_segments.empty()) {
            // Process new vertices - find all segments connected to them
            if (!new_vertices.empty()) {
                VertexPtr vtx = new_vertices.back();
                new_vertices.pop_back();
                
                // Find all segments connected to this vertex
                if (vtx->descriptor_valid()) {
                    for (auto edesc : sorted_out_edges(vtx->get_descriptor(), m_full_graph)) {
                        SegmentPtr seg = m_full_graph[edesc].segment;
                        if (seg && seg->descriptor_valid() && used_segments.find(seg) == used_segments.end()) {
                            // F12 (doc pr/40 round 6): see guard_excludes above.
                            if (guard_excludes(seg)) {
                                if (absorb_dbg()) {
                                    std::fprintf(stderr, "SHOWER_ABSORB EXCLUDE shower_start_seg=%d seg=%d pdg=%d\n",
                                                 seg_display_id(m_start_segment), seg_display_id(seg),
                                                 seg->has_particle_info() && seg->particle_info() ? seg->particle_info()->pdg() : 0);
                                }
                                continue;
                            }
                            if (absorb_dbg()) {
                                std::fprintf(stderr, "SHOWER_ABSORB ADD shower_start_seg=%d seg=%d pdg=%d len_cm=%.2f straight=%d\n",
                                             seg_display_id(m_start_segment), seg_display_id(seg),
                                             seg->has_particle_info() && seg->particle_info() ? seg->particle_info()->pdg() : 0,
                                             segment_track_length(seg)/units::cm,
                                             (int)segment_is_straight_long_track(seg));
                            }
                            // add_segment() already performs the fit/associate merge
                            // (and merge_wpid_params, which the block that used to sit
                            // here omitted).  It previously ran with the DEFAULT cloud
                            // names while a second, hand-rolled copy of the same merge
                            // followed for cloud_name_fit/cloud_name_associate -- so
                            // every segment's points were added to the shower TWICE
                            // (instrumented on SBND nueCC evt 172230: 22/22 second
                            // merges duplicated a first one).  Thread the names through
                            // and merge exactly once.
                            this->add_segment(seg, false, cloud_name_fit, cloud_name_associate);
                            new_segments.push_back(seg);
                            used_segments.insert(seg);
                        }
                    }
                }
            }
            
            // Process new segments - find all vertices connected to them
            if (!new_segments.empty()) {
                SegmentPtr seg = new_segments.back();
                new_segments.pop_back();
                
                // Find vertices connected to this segment (excluding start_vertex)
                auto vertices = find_vertices(m_full_graph, seg);
                // doc pr/91 round 3: legacy frontier test is pure view
                // MEMBERSHIP (!has_node) -- a vertex added by
                // set_start_vertex()/set_start_segment()/add_segment()/
                // add_shower() but never actually scanned by this worklist is
                // permanently skipped.  walk_visited_parity switches the test
                // to m_walked_nodes, the prototype's map_vtx_segs equivalent
                // (WCShower.cxx:735-742): a present-but-unwalked vertex is
                // re-enqueued.  See PRShower.h for the full mechanism.
                if (vertices.first && vertices.first != m_start_vertex) {
                    const auto vd = vertices.first->get_descriptor();
                    const bool unseen = walk_visited_parity ? !m_walked_nodes.count(vd) : !this->has_node(vd);
                    if (unseen) {
                        if (walk_dbg() && this->has_node(vd)) {
                            std::fprintf(stderr,
                                         "SHOWER_WALK rewalk shower_id=%d via_seg=%d vtx=%d\n",
                                         m_shower_id, seg_display_id(seg), vtx_display_id(vertices.first));
                        }
                        this->add_vertex(vertices.first);
                        new_vertices.push_back(vertices.first);
                        m_walked_nodes.insert(vd);  // about to be walked
                    }
                }
                if (vertices.second && vertices.second != m_start_vertex) {
                    const auto vd = vertices.second->get_descriptor();
                    const bool unseen = walk_visited_parity ? !m_walked_nodes.count(vd) : !this->has_node(vd);
                    if (unseen) {
                        if (walk_dbg() && this->has_node(vd)) {
                            std::fprintf(stderr,
                                         "SHOWER_WALK rewalk shower_id=%d via_seg=%d vtx=%d\n",
                                         m_shower_id, seg_display_id(seg), vtx_display_id(vertices.second));
                        }
                        this->add_vertex(vertices.second);
                        new_vertices.push_back(vertices.second);
                        m_walked_nodes.insert(vd);  // about to be walked
                    }
                }
            }
        }
    }

    void Shower::fill_sets(IndexedVertexSet& used_vertices, IndexedSegmentSet& used_segments, bool flag_exclude_start_segment,
                           bool exclude_start_vertex){
        // Fill used_vertices with all vertices in this shower's view (index-stable order).
        // exclude_start_vertex: prototype map_vtx_segs parity -- the start
        // vertex (a main-track attachment point, present in the toolkit view
        // via set_start_vertex) is omitted so consumers' BFS barriers block
        // only shower-INTERIOR vertices (see header comment, doc pr/38).
        for (auto vdesc : ordered_nodes(*this, m_full_graph)) {
            VertexPtr vtx = m_full_graph[vdesc].vertex;
            if (vtx) {
                if (exclude_start_vertex && vtx == m_start_vertex) continue;
                used_vertices.insert(vtx);
            }
        }

        // Fill used_segments with all segments in this shower's view (index-stable order)
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (seg) {
                // Skip start_segment if flag is set
                if (flag_exclude_start_segment && seg == m_start_segment) {
                    continue;
                }
                used_segments.insert(seg);
            }
        }
    }

    void Shower::fill_point_vector(std::vector<WireCell::Point>& points, bool flag_main){
        // Get the main cluster ID if flag_main is true
        const Facade::Cluster* main_cluster = nullptr;
        if (flag_main && m_start_segment && m_start_segment->cluster()) {
            main_cluster = m_start_segment->cluster();
        }

        // Pre-count to avoid repeated reallocations
        size_t n_pts = 0;
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (!seg) continue;
            if (flag_main && main_cluster && seg->cluster() != main_cluster) continue;
            const auto& sz = seg->fits().size();
            if (sz >= 2) n_pts += sz - 2;
        }
        for (auto vdesc : ordered_nodes(*this, m_full_graph)) {
            VertexPtr vtx = m_full_graph[vdesc].vertex;
            if (!vtx) continue;
            if (flag_main && main_cluster && vtx->cluster() != main_cluster) continue;
            ++n_pts;
        }
        points.reserve(points.size() + n_pts);

        // Fill points from all segments in the shower's view (index-stable order)
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = m_full_graph[edesc].segment;
            if (seg) {
                // Skip if flag_main is set and segment is not in the main cluster
                if (flag_main && main_cluster && seg->cluster() != main_cluster) {
                    continue;
                }

                // Get segment fit points and add all except first and last
                const auto& fits = seg->fits();
                for (size_t i = 1; i + 1 < fits.size(); i++) {
                    points.push_back(fits[i].point);
                }
            }
        }

        // Fill points from all vertices in the shower's view (index-stable order)
        for (auto vdesc : ordered_nodes(*this, m_full_graph)) {
            VertexPtr vtx = m_full_graph[vdesc].vertex;
            if (vtx) {
                // Skip if flag_main is set and vertex is not in the main cluster
                if (flag_main && main_cluster && vtx->cluster() != main_cluster) {
                    continue;
                }

                // Add the vertex fit point
                points.push_back(vtx->fit().point);
            }
        }
    }

    TrajectoryView& Shower::fill_maps() {
        return *this;
    }

    int Shower::count_connected_segments(SegmentPtr seg) {
        if (!seg || !seg->descriptor_valid() || !this->has_edge(seg->get_descriptor())) return 0;

        std::unordered_set<SegmentPtr> used_segments;
        std::unordered_set<VertexPtr> used_vertices;
        std::vector<SegmentPtr> new_segments;
        std::vector<VertexPtr> new_vertices;

        new_segments.push_back(seg);
        used_segments.insert(seg);

        while (!new_vertices.empty() || !new_segments.empty()) {
            while (!new_vertices.empty()) {
                VertexPtr vtx = new_vertices.back(); new_vertices.pop_back();
                if (!vtx->descriptor_valid()) continue;
                for (auto edesc : sorted_out_edges(vtx->get_descriptor(), m_full_graph)) {
                    if (!this->has_edge(edesc)) continue;
                    SegmentPtr s = m_full_graph[edesc].segment;
                    if (s && used_segments.insert(s).second) new_segments.push_back(s);
                }
            }
            while (!new_segments.empty()) {
                SegmentPtr s = new_segments.back(); new_segments.pop_back();
                auto [va, vb] = find_vertices(m_full_graph, s);
                if (va && this->has_node(va->get_descriptor()) && used_vertices.insert(va).second) new_vertices.push_back(va);
                if (vb && this->has_node(vb->get_descriptor()) && used_vertices.insert(vb).second) new_vertices.push_back(vb);
            }
        }
        return (int)used_segments.size();
    }

    std::pair<IndexedVertexSet, IndexedSegmentSet> Shower::get_connected_pieces(SegmentPtr seg){
        IndexedVertexSet  result_vertices;
        IndexedSegmentSet result_segments;

        if (!seg || !seg->descriptor_valid() || !this->has_edge(seg->get_descriptor())) {
            return std::make_pair(result_vertices, result_segments);
        }

        // Use unordered_set internally for O(1) membership checks; copy to ordered sets at return.
        std::unordered_set<SegmentPtr> used_segments;
        std::unordered_set<VertexPtr> used_vertices;

        std::vector<SegmentPtr> new_segments;
        std::vector<VertexPtr> new_vertices;

        new_segments.push_back(seg);
        used_segments.insert(seg);

        while (!new_vertices.empty() || !new_segments.empty()) {
            // Drain the vertex queue fully before switching to segments
            while (!new_vertices.empty()) {
                VertexPtr vtx = new_vertices.back();
                new_vertices.pop_back();

                if (!vtx->descriptor_valid()) continue;
                for (auto edesc : sorted_out_edges(vtx->get_descriptor(), m_full_graph)) {
                    if (!this->has_edge(edesc)) continue;
                    SegmentPtr seg1 = m_full_graph[edesc].segment;
                    if (seg1 && used_segments.insert(seg1).second) {
                        new_segments.push_back(seg1);
                    }
                }
            }

            // Drain the segment queue fully before switching to vertices
            while (!new_segments.empty()) {
                SegmentPtr seg1 = new_segments.back();
                new_segments.pop_back();

                auto [va, vb] = find_vertices(m_full_graph, seg1);
                if (va && this->has_node(va->get_descriptor()) && used_vertices.insert(va).second) {
                    new_vertices.push_back(va);
                }
                if (vb && this->has_node(vb->get_descriptor()) && used_vertices.insert(vb).second) {
                    new_vertices.push_back(vb);
                }
            }
        }

        result_vertices.insert(used_vertices.begin(), used_vertices.end());
        result_segments.insert(used_segments.begin(), used_segments.end());
        return std::make_pair(result_vertices, result_segments);
    }

    std::pair<SegmentPtr, VertexPtr> Shower::get_last_segment_vertex_long_muon(IndexedSegmentSet& segments_in_muons) {
        VertexPtr s_vtx = m_start_vertex;
        SegmentPtr s_seg = m_start_segment;
        
        if (!s_vtx || !s_seg) {
            return std::make_pair(s_seg, s_vtx);
        }
        
        std::unordered_set<SegmentPtr> used_segments;
        used_segments.insert(s_seg);
        
        bool flag_continue = true;
        while (flag_continue) {
            flag_continue = false;
            // Progress guard: the s_vtx == m_start_vertex branch requests
            // another pass without touching any state, and the only escape is
            // find_vertices() moving s_vtx below.  find_vertices() returns
            // {nullptr, nullptr} for an invalid segment descriptor
            // (PRGraph.cxx), and a self-loop segment leaves s_vtx unchanged
            // too -- either way the pass is an exact no-op and the loop never
            // terminates (the same silent-no-op shape as the fixed class-A
            // hang, doc pr/11 sec 6.1).  Detect "nothing moved" at the end of
            // the pass and stop.
            const SegmentPtr pass_start_seg = s_seg;
            const VertexPtr pass_start_vtx = s_vtx;

            // If current vertex is start_vertex, continue
            if (s_vtx == m_start_vertex) {
                flag_continue = true;
            } else {
                // Look for a new segment connected to s_vtx that is in segments_in_muons and not used
                if (s_vtx->descriptor_valid()) {
                    auto vdesc = s_vtx->get_descriptor();
                    // sorted_out_edges gives index-stable selection when multiple muon segments are available
                    for (auto edesc : sorted_out_edges(vdesc, m_full_graph)) {
                        if (this->has_edge(edesc)) {
                            SegmentPtr sg = m_full_graph[edesc].segment;
                            if (sg && segments_in_muons.find(sg) != segments_in_muons.end() 
                                && used_segments.find(sg) == used_segments.end()) {
                                s_seg = sg;
                                used_segments.insert(s_seg);
                                flag_continue = true;
                                break;
                            }
                        }
                    }
                }
            }
            
            // If we found a new segment, find the other vertex connected to it
            if (flag_continue) {
                auto vertices = find_vertices(m_full_graph, s_seg);
                if (vertices.first && vertices.first != s_vtx) {
                    s_vtx = vertices.first;
                } else if (vertices.second && vertices.second != s_vtx) {
                    s_vtx = vertices.second;
                }
                if (s_seg == pass_start_seg && s_vtx == pass_start_vtx) {
                    SPDLOG_LOGGER_WARN(s_log, "get_last_segment_vertex_long_muon: a pass requested continuation but moved neither segment nor vertex (shower {}); stopping to avoid an unbounded loop", m_shower_id);
                    break;
                }
            }
        }
        
        return std::make_pair(s_seg, s_vtx);
    }

    int Shower::get_num_main_segments() {
        if (m_num_main_segs_cache >= 0) return m_num_main_segs_cache;

        int num = 0;
        if (!m_start_segment) {
            m_num_main_segs_cache = 0;
            return 0;
        }
        auto start_cluster = m_start_segment->cluster();
        if (!start_cluster) {
            m_num_main_segs_cache = 0;
            return 0;
        }
        const auto& view = this->view_graph();
        // Left on the raw (address-hashed) edge set deliberately: this body is
        // an unconditional integer count, so no walk order is observable and
        // ordered_edges() would only buy a vector allocation.  Every loop in
        // this file that ACCUMULATES or picks does use ordered_edges.
        for (auto edesc : this->edges()) {
            SegmentPtr seg = view[edesc].segment;
            if (seg && seg->cluster() == start_cluster) ++num;
        }
        m_num_main_segs_cache = num;
        return num;
    }

    int Shower::get_num_segments() {
        return this->edges().size();
    }

    double Shower::get_total_length(){
        if (m_total_length_cache >= 0.0) return m_total_length_cache;

        double total_length = 0;
        const auto& view = this->view_graph();
        // ordered_edges, not this->edges(): edges() is an unordered_set whose
        // hash is over two node_descriptors, and with boost::setS vertices a
        // node_descriptor is a HEAP ADDRESS -- so the bucket walk differs
        // between runs of the identical binary.  This += is an FP accumulation,
        // so that walk order is observable at the ULP: one leaf of SBND evt 388
        // (doc pr/28 sec 15).
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = view[edesc].segment;
            if (seg) total_length += segment_track_length(seg);
        }
        m_total_length_cache = total_length;
        return total_length;
    }
    double Shower::get_total_length(Facade::Cluster* cluster){
        double total_length = 0;
        
        if (!cluster) {
            return 0;
        }
        
        // Get the view graph to access segments
        const auto& view = this->view_graph();
        
        // Iterate through all segments in the shower (ordered_edges: FP
        // accumulation, see get_total_length() above)
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = view[edesc].segment;
            if (!seg) continue;

            // Check if segment's cluster matches the input cluster
            if (seg->cluster() == cluster) {
                total_length += segment_track_length(seg);
            }
        }
        
        return total_length;
    }
    double Shower::get_total_track_length(){
        double total_length = 0;
        
        // Get the view graph to access segments
        const auto& view = this->view_graph();
        
        // Iterate through all segments in the shower (ordered_edges: FP
        // accumulation, see get_total_length() above)
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = view[edesc].segment;
            if (!seg) continue;

            // Mirrors prototype: only count if !get_flag_shower()
            // = !(trajectory || topology || |pdg|==11)
            bool is_shower_seg = seg->flags_any(SegmentFlags::kShowerTrajectory) ||
                                 seg->flags_any(SegmentFlags::kShowerTopology) ||
                                 (seg->has_particle_info() && std::abs(seg->particle_info()->pdg()) == 11);
            if (!is_shower_seg) {
                total_length += segment_track_length(seg);
            }
        }
        
        return total_length;
    }

    void Shower::update_particle_type(const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, double mip_dqdx, VertexPtr main_vertex, bool protect_proton_daughter_pion, double proton_daughter_mip_dqdx, bool vote_track_pid_counts, bool accept_pid_guard, double accept_pid_min_len){
        double track_length = 0;
        double shower_length = 0;
        
        // Only process if there's more than one segment
        if (this->edges().size() <= 1) {
            return;
        }
        
        // Get the view graph to access segments
        const auto& view = this->view_graph();
        
        // Iterate through all segments in the shower.  ordered_edges, NOT
        // this->edges(): the two accumulators below feed the
        // `shower_length > track_length` branch, so an ULP move here can flip
        // the branch and retype the start segment to an electron.  This is the
        // one converted site whose effect is not confined to output rounding.
        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            SegmentPtr seg = view[edesc].segment;
            if (!seg) continue;

            double length = segment_track_length(seg);
            
            // Check if segment is a shower segment OR not a proton (PDG 2212)
            bool is_shower = seg->flags_any(SegmentFlags::kShowerTrajectory) ||
                           seg->flags_any(SegmentFlags::kShowerTopology);

            bool is_not_proton = true;
            if (seg->has_particle_info()) {
                int pdg = seg->particle_info()->pdg();
                is_not_proton = (std::abs(pdg) != 2212);
            }

            // doc sbnd_xin/docs/pr/93 Cause C (SBND 18255-292643): legacy
            // vote counts ONLY confirmed protons as track, so an unflagged
            // muon/pion chain always votes electron.  When enabled, any
            // unflagged member carrying a track pdg {13,211,2212} counts as
            // track.  No score conjunct: score 100 is the median-fallback
            // population itself, not a confidence signal.
            // C++ default false => byte-identical legacy vote.
            bool counts_as_track = !is_not_proton;
            if (vote_track_pid_counts && seg->has_particle_info() && seg->particle_info()) {
                const int apdg = std::abs(seg->particle_info()->pdg());
                counts_as_track = (apdg == 13 || apdg == 211 || apdg == 2212);
            }

            if (is_shower || !counts_as_track) {
                shower_length += length;
            } else {
                track_length += length;
            }
        }
        
        // If shower_length dominates, update start_segment to electron
        if (shower_length > track_length && m_start_segment) {
            // doc sbnd_xin/docs/pr/40 round 3: an electron cannot father a
            // proton -- if set_default_shower_particle_info already relabelled
            // m_start_segment pion via this same topology test, don't
            // silently revert it back to electron here.  is_not_proton above
            // only exempts a confirmed PROTON from the shower-length bucket,
            // so a pion-labelled start_segment still landed in shower_length
            // and still trips this branch; C++ default false/nullptr = legacy
            // = byte-identical.
            const bool protected_pion = protect_proton_daughter_pion && main_vertex &&
                segment_has_proton_daughter(m_full_graph, m_start_segment, main_vertex, proton_daughter_mip_dqdx);
            // doc sbnd_xin/docs/pr/93 Cause B: shower_accept_pid_guard also
            // guards THIS overwrite -- otherwise the vote re-flips a
            // confidently-PID'd non-electron start segment one call after
            // the acceptance-site guard spared it (348471's 0.23-score
            // proton).  Same predicate as the acceptance sites, incl. the
            // min-length floor.  C++ default false = legacy = byte-identical.
            const bool protected_confident_pid = accept_pid_guard &&
                segment_track_length(m_start_segment) > accept_pid_min_len &&
                segment_confident_nonelectron_pid(m_start_segment);
            // if-guarded rather than an early return: keeps any code appended
            // to this function later from being silently skipped for a
            // protected shower.
            if (!protected_pion && !protected_confident_pid) {
                // Calculate 4-momentum for electron (PDG = 11)
                auto four_momentum = segment_cal_4mom(m_start_segment, 11, particle_data, recomb_model, mip_dqdx);

                // Create ParticleInfo for electron
                auto pinfo = std::make_shared<Aux::ParticleInfo>(
                    11,                                          // electron PDG
                    particle_data->get_particle_mass(11),       // electron mass
                    particle_data->pdg_to_name(11),             // "electron"
                    four_momentum                                // 4-momentum
                );

                // Store particle info in start_segment
                m_start_segment->particle_info(pinfo);
            }
        }
    }

    std::vector<double> Shower::get_stem_dQ_dx(VertexPtr vertex, SegmentPtr segment, int limit /*=20*/, double mip_dqdx_median /*=43000/units::cm*/){
        std::vector<double> vec_dQ_dx;
        const double MIP_dQdx = mip_dqdx_median;
        
        if (!vertex || !segment) {
            return vec_dQ_dx;
        }
        
        // Get dQ and dx from segment's fits
        const auto& fits = segment->fits();
        if (fits.empty()) {
            return vec_dQ_dx;
        }
        
        // Determine direction based on vertex position relative to segment
        // Use squared distances to avoid sqrt
        bool vertex_at_front = false;
        if (!segment->wcpts().empty()) {
            const auto& vp = vertex->wcpt().point;
            const auto& fp = segment->wcpts().front().point;
            const auto& bp = segment->wcpts().back().point;
            double d1sq = (vp - fp).magnitude2();
            double d2sq = (vp - bp).magnitude2();
            vertex_at_front = (d1sq < d2sq);
        }
        
        // Fill vec_dQ_dx based on direction
        if (vertex_at_front) {
            for (size_t i = 0; i < fits.size(); i++) {
                double dQ_dx_normalized = fits[i].dQ / (fits[i].dx + 1e-9) / MIP_dQdx;
                vec_dQ_dx.push_back(dQ_dx_normalized);
                if (vec_dQ_dx.size() >= (size_t)limit) break;
            }
        } else {
            for (int i = (int)fits.size() - 1; i >= 0; i--) {
                double dQ_dx_normalized = fits[i].dQ / (fits[i].dx + 1e-9) / MIP_dQdx;
                vec_dQ_dx.push_back(dQ_dx_normalized);
                if (vec_dQ_dx.size() >= (size_t)limit) break;
            }
        }
        
        // If this is the start_segment and we don't have enough points, continue to next segments
        if (segment == m_start_segment && vec_dQ_dx.size() < (size_t)limit) {
            VertexPtr curr_vertex = vertex;
            SegmentPtr curr_segment = segment;
            int count = 0;
            
            while (vec_dQ_dx.size() < (size_t)limit && count < 3) {
                // Find next vertex (the other end of current segment)
                VertexPtr next_vertex = find_other_vertex(m_full_graph, curr_segment, curr_vertex);
                if (!next_vertex) break;
                
                // Direction from current vertex to next vertex
                WireCell::Vector dir1 = curr_vertex->fit().point - next_vertex->fit().point;
                const double dir1_mag = dir1.magnitude();

                // Single pass over out-edges: find best next segment AND check flag_bad,
                // caching computed directions to avoid recomputing them in a second pass.
                SegmentPtr next_segment = nullptr;
                WireCell::Vector dir2;
                double max_angle = 0;

                // Collect candidate (other) segments with their cached directions and lengths
                struct SegCandidate { SegmentPtr seg; WireCell::Vector dir; double length; };
                std::vector<SegCandidate> other_candidates;

                auto next_vdesc = next_vertex->get_descriptor();
                for (auto edesc : sorted_out_edges(next_vdesc, m_full_graph)) {
                    if (!has_edge(edesc)) continue;

                    SegmentPtr seg = m_full_graph[edesc].segment;
                    if (seg == curr_segment) continue;

                    WireCell::Vector tmp_dir = segment_cal_dir_3vector(seg, next_vertex->fit().point, 10 * units::cm);
                    double denom = dir1_mag * tmp_dir.magnitude() + 1e-9;
                    double angle = std::acos(std::clamp(dir1.dot(tmp_dir) / denom, -1.0, 1.0)) * 180.0 / M_PI;

                    if (angle > max_angle) {
                        // Demote previous best to other_candidates before replacing
                        if (next_segment) {
                            other_candidates.push_back({next_segment, dir2, segment_track_length(next_segment)});
                        }
                        max_angle = angle;
                        next_segment = seg;
                        dir2 = tmp_dir;
                    } else {
                        other_candidates.push_back({seg, tmp_dir, segment_track_length(seg)});
                    }
                }

                if (!next_segment) break;

                // Check flag_bad using cached directions — no second edge traversal needed
                bool flag_bad = false;
                const double dir2_mag = dir2.magnitude();
                for (const auto& cand : other_candidates) {
                    if (cand.length > 3 * units::cm) {
                        double denom = dir2_mag * cand.dir.magnitude() + 1e-9;
                        double angle = std::acos(std::clamp(dir2.dot(cand.dir) / denom, -1.0, 1.0)) * 180.0 / M_PI;
                        if (angle < 25) {
                            flag_bad = true;
                            break;
                        }
                    }
                }

                if (flag_bad) break;
                
                // Remove last element and add points from next segment
                if (!vec_dQ_dx.empty()) {
                    vec_dQ_dx.pop_back();
                }
                
                const auto& next_fits = next_segment->fits();
                if (next_fits.empty()) break;
                
                // Determine direction for next segment (use squared distances to avoid sqrt)
                bool next_vertex_at_front = false;
                if (!next_segment->wcpts().empty()) {
                    const auto& nvp = next_vertex->wcpt().point;
                    const auto& nfp = next_segment->wcpts().front().point;
                    const auto& nbp = next_segment->wcpts().back().point;
                    double d1sq = (nvp - nfp).magnitude2();
                    double d2sq = (nvp - nbp).magnitude2();
                    next_vertex_at_front = (d1sq < d2sq);
                }
                
                // Add dQ/dx from next segment
                if (next_vertex_at_front) {
                    for (size_t i = 0; i < next_fits.size(); i++) {
                        double dQ_dx_normalized = next_fits[i].dQ / (next_fits[i].dx + 1e-9) / MIP_dQdx;
                        vec_dQ_dx.push_back(dQ_dx_normalized);
                        if (vec_dQ_dx.size() >= (size_t)limit) break;
                    }
                } else {
                    for (int i = (int)next_fits.size() - 1; i >= 0; i--) {
                        double dQ_dx_normalized = next_fits[i].dQ / (next_fits[i].dx + 1e-9) / MIP_dQdx;
                        vec_dQ_dx.push_back(dQ_dx_normalized);
                        if (vec_dQ_dx.size() >= (size_t)limit) break;
                    }
                }
                
                if (vec_dQ_dx.size() >= (size_t)limit) break;
                
                // Prepare for next iteration
                curr_vertex = next_vertex;
                curr_segment = next_segment;
                count++;
            }
        }
        
        return vec_dQ_dx;
    }

    void Shower::calculate_kinematics(const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool exclude_start_vertex_from_endpoint, bool endpoint_skip_orphan_vertices){
        // doc pr/91 round 1 -- the set of view nodes an actual member segment
        // touches.  Built once here and consulted by both farthest-vertex
        // searches below.  Empty (and never consulted) when the knob is off, so
        // the legacy path does not even pay for the walk.
        std::set<node_descriptor> pr91_touched;
        if (endpoint_skip_orphan_vertices) {
            for (auto edesc : ordered_edges(*this, m_full_graph)) {
                SegmentPtr seg = m_full_graph[edesc].segment;
                if (!seg || !seg->descriptor_valid()) continue;
                pr91_touched.insert(boost::source(edesc, m_full_graph));
                pr91_touched.insert(boost::target(edesc, m_full_graph));
            }
        }
        int nsegments = this->edges().size();
        
        if (nsegments == 1) {
            // Single segment case
            if (!m_start_segment) return;
            
            // Set particle type from start segment
            if (m_start_segment->has_particle_info()) {
                data.particle_type = m_start_segment->particle_info()->pdg();
            }
            
            // Mirrors prototype: ProtoSegment::get_flag_shower() = flag_shower_trajectory || flag_shower_topology || get_flag_shower_dQdx()
            // where get_flag_shower_dQdx() = (|particle_type| == 11)
            bool flag_shower = m_start_segment->flags_any(SegmentFlags::kShowerTrajectory) ||
                             m_start_segment->flags_any(SegmentFlags::kShowerTopology) ||
                             (m_start_segment->has_particle_info() && std::abs(m_start_segment->particle_info()->pdg()) == 11);
            if (flag_shower) set_flags(ShowerFlags::kShower);
            else unset_flags(ShowerFlags::kShower);
            
            // Calculate energies
            double seg_length = segment_track_length(m_start_segment);
            data.kenergy_range = cal_kine_range(seg_length, data.particle_type, particle_data);
            data.kenergy_dQdx = segment_cal_kine_dQdx(m_start_segment, recomb_model);
            
            // Calculate kenergy_best
            if (data.start_connection_type == 1) {
                data.kenergy_best = (seg_length < 4 * units::cm) ? data.kenergy_dQdx : data.kenergy_range;
            } else {
                if (flag_shower) {
                    data.kenergy_best = 0;
                } else {
                    data.kenergy_best = (seg_length < 4 * units::cm) ? data.kenergy_dQdx : data.kenergy_range;
                }
            }
            
            // Calculate start_point and end_point
            const auto& fits = m_start_segment->fits();
            bool has_fit_pcloud = (this->dpcloud("fit") != nullptr);
            SPDLOG_LOGGER_TRACE(s_log,
                "calculate_kinematics start_point: nseg=1 conn_type={} nfits={} dirsign={}"
                " has_fit_pcloud={} has_start_vertex={} vtx_fit_valid={}",
                data.start_connection_type, fits.size(), m_start_segment->dirsign(),
                has_fit_pcloud ? 1 : 0,
                m_start_vertex ? 1 : 0,
                (m_start_vertex && m_start_vertex->fit().valid()) ? 1 : 0);
            if (data.start_connection_type == 1 || !has_fit_pcloud) {
                if (!fits.empty()) {
                    if (m_start_segment->dirsign() == 1) {
                        data.start_point = fits.front().point;
                        data.end_point = fits.back().point;
                    } else if (m_start_segment->dirsign() == -1) {
                        data.start_point = fits.back().point;
                        data.end_point = fits.front().point;
                    }
                    // std::cout << "cluster: " << m_start_segment->cluster()->get_cluster_id() << " segment: " << m_start_segment->get_graph_index() << " " << data.end_point << std::endl;
                }
                SPDLOG_LOGGER_TRACE(s_log,
                    "calculate_kinematics start_point:   branch=fits start=({:.1f},{:.1f},{:.1f})cm",
                    data.start_point.x()/units::cm, data.start_point.y()/units::cm, data.start_point.z()/units::cm);
            } else {
                if (m_start_vertex) {
                    auto [sgcp_dis, sgcp_pt] = shower_get_closest_point(*this, m_start_vertex->fit().point, "fit");
                    SPDLOG_LOGGER_TRACE(s_log,
                        "calculate_kinematics start_point:   branch=shower_get_closest vtx=({:.1f},{:.1f},{:.1f})cm"
                        " closest_dis={:.3f}cm closest_pt=({:.1f},{:.1f},{:.1f})cm fit_pcloud_npts={}",
                        m_start_vertex->fit().point.x()/units::cm,
                        m_start_vertex->fit().point.y()/units::cm,
                        m_start_vertex->fit().point.z()/units::cm,
                        sgcp_dis/units::cm,
                        sgcp_pt.x()/units::cm, sgcp_pt.y()/units::cm, sgcp_pt.z()/units::cm,
                        this->dpcloud("fit") ? (int)this->dpcloud("fit")->npoints() : -1);
                    data.start_point = sgcp_pt;
                    // Fallback: if "fit" pcloud is absent or empty, use fits directly (same as conn_type==1)
                    if (data.start_point.x() == 0 && data.start_point.y() == 0 && data.start_point.z() == 0) {
                        SPDLOG_LOGGER_TRACE(s_log, "calculate_kinematics start_point:   shower_get_closest returned (0,0,0), applying fits fallback");
                        if (!fits.empty()) {
                            data.start_point = (m_start_segment->dirsign() == -1) ? fits.back().point : fits.front().point;
                        }
                    }

                    // Find farthest vertex — ordered_nodes gives index-stable tie-breaking
                    probe_endpoint(*this, m_full_graph, data.start_point, m_start_vertex,
                                   data.start_connection_type,
                                   exclude_start_vertex_from_endpoint,
                                   endpoint_skip_orphan_vertices, "single_seg");
                    double max_dis = 0;
                    const auto& view = this->view_graph();
                    for (auto vdesc : ordered_nodes(*this, m_full_graph)) {
                        VertexPtr vtx = view[vdesc].vertex;
                        if (!vtx) continue;
                        if (exclude_start_vertex_from_endpoint && vtx == m_start_vertex) continue;
                        // doc pr/91 F1: never end on a vertex no member segment reaches.
                        if (endpoint_skip_orphan_vertices && !pr91_touched.count(vdesc)) continue;
                        double dis = (data.start_point - vtx->fit().point).magnitude();
                        if (dis > max_dis) {
                            max_dis = dis;
                            data.end_point = vtx->fit().point;
                        }
                    }
                    // std::cout << "Vertex 1: " <<  data.end_point << std::endl;

                }
            }
            
            // Calculate init_dir
            if (data.start_connection_type == 1) {
                data.init_dir = segment_cal_dir_3vector(m_start_segment);
            } else if (data.start_connection_type == 2 || data.start_connection_type == 3) {
                if (m_start_vertex) {
                    data.init_dir = (data.start_point - m_start_vertex->fit().point).norm();
                }
            }
            // Fallback: if direction is still zero, use shower_cal_dir_3vector from start vertex
            if (data.init_dir.magnitude() == 0 && m_start_vertex) {
                data.init_dir = shower_cal_dir_3vector(*this, m_start_vertex->fit().point, 12 * units::cm);
            }
            
        } else {
            // Multiple segments case
            if (!m_start_segment) return;
            
            // Count connected segments via BFS without building the full result sets
            int nconnected_segs = count_connected_segments(m_start_segment);
            
            // Set particle type
            if (m_start_segment->has_particle_info()) {
                data.particle_type = m_start_segment->particle_info()->pdg();
            }
            
            // Mirrors prototype: ProtoSegment::get_flag_shower() = flag_shower_trajectory || flag_shower_topology || get_flag_shower_dQdx()
            // where get_flag_shower_dQdx() = (|particle_type| == 11)
            bool flag_shower = m_start_segment->flags_any(SegmentFlags::kShowerTrajectory) ||
                             m_start_segment->flags_any(SegmentFlags::kShowerTopology) ||
                             (m_start_segment->has_particle_info() && std::abs(m_start_segment->particle_info()->pdg()) == 11);
            if (flag_shower) set_flags(ShowerFlags::kShower);
            else unset_flags(ShowerFlags::kShower);

            // Common preamble for both single- and multi-track cases — start_point,
            // init_dir, end_point, and dQ/dx collection are identical; only the
            // final energy quantities differ.
            const auto& fits = m_start_segment->fits();
            double seg_length = segment_track_length(m_start_segment);
            const auto& view = this->view_graph();

            // Calculate start_point
            if (data.start_connection_type == 1 || !this->dpcloud("fit")) {
                if (!fits.empty()) {
                    if (m_start_segment->dirsign() == 1) {
                        data.start_point = fits.front().point;
                    } else if (m_start_segment->dirsign() == -1) {
                        data.start_point = fits.back().point;
                    }
                }
            } else {
                if (m_start_vertex) {
                    auto [sgcp_dist, sgcp_pt] = shower_get_closest_point(*this, m_start_vertex->fit().point, "fit");
                    if (sgcp_dist >= 0) {
                        // Valid closest-point found in "fit" pcloud.
                        data.start_point = sgcp_pt;
                    } else {
                        // Fallback: "fit" pcloud absent or empty — use segment endpoints directly.
                        // NOTE: do NOT test sgcp_pt == (0,0,0) as sentinel; a legitimate hit at
                        // the origin would falsely trigger the fallback (B16.1 in review).
                        if (!fits.empty()) {
                            data.start_point = (m_start_segment->dirsign() == -1) ? fits.back().point : fits.front().point;
                        }
                    }
                }
            }

            // Calculate init_dir
            if (data.start_connection_type == 1) {
                if (seg_length > 8 * units::cm) {
                    data.init_dir = segment_cal_dir_3vector(m_start_segment);
                } else if (m_start_vertex) {
                    data.init_dir = shower_cal_dir_3vector(*this, m_start_vertex->fit().point, 12 * units::cm);
                }
            } else if (data.start_connection_type == 2 || data.start_connection_type == 3) {
                if (m_start_vertex) {
                    data.init_dir = (data.start_point - m_start_vertex->fit().point).norm();
                }
            }
            // Fallback: if direction is still zero, use shower_cal_dir_3vector from start vertex
            if (data.init_dir.magnitude() == 0) {
                SPDLOG_LOGGER_TRACE(s_log,
                    "calculate_kinematics: nseg={} conn_type={} seg_length={:.1f}cm"
                    " seg_nfits={} seg_dirsign={} — init_dir is zero, applying fallback",
                    get_num_segments(), data.start_connection_type,
                    seg_length / units::cm,
                    m_start_segment->fits().size(), m_start_segment->dirsign());
                if (m_start_vertex) {
                    data.init_dir = shower_cal_dir_3vector(*this, m_start_vertex->fit().point, 12 * units::cm);
                }
            }

            // Find farthest vertex for end_point — ordered_nodes gives index-stable tie-breaking
            probe_endpoint(*this, m_full_graph, data.start_point, m_start_vertex,
                           data.start_connection_type,
                           exclude_start_vertex_from_endpoint,
                           endpoint_skip_orphan_vertices, "multi_seg");
            double max_dis = 0;
            for (auto vdesc : ordered_nodes(*this, m_full_graph)) {
                VertexPtr vtx = view[vdesc].vertex;
                if (!vtx) continue;
                if (exclude_start_vertex_from_endpoint && vtx == m_start_vertex) continue;
                // doc pr/91 F1: never end on a vertex no member segment reaches.
                if (endpoint_skip_orphan_vertices && !pr91_touched.count(vdesc)) continue;
                double dis = (data.start_point - vtx->fit().point).magnitude();
                if (dis > max_dis) {
                    max_dis = dis;
                    data.end_point = vtx->fit().point;
                }
            }

            // Collect dQ and dx from all segments; accumulate total_length for range-based energy
            double total_length = 0;
            std::vector<double> vec_dQ, vec_dx;
            // ordered_edges: this is the site that was actually caught.  Both
            // the `total_length +=` and the push order of vec_dQ/vec_dx are
            // observable -- cal_kine_dQdx() is a plain `kine_energy += dE` over
            // the vector (PRSegmentFunctions.cxx:1281), so the summation order
            // sets the last bit.  showers[1]/kine_dQdx moved 1314.124434586102
            // -> ...103 (rel 6.9e-16) between two `setarch -R` runs of the same
            // binary on SBND evt 388 (doc pr/28 sec 15).
            for (auto edesc : ordered_edges(*this, m_full_graph)) {
                SegmentPtr seg = view[edesc].segment;
                if (!seg) continue;
                total_length += segment_track_length(seg);
                const auto& seg_fits = seg->fits();
                for (const auto& fit : seg_fits) {
                    vec_dQ.push_back(fit.dQ);
                    vec_dx.push_back(fit.dx);
                }
            }

            // Calculate energies — only final quantities differ between single/multi-track
            data.kenergy_dQdx = cal_kine_dQdx(vec_dQ, vec_dx, recomb_model);
            if (nsegments == nconnected_segs) {
                // Single track: range-based energy is meaningful
                data.kenergy_range = cal_kine_range(total_length, data.particle_type, particle_data);
                if (data.start_connection_type == 1) {
                    data.kenergy_best = (seg_length < 4 * units::cm) ? data.kenergy_dQdx : data.kenergy_range;
                } else {
                    data.kenergy_best = flag_shower ? 0 : ((seg_length < 4 * units::cm) ? data.kenergy_dQdx : data.kenergy_range);
                }
            } else {
                // Multiple tracks: range energy not meaningful
                data.kenergy_range = 0;
                data.kenergy_best = 0;
            }
        }

        // size_t nclusters = 0;
        // {
        //     std::unordered_set<const Facade::Cluster*> cluster_set;
        //     const auto& view = this->view_graph();
        //     for (auto edesc : this->edges()) {
        //         SegmentPtr seg = view[edesc].segment;
        //         if (seg && seg->cluster()) {
        //             cluster_set.insert(seg->cluster());
        //         }
        //     }
        //     nclusters = cluster_set.size();
        // }

        // std::cout << "Shower::calculate_kinematics: nsegments=" << nsegments << " nvertices=" << this->nodes().size() << " " << this->edges().size()
        //           << " nclusters:" << nclusters
        //           << " particle_type=" << data.particle_type
        //           << " kenergy_range=" << data.kenergy_range / units::MeV << " MeV"
        //           << " kenergy_dQdx=" << data.kenergy_dQdx / units::MeV << " MeV"
        //           << " kenergy_best=" << data.kenergy_best / units::MeV << " MeV"
        //           << " kenergy_charge=" << data.kenergy_charge / units::MeV << " MeV"
        //           << " end_point=(" << data.end_point.x() / units::cm << "," << data.end_point.y() / units::cm << "," << data.end_point.z() / units::cm << ") cm"
        //           << std::endl;
    }

    void Shower::calculate_kinematics_long_muon(IndexedSegmentSet& segments_in_muons, const Clus::ParticleDataSet::pointer& particle_data, const IRecombinationModel::pointer& recomb_model, bool exclude_start_vertex_from_endpoint, int best_mode, double ratio_lo, double ratio_hi){
        // Invariant: this function is only called when shower->get_particle_type() == 13
        // (NeutrinoEnergyReco.cxx), which requires shower->set_particle_type(13) to have been
        // called (NeutrinoShowerClustering.cxx:118), which in turn requires m_start_segment to
        // already have particle_info() with pdg ±13.  The guard below is therefore logically
        // unreachable in normal execution; it is kept as a defensive check so that a caller
        // mistake produces a silent no-op rather than a null-dereference.
        if (!m_start_segment || !m_start_segment->has_particle_info()) return;
        int particle_type = abs(m_start_segment->particle_info()->pdg());
        // double particle_mass = m_start_segment->particle_info()->mass();
        
        unset_flags(ShowerFlags::kKinematics);
        unset_flags(ShowerFlags::kShower);  // long muon is not a shower (mirrors prototype's flag_shower = false)

        // Single pass over edges: accumulate length for muon segments, collect all dQ/dx,
        // and record which vertices touch a muon segment (avoids a second nested loop later).
        double total_length = 0;
        std::vector<double> vec_dQ;
        std::vector<double> vec_dx;
        // Use a set keyed by vertex index (stable across runs) to deduplicate muon vertices.
        // Ordered by index so subsequent max-distance search is deterministic on ties.
        std::map<size_t, VertexPtr> muon_vertices_by_index;

        for (auto edesc : ordered_edges(*this, m_full_graph)) {
            if (!has_edge(edesc)) continue;
            auto seg = view_graph()[edesc].segment;
            if (!seg) continue;

            bool in_muons = (segments_in_muons.find(seg) != segments_in_muons.end());
            if (in_muons) {
                total_length += segment_track_length(seg);
                auto [va, vb] = find_vertices(m_full_graph, seg);
                if (va) muon_vertices_by_index[m_full_graph[va->get_descriptor()].index] = va;
                if (vb) muon_vertices_by_index[m_full_graph[vb->get_descriptor()].index] = vb;
            }

            const auto& seg_fits = seg->fits();
            for (const auto& fit : seg_fits) {
                vec_dQ.push_back(fit.dQ);
                vec_dx.push_back(fit.dx);
            }
        }

        // Calculate kinetic energies
        data.kenergy_range = cal_kine_range(total_length, particle_type, particle_data);
        data.kenergy_dQdx = cal_kine_dQdx(vec_dQ, vec_dx, recomb_model);

        // For long muon, use dQdx as best energy
        data.kenergy_best = data.kenergy_dQdx;

        // Calculate initial direction from start segment
        data.init_dir = segment_cal_dir_3vector(m_start_segment);

        // Set start point based on direction
        auto& fits = m_start_segment->fits();
        int dirsign_val = m_start_segment->dirsign();
        if (dirsign_val == 1) {
            data.start_point = fits.front().point;
        } else {
            data.start_point = fits.back().point;
        }

        // Find farthest muon-connected vertex. Iteration is in ascending index order,
        // so on a distance tie the vertex with the smaller index wins — deterministic.
        double max_dis = 0;
        VertexPtr farthest_vertex = nullptr;
        for (auto& [idx, vtx] : muon_vertices_by_index) {
            if (!vtx) continue;
            if (exclude_start_vertex_from_endpoint && vtx == m_start_vertex) continue;
            double dis = (vtx->fit().point - data.start_point).magnitude();
            if (dis > max_dis) {
                max_dis = dis;
                farthest_vertex = vtx;
            }
        }
        
        // doc pr/101 (K4): range over the muon chain as best energy.
        // best_mode 0 keeps the legacy dQdx assignment above untouched.
        if (best_mode != 0) {
            int end_degree = -1;
            if (farthest_vertex && farthest_vertex->descriptor_valid()) {
                end_degree = static_cast<int>(boost::out_degree(farthest_vertex->get_descriptor(), m_full_graph));
            }
            const double ratio = data.kenergy_range > 0 ? data.kenergy_dQdx / data.kenergy_range : -1.0;
            bool use_range = data.kenergy_range > 0;
            if (best_mode == 2) {
                use_range = use_range && end_degree == 1
                    && ratio >= 1.0 - ratio_lo && ratio <= 1.0 + ratio_hi;
            }
            if (use_range) data.kenergy_best = data.kenergy_range;
            SPDLOG_LOGGER_DEBUG(s_log,
                "kine_long_muon: shower id={} nseg_chain={} L_cm={:.1f} range={:.1f} dqdx={:.1f} ratio={:.2f} "
                "end_degree={} mode={} used={}",
                m_shower_id, muon_vertices_by_index.size() ? muon_vertices_by_index.size() - 1 : 0,
                total_length / units::cm, data.kenergy_range / units::MeV, data.kenergy_dQdx / units::MeV,
                ratio, end_degree, best_mode, use_range ? "range" : "dqdx");
        }

        // Set end point to the farthest vertex
        if (farthest_vertex) {
            data.end_point = farthest_vertex->fit().point;
        } else {
            // Fallback: use the other end of start segment
            auto& fits = m_start_segment->fits();
            if (m_start_segment->dirsign() == 1) {
                data.end_point = fits.back().point;
            } else {
                data.end_point = fits.front().point;
            }
        }
    }
}