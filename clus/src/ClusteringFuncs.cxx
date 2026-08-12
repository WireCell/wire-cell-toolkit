#include <WireCellClus/ClusteringFuncs.h>
#include "WireCellUtil/Array.h"
#include "WireCellUtil/Logging.h"

#include <unordered_map>


#include <iostream>              // temp debug

using namespace WireCell::Clus::Facade;

namespace {
    // Named logger for the iso-band veto (doc pr/66).  Emitted from a hot,
    // generic function (merge_clusters), so it goes on its own sink rather
    // than std::cout -- unlike ClusteringNeutrino's own refusal marker, which
    // stays on std::cout to leave today's SBND stdout untouched.
    WireCell::Log::logptr_t band_veto_log() {
        static WireCell::Log::logptr_t l = WireCell::Log::logger("clus.bandveto");
        return l;
    }

    // Blob-center drift-x extent, the same measure ClusteringNeutrino's
    // iso_band_like()/blob_center_xext() use (clustering_neutrino.cxx).  Kept
    // as a separate copy here (M10) so the veto's log line can report the same
    // quantity ClusteringNeutrino's own marker line reports, without adding a
    // cross-file dependency for one helper.
    double band_veto_xext(const WireCell::Clus::Facade::Cluster* c) {
        double xmin = 1e300, xmax = -1e300;
        for (const auto* b : c->children()) {
            const double x = b->center_pos().x();
            xmin = std::min(xmin, x);
            xmax = std::max(xmax, x);
        }
        return xmax - xmin;
    }

    // True if cluster `c` carries at least one blob whose "nu_band_veto_role"
    // value equals `role`.  Fail-open (false) on any absent/mismatched array,
    // per band_veto_forbids()'s contract in the header.
    bool band_veto_has_role(const WireCell::Clus::Facade::Cluster* c, int role,
                            const std::string& pcname) {
        if (!c->has_pcarray<int>("nu_band_veto_role", pcname)) return false;
        auto sr = c->get_pcarray<int>("nu_band_veto_role", pcname);
        if (sr.size() != (size_t) c->nchildren()) return false;
        for (size_t i = 0; i < sr.size(); ++i) {
            if (sr[i] == role) return true;
        }
        return false;
    }
}

bool WireCell::Clus::Facade::band_veto_forbids(const Cluster* a, const Cluster* b,
                                               const std::string& pcname)
{
    if (!a || !b) return false;
    return (band_veto_has_role(a, band_veto_band, pcname) &&
            band_veto_has_role(b, band_veto_nonband, pcname))
        || (band_veto_has_role(a, band_veto_nonband, pcname) &&
            band_veto_has_role(b, band_veto_band, pcname));
}

bool WireCell::Clus::Facade::cluster_has_band_veto_role(const Cluster* c, int role,
                                                        const std::string& pcname)
{
    if (!c) return false;
    return band_veto_has_role(c, role, pcname);
}


// Add this to your clustering_util.cxx file

std::tuple<geo_point_t, double, double, double> 
WireCell::Clus::Facade::extract_geometry_params(
    const Grouping& grouping,
    const IDetectorVolumes::pointer dv)
{
    geo_point_t drift_dir(1, 0, 0);  // initialize drift direction
    double angle_u = 0, angle_v = 0, angle_w = 0;  // initialize angles

    // Find the first valid WirePlaneId in the grouping
    for (const auto& gwpid : grouping.wpids()) {
        // Update drift direction based on face orientation
        int face_dirx = dv->face_dirx(gwpid);
        drift_dir.x(face_dirx);
        
        // Create wpids for all three planes with the same APA and face
        WirePlaneId wpid_u(kUlayer, gwpid.face(), gwpid.apa());
        WirePlaneId wpid_v(kVlayer, gwpid.face(), gwpid.apa());
        WirePlaneId wpid_w(kWlayer, gwpid.face(), gwpid.apa());
        
        // Get wire directions for all planes
        Vector wire_dir_u = dv->wire_direction(wpid_u);
        Vector wire_dir_v = dv->wire_direction(wpid_v);
        Vector wire_dir_w = dv->wire_direction(wpid_w);
        
        // Calculate angles
        angle_u = std::atan2(wire_dir_u.z(), wire_dir_u.y());
        angle_v = std::atan2(wire_dir_v.z(), wire_dir_v.y());
        angle_w = std::atan2(wire_dir_w.z(), wire_dir_w.y());
        
        // Only need to process the first valid WirePlaneId
        break;
    }
    
    return std::make_tuple(drift_dir, angle_u, angle_v, angle_w);
}

void WireCell::Clus::Facade::validate_drift_group(
    const std::set<WirePlaneId>& wpids,
    const IDetectorVolumes::pointer dv,
    bool allow_mixed_faces,
    const std::string& who)
{
    if (wpids.empty()) return;

    auto wit = wpids.begin();
    const int common_face = wit->face();
    double FV_xmin = dv->metadata(*wit)["FV_xmin"].asDouble();
    double FV_xmax = dv->metadata(*wit)["FV_xmax"].asDouble();
    double FV_xmin_margin = dv->metadata(*wit)["FV_xmin_margin"].asDouble();
    double FV_xmax_margin = dv->metadata(*wit)["FV_xmax_margin"].asDouble();

    for (++wit; wit != wpids.end(); ++wit) {
        const auto md = dv->metadata(*wit);
        if ((!allow_mixed_faces && wit->face() != common_face) ||
            md["FV_xmin"].asDouble() != FV_xmin || md["FV_xmax"].asDouble() != FV_xmax ||
            md["FV_xmin_margin"].asDouble() != FV_xmin_margin ||
            md["FV_xmax_margin"].asDouble() != FV_xmax_margin) {
            for (const auto& wpid : wpids) {
                std::cout << who << " wpid: " << wpid.name() << std::endl;
            }
            raise<ValueError>("%s: grouping has %d wpids with mixed faces or differing FV x metadata",
                              who, wpids.size());
        }
    }
}

std::vector<Cluster*> WireCell::Clus::Facade::merge_clusters(
    cluster_connectivity_graph_t& g,
    Grouping& grouping,
    const std::string& aname, const std::string& pcname,
    const std::string& orig_id_aname, bool flags_from_longest,
    const std::string& orig_main_aname,
    const std::string& orig_wasmain_aname)
{
    // Note, here we do an unusual thing and COPY the vector of children
    // facades.  In most simple access we would get the reference to the child
    // vector to save a little copy time.  We explicitly copy here as we must
    // preserve the original order of children facades even as we remove them
    // from the grouping.  As each child facade is removed, it's
    // unique_ptr<node> is returned which we ignore/drop and thus the child
    // facade dies along with its node.  This leaves the orig_clusters element
    // that was just holding the pointer to the doomed facade now holding
    // invalid memory.  But, it is okay as we never revisit the same cluster in
    // the grouping.  All that to explain a missing "&"! :)
    //
    // Hoisted above connected_components (was below, unchanged after that):
    // the iso-band veto below needs the vertex-index -> Cluster* mapping
    // BEFORE components are computed, to filter edges first.  Reading it here
    // is still safe -- nothing has destroyed a child yet at this point in the
    // function.
    auto orig_clusters = grouping.children();

    // Iso-band veto (ClusteringNeutrino's record_band_veto, doc pr/66).  A
    // per-APA clustering_neutrino pass may have REFUSED to merge an
    // isochronous band with a drift-spanning partner and recorded that
    // refusal as per-blob "nu_band_veto_role" provenance; the all-APA chain is
    // a SEPARATE MultiAlgBlobClustering pnode with no iso-band guard of its
    // own and would otherwise re-merge them.  Honor the refusal here, once,
    // for every merge stage in the codebase -- ClusteringExtendLoop's inner
    // clustering_extend() calls have no configuration surface of their own,
    // so a per-stage veto could not reach them at all.
    //
    // Purely PRESENCE gated: with the writer knob off "nu_band_veto_role" is
    // never created, so band_veto_forbids() is constant-false and not one
    // edge is ever removed => byte-identical.  Drops the EDGE only (not the
    // whole component): a third, unmarked cluster may still bridge the two,
    // exactly as at the per-APA refusal site (ClusteringNeutrino only skips
    // add_edge for the one pair it tested).
    {
        std::vector<std::pair<cluster_connectivity_graph_t::vertex_descriptor,
                              cluster_connectivity_graph_t::vertex_descriptor>> drop;
        for (auto er = boost::edges(g); er.first != er.second; ++er.first) {
            const auto u = boost::source(*er.first, g);
            const auto v = boost::target(*er.first, g);
            const int i1 = g[u], i2 = g[v];
            if (i1 < 0 || i2 < 0) continue;
            if ((size_t) i1 >= orig_clusters.size() || (size_t) i2 >= orig_clusters.size()) continue;
            const Cluster* c1 = orig_clusters[i1];
            const Cluster* c2 = orig_clusters[i2];
            if (!band_veto_forbids(c1, c2, pcname)) continue;
            band_veto_log()->info(
                "nu_band_veto: dropped edge, lens {:.1f}/{:.1f} cm, xext {:.1f}/{:.1f} cm",
                c1->get_length() / units::cm, c2->get_length() / units::cm,
                band_veto_xext(c1) / units::cm, band_veto_xext(c2) / units::cm);
            drop.emplace_back(u, v);
        }
        // Remove by (u,v) pair after the scan completes: edge descriptors of a
        // vecS/vecS adjacency_list are invalidated by remove_edge, so a second
        // pass over descriptors collected during the scan would be UB.
        for (const auto& [u, v] : drop) {
            boost::remove_edge(u, v, g);
        }
    }

    std::unordered_map<int, int> desc2id;
    std::map<int, std::set<int> > id2desc;
    /*int num_components =*/ boost::connected_components(g, boost::make_assoc_property_map(desc2id));
    for (const auto& [desc, id] : desc2id) {
        id2desc[id].insert(desc);
    }

    std::vector<Cluster*> fresh;

    const bool savecc = aname.size() > 0 && pcname.size() > 0;
    const bool save_origid = orig_id_aname.size() > 0 && pcname.size() > 0;
    const bool save_origmain = orig_main_aname.size() > 0 && pcname.size() > 0;
    // Per-blob "the member this blob came from was a bundle main before this
    // merge" (doc pr/20 Part I P1).  Independent of save_origid: the value is
    // read off each member's own flag, not derived from the recorded idents.
    const bool save_wasmain = orig_wasmain_aname.size() > 0 && pcname.size() > 0;

    // Provenance pairs CARRIED across this merge (doc 52 Stage 2, defect D4).
    // A merge can WRITE provenance (orig_id_aname/orig_main_aname) and a split
    // can CARVE it (clustering_switch_scope), but until now nothing could CARRY
    // an existing per-blob array THROUGH a merge -- so any provenance array died
    // at the next of the 12 merge_clusters() call sites.  Handled here, for the
    // whole codebase at once, rather than per call site: cathode_connect and the
    // generic all-APA passes merge with no provenance arguments at all, and
    // future passes would silently reintroduce the same bug.
    //
    // Semantics, per doc 52 4b:
    //  * the MAIN array is concatenated VERBATIM.  A "connection" merge (extend,
    //    regular, parallel_prolong, close, connect, deghost, live_dead,
    //    neutrino, cathode_connect) asserts that the charge is one particle; it
    //    must not rewrite the members' roles.  Promoting every member to main
    //    was considered and is WRONG: it would sweep up the associated fragments
    //    already grouped onto either half, and the isolated un-merge could then
    //    never split them off -- the original bug.  Carrying roles unchanged
    //    keeps a cathode crosser whole (both halves are mains, retained by
    //    ClusteringUnmergeBundle's union rule) while its delta rays still leave.
    //  * the ID array is REBASED per member into a fresh dense range, in
    //    ascending original-id order.  Ids are only unique within the scope that
    //    wrote them: the per-APA (SBND) / per-drift-group (PDHD, PDVD) stages
    //    renumber independently, so member A's group 7 and member B's group 7 are
    //    different objects and must not be pooled into one split-off cluster.
    //  * a member with NO array contributes one fresh id and main=1.  An
    //    ungrouped cluster is a main, not an associated fragment; defaulting it
    //    the other way would let the un-merge discard the far half of any
    //    crosser that never went through clustering_isolated.
    //
    // Fail open, exactly as clustering_switch_scope's carve() does: any size
    // disagreement drops the pair for that merged cluster rather than raising.
    // When no member carries the arrays nothing is written, so a build with the
    // writer knob off is byte-identical.
    static const std::vector<std::pair<std::string, std::string>> carry_pairs = {
        {"assoc_cluster_id", "assoc_cluster_main"},
    };

    // SINGLE per-blob arrays carried VERBATIM across this merge (doc pr/66).
    // Distinct from carry_pairs above: those hold ids that are only unique
    // within the scope that wrote them and so must be REBASED per member
    // (see the long comment above); these hold a value that is meaningful
    // across the whole event -- e.g. "nu_band_veto_role" is an intrinsic
    // per-cluster classification (band vs non-band), so rebasing it would
    // destroy the very fact being carried.  A member with no array
    // contributes zeros -- the documented "unmarked" value, not a sentinel.
    static const std::vector<std::string> carry_singles = {
        "nu_band_veto_role",   // ClusteringNeutrino record_band_veto
    };
    const bool do_carry = pcname.size() > 0;

    for (const auto& [id, descs] : id2desc) {
        if (descs.size() < 2) {
            continue;
        }

        // it starts with no cluster facade
        Cluster& fresh_cluster = grouping.make_child();

        std::vector<int> cc;
        int parent_id = 0;

        // Per-blob original cluster ident (the pre-merge ident() of each
        // sub-cluster), parallel to cc but keyed by the member's own ident.
        std::vector<int> orig_id;
        // Per-blob marker of the merged cluster's REPRESENTATIVE member (the
        // flash/flags donor chosen below), so ITS points stay identifiable
        // after the merge.  NB the "main_cluster" flag cannot serve here:
        // with flag_matched_mains every bundle's main carries it, and a flash
        // group can merge several bundles -- all mains of their own bundle.
        // The representative is the one main the merged cluster reports.
        std::vector<int> orig_main;
        // Per-blob record of flag_main_cluster as each member carried it ON THE
        // WAY IN.  Distinct from orig_main above: that one marks the single
        // representative, this one marks every member that was a matched bundle
        // main -- i.e. exactly the ones this merge is about to demote, which is
        // the fact that has no other surviving witness afterwards.
        std::vector<int> orig_wasmain;

        // Accumulators for the carried provenance pairs (see carry_pairs above).
        //
        // Keyed by BLOB POINTER, not by row position.  Positional copying is
        // wrong and was measurably wrong: a member's blobs do NOT necessarily
        // end up as a contiguous, order-preserving run of the merged cluster's
        // children, so the k-th row of a member's array is not the k-th of its
        // blobs afterwards.  The orig_id/orig_main arrays above cannot see this
        // because their value is CONSTANT within a member -- any permutation of
        // a member's rows is a no-op there -- which is why the flash pair looked
        // healthy while the carried pair came out scrambled (doc 52 §11:
        // evt284349 cluster 10 had the right group sizes, 4+5+4+311, attached to
        // the wrong blobs, so the un-merge split the trunk's tip and kept the
        // detached piece).  Blob facades survive take_children (only the nodes
        // move), so the pointer is a stable key.
        struct CarryAcc {
            std::vector<std::pair<Blob*, std::pair<int, int>>> rows;  // blob -> (id, main)
            int next_id{0};      // next fresh dense id to hand out
            bool any{false};     // did at least one member actually carry it?
        };
        std::vector<CarryAcc> carried(do_carry ? carry_pairs.size() : 0);

        // Same blob-pointer keying as CarryAcc above, and for the same
        // reason, for the single-array registry (carry_singles).
        struct CarrySingleAcc {
            std::vector<std::pair<Blob*, int>> rows;
            bool any{false};
        };
        std::vector<CarrySingleAcc> carried1(do_carry ? carry_singles.size() : 0);

        // Flash bookkeeping: Cluster::from() copies the first-encountered
        // member's cluster_t0/flash/matched_flash_gid (arbitrary std::set order).
        // Track the longest contributing member that carries a valid matched
        // flash so we can set a deterministic, physically-meaningful
        // representative flash on the merged cluster after the merge.  When no
        // member carries a flash (e.g. the entire per-APA stage, where flash
        // matching has not yet run) this stays a no-op, preserving prior output.
        bool have_flash = false;
        double best_len = -1;
        double best_t0 = 0;
        int best_flash = -1;
        int best_gid = -1;
        int best_flash_ident = -1;   // ident of the flash donor
        int best_any_ident = -1;     // ident of the longest member overall
        double best_any_len2 = -1;

        // Flag bookkeeping (flags_from_longest).  from() copies each member's
        // flag values in turn, so the last member visited wins -- a matched main
        // can lose flag_main_cluster to a tiny co-merged fragment.  Capture the
        // representative member's flags by value here (the member is destroyed
        // inside the loop) and re-apply them after the merge.  Representative =
        // the flash donor when one exists, else the longest member overall, so
        // flags and flash always describe the same member.
        std::vector<std::pair<std::string, int>> flash_flags, longest_flags;
        double best_any_len = -1;
        auto snapshot_flags = [](const Cluster* c) {
            std::vector<std::pair<std::string, int>> kv;
            for (const auto& fname : c->flag_names()) kv.emplace_back(fname, c->get_flag(fname));
            return kv;
        };

        for (const auto& desc : descs) {
            const int idx = g[desc];
            if (idx < 0) {  // no need anymore ...
                continue;
            }

            auto live = orig_clusters[idx];

            const int live_flash = live->get_scalar<int>("flash", -1);
            if (live_flash >= 0) {
                const double live_len = live->get_length();
                if (live_len > best_len) {
                    best_len = live_len;
                    best_t0 = live->get_cluster_t0();
                    best_flash = live_flash;
                    best_gid = live->get_scalar<int>("matched_flash_gid", -1);
                    best_flash_ident = live->ident();
                    have_flash = true;
                    if (flags_from_longest) flash_flags = snapshot_flags(live);
                }
            }
            if (save_origmain) {
                const double any_len2 = live->get_length();
                if (any_len2 > best_any_len2) {
                    best_any_len2 = any_len2;
                    best_any_ident = live->ident();
                }
            }
            if (flags_from_longest) {
                const double any_len = live->get_length();
                if (any_len > best_any_len) {
                    best_any_len = any_len;
                    longest_flags = snapshot_flags(live);
                }
            }

            // Snapshot the carried arrays while `live` is still alive -- it is
            // destroyed at the bottom of this loop body.
            std::vector<std::pair<std::vector<int>, std::vector<int>>> snap(carried.size());
            for (size_t ip = 0; ip < carried.size(); ++ip) {
                const auto& [id_aname, main_aname] = carry_pairs[ip];
                if (!live->has_pcarray<int>(id_aname, pcname)) continue;
                if (!live->has_pcarray<int>(main_aname, pcname)) continue;
                auto sid = live->get_pcarray<int>(id_aname, pcname);
                auto smain = live->get_pcarray<int>(main_aname, pcname);
                if (sid.size() != smain.size()) continue;
                if (sid.size() != (size_t) live->nchildren()) continue;
                snap[ip].first.assign(sid.begin(), sid.end());
                snap[ip].second.assign(smain.begin(), smain.end());
            }
            // Same snapshot, for the single-array (carry_singles) registry.
            std::vector<std::vector<int>> snap1(carried1.size());
            for (size_t is = 0; is < carried1.size(); ++is) {
                if (!live->has_pcarray<int>(carry_singles[is], pcname)) continue;
                auto sv = live->get_pcarray<int>(carry_singles[is], pcname);
                if (sv.size() != (size_t) live->nchildren()) continue;
                snap1[is].assign(sv.begin(), sv.end());
            }
            // The member's children in the SAME order as its per-blob array
            // rows (both are "children order" while the member is still intact).
            std::vector<Blob*> member_kids;
            if (!carried.empty() || !carried1.empty()) member_kids = live->children();
            const size_t nb_before = fresh_cluster.nchildren();

            fresh_cluster.from(*live);
            fresh_cluster.take_children(*live, true);

            // Record this member's rows against the blobs they belong to.
            const size_t nb_added = (size_t) fresh_cluster.nchildren() - nb_before;
            for (size_t ip = 0; ip < carried.size(); ++ip) {
                auto& acc = carried[ip];
                const auto& sid = snap[ip].first;
                if (sid.size() == member_kids.size() && sid.size() == nb_added
                    && nb_added > 0) {
                    // Rebase into a fresh dense range, ascending by original id
                    // => deterministic and collision-free across members.
                    std::map<int, int> remap;
                    for (int v : sid) remap.emplace(v, 0);
                    for (auto& [orig, fresh_id] : remap) { (void) orig; fresh_id = acc.next_id++; }
                    for (size_t k = 0; k < sid.size(); ++k) {
                        acc.rows.emplace_back(member_kids[k],
                                              std::make_pair(remap.at(sid[k]),
                                                             snap[ip].second[k]));
                    }
                    acc.any = true;
                }
                else {
                    // No usable array on this member: one fresh id, all main.
                    const int fresh_id = acc.next_id++;
                    for (Blob* b : member_kids) {
                        acc.rows.emplace_back(b, std::make_pair(fresh_id, 1));
                    }
                }
            }
            for (size_t is = 0; is < carried1.size(); ++is) {
                auto& acc = carried1[is];
                const auto& sv = snap1[is];
                if (sv.size() == member_kids.size() && sv.size() == nb_added && nb_added > 0) {
                    for (size_t k = 0; k < sv.size(); ++k) acc.rows.emplace_back(member_kids[k], sv[k]);
                    acc.any = true;
                }
                else {
                    // No usable array on this member: every blob is unmarked
                    // (0), the array's own documented default -- unlike
                    // carry_pairs, there is no "one fresh id" concept here.
                    for (Blob* b : member_kids) acc.rows.emplace_back(b, 0);
                }
            }

            if (fresh_cluster.ident() < 0) {
                fresh_cluster.set_ident(live->ident());
            }

            if (savecc) {
                cc.resize(fresh_cluster.nchildren(), parent_id);
                ++parent_id;
            }

            if (save_origid) {
                orig_id.resize(fresh_cluster.nchildren(), live->ident());
            }

            if (save_wasmain) {
                // resize-fills-the-new-rows, same idiom as orig_id above; the
                // value is constant within a member.
                orig_wasmain.resize(fresh_cluster.nchildren(),
                                    live->get_flag(Flags::main_cluster) ? 1 : 0);
            }

            grouping.destroy_child(live);
            assert(live == nullptr);
        }
        if (savecc) {
            fresh_cluster.put_pcarray(cc, aname, pcname);
        }
        if (save_origid) {
            fresh_cluster.put_pcarray(orig_id, orig_id_aname, pcname);
        }
        if (save_wasmain) {
            fresh_cluster.put_pcarray(orig_wasmain, orig_wasmain_aname, pcname);
        }
        // Re-attach the carried provenance, placing every row at the position
        // its own blob ended up in.  Only when some member really had the
        // array (else this is a no-op and output is unchanged); fail open on
        // any blob we cannot place or any row left unfilled.  Shared between
        // carry_pairs and carry_singles below -- both need the same
        // blob-pointer -> merged-row-index map.
        if (!carried.empty() || !carried1.empty()) {
            const auto merged_kids = fresh_cluster.children();
            std::unordered_map<const Blob*, size_t> where;
            where.reserve(merged_kids.size());
            for (size_t i = 0; i < merged_kids.size(); ++i) where[merged_kids[i]] = i;
            for (size_t ip = 0; ip < carried.size(); ++ip) {
                const auto& acc = carried[ip];
                if (!acc.any) continue;
                if (acc.rows.size() != merged_kids.size()) continue;
                std::vector<int> ids(merged_kids.size(), -1), mains(merged_kids.size(), -1);
                bool ok = true;
                for (const auto& [blob, val] : acc.rows) {
                    auto it = where.find(blob);
                    if (it == where.end()) { ok = false; break; }
                    ids[it->second] = val.first;
                    mains[it->second] = val.second;
                }
                if (!ok) continue;
                for (int v : mains) if (v < 0) { ok = false; break; }
                if (!ok) continue;
                fresh_cluster.put_pcarray(ids, carry_pairs[ip].first, pcname);
                fresh_cluster.put_pcarray(mains, carry_pairs[ip].second, pcname);
            }
            for (size_t is = 0; is < carried1.size(); ++is) {
                const auto& acc = carried1[is];
                if (!acc.any) continue;                                // nobody carried => no key written
                if (acc.rows.size() != merged_kids.size()) continue;   // fail open
                std::vector<int> vals(merged_kids.size(), -1);         // -1 = unfilled sentinel;
                                                                        // every legal value is >= 0
                bool ok = true;
                for (const auto& [blob, v] : acc.rows) {
                    auto it = where.find(blob);
                    if (it == where.end()) { ok = false; break; }
                    vals[it->second] = v;
                }
                if (!ok) continue;
                for (int v : vals) if (v < 0) { ok = false; break; }
                if (!ok) continue;
                fresh_cluster.put_pcarray(vals, carry_singles[is], pcname);
            }
        }
        if (save_origmain && save_origid) {
            // Representative = the same member whose flash/flags the merged
            // cluster carries (flash donor when one exists, else the longest).
            const int rep_ident = have_flash ? best_flash_ident : best_any_ident;
            orig_main.resize(orig_id.size());
            for (size_t i = 0; i < orig_id.size(); ++i) {
                orig_main[i] = (orig_id[i] == rep_ident) ? 1 : 0;
            }
            fresh_cluster.put_pcarray(orig_main, orig_main_aname, pcname);
        }

        // Override from()'s arbitrary first-wins flash with the longest
        // flash-bearing member's flash, so the merged cluster's cluster_t0 is a
        // real, deterministic member value (never 0 by accident).
        if (have_flash) {
            fresh_cluster.set_cluster_t0(best_t0);
            fresh_cluster.set_scalar<int>("flash", best_flash);
            fresh_cluster.set_scalar<int>("matched_flash_gid", best_gid);
        }

        // Same treatment for the flags: replace from()'s last-member-wins values
        // with the representative member's, so e.g. flag_main_cluster follows the
        // cluster the bundle was actually matched on instead of an arbitrary
        // co-merged fragment.  Flag names carried by other members but not by the
        // representative keep their union value (QLMatching materializes its flags
        // on every cluster, so in practice the name sets are identical).
        if (flags_from_longest) {
            const auto& rep = have_flash ? flash_flags : longest_flags;
            for (const auto& [fname, fval] : rep) fresh_cluster.set_flag(fname, fval);
        }

        // Normally, it would be weird/wrong to store an address of a reference.
        // But, we know the Cluster facade is held by the pc tree node that we
        // just added to the grouping node.
        fresh.push_back(&fresh_cluster);
    }

    return fresh;
}


std::map<const Cluster*, int> WireCell::Clus::Facade::assign_flash_t0_groups(
    const std::vector<Cluster*>& clusters, double window)
{
    std::map<const Cluster*, int> group_of;

    // Collect the matched clusters (valid "flash" scalar) and remember the rest
    // so they can be given unique singleton groups.
    std::vector<const Cluster*> matched;
    for (const Cluster* c : clusters) {
        if (c->get_scalar<int>("flash", -1) >= 0) {
            matched.push_back(c);
        }
    }

    // Sort matched clusters by their matched flash time.  Break ties on ident to
    // keep the assignment deterministic across rebuilds.
    std::sort(matched.begin(), matched.end(), [](const Cluster* a, const Cluster* b) {
        const double ta = a->get_cluster_t0();
        const double tb = b->get_cluster_t0();
        if (ta != tb) return ta < tb;
        return a->ident() < b->ident();
    });

    // Greedily start a new group whenever the gap to the previous flash time
    // exceeds the window.
    int next_group = 0;
    double prev_t0 = 0;
    for (size_t i = 0; i < matched.size(); ++i) {
        const double t0 = matched[i]->get_cluster_t0();
        if (i == 0 || (t0 - prev_t0) > window) {
            ++next_group;  // group ids for matched clusters are >= 1
        }
        group_of[matched[i]] = next_group;
        prev_t0 = t0;
    }

    // Every remaining (unmatched) cluster gets a unique singleton id, distinct
    // from all matched-group ids, so it can never share a group with anyone.
    for (const Cluster* c : clusters) {
        if (group_of.find(c) == group_of.end()) {
            group_of[c] = ++next_group;
        }
    }

    return group_of;
}


geo_vector_t WireCell::Clus::Facade::calc_pca_dir(const geo_point_t& center, const std::vector<geo_point_t>& points)
{
    // Create covariance matrix
    Eigen::MatrixXd cov_matrix(3, 3);

    // Calculate covariance matrix elements
    for (int i = 0; i != 3; i++) {
        for (int j = i; j != 3; j++) {
            cov_matrix(i, j) = 0;
            for (const auto& p : points) {
                if (i == 0 && j == 0) {
                    cov_matrix(i, j) += (p.x() - center.x()) * (p.x() - center.x());
                }
                else if (i == 0 && j == 1) {
                    cov_matrix(i, j) += (p.x() - center.x()) * (p.y() - center.y());
                }
                else if (i == 0 && j == 2) {
                    cov_matrix(i, j) += (p.x() - center.x()) * (p.z() - center.z());
                }
                else if (i == 1 && j == 1) {
                    cov_matrix(i, j) += (p.y() - center.y()) * (p.y() - center.y());
                }
                else if (i == 1 && j == 2) {
                    cov_matrix(i, j) += (p.y() - center.y()) * (p.z() - center.z());
                }
                else if (i == 2 && j == 2) {
                    cov_matrix(i, j) += (p.z() - center.z()) * (p.z() - center.z());
                }
            }
        }
    }

    // std::cout << "Test: " << center << " " << points.at(0) << std::endl;
    // std::cout << "Test: " << center << " " << points.at(1) << std::endl;
    // std::cout << "Test: " << center << " " << points.at(2) << std::endl;

    // Fill symmetric parts
    cov_matrix(1, 0) = cov_matrix(0, 1);
    cov_matrix(2, 0) = cov_matrix(0, 2);
    cov_matrix(2, 1) = cov_matrix(1, 2);

    // Calculate eigenvalues/eigenvectors using Eigen
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(cov_matrix);
    auto eigen_vectors = eigenSolver.eigenvectors();

    // std::cout << "Test: " << eigen_vectors(0,0) << " " << eigen_vectors(1,0) << " " << eigen_vectors(2,0) << std::endl;

    // Get primary direction (first eigenvector)
    double norm = sqrt(eigen_vectors(0, 2) * eigen_vectors(0, 2) + 
                      eigen_vectors(1, 2) * eigen_vectors(1, 2) + 
                      eigen_vectors(2, 2) * eigen_vectors(2, 2));

    return geo_vector_t(eigen_vectors(0, 2) / norm,
                       eigen_vectors(1, 2) / norm, 
                       eigen_vectors(2, 2) / norm);
}
