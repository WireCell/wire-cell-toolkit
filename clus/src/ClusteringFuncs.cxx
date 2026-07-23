#include <WireCellClus/ClusteringFuncs.h>
#include "WireCellUtil/Array.h"


#include <iostream>              // temp debug

using namespace WireCell::Clus::Facade;


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
    const std::string& orig_id_aname, bool flags_from_longest)
{
    std::unordered_map<int, int> desc2id;
    std::map<int, std::set<int> > id2desc;
    /*int num_components =*/ boost::connected_components(g, boost::make_assoc_property_map(desc2id));
    for (const auto& [desc, id] : desc2id) {
        id2desc[id].insert(desc);
    }

    std::vector<Cluster*> fresh;

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
    auto orig_clusters = grouping.children();

    const bool savecc = aname.size() > 0 && pcname.size() > 0;
    const bool save_origid = orig_id_aname.size() > 0 && pcname.size() > 0;

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
                    have_flash = true;
                    if (flags_from_longest) flash_flags = snapshot_flags(live);
                }
            }
            if (flags_from_longest) {
                const double any_len = live->get_length();
                if (any_len > best_any_len) {
                    best_any_len = any_len;
                    longest_flags = snapshot_flags(live);
                }
            }

            fresh_cluster.from(*live);
            fresh_cluster.take_children(*live, true);

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

            grouping.destroy_child(live);
            assert(live == nullptr);
        }
        if (savecc) {
            fresh_cluster.put_pcarray(cc, aname, pcname);
        }
        if (save_origid) {
            fresh_cluster.put_pcarray(orig_id, orig_id_aname, pcname);
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
