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

std::vector<Cluster*> WireCell::Clus::Facade::merge_clusters(
    cluster_connectivity_graph_t& g,
    Grouping& grouping,
    const std::string& aname, const std::string& pcname)
{
    std::unordered_map<int, int> desc2id;
    std::unordered_map<int, std::set<int> > id2desc;
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

    for (const auto& [id, descs] : id2desc) {
        if (descs.size() < 2) {
            continue;
        }

        // it starts with no cluster facade
        Cluster& fresh_cluster = grouping.make_child();

        std::vector<int> cc;
        int parent_id = 0;
        
        for (const auto& desc : descs) {
            const int idx = g[desc];
            if (idx < 0) {  // no need anymore ...
                continue;
            }

            auto live = orig_clusters[idx];
            fresh_cluster.from(*live);
            fresh_cluster.take_children(*live, true);

            if (fresh_cluster.ident() < 0) {
                fresh_cluster.set_ident(live->ident());
            }

            if (savecc) {
                cc.resize(fresh_cluster.nchildren(), parent_id);
                ++parent_id;
            }

            grouping.destroy_child(live);
            assert(live == nullptr);
        }
        if (savecc) {
            fresh_cluster.put_pcarray(cc, aname, pcname);
        }

        // Normally, it would be weird/wrong to store an address of a reference.
        // But, we know the Cluster facade is held by the pc tree node that we
        // just added to the grouping node.
        fresh.push_back(&fresh_cluster);
    }

    return fresh;
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
