#include <WireCellClus/ClusteringFuncs.h>

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Aux;
using namespace WireCell::Aux::TensorDM;
using namespace WireCell::PointCloud::Facade;
using namespace WireCell::PointCloud::Tree;

// The original developers do not care.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"

// #define __DEBUG__
#ifdef __DEBUG__
#define LogDebug(x) std::cout << "[yuhw]: " << __LINE__ << " : " << x << std::endl
#else
#define LogDebug(x)
#endif

// All APA Faces 
void WireCell::PointCloud::Facade::clustering_examine_bundles(Grouping& live_grouping, const IDetectorVolumes::pointer dv ,
    const std::string& pc_name,                        // point cloud name
    const std::vector<std::string>& coords,            // coordinate names
    const bool use_ctpc)
{
    // std::cout << "Test Examine Bundles" << std::endl;

    std::vector<Cluster *> live_clusters = live_grouping.children();
    Tree::Scope scope{pc_name, coords};
    for (size_t i=0;i!=live_clusters.size();i++){
        if (!live_clusters.at(i)->get_scope_filter(scope)) continue; // move on if the cluster is not in the scope filter ...
        if (live_clusters.at(i)->get_default_scope().hash() != scope.hash()) {
            live_clusters.at(i)->set_default_scope(scope);
            // std::cout << "Test: Set default scope: " << pc_name << " " << coords[0] << " " << coords[1] << " " << coords[2] << " " << cluster->get_default_scope().hash() << " " << scope.hash() << std::endl;
        }

        // if there is a cc component, record the main cluster as id of the blobs???
        auto old_cc_array = live_clusters.at(i)->get_pcarray("isolated", "perblob");
        
        // currently reset the cc component (todo: find the main component)

        // do the examine graph
        auto b2groupid = live_clusters.at(i)->examine_graph(dv, true);
        
        bool flag_largest = false;
        // Compare old and new cluster groupings
        if (old_cc_array.size() == b2groupid.size()) {
            // Check if main cluster exists in old array
            bool has_main = std::find(old_cc_array.begin(), old_cc_array.end(), -1) != old_cc_array.end();
            if (has_main) {
                // Count overlap between old main cluster and each new cluster
                std::map<int, int> overlap_counts;
                for (size_t j = 0; j < old_cc_array.size(); j++) {
                    if (old_cc_array[j] == -1) {  // this blob was in main cluster
                        overlap_counts[b2groupid[j]]++;
                    }
                }
                
                // Find cluster with maximum overlap
                int max_overlap = 0;
                int main_cluster_id = -1;
                for (const auto& pair : overlap_counts) {
                    if (pair.second > max_overlap) {
                        max_overlap = pair.second;
                        main_cluster_id = pair.first;
                    }
                }
                
                // Mark new main cluster
                for (auto& id : b2groupid) {
                    if (id == main_cluster_id) {
                        id = -1;
                    }
                }
            } else{
                flag_largest = true;
            }
        }else{
            flag_largest = true;
        }

        //run the longest ...
        if (flag_largest){
            // No main cluster existed, find longest cluster
            std::map<int, double> cluster_lengths;
            std::set<int> unique_ids;
            for (const auto& id : b2groupid) {
                if (id >= 0) {  // skip any existing -1 values
                    unique_ids.insert(id);
                }
            }
            
            // Calculate length for each unique cluster ID
            for (const auto& id : unique_ids) {
                double length = get_length(live_clusters.at(i) , b2groupid, id);
                cluster_lengths[id] = length;
            }
            
            // Find cluster with maximum length
            double max_length = 0;
            int longest_cluster_id = -1;
            for (const auto& pair : cluster_lengths) {
                if (pair.second > max_length) {
                    max_length = pair.second;
                    longest_cluster_id = pair.first;
                }
            }
           // std::cout << "Test: " << longest_cluster_id << " " << max_length << std::endl;
            
            // Mark longest cluster as main
            for (auto& id : b2groupid) {
                if (id == longest_cluster_id) {
                    id = -1;
                }
            }
        }

        live_clusters.at(i)->put_pcarray(b2groupid, "isolated", "perblob");

    }

    // set cluster id ... 
    int cluster_id = 1;
    for (auto* cluster : live_grouping.children()) {
        cluster->set_cluster_id(cluster_id++);
    }
}