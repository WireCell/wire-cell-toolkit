#include "WireCellUtil/Exceptions.h"
#include "SteinerGrapher.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

Steiner::Grapher::graph_type Steiner::Grapher::create_steiner_graph()
{
    return fake_steiner_graph(); // see SteinerGrapher_face.cxx
}


void Steiner::Grapher::calc_sampling_points(/*, ...*/)
{
    raise<LogicError>("not implemented yet");

    // Xin, this either needs access to IBlobs in order to sample from scratch
    // or it needs to take points from previously sampled clusters.
}

Steiner::Grapher::vertex_set Steiner::Grapher::find_peak_point_indices(bool disable_dead_mix_cell)
{
    Steiner::Grapher::vertex_set peak_points;
    raise<LogicError>("not implemented yet");
    return peak_points;
}


Steiner::Grapher::vertex_set Steiner::Grapher::find_steiner_terminals(bool disable_dead_mix_cell)
{
    Steiner::Grapher::vertex_set steiner_terminals;
    raise<LogicError>("not implemented yet");
    return steiner_terminals;
}

void Steiner::Grapher::establish_same_blob_steiner_edges(graph_type& graph, bool disable_dead_mix_cell, int flag)
{
    raise<LogicError>("not implemented yet");

}

Steiner::Grapher::graph_type Steiner::Grapher::create_steiner_tree(/*what type for point_cloud_steiner?*/)
{
    Steiner::Grapher::graph_type graph_steiner;
    raise<LogicError>("not implemented yet");
    return graph_steiner;
}




Steiner::Grapher::blob_vertex_map Steiner::Grapher::form_cell_points_map()
{
    blob_vertex_map cell_points;
    
    // Get the 3D scoped view with x,y,z coordinates
    const auto& sv = m_cluster.sv3d();   // This is default scoped view with 3D coordinates ... 
    const auto& nodes = sv.nodes();   // These are the blob nodes
    const auto& skd = sv.kd();        // K-d tree with point data
    
    // Check if we have valid data
    if (nodes.empty() || skd.npoints() == 0) {
        return cell_points;  // Return empty map if no data
    }
    
    // The major indices tell you which blob each point belongs to
    const auto& majs = skd.major_indices();  // blob index for each point
    
    // Iterate through all points and assign them to their respective blobs
    for (size_t point_idx = 0; point_idx < skd.npoints(); ++point_idx) {
        size_t blob_idx = majs[point_idx];
        
        // Bounds check
        if (blob_idx >= nodes.size()) {
            continue;  // Skip invalid blob indices
        }
        
        // Get the blob facade from the node
        const auto* blob = nodes[blob_idx]->value.facade<Blob>();
        if (!blob) {
            continue;  // Skip if facade creation failed
        }
        
        // Initialize the set for this blob if it doesn't exist
        if (cell_points.find(blob) == cell_points.end()) {
            cell_points[blob] = vertex_set();
        }
        
        // Add this point index to the blob's set
        cell_points[blob].insert(point_idx);
    }
    
    return cell_points;
}