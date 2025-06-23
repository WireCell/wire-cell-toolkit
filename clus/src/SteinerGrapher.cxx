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

Steiner::Grapher::blob_vertex_map Steiner::Grapher::form_cell_points_map()
{
    Steiner::Grapher::blob_vertex_map cell_points;
    raise<LogicError>("not implemented yet");
    return cell_points;
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




