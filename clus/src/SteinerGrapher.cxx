#include "WireCellUtil/Exceptions.h"
#include "SteinerGrapher.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;

Steiner::Grapher::Grapher(Cluster& cluster, const Steiner::Grapher::Config& cfg)
    : m_cluster(cluster), m_config(cfg)
{

}


Steiner::Grapher::graph_type Steiner::Grapher::create_steiner_graph()
{
    raise<LogicError>("not implemented yet");

    // Xin, the start of WCP's version is to call Improve_Cluster_2 which
    // returns a new PR3DCluster.  I think should start out to separate
    // "creation" and "improvement".  To "create" I think we need some new
    // support in NaryFacade and Facade_*.h to provide a "clone" method.
    // I think this has implications for calc_sampling_points().


    Steiner::Grapher::graph_type graph_steiner;  // our purpose.


    // Finally, return the results
    return graph_steiner;       // no std::move
}


void Steiner::Grapher::calc_sampling_points(/*, ...*/)
{
    raise<LogicError>("not implemented yet");

    // Xin, this either needs access to IBlobs in order to sample from scratch
    // or it needs to take points from previously sampled clusters.
}

Steiner::Grapher::graph_vertex_set Steiner::Grapher::find_peak_point_indices(bool disable_dead_mix_cell)
{
    Steiner::Grapher::graph_vertex_set peak_points;
    raise<LogicError>("not implemented yet");
    return peak_points;
}

Steiner::Grapher::blob_vertex_map Steiner::Grapher::form_cell_points_map()
{
    Steiner::Grapher::blob_vertex_map cell_points;
    raise<LogicError>("not implemented yet");
    return cell_points;
}

Steiner::Grapher::graph_vertex_set Steiner::Grapher::find_steiner_terminals(bool disable_dead_mix_cell)
{
    Steiner::Grapher::graph_vertex_set steiner_terminals;
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

