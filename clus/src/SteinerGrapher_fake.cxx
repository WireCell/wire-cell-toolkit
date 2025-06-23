#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/GraphTools.h"
#include "SteinerGrapher.h"
#include "PAAL.h"               // eventually, subsume int GraphAlgorithms! 

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::Facade;
using WireCell::GraphTools::vertex_range;


static
Steiner::Grapher::vertex_set fake_find_steiner_vertices(const Steiner::Grapher::graph_type& graph)
{
    Steiner::Grapher::vertex_set ret;
    for (auto vtx : vertex_range(graph)) {
        ret.insert(vtx);
    }
    return ret;    
}

Steiner::Grapher::graph_type Steiner::Grapher::fake_steiner_graph()
{
    log->warn("running fake_steiner_graph()");

    /// This is not the real algorithm yet, just to show some ways of doing
    /// various things.

    /// Here is an example of filtering a graph.
    {
        // We first make one giving its number of vertices.  This gets stored by
        // the name on the cluster for later recall.  We can also get a
        // reference immediately from return..
        graph_type& dummy = m_cluster.make_graph("dummy", 3);
        assert(boost::num_vertices(dummy) == 3);

        // Add two edges.  Here we can use literal vertex descriptors which are
        // coincident with being simple integer (size_t) indices.
        auto [e0,ok0] = boost::add_edge(0, 1, dummy);
        auto [e1,ok1] = boost::add_edge(0, 2, dummy);
        assert (ok0 && ok1);    
        // We don't usually need to catch both edge and the "added ok" bool,
        // just showing how to do it.

        // Graph is now: 2--0--1.
        log->debug("dummy graph with {} vertices and {} edges",
                   boost::num_vertices(dummy), boost::num_edges(dummy));

        // The graph is also available for use with a graph algorithms class.
        auto& dummy_ga = m_cluster.graph_algorithms("dummy");

        // Form a subset of edges
        edge_set to_filter = {e0};

        // We can get a filtered graph with a subset of edges:
        auto filtered_with_e0 = dummy_ga.reduce(to_filter);
        log->debug("filtered with e0 graph with {} vertices and {} edges",
                   boost::num_vertices(filtered_with_e0), boost::num_edges(filtered_with_e0));


        // Or which excludes the set:
        auto filtered_without_e0 = dummy_ga.reduce(to_filter, false);
        log->debug("filtered without e0 graph with {} vertices and {} edges",
                   boost::num_vertices(filtered_without_e0), boost::num_edges(filtered_without_e0));


        // The filtered graph relies on the underlying graph to be kept around.
        // It also keeps the vertex descriptors / indices of the original graph.
        // You can make a "real" graph with new, contiguously numbered vertices.
        graph_type real_with_e0;
        boost::copy_graph(filtered_with_e0, real_with_e0);
    }

    // Exercise the PAAL excerpt to assure it compiles.  The "real" usage should
    // be via GraphAlgorithms as an expensive dijkstra's shortest paths is
    // needed and so results should be cached.  But, I start here as I don't
    // really know how it is to be used.
    {
        // For this example, we make the graph on the stack instead of holding
        // it on the cluster.  We can always "give" the graph to the cluster
        // later.
        const size_t npoints = 3;
        graph_type graph(npoints);
        boost::add_edge(0, 1, graph);
        boost::add_edge(0, 2, graph);

        // This is "ported" from void PR3DCluster::recover_steiner_graph() with
        // very little understanding of what it is doing.  Again, just to check
        // compilation.

        // It is not clear that this is needed since vertex descriptor == vertex index.
        auto index = get(boost::vertex_index, graph);

        std::vector<vertex_type> terminals = {1,2};
        std::vector<vertex_type> nearest_terminal(npoints);
        auto nearest_terminal_map = boost::make_iterator_property_map(nearest_terminal.begin(), index);
        for (auto terminal : terminals) {
            nearest_terminal_map[terminal] = terminal;
        }

        /// This is a highly suspect line of code from WCP:
        // auto edge_weight = boost::choose_pmap(boost::get_param(boost::no_named_parameters(), boost::edge_weight), graph, boost::edge_weight);
        // Let's try something a bit less baroque.
        auto edge_weight = get(boost::edge_weight, graph);

        std::vector<edge_weight_type> distance(npoints);
        auto distance_map = boost::make_iterator_property_map(distance.begin(), index);

        std::vector<edge_type> vpred(npoints);
        auto last_edge = boost::make_iterator_property_map(vpred.begin(), index);

        boost::dijkstra_shortest_paths(
            graph, terminals.begin(), terminals.end(), boost::dummy_property_map(),
            distance_map, edge_weight, index, PAAL::less(),
            boost::closed_plus<edge_weight_type>(),
            std::numeric_limits<edge_weight_type>::max(), 0,
            boost::make_dijkstra_visitor(
                PAAL::make_nearest_recorder(
                    nearest_terminal_map, last_edge, boost::on_edge_relaxed{})));
    }


    /// Here is how you'd get an existing graph held by the cluster.
    const std::string flavor = "basic"; // or "ctpc" or "relaxed".
    const auto& graph = get_graph(flavor);

    // And a starting point cloud.
    const auto& point_cloud = get_point_cloud(flavor);

    // We can also get graph algorithms on the graph.
    auto ga = m_cluster.graph_algorithms(flavor);

    // Now, assume you select a subset of vertices in some way.
    vertex_set selected_vertices = fake_find_steiner_vertices(graph);

    // We can get a filtered graph with that subset.
    auto filtered_steiner = ga.reduce(selected_vertices);

    // The vertex descriptors of this filtered graph can still be used to
    // identify points in the point_cloud.  We can also reduce both.  First the
    // graph:
    graph_type graph_steiner;
    boost::copy_graph(filtered_steiner, graph_steiner);
    
    // Then the PC.  This needs to convert from set to vector.
    std::vector<size_t> selected_indices(selected_vertices.begin(), selected_vertices.end());
    auto point_cloud_steiner = point_cloud.subset(selected_indices);


    // Finally, return the results
    return graph_steiner;       // no std::move
}

