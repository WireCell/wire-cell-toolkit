#include "WireCellAux/ClusterShadow.h"
#include "WireCellUtil/GraphTools.h"

#include "WireCellUtil/Graph.h"

using namespace WireCell;
using namespace WireCell::GraphTools;
using namespace WireCell::Aux;
using namespace WireCell::Aux::ClusterShadow;


namespace {
    // Shared "vertex" phase: find geometric clusters and summarize each as a
    // cs_graph node.
    ClusterShadow::graph_t cluster_vertices(const cluster_graph_t& cgraph,
                                            ClusterShadow::blob_cluster_map_t& clusters)
    {
        using coverage_t = ClusterShadow::coverage_t;
        using Node = ClusterShadow::Node;

        // First get just the b-b edged subgraph from cgraph.
        using cvertex_t = typename boost::graph_traits<cluster_graph_t>::vertex_descriptor;
        using Filtered = typename boost::filtered_graph<cluster_graph_t, boost::keep_all,
                                                        std::function<bool(cvertex_t)> >;
        Filtered bcg(cgraph, {}, [&](auto vtx) { return cgraph[vtx].code() == 'b'; });

        // Find the clusters
        // std::unordered_map<cvdesc_t, int> desc2id;
        // auto nclust = boost::connected_components(bcg, boost::make_assoc_property_map(desc2id));
        auto nclust = boost::connected_components(bcg, boost::make_assoc_property_map(clusters));

        // Add the cluster vertices.
        ClusterShadow::graph_t cs_graph(nclust);

        for (auto& [bs_vtx, cs_vtx] : clusters) { // each blob with its geom cluster ID
            auto& cs_node = cs_graph[cs_vtx];
            ++cs_node.nblobs;

            const auto& iblob = std::get<cluster_node_t::blob_t>(cgraph[bs_vtx].ptr);
            cs_node.value += Node::value_t(iblob->value(), iblob->uncertainty());

            const auto& islice = iblob->slice();
            coverage_t::xinterval_t xi(islice->start(), islice->start()+islice->span());

            // Each view/strip of the blob shape.
            for (const auto& strip : iblob->shape().strips()) {

	  // WCT convention ...
	  if (strip.layer < 2) {
	    // skip layers defining anode sensitive area
	    continue;
	  }
	  const size_t view = strip.layer - 2;

                const auto b = strip.bounds;
                coverage_t::yinterval_t yi(b.first, b.second);
                if (cs_node.coverage.size() <= view) {
                    cs_node.coverage.resize(view+1);
                }
                cs_node.coverage[view].add(xi, yi, bs_vtx);
            }
        }
        return cs_graph;
    }

    // Shared "edge" phase: add the cs edge for one blob shadow, dedup'ed by
    // (cluster pair, layer).
    void add_cluster_edge(ClusterShadow::graph_t& cs_graph,
                          ClusterShadow::vdesc_t cs_tail,
                          ClusterShadow::vdesc_t cs_head,
                          const WirePlaneId& wpid)
    {
        // check if cs edge exists for this layer
        bool c_layer_edge_exist = false;
        for (const auto& cedge : mir(boost::edge_range(cs_tail, cs_head, cs_graph))) {
            const auto& eobj = cs_graph[cedge];
            if (eobj.wpid.layer() == wpid.layer() ) {
                c_layer_edge_exist = true;
            }
        }
        if (c_layer_edge_exist) {
            return;
        }

        auto [cs_edge, added] = boost::add_edge(cs_tail, cs_head, cs_graph);
        if (added) {
            cs_graph[cs_edge].wpid = wpid;
        }
    }
}

ClusterShadow::graph_t ClusterShadow::shadow(const cluster_graph_t& cgraph,
                                             const BlobShadow::graph_t& bs_graph,
                                             ClusterShadow::blob_cluster_map_t& clusters)
{
    auto cs_graph = cluster_vertices(cgraph, clusters);
    cs_graph[boost::graph_bundle].stype = bs_graph[boost::graph_bundle].stype;

    // Add edges based on mapping b's from blob shadow to b's from
    // clusters.
    for (const auto& bs_edge : mir(boost::edges(bs_graph))) {
        auto bs_tail = boost::source(bs_edge, bs_graph);
        auto bs_head = boost::target(bs_edge, bs_graph);
        auto c_tail = bs_graph[bs_tail].desc;
        auto c_head = bs_graph[bs_head].desc;
        add_cluster_edge(cs_graph, clusters[c_tail], clusters[c_head],
                         bs_graph[bs_edge].wpid);
    }

    return cs_graph;
}

ClusterShadow::graph_t ClusterShadow::shadow(const cluster_graph_t& cgraph,
                                             const BlobShadow::Shadows& shadows,
                                             ClusterShadow::blob_cluster_map_t& clusters)
{
    auto cs_graph = cluster_vertices(cgraph, clusters);
    cs_graph[boost::graph_bundle].stype = shadows.stype;

    // Add edges based on mapping b's from blob shadows to b's from
    // clusters.  shadows.edges is in the same order as edges() of the
    // equivalent blob shadow graph.
    for (const auto& er : shadows.edges) {
        auto c_tail = shadows.nodes[er.v1].desc;
        auto c_head = shadows.nodes[er.v2].desc;
        add_cluster_edge(cs_graph, clusters[c_tail], clusters[c_head], er.prop.wpid);
    }

    return cs_graph;
}

ClusterShadow::graph_t ClusterShadow::shadow(const cluster_graph_t& cgraph,
                                             const BlobShadow::graph_t& bs_graph)
{
    ClusterShadow::blob_cluster_map_t clusters;
    return ClusterShadow::shadow(cgraph, bs_graph, clusters);
}
