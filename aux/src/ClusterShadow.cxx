#include "WireCellAux/ClusterShadow.h"
#include "WireCellUtil/GraphTools.h"

#include "WireCellUtil/Graph.h"

using namespace WireCell;
using namespace WireCell::GraphTools;
using namespace WireCell::Aux;
using namespace WireCell::Aux::ClusterShadow;


ClusterShadow::graph_t ClusterShadow::shadow(const cluster_graph_t& cgraph,
                                             const BlobShadow::graph_t& bs_graph,
                                             ClusterShadow::blob_cluster_map_t& clusters)
{
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
    cs_graph[boost::graph_bundle].stype = bs_graph[boost::graph_bundle].stype;

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

    // Add edges based on mapping b's from blob shadow to b's from
    // clusters.
    for (const auto& bs_edge : mir(boost::edges(bs_graph))) {
        auto bs_tail = boost::source(bs_edge, bs_graph);
        auto bs_head = boost::target(bs_edge, bs_graph);
        auto c_tail = bs_graph[bs_tail].desc;
        auto c_head = bs_graph[bs_head].desc;
        auto cs_tail = clusters[c_tail];
        auto cs_head = clusters[c_head];

        // check if cs edge exists for this layer
        bool c_layer_edge_exist = false;
        for (const auto& cedge : mir(boost::edge_range(cs_tail, cs_head, cs_graph))) {
            const auto& eobj = cs_graph[cedge];
            if (eobj.wpid.layer() == bs_graph[bs_edge].wpid.layer() ) {
                c_layer_edge_exist = true;
            }
        }
        if (c_layer_edge_exist) {
            continue;
        }

        auto [cs_edge, added] = boost::add_edge(cs_tail, cs_head, cs_graph);
        if (added) {
            cs_graph[cs_edge].wpid = bs_graph[bs_edge].wpid;
        }
    }

    return cs_graph;
}

ClusterShadow::graph_t ClusterShadow::shadow(const cluster_graph_t& cgraph,
                                             const BlobShadow::graph_t& bs_graph)
{
    ClusterShadow::blob_cluster_map_t clusters;
    return ClusterShadow::shadow(cgraph, bs_graph, clusters);
}
