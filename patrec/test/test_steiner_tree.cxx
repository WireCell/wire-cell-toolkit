#include "sample_graph.h"

#include <iostream>
#include <vector>

const char* name[] = {"A", "B", "C", "D", "E",
                     "F", "G", "H", "I", "J",
                     "K", "L"};

int main() {
    using SGM = sample_graphs_metrics;

    //
    

    auto metric = SGM::get_graph_metric_steiner();
    auto g = SGM::get_graph_steiner();

    boost::dynamic_properties dp;
    dp.property("node_id", boost::get(boost::vertex_index,g));
    dp.property("color", boost::get(boost::vertex_color,g));
    dp.property("weight", boost::get(boost::edge_weight,g));

    //   dp.property("id", name);
    // dp.property("node_id", boost::make_label_writer(name));
    //boost::make_label_writer(boost::get(boost::vertex_color,g)),
    
    boost::write_graphviz(std::cout, g, boost::make_label_writer(name),
			  boost::make_label_writer(boost::get(boost::edge_weight, g))
			  );
    
    
    std::vector<int> terminals = {SGM::A, SGM::B, SGM::C, SGM::D};
    std::vector<int> nonterminals = {SGM::E, SGM::F};
    std::vector<int> selected_nonterminals;
    
    // solve it
    paal::ir::steiner_tree_iterative_rounding(metric, terminals,
            nonterminals, std::back_inserter(selected_nonterminals));
    // print result
    std::cout << "Selected vertices:" << std::endl;
    for (auto v : terminals) {
        std::cout << v << std::endl;
    }
    for (auto v : selected_nonterminals) {
      std::cout << v << " N" <<  std::endl;
    }
    auto cost = paal::ir::steiner_utils::count_cost(selected_nonterminals,
            terminals, metric);
    std::cout << "Cost of the solution: " << cost << std::endl;

    auto all_elements=boost::range::join(terminals,selected_nonterminals);
    using Vertex = typename paal::data_structures::metric_traits<sample_graphs_metrics::GraphMT>::VertexType;
    paal::data_structures::bimap<Vertex> idx;
    auto g_sub = paal::data_structures::metric_to_bgl_with_index(metric,
							     all_elements, idx);
    std::vector<std::size_t> pm(all_elements.size());
    boost::prim_minimum_spanning_tree(g_sub, &pm[0]);
    auto idx_m = paal::data_structures::make_metric_on_idx(metric, idx);
    for (size_t i=0; i!=pm.size(); i++) {
      std::cout << i << " " << pm[i] << " " << idx_m(i, pm[i]) << std::endl;
    }
    
      
    paal::lp::glp::free_env();
}
