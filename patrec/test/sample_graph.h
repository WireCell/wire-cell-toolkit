//=======================================================================
// Copyright (c)
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
/**
 * @file sample_graph.hpp
 * @brief
 * @author Piotr Wygocki, Piotr Godlewski
 * @version 1.0
 * @date 2013-08-04
 */
#ifndef PAAL_SAMPLE_GRAPH_HPP
#define PAAL_SAMPLE_GRAPH_HPP

#include "WireCellPatRec/paal.h"

struct sample_graphs_metrics {
    using EdgeProp = boost::property<boost::edge_weight_t, int>;
    using Graph = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::undirectedS,
        boost::property<boost::vertex_color_t, int>, EdgeProp>;
    using GraphWithoutEdgeWeight = boost::adjacency_list<
        boost::vecS, boost::vecS, boost::undirectedS,
        boost::property<boost::vertex_color_t, int>, boost::no_property>;
    using ListGraph = boost::adjacency_list<
        boost::listS, boost::listS, boost::undirectedS,
        boost::property<boost::vertex_color_t, int>, EdgeProp>;
    using Edge = std::pair<int, int>;
    using GraphMT = paal::data_structures::graph_metric<Graph, int>;
    using Terminals = std::vector<int>;

    enum nodes { A, B, C, D, E, F, G, H, I, J, K, L };

    // graph small
    static Graph get_graph_small() {
        const int num_nodes = 5;
        Edge edge_array[] = { Edge(A, C), Edge(B, B), Edge(B, D),
                              Edge(B, E), Edge(C, B), Edge(C, D),
                              Edge(D, E), Edge(E, A), Edge(E, B) };
        int weights[] = { 1, 2, 1, 2, 7, 3, 1, 1, 1 };
        int num_arcs = sizeof(edge_array) / sizeof(Edge);

        Graph g(edge_array, edge_array + num_arcs, weights, num_nodes);

        return g;
    }

    static Graph get_graph_medium() {
        const int num_nodes = 9;
        Edge edge_array[] = { Edge(A, B), Edge(B, C), Edge(B, D), Edge(B, E),
                              Edge(D, E), Edge(E, F), Edge(F, G), Edge(G, H),
                              Edge(E, H), Edge(H, I) };

        int weights[] = {64, 23, 54, 63, 25, 49, 32, 15, 74, 31};
        int num_arcs = sizeof(edge_array) / sizeof(Edge);

        Graph g(edge_array, edge_array + num_arcs, weights, num_nodes);

        return g;
    }

    static ListGraph get_list_graph_medium() {
        const int num_nodes = 9;
        Edge edge_array[] = { Edge(A, B), Edge(B, C), Edge(B, D), Edge(B, E),
                              Edge(D, E), Edge(E, F), Edge(F, G), Edge(G, H),
                              Edge(E, H), Edge(H, I) };

        int weights[] = {64, 23, 54, 63, 25, 49, 32, 15, 74, 31};
        int num_arcs = sizeof(edge_array) / sizeof(Edge);

        ListGraph g(edge_array, edge_array + num_arcs, weights, num_nodes);

        return g;
    }

    static Graph get_star_medium() {
        const int num_nodes = 11;
        Edge edge_array[] = { Edge(A, B), Edge(A, C), Edge(A, D), Edge(A, E),
                              Edge(A, F), Edge(A, G), Edge(A, H), Edge(A, I),
                              Edge(A, J), Edge(A, K) };
        int weights[] = { 64, 23, 54, 63, 25, 49, 32, 15, 74, 31};
        int num_arcs = sizeof(edge_array) / sizeof(Edge);

        Graph g(edge_array, edge_array + num_arcs, weights, num_nodes);

        return g;
    }

    static Graph get_star_random(int seed, int num_nodes, int minW, int maxW) {
        std::srand(seed);

        int center = std::rand() % num_nodes;

        std::vector< std::pair<int, int> > edges;
        std::vector<int> weights(num_nodes-1);

        for (int v: paal::irange(num_nodes)) if (v != center) edges.push_back(std::make_pair(center, v));
        for (int i: paal::irange(num_nodes-1)) weights[i] = rand() % (maxW - minW + 1) + minW;

        Graph g(edges.begin(), edges.end(), weights.begin(), num_nodes);

        return g;
    }

    static GraphMT get_graph_metric_small() {
        return GraphMT(get_graph_small());
    }

    // graph steiner
    static Graph get_graph_steiner() {
        const int num_nodes = 6;
        Edge edge_array[] = { Edge(A, B), Edge(B, C), Edge(C, D),
                              Edge(D, A), Edge(A, E), Edge(B, E),
                              Edge(C, E), Edge(D, E), Edge(A, F) };
        int weights[] = { 2, 2, 2, 2, 1, 1, 1, 1, 1 };
        int num_arcs = sizeof(edge_array) / sizeof(Edge);

        Graph g(edge_array, edge_array + num_arcs, weights, num_nodes);
        auto color = get(boost::vertex_color, g);
        put(color, A, 1);
        put(color, B, 1);
        put(color, C, 1);
        put(color, D, 1);

        return g;
    }

    static Graph two_points_steiner() {
        int const  num_nodes = 2;
        Edge * edge_array = nullptr;
        int * weights = nullptr;
        Graph g(edge_array, edge_array, weights, num_nodes);
        auto color = get(boost::vertex_color, g);
        put(color, A, 1);
        put(color, B, 1);
        return g;
    }

    static Graph get_graph_steiner_multi_edges() {
        const int num_nodes = 3;
        Edge edge_array[] = { Edge(A, B), Edge(A, C), Edge(B, C), Edge(B, C),
                              Edge(B, C) };
        int weights[] = { 20, 30, 3, 5, 9 };
        int num_arcs = sizeof(edge_array) / sizeof(Edge);

        Graph g(edge_array, edge_array + num_arcs, weights, num_nodes);
        auto color = get(boost::vertex_color, g);
        put(color, A, 1);
        put(color, B, 1);
        put(color, C, 1);

        return g;
    }
    static GraphWithoutEdgeWeight get_graph_steiner_edge() {

        const int num_nodes = 2;
        Edge edge_array[] = { Edge(A, B) };
        int num_arcs = sizeof(edge_array) / sizeof(Edge);

        GraphWithoutEdgeWeight g(edge_array, edge_array + num_arcs, num_nodes);
        auto color = get(boost::vertex_color, g);
        put(color, A, 1);
        put(color, B, 1);

        return g;
    }

    static Graph get_graph_stainer_tree_cycle() {
        const int num_nodes = 5;
        Edge edge_array[] = { Edge(A, B), Edge(A, C), Edge(B, C), Edge(C, D),
                              Edge(C, E) };
        int weights[] = { 4, 2, 2, 2, 2 };
        int num_arcs = sizeof(edge_array) / sizeof(Edge);

        Graph g(edge_array, edge_array + num_arcs, weights, num_nodes);
        auto color = get(boost::vertex_color, g);
        put(color, A, 1);
        put(color, B, 1);
        put(color, D, 1);
        put(color, E, 1);

        return g;
    }

    static GraphMT get_graph_metric_steiner() {
        return GraphMT(get_graph_steiner());
    }

    static std::pair<Terminals, Terminals> get_graph_steiner_vertices() {
        Terminals terminals = { A, B, C, D }, non_terminals = { E, F };
        return std::make_pair(terminals, non_terminals);
    }

    // graph medium
    static GraphMT get_graph_metric_medium() {
        const int num_nodes = 8;
        Edge edge_array[] = { Edge(A, C), Edge(A, F), Edge(B, E), Edge(B, G),
                              Edge(B, H), Edge(C, E), Edge(C, G), Edge(C, H),
                              Edge(D, E), Edge(D, G), Edge(D, H), Edge(F, H) };
        int weights[] = { 4, 2, 1, 2, 7, 3, 1, 8, 1, 3, 4, 10 };
        int num_arcs = sizeof(edge_array) / sizeof(Edge);

        Graph g(edge_array, edge_array + num_arcs, weights, num_nodes);

        return GraphMT(g);
    }

    // graph steiner bigger
    static Graph get_graph_steiner_bigger(int p = 3, int q = 2) {
        bool b;
        int n = p + p * q;
        Graph g(n);
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                int cost = 3;
                if (i < p && j < p) {
                    cost = 1;
                } else if ((j - p) / q == i) {
                    cost = 2;
                }
                b = add_edge(i, j, EdgeProp(cost), g).second;
                assert(b);
            }
        }
        return g;
    }

    static GraphMT get_graph_metric_steiner_bigger() {
        return GraphMT(get_graph_steiner_bigger());
    }

    static std::pair<Terminals, Terminals>
    get_graph_steiner_bigger_vertices(int p = 3, int q = 2) {
        int n = p + p * q;
        Terminals terminals, non_terminals;
        for (int i = 0; i < n; i++) {
            if (i >= p)
                terminals.push_back(i);
            else
                non_terminals.push_back(i);
        }
        return make_pair(terminals, non_terminals);
    }

    // eucildean steiner
    template <typename Points = std::vector<std::pair<int, int>>>
    static std::tuple<paal::data_structures::euclidean_metric<int>, Points,
                      Points>
    get_euclidean_steiner_sample() {
        return std::make_tuple(paal::data_structures::euclidean_metric<int>{},
                               Points{ { 0, 0 }, { 0, 2 }, { 2, 0 }, { 2, 2 } },
                               Points{ { 1, 1 } });
    }
};

#endif // PAAL_SAMPLE_GRAPH_HPP
