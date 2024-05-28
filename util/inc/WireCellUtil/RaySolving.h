/** Solving attempts to invert the matrix equation m = G*b where m is
 * a vector of measurements on channels, b is a vector of blobs and G
 * is a matrix that relates the two.
 */
#ifndef WIRECELL_RAYSOLVING_H
#define WIRECELL_RAYSOLVING_H

#include "WireCellUtil/Graph.h"



#include "WireCellUtil/Eigen.h"

#include <vector>
#include <unordered_map>

namespace WireCell {

    namespace RayGrid {

        class Grouping {
           public:
            Grouping() = default;
            virtual ~Grouping() = default;

            typedef size_t ident_t;
            struct node_t {
                char ntype;     // 'm', 'w', or 's'
                ident_t ident;  // an ident for a 'm' or an 's'
                float value;
                float weight;
            };
            typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, node_t> graph_t;
            typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_t;

            // Add an 'm' or 's' type node.
            virtual void add(char ntype, ident_t chid, std::vector<ident_t> wids, float value, float weight = 1.0);

            typedef std::vector<node_t> cluster_t;
            typedef std::unordered_map<int, cluster_t> clusterset_t;

            // Return nodes in connected subgraphs
            clusterset_t clusters();

            graph_t& graph() { return m_graph; }

           private:
            graph_t m_graph;
            std::unordered_map<ident_t, vertex_t> m_wid2vtx;

            // helper
            vertex_t wire_node(ident_t wid);
        };

        class Solving {
           public:
            Solving() = default;
            ~Solving() = default;

            typedef size_t ident_t;
            struct node_t {
                char ntype;     // 'm' or 's'
                ident_t ident;  // an ident number for a 'm' or an 's'
                float value;
                float weight;
            };
            typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, node_t> graph_t;
            typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_t;

            //
            void add(const Grouping::clusterset_t& cset);

            // Solve for "S", return map from original source ident to solved value.
            typedef std::unordered_map<ident_t, float> solution_t;
            solution_t solve();

            graph_t& graph() { return m_graph; }

           private:
            graph_t m_graph;
            std::unordered_map<ident_t, vertex_t> m_sid2vtx;

            vertex_t measurement_node(float value, float weight);
            vertex_t source_node(ident_t sid, float value, float weight);

            void solve_one(solution_t& answer, const std::vector<vertex_t>& sources,
                           const std::vector<vertex_t>& measures);
        };

    }  // namespace RayGrid

}  // namespace WireCell

#endif /* WIRECELL_RAYSOLVING_H */
