// pr_ml doc 44 sec 8 (3e4 FD-HD nc 20022 e14, SDCC 5157.7192): connect_graph_ctpc_with_reference sizes its five
// pair tables and its two MST graphs by the components that kept points after the reference filter, not by every
// component.  The emptied components only ever were isolated vertices of the MST graphs.  What makes that change
// output-neutral is that process_mst_deterministically gives the same tree on the graph without those isolated
// vertices, relabelled in the same increasing order.  These tests pin exactly that, on random graphs with many tied
// weights (ties are where a relabelling could leak into Prim's choice).
#include "WireCellUtil/doctest.h"
#include "WireCellClus/Facade_Util.h"
#include "WireCellClus/Graphs.h"

#include <random>
#include <tuple>
#include <vector>

using namespace WireCell::Clus;
using Table = std::vector<std::vector<std::tuple<int, int, double>>>;

namespace {
    const std::tuple<int, int, double> invalid{-1, -1, 1e9};

    Table make_table(size_t n) { return Table(n, std::vector<std::tuple<int, int, double>>(n, invalid)); }

    // the MST of the valid upper-triangle entries of `dis`, edges added in j < k order (as the caller does)
    Table mst_of(const Table& dis)
    {
        const size_t n = dis.size();
        Graphs::Weighted::Graph g(n);
        for (size_t j = 0; j != n; ++j)
            for (size_t k = j + 1; k != n; ++k)
                if (std::get<0>(dis[j][k]) >= 0 && !boost::edge(j, k, g).second) add_edge(j, k, std::get<2>(dis[j][k]), g);
        Table copy = dis, mst = make_table(n);
        Facade::process_mst_deterministically(g, copy, mst);
        return mst;
    }
}

TEST_CASE("ctpc_with_reference compaction: MST on live components == MST on all components (random, tied weights)")
{
    size_t n_ties = 0, n_dead = 0, n_mst_edges = 0;
    for (unsigned seed = 0; seed != 300; ++seed) {
        std::mt19937 rng(seed);
        const size_t num = 2 + rng() % 60;
        const double p_dead = (rng() % 80) / 100.0, p_edge = 0.05 + (rng() % 60) / 100.0;
        std::vector<char> dead(num);
        std::vector<size_t> live;
        for (size_t c = 0; c != num; ++c) {
            dead[c] = std::uniform_real_distribution<double>(0, 1)(rng) < p_dead;
            if (!dead[c]) live.push_back(c);
            n_dead += dead[c];
        }
        Table full = make_table(num);
        for (size_t j = 0; j != num; ++j) {
            for (size_t k = j + 1; k != num; ++k) {
                if (dead[j] || dead[k]) continue;   // the caller skips these pairs: entries stay invalid
                if (std::uniform_real_distribution<double>(0, 1)(rng) >= p_edge) continue;
                const double w = 1.0 + rng() % 4;   // four weights only: many ties
                full[j][k] = std::make_tuple(int(rng() % 50), int(rng() % 50), w);
            }
        }
        const size_t nl = live.size();
        Table comp = make_table(nl);
        for (size_t a = 0; a != nl; ++a)
            for (size_t b = 0; b != nl; ++b) comp[a][b] = full[live[a]][live[b]];

        const Table mst_full = mst_of(full), mst_comp = mst_of(comp);

        for (size_t j = 0; j != num; ++j)
            for (size_t k = 0; k != num; ++k)
                if (dead[j] || dead[k]) CHECK(mst_full[j][k] == invalid);
        for (size_t a = 0; a != nl; ++a) {
            for (size_t b = 0; b != nl; ++b) {
                CHECK(mst_comp[a][b] == mst_full[live[a]][live[b]]);
                n_mst_edges += std::get<0>(mst_comp[a][b]) >= 0;
            }
        }
        // count tied weights among the valid entries, to show the test exercises ties
        std::vector<size_t> per_w(5, 0);
        for (size_t j = 0; j != num; ++j)
            for (size_t k = j + 1; k != num; ++k)
                if (std::get<0>(full[j][k]) >= 0) ++per_w[size_t(std::get<2>(full[j][k]))];
        for (size_t w = 1; w != 5; ++w) n_ties += per_w[w] > 1;
    }
    CHECK(n_dead > 1000);
    CHECK(n_ties > 300);
    CHECK(n_mst_edges > 1000);
}

TEST_CASE("ctpc_with_reference compaction: the relabelled graph is not a no-op (a dead vertex between live ones)")
{
    // 0 live, 1 dead, 2 live, 3 live; edges 0-2 (w 1), 2-3 (w 1), 0-3 (w 1): a three-way tie
    Table full = make_table(4);
    full[0][2] = std::make_tuple(1, 2, 1.0);
    full[2][3] = std::make_tuple(3, 4, 1.0);
    full[0][3] = std::make_tuple(5, 6, 1.0);
    const std::vector<size_t> live{0, 2, 3};
    Table comp = make_table(3);
    for (size_t a = 0; a != 3; ++a)
        for (size_t b = 0; b != 3; ++b) comp[a][b] = full[live[a]][live[b]];
    const Table mf = mst_of(full), mc = mst_of(comp);
    size_t nedges = 0;
    for (size_t a = 0; a != 3; ++a)
        for (size_t b = 0; b != 3; ++b) {
            CHECK(mc[a][b] == mf[live[a]][live[b]]);
            nedges += std::get<0>(mc[a][b]) >= 0;
        }
    CHECK(nedges == 2);
}
