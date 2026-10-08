// Steiner-style repair of the kept cells (wcfm doc 14).  See WireCellAux/CellSteiner.h.
// Python reference: wcp-porting-img wcfm/scripts/d11_repair.py repair_v L77-132 (variant V2),
// d10_connectivity.py repair L114-158.

#include "WireCellAux/CellSteiner.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <queue>
#include <tuple>

using namespace WireCell::Aux;

namespace {
    struct DSU {
        std::vector<int> p;
        explicit DSU(size_t n) : p(n) { std::iota(p.begin(), p.end(), 0); }
        int find(int x)
        {
            while (p[x] != x) {
                p[x] = p[p[x]];
                x = p[x];
            }
            return x;
        }
    };

    // CSR adjacency, neighbours in ascending order
    struct Adj {
        std::vector<size_t> ptr;
        std::vector<int64_t> nbr;
        std::vector<size_t> eid;
    };
    Adj make_adj(size_t n, const std::vector<std::array<int64_t, 2>>& edges)
    {
        Adj a;
        a.ptr.assign(n + 1, 0);
        for (const auto& e : edges) {
            a.ptr[e[0] + 1] += 1;
            a.ptr[e[1] + 1] += 1;
        }
        for (size_t i = 0; i < n; ++i) a.ptr[i + 1] += a.ptr[i];
        a.nbr.resize(a.ptr[n]);
        a.eid.resize(a.ptr[n]);
        std::vector<size_t> fill(a.ptr.begin(), a.ptr.end() - 1);
        for (size_t k = 0; k < edges.size(); ++k) {
            const auto& e = edges[k];
            a.nbr[fill[e[0]]] = e[1]; a.eid[fill[e[0]]++] = k;
            a.nbr[fill[e[1]]] = e[0]; a.eid[fill[e[1]]++] = k;
        }
        return a;
    }
}  // namespace

std::vector<int> Cascade::components(size_t n, const std::vector<std::array<int64_t, 2>>& edges,
                                     const std::vector<bool>& mask)
{
    DSU d(n);
    for (const auto& e : edges) {
        if (!mask[e[0]] || !mask[e[1]]) continue;
        const int a = d.find((int) e[0]), b = d.find((int) e[1]);
        if (a != b) d.p[std::max(a, b)] = std::min(a, b);
    }
    std::vector<int> lab(n, -1), root_lab(n, -1);
    int nl = 0;
    for (size_t i = 0; i < n; ++i) {
        if (!mask[i]) continue;
        const int r = d.find((int) i);
        if (root_lab[r] < 0) root_lab[r] = nl++;
        lab[i] = root_lab[r];
    }
    return lab;
}

Cascade::SteinerResult Cascade::steiner_repair(size_t n, const std::vector<std::array<int64_t, 2>>& edges,
                                               const std::vector<float>& logit, const std::vector<float>& qhat,
                                               const std::vector<bool>& keep_in, const SteinerParams& par)
{
    SteinerResult res;
    std::vector<bool> K = keep_in;
    std::vector<double> P(n);
    for (size_t i = 0; i < n; ++i) P[i] = 1.0 / (1.0 + std::exp(-(double) logit[i]));

    // 1. terminal filter (V2)
    {
        auto kl = components(n, edges, K);
        int nf = 0;
        for (int l : kl) nf = std::max(nf, l + 1);
        std::vector<double> maxp(nf, 0.0), qs(nf, 0.0);
        for (size_t i = 0; i < n; ++i) {
            if (!K[i]) continue;
            maxp[kl[i]] = std::max(maxp[kl[i]], P[i]);
            qs[kl[i]] += (double) qhat[i] * par.q_unit;
        }
        std::vector<bool> term(nf);
        for (int f = 0; f < nf; ++f) {
            term[f] = (maxp[f] >= par.p_term) || (qs[f] >= par.q_floor);
            if (!term[f]) ++res.nweak_fragments;
        }
        for (size_t i = 0; i < n; ++i) {
            if (K[i] && !term[kl[i]]) {
                K[i] = false;
                ++res.nweak_cells;
            }
        }
    }
    const auto kl = components(n, edges, K);
    std::vector<int64_t> kidx;
    for (size_t i = 0; i < n; ++i) if (K[i]) kidx.push_back((int64_t) i);
    res.keep = K;
    if (kidx.size() < 2) return res;

    // 2. costs
    std::vector<double> c(n);
    for (size_t i = 0; i < n; ++i) {
        c[i] = std::clamp(-std::log(std::clamp(P[i], 1e-12, 1.0)), 0.0, par.cost_max);
    }
    std::vector<double> w(edges.size());
    for (size_t k = 0; k < edges.size(); ++k) w[k] = 0.5 * (c[edges[k][0]] + c[edges[k][1]]) + 1e-6;

    // 3. multi-source Dijkstra, limited to the budget; ties by (distance, cell index)
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double> dist(n, inf);
    std::vector<int64_t> pred(n, -1), src(n, -1);
    std::vector<bool> done(n, false);
    using item = std::pair<double, int64_t>;
    std::priority_queue<item, std::vector<item>, std::greater<item>> pq;
    for (auto s : kidx) {
        dist[s] = 0.0;
        src[s] = s;
        pq.push({0.0, s});
    }
    const auto adj = make_adj(n, edges);
    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (done[u]) continue;
        done[u] = true;
        for (size_t a = adj.ptr[u]; a < adj.ptr[u + 1]; ++a) {
            const int64_t v = adj.nbr[a];
            const double nd = d + w[adj.eid[a]];
            if (nd > par.budget) continue;
            if (nd < dist[v]) {
                dist[v] = nd;
                pred[v] = u;
                src[v] = src[u];
                pq.push({nd, v});
            }
        }
    }

    // 4. bridges: cheapest per fragment pair
    struct Bridge {
        int a, b;
        double cost;
        int64_t i, j;
    };
    std::vector<Bridge> br;
    for (size_t k = 0; k < edges.size(); ++k) {
        const int64_t i = edges[k][0], j = edges[k][1];
        if (!std::isfinite(dist[i]) || !std::isfinite(dist[j])) continue;
        const int fi = kl[src[i]], fj = kl[src[j]];
        if (fi == fj) continue;
        const double cost = dist[i] + w[k] + dist[j];
        if (cost > par.budget) continue;
        br.push_back({std::min(fi, fj), std::max(fi, fj), cost, i, j});
    }
    // np.lexsort((cost, b, a)) then first per (a, b)
    std::stable_sort(br.begin(), br.end(), [](const Bridge& x, const Bridge& y) {
        return std::tie(x.a, x.b, x.cost) < std::tie(y.a, y.b, y.cost);
    });
    std::vector<Bridge> best;
    for (size_t k = 0; k < br.size(); ++k) {
        if (k == 0 || br[k].a != br[k - 1].a || br[k].b != br[k - 1].b) best.push_back(br[k]);
    }
    // 5. Kruskal in ascending cost (stable)
    std::stable_sort(best.begin(), best.end(), [](const Bridge& x, const Bridge& y) { return x.cost < y.cost; });
    int nf = 0;
    for (int l : kl) nf = std::max(nf, l + 1);
    DSU dsu(nf);
    std::vector<bool> out = K;
    for (const auto& b : best) {
        const int ra = dsu.find(b.a), rb = dsu.find(b.b);
        if (ra == rb) continue;
        dsu.p[ra] = rb;
        ++res.nbridges;
        for (int64_t v : {b.i, b.j}) {
            while (v >= 0 && !out[v]) {
                out[v] = true;
                ++res.nadded;
                v = pred[v];
            }
        }
    }
    res.keep = out;
    return res;
}
