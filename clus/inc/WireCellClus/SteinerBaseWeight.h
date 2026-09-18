/** Charge-aware pricing of the Steiner BASE graph before the Voronoi step.

    doc pdvd/115.  This is the pure core of the `base_weight_blank_alpha` /
    `base_weight_scope` knobs of CreateSteinerGraph / Steiner::Grapher, split
    out of the (src-private) SteinerGrapher.h so a doctest can exercise it
    without a cluster -- as SteinerBlankPlane.h, SteinerThinning.h and
    CtpcAnisoMetric.h are.

    Background.  create_enhanced_steiner_graph builds the Steiner tree in two
    stages on the cluster's base graph ("ctpc_ref_pid"): a multi-source
    Dijkstra from the terminals (the Voronoi tessellation) and, per terminal
    pair, the cheapest bridging edge; the vertices on the Dijkstra paths that
    the selected bridges walk back along become the tree's INTERIORS.  The
    base graph is priced by geometry alone (every edge weight is the Euclidean
    length of the edge, on both the prototype and the toolkit), so an interior
    is admitted without any charge test: doc 114 measured that the Steiner
    seed's remaining off-image detours run through interiors that see charge
    on two planes or on none (the retiler's painted halo and empty cells),
    which no terminal rule reaches.  The reduced graph's own charge factor
    (0.8 + 0.4 * mean Q0/(Q+Q0), endpoint charges) is applied AFTER the tree
    exists, so it can re-route among the admitted edges but cannot admit the
    on-image vertices the geometric Voronoi left out.

    The lever prices the base graph edge by the number of planes at charge
    exactly 0 at its two endpoints, BEFORE the Voronoi step:

        w' = w * (1 + alpha * 0.5 * (nz(s) + nz(t)))

    with nz(v) in 0..3 the number of planes whose charge value at vertex v is
    exactly 0.  alpha = 0 (the C++ default) leaves every weight unchanged and
    the code path is not taken at all.  A dead plane counts as a zero plane:
    a pricing lever needs no dead-channel exemption to keep gap jumping,
    because inside a dead or empty region every route is priced alike and
    Dijkstra still takes the only one -- the penalty reorders routes only where
    an on-image alternative exists.

    Scope.  "tree" (the default when the knob is on) uses the priced graph for
    the Voronoi and the bridge selection only; the reduced graph's edge weight
    is still the GEOMETRIC length times the production charge factor, so only
    the tree's topology moves and every downstream consumer sees the same
    pricing model.  "tree+path" also carries the priced length into the
    reduced graph, so the STM rough path (a Dijkstra on the reduced graph)
    avoids blank interiors too.
 */

#ifndef WIRECELLCLUS_STEINERBASEWEIGHT
#define WIRECELLCLUS_STEINERBASEWEIGHT

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <string>

namespace WireCell::Clus::Steiner {

    enum class BaseWeightScope { tree = 0, tree_path };

    /// Parse the config string.  Unknown strings return false and leave
    /// `scope` untouched, so a caller can refuse a typo instead of silently
    /// running a different level.
    inline bool parse_base_weight_scope(const std::string& s, BaseWeightScope& scope)
    {
        if (s == "tree") { scope = BaseWeightScope::tree; return true; }
        if (s == "tree+path") { scope = BaseWeightScope::tree_path; return true; }
        return false;
    }

    /// The pricing factor of one edge from the zero-plane counts of its two
    /// endpoints.  alpha <= 0 => 1 exactly.
    inline double blank_weight_factor(int nz_source, int nz_target, double alpha)
    {
        if (alpha <= 0) return 1.0;
        return 1.0 + alpha * 0.5 * (double(nz_source) + double(nz_target));
    }

    /// A priced COPY of `base`: same vertices (same indices), same edges,
    /// each weight multiplied by blank_weight_factor of its endpoints.
    /// `base` is not touched.  With alpha <= 0 the copy has the base weights.
    ///
    /// @param nzero  nzero(v) -> number of planes at charge exactly 0 at
    ///               vertex v (0..3).
    template <class Graph, class NZero>
    Graph reweight_base_graph(const Graph& base, NZero nzero, double alpha)
    {
        Graph out = base;
        if (alpha <= 0) return out;
        auto wmap = boost::get(boost::edge_weight, out);
        for (auto [ei, ee] = boost::edges(out); ei != ee; ++ei) {
            const auto s = boost::source(*ei, out);
            const auto t = boost::target(*ei, out);
            boost::put(wmap, *ei, boost::get(wmap, *ei) * blank_weight_factor(nzero(s), nzero(t), alpha));
        }
        return out;
    }

}  // namespace WireCell::Clus::Steiner

#endif
