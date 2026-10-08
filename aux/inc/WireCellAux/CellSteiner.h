/** A Steiner-style connectivity repair of a kept cell set (CascadeDeghosting, wcfm doc 14).

    Port of wcp-porting-img wcfm/scripts/d11_repair.py repair_v, variant V2 (doc 11 sec 2, chosen on dev),
    which extends d10_connectivity.py repair (doc 10 sec 2.3):

      1. terminal filter: connected components of the kept set over the edges; a component is a terminal if
         its max P >= p_term or its summed qhat (e) >= q_floor; the other ("weak") components leave the keep
         set but stay crossable candidates;
      2. node cost c = clip(-log P, 0, cost_max); edge weight w_ij = (c_i + c_j)/2 + 1e-6;
      3. one multi-source Dijkstra from every kept cell, limited to the budget (distance, predecessor, source);
      4. for every edge whose ends reach different kept fragments: bridge cost d_i + w_ij + d_j <= budget;
         the cheapest bridge per fragment pair;
      5. Kruskal over the fragment graph in ascending cost; each accepted bridge adds the cells on the two
         predecessor walks back to kept cells.

    It is the Voronoi-bridge rule of clus SteinerGrapher (create_enhanced_steiner_graph) on cells instead of
    points, with the Kruskal step of the prototype's recover_steiner_graph.  Ties are broken by cell index
    (scipy's Dijkstra may break them differently).
 */
#ifndef WIRECELLAUX_CELLSTEINER
#define WIRECELLAUX_CELLSTEINER

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace WireCell::Aux::Cascade {

    struct SteinerParams {
        double p_term{0.8};
        double q_floor{1e4};       // electrons
        double budget{5.5};
        double cost_max{10.0};
        double q_unit{1e4};        // qhat is in these units (Q0)
    };

    struct SteinerResult {
        std::vector<bool> keep;    // the repaired keep set
        size_t nweak_fragments{0}, nweak_cells{0};
        size_t nbridges{0}, nadded{0};
    };

    /// edges: undirected pairs (each once) over the cells, e.g. bb + bb_in of the final level.
    SteinerResult steiner_repair(size_t n, const std::vector<std::array<int64_t, 2>>& edges,
                                 const std::vector<float>& logit, const std::vector<float>& qhat,
                                 const std::vector<bool>& keep, const SteinerParams& par);

    /// Connected-component labels of the cells in `mask` over the edges (-1 elsewhere), labelled in order of
    /// the lowest cell index.
    std::vector<int> components(size_t n, const std::vector<std::array<int64_t, 2>>& edges, const std::vector<bool>& mask);

}  // namespace WireCell::Aux::Cascade

#endif
