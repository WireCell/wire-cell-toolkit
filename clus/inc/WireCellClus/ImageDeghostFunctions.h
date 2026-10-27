/** Free function behind the ClusteringImageDeghost visitor (pr_ml doc 28,
    "combination 3": traditional imaging, clustering and Q/L, then the
    deghosting cascade's cells replace the blobs of every live cluster).

    The cells of one cascade cluster are T0-corrected together (the correction
    is applied while they still sit in their cascade cluster, doc pdvd/125), so
    all cells taken out of one cascade cluster must go to live clusters with the
    same cluster_t0.  ClusteringBeamDeghost had one bundle = one t0; with every
    cluster in scope the owners can differ.  This function decides, per cascade
    cluster, the t0 its taken cells get and which cells must stay behind; plain
    integers and doubles, so it is unit-tested without a Grouping
    (clus/test/doctest_image_deghost.cxx).

    Deterministic: cascade clusters and owners are indexed by position, ties
    break on the smaller owner key.
 */
#ifndef WIRECELL_CLUS_IMAGEDEGHOSTFUNCTIONS
#define WIRECELL_CLUS_IMAGEDEGHOSTFUNCTIONS

#include <vector>

namespace WireCell::Clus::PR {

    struct ImageDeghostT0 {
        bool touched{false};   ///< at least one cell of the cascade cluster is taken
        double t0{0};          ///< the cluster_t0 the cascade cluster borrows
        bool conflict{false};  ///< its takers carry more than one t0 (beyond tolerance)
    };

    /// cell_dgc[i]: the cascade cluster (0 .. n_dgc-1) of cell i.  owner_of[i]:
    /// the live-cluster key taking cell i, -1 = none.  selected[k] != 0: live
    /// cluster k is replaced.  owner_t0[k]: its cluster_t0.  Per cascade
    /// cluster: touched when a selected owner takes one of its cells; t0 = the
    /// t0 of the selected owner taking most of its cells (tie: smaller key);
    /// conflict when another taker's t0 differs by more than tol.  take[i] (out)
    /// = 1 when cell i goes to its owner: the owner is selected and its t0 is
    /// within tol of the cascade cluster's t0.
    std::vector<ImageDeghostT0> image_deghost_cluster_t0(const std::vector<int>& cell_dgc,
                                                         const std::vector<int>& owner_of,
                                                         const std::vector<char>& selected,
                                                         const std::vector<double>& owner_t0,
                                                         int n_dgc, double tol,
                                                         std::vector<char>& take);

}  // namespace WireCell::Clus::PR

#endif
