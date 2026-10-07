/** Free function behind the ClusteringBeamDeghost visitor (doc pdvd/125).

    The deghosting cascade does not start from the production blobs: it cuts
    the uncut tiling of an anode to 4-strip cells, so a deghosted cell has no
    parent production blob.  The link between the two images is geometric --
    same anode face, overlapping time slices, overlapping wire ranges in all
    three views -- and this function computes it on plain integers so it can be
    unit-tested (clus/test/doctest_beam_deghost.cxx) without a Grouping.

    Deterministic by construction: owners are bucketed by their int wpid and
    visited in a (tick, input index) order, and every tie breaks on the int
    cluster id, never on an address.
 */
#ifndef WIRECELL_CLUS_BEAMDEGHOSTFUNCTIONS
#define WIRECELL_CLUS_BEAMDEGHOSTFUNCTIONS

#include <vector>

namespace WireCell::Clus::PR {

    /// One blob as a box in (face, tick, wire) space.  All ranges half-open
    /// [min, max), the Facade::Blob convention (slice_index_min/max in ticks,
    /// {u,v,w}_wire_index_min/max face-local wire-in-plane).
    struct DeghostBox {
        int wpid{0};          ///< WirePlaneId::ident() of the face (anode + face)
        int t0{0}, t1{0};
        int u0{0}, u1{0}, v0{0}, v1{0}, w0{0}, w1{0};
        int cluster_id{-1};   ///< owners only: the cluster holding this blob
    };

    /// Overlap volume of two boxes, ticks x U wires x V wires x W wires; 0 when
    /// the faces differ or any of the four ranges does not intersect.
    long long beam_deghost_overlap(const DeghostBox& a, const DeghostBox& b);

    /// For each cell the cluster_id of the owner cluster it overlaps most: the
    /// overlap volumes with all blobs of one cluster are summed, the largest
    /// sum wins, a tie goes to the smaller cluster_id.  -1 when the cell
    /// overlaps no owner blob.
    std::vector<int> beam_deghost_assign(const std::vector<DeghostBox>& cells,
                                         const std::vector<DeghostBox>& owners);

}  // namespace WireCell::Clus::PR

#endif
