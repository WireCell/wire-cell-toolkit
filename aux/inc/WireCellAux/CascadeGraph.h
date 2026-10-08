/** The level graph of the coarse-to-fine learned deghosting (CascadeDeghosting, wcfm doc 14).

    One level = the blobs alive at one cut width, their wire (channel) nodes -- or super-wire nodes, channels
    grouped by channel // k, on the coarse levels -- and the four edge sets the e2-projq model reads:

      bw     blob -> (super-)wire, one per (blob, bin); bw_w = the blob's channels in the bin
      bb     cross-slice blob-blob, from Img::geom_clustering on the level's blobs (the BlobClustering rule)
      bb_in  in-slice blob-blob: same face and slice, half-open wire ranges overlap or abut on all three planes
      ww     wire-wire: consecutive covered wires of one (face, plane, slice) on different channels, binned

    It reproduces, from the toolkit's own objects, the arrays of the python that trained the models:
    wcp-porting-img wcfm/scripts/gnn_dataset.py build_graph (the 4-wire level, k = 1) and
    wcfm/scripts/d12_cascade.py contract_level (the super-wire levels, k > 1).  Documented differences:
    the charge is the slice activity (x charge_scale) instead of the SP frame; bb and bb_in are computed
    geometrically on the level's blobs instead of contracted from the 4-wire cells; the wire universe (every
    wire node of the level, including ones no alive blob touches) is the coverage of the uncut input blobs.

    Every container is ordered by stable indices (vertex order, slice start, channel ident).
 */
#ifndef WIRECELLAUX_CASCADEGRAPH
#define WIRECELLAUX_CASCADEGRAPH

#include "WireCellIface/IBlob.h"
#include "WireCellIface/ICluster.h"
#include "WireCellIface/IFrame.h"
#include "WireCellIface/ISlice.h"
#include "WireCellIface/ITrace.h"
#include "WireCellUtil/RayGrid.h"

#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace WireCell::Aux::Cascade {

    /// Canonical charge per slice time and channel, merged over the tiling passes' slices.
    struct SliceCharge {
        std::unordered_map<const ISlice*, int> sidx;   // lookup only: slice -> slice-time index
        std::vector<ISlice::pointer> slice_of;         // one representative slice per index, by start time
        std::vector<std::unordered_map<int, float>> q;  // activity mode: per index, channel ident -> charge x scale

        // frame mode (the training charge, crossview_probe.load_charge): scale x the sum of the frame's
        // tagged trace over the slice's ticks, summed in double and stored as float
        bool from_frame{false};
        double scale{0.25};
        std::vector<int> t0, nt;                                     // per index: first tick, ticks
        std::unordered_map<int, ITrace::pointer> traces;            // channel ident -> trace (lookup only)

        int index(const ISlice::pointer& s) const;
        float charge(int s, int ch) const;
    };

    /// From the s-nodes of a cluster graph.  Activity entries with uncertainty >= uncer_cut (dummy or
    /// masked channels) are ignored; a channel seen in several passes' slices of one time keeps the max.
    SliceCharge make_slice_charge(const cluster_graph_t& gr, double scale, double uncer_cut);
    /// The same from a list of slices (pdvd doc 130: the clustering retiler's slices); the graph version is this on
    /// the graph's s-nodes in vertex order.
    SliceCharge make_slice_charge(const std::vector<ISlice::pointer>& slices_in, double scale, double uncer_cut);

    /// The same slice indexing, with the charge taken from the frame's traces of `tag` ("" = all traces):
    /// charge(s, ch) = float(scale x sum over the slice's ticks of the trace), 0 for a channel without a trace.
    /// The slice's first tick is (slice start - frame time) / tick: right for slicers whose slice start is absolute
    /// (SumSlice).  MaskSlice creates its slices with a start relative to the frame, so for a frame with a non-zero
    /// time that is off by frame time / tick (pdvd doc 122 sec 3); slice_start_relative = true takes slice start / tick.
    /// The two are the same for a frame with time 0.
    SliceCharge make_slice_charge_frame(const cluster_graph_t& gr, const IFrame::pointer& frame,
                                        const std::string& tag, double scale, bool slice_start_relative = false);

    /// Channel ident of every wire of every blob strip, per plane: chans[p] sorted unique.
    std::array<std::vector<int>, 3> blob_channels(const IBlob::pointer& blob);

    /// The fine wire nodes (plane, slice index, channel) covered by the given (uncut) blobs, their charge,
    /// and the wire-wire pairs (gnn_dataset.wire_adjacency on this coverage).
    struct Universe {
        std::vector<std::array<int, 3>> wires;          // (plane, sidx, channel), sorted
        std::vector<float> q;
        std::vector<std::array<int64_t, 2>> ww;         // i < j, sorted unique
        std::unordered_map<int64_t, int> lookup;
        static int64_t key(int plane, int sidx, int ch) { return ((int64_t) plane * 100000000LL + sidx) * 1000000LL + ch; }
        int index(int plane, int sidx, int ch) const;
    };
    Universe make_universe(const std::vector<IBlob::pointer>& blobs, const SliceCharge& sc);

    struct Level {
        int k{1};
        std::vector<IBlob::pointer> blobs;
        std::vector<int> sidx;
        std::vector<float> xb;                 // N x 15 (b_charge)
        std::vector<float> wq;                 // per (super-)wire node
        std::vector<int64_t> wplane;
        std::vector<int> wchan, wsidx;         // per (super-)wire node: first channel of the bin, slice index
        std::vector<int64_t> bw_src, bw_dst;
        std::vector<float> bw_w;
        std::vector<int64_t> bb, bb_in, ww;    // flattened (i, j) pairs, each pair once
        std::array<double, 4> tpart{};         // build wall seconds: wires + features, bb, bb_in, ww (logging only)
        size_t nnodes() const { return blobs.size(); }
        size_t nedges() const { return bw_src.size() + bb.size() / 2 + bb_in.size() / 2; }
    };

    /// Build one level.  k = super-wire size (1 = real wires).  policy = the geom_clustering policy.
    /// nthreads: threads for the node features / bw edges and the bb search (the result does not depend on it).
    Level build_level(const std::vector<IBlob::pointer>& blobs, const SliceCharge& sc, const Universe& uni,
                      int k, const std::string& policy, int nthreads = 1);

    /// Run f(begin, end) over disjoint chunks covering [0, n), handed out dynamically to at most nthreads
    /// std::threads (nthreads <= 1 or n <= 1: f(0, n) in the calling thread); the first exception is rethrown.
    /// f must write only slots of its own items (the result must not depend on the schedule).
    void parallel_blocks(size_t n, int nthreads, const std::function<void(size_t, size_t)>& f);
    /// Run f(worker) on nthreads std::threads (1: in the calling thread); the first exception is rethrown.
    void parallel_workers(int nthreads, const std::function<void(int)>& f);

    /// Cross-slice blob pairs (i < j, sorted unique) that Img::geom_clustering with `policy` makes between the
    /// given blobs grouped into one blob set per slice index (sidx[i] indexes slice_of, time order).
    /// geom_pairs computes them directly: the same slice-pair loop and tolerances, and the RayGrid::overlap rule
    /// written as an interval test per wire layer (wcfm doc 15 round 1).  geom_pairs_graph is the reference, the
    /// same pairs read back from a cluster graph filled by Img::geom_clustering itself (doc 14's code path).
    /// Moved from img to aux (pdvd doc 130 phase 2) so that clus can run the cascade on re-tiled blobs.
    std::vector<std::array<int64_t, 2>> geom_pairs(const std::vector<IBlob::pointer>& blobs, const std::vector<int>& sidx,
                                                   const std::vector<ISlice::pointer>& slice_of,
                                                   const std::string& policy, int nthreads = 1);
    /// (geom_pairs_graph, the reference built with Img::geom_clustering itself, is declared in
    /// WireCellImg/CascadeGraph.h and lives in the img plugin, which owns geom_clustering.)

    /// Per slice index, the node indices (ascending); the non-empty ones in slice-index (time) order, with
    /// set_slice[k] = the slice index of set k.
    std::vector<std::vector<int64_t>> slice_sets(const std::vector<int>& sidx, size_t nslices, std::vector<int>& set_slice);

    /// In-slice adjacency (gnn_dataset.inslice_adjacency) of blobs given (group key, u0,u1,v0,v1,w0,w1).
    std::vector<std::array<int64_t, 2>> inslice_pairs(const std::vector<int64_t>& group,
                                                      const std::vector<std::array<int, 6>>& rng);

    /// The refinement operator: BlobCutting's bisection (duplicated in CascadeCut.cxx) continuing a node's
    /// bisection depth from the uncut blob, so that cutting level by level gives exactly one direct cut.
    struct CutParams {
        double nudge{0.01};
        int max_depth{10};    // bisections from the uncut blob (BlobCutting max_depth of the direct cut)
        int min_length{2};
    };
    bool needs_cutting(const RayGrid::Blob& blob, int length_threshold);
    /// Leaves (shape, depth) of the bisection of `blob` (at `depth`) until no strip is wider than `width`.
    std::vector<std::pair<RayGrid::Blob, int>> cut_shape(const RayGrid::Coordinates& coords, const RayGrid::Blob& blob,
                                                         int depth, int width, const CutParams& par);

    /// The coverage guard (d12_cascade.guard_prune): candidates in the given order are pruned only while
    /// every (super-)wire they touch keeps another surviving node.  Returns the pruned mask.
    std::vector<bool> guard_prune(const Level& lev, const std::vector<int>& cand_order, size_t& nguarded);

    /// The coverage guard of the FINAL keep set (wcfm doc 22, default off in CascadeDeghosting): every wire node with
    /// measured charge (wq > 0) that no kept node touches gets back its highest-logit node (ties: the lowest node
    /// index), wire nodes visited in index order, so every 2-D measurement keeps at least one 3-D explanation.
    /// keep is extended in place; returns the number of nodes added.
    size_t final_guard(const Level& lev, const std::vector<float>& logit, std::vector<bool>& keep);

    /// The dense-ambiguous-slice fallback of the final keep set (wcfm doc 17, default off in CascadeDeghosting).
    /// A (face, slice) group of level nodes is an ambiguous dense slab when it has >= nmin nodes, >= mmin nodes per
    /// charged wire node of its slice (wq > 0, all planes, both faces: wire nodes carry no face) and a fraction
    /// >= amin of its nodes with amb_lo < logit < amb_hi; there every node with logit >= t_keep is also kept.
    /// Truth-free, order-independent (wcfm scripts/d17_shower.py iso_fallback is the offline twin).
    struct IsoParams {
        int nmin{500};
        double mmin{4};
        double amin{0.8};
        double t_keep{-1.5};
        double amb_lo{-3.0};
        double amb_hi{0.5};
    };
    struct IsoResult {
        size_t nslices{0};   // triggered (face, slice) groups
        size_t nadded{0};    // nodes kept by the fallback that the keep set did not hold
    };
    /// face, sidx: per node; wq, wsidx: per wire node (slice index as sidx); keep is extended in place.
    IsoResult iso_fallback(const std::vector<int>& face, const std::vector<int>& sidx, const std::vector<float>& wq,
                           const std::vector<int>& wsidx, const std::vector<float>& logit, std::vector<bool>& keep,
                           const IsoParams& par);

}  // namespace WireCell::Aux::Cascade

#endif
