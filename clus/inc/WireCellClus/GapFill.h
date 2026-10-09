/** GapFill -- the pure decision core of the retiler's "gapfill" mode (pdvd doc 131).

    The retile of a beam-bundle cluster, redesigned so that the original image decides who owns which charge, the
    retile only bridges gaps of the cluster's own image along its path, and the bridging is done at the resolution
    of the deghosting cells.  Everything here is index arithmetic on slices and wire indices, so it is unit-tested
    without a detector (clus/test/doctest_gap_fill.cxx).  ImproveCluster_2 does the geometry and the tiling.
 */
#ifndef WIRECELLCLUS_GAPFILL
#define WIRECELLCLUS_GAPFILL

#include <array>
#include <cstddef>
#include <map>
#include <utility>
#include <vector>

namespace WireCell::Clus::GapFill {

    /// One cell / blob as a box of indices: slices [smin, smax) in ticks, wires [lo, hi) per plane.
    struct Box {
        int smin{0}, smax{0};
        std::array<std::pair<int, int>, 3> w{};
    };

    /// True when the box contains (slice, wire) within the margins: the slice within `slice_margin` ticks of
    /// [smin, smax) and every plane's wire within `wire_margin` wires of [lo, hi).
    bool covers(const Box& b, int slice, const std::array<int, 3>& wire, int wire_margin, int slice_margin);

    /// A sample of the cluster's path, converted to the face's slice (start tick, aligned) and wire indices.
    struct PathPt {
        int slice{0};
        std::array<int, 3> wire{};
        double s{0};          // arc length along the path
        bool on_face{true};   // false: the sample lies on another face / anode (never covered, never a gap)
    };

    /// A gap: path samples [first, last] (inclusive), all uncovered and on the face, bounded on both sides by a
    /// covered sample on the face.  A run that reaches an end of the path is not a gap (the track leaves the
    /// base there; nothing to bridge to).
    struct Gap {
        size_t first{0}, last{0};
        double length{0};   // s[last] - s[first]
        int nslices{0};     // distinct slices of the run
    };

    /// The gaps of a path given per-sample coverage; a run counts when it spans at least `min_slices` distinct
    /// slices and, if `max_length` > 0, is no longer than `max_length`.
    std::vector<Gap> find_gaps(const std::vector<PathPt>& pts, const std::vector<bool>& covered, int min_slices,
                               double max_length);

    /// Who owns a (plane, slice, wire) of the region besides the cluster: the charge other clusters' original
    /// blobs explain there, a blob's charge spread uniformly over its wires in each plane.  Sparse: only the
    /// region's slices are entered.
    class OwnerMap {
       public:
        void add(int plane, int slice, int lo, int hi, double charge);   // wires [lo, hi)
        bool is_other(int plane, int slice, int wire) const;
        double other_share(int plane, int slice, int wire) const;       // 0 when nobody else owns it
        size_t size() const;

       private:
        std::array<std::map<std::pair<int, int>, double>, 3> m_q;
    };

    /// The residual charge left for the cluster: measured minus what others explain, never negative.
    inline double residual(double measured, double other_share) { return measured > other_share ? measured - other_share : 0.0; }

    /// The corridor rule for one path sample: `nlive[p]` = live wires with charge (not owned by another cluster)
    /// in the corridor of plane p.  Candidates are made only when at most one plane has none (owner ruling 3 of
    /// doc 131: bridge wires in at most one plane).
    inline int missing_planes(const std::array<int, 3>& nlive)
    {
        return (nlive[0] == 0) + (nlive[1] == 0) + (nlive[2] == 0);
    }

}  // namespace WireCell::Clus::GapFill

#endif
