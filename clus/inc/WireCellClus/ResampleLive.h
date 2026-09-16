/** Pure helpers of ClusteringResampleLive (clus/src/clustering_resample_live.cxx).

    doc sbnd_xin/pr/149 round 2.  The visitor re-samples every live blob's "3d"
    point cloud inside the pattern-recognition job, which holds only the saved
    point-cloud tree (no IBlob, no ISlice).  It therefore has to rebuild, from a
    blob's "scalar" PC and the grouping's ctpc / dead-wind records, the two
    objects a BlobSampler reads: the blob SHAPE and the slice ACTIVITY.  Both
    rebuilds are split out here so a doctest can pin them directly -- the same
    reason SteinerThinning.h exists.

    prototype (wire-cell-prod-nue.cxx:1289-1294, wire-cell-prod-stm.cxx:734):
        // replace by the new sampling points ...
        for (size_t i=0; i!=live_clusters.size();i++){
          WCPPID::calc_sampling_points(gds,live_clusters.at(i),nrebin, frame_length, unit_dis);
          live_clusters.at(i)->Create_point_cloud();
        }
    i.e. the PR executables replace the clustering job's stepped cloud
    (WCP2dToy::calc_sampling_points) with the charge_stepped rule,
    disable_mix_dead_cell = true (pid/inc/WCPPID/CalcPoints.h:9-10).
 */

#ifndef WIRECELLCLUS_RESAMPLELIVE
#define WIRECELLCLUS_RESAMPLELIVE

#include "WireCellUtil/RayGrid.h"
#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/PointCloudDataset.h"

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace WireCell::Clus::ResampleLive {

    /// Half-open [min, max) wire-index bounds per plane (U, V, W), exactly as a
    /// blob's "scalar" PC stores them (aux/src/SamplingHelpers.cxx:100-106).
    using plane_bounds_t = std::array<std::pair<int, int>, 3>;

    /// The value MaskSlice writes for a dead / masked channel
    /// (img/inc/WireCellImg/MaskSlice.h: m_dummy_charge / m_dummy_error,
    /// m_masked_charge / m_masked_error).  The ctpc keeps only channels whose
    /// uncertainty is <= the dead threshold (aux/src/SamplingHelpers.cxx add_ctpc),
    /// so a dead channel survives only as a dead-wind x range and is put back with
    /// these values.
    constexpr double dead_charge = 0.0;
    constexpr double dead_error = 1e12;

    /// The five strips an imaging blob carries: the two dummy bounding layers
    /// {0, 1} with bounds {0, 1}, then U, V, W at layers 2, 3, 4 -- the layout of
    /// aux/src/ClusterArrays.cxx:777-815.  The dummy layers are load-bearing:
    /// BlobSampler derives its plane offset from strips.size() == 5
    /// (BlobSampler.cxx Stepped/ChargeStepped sample()).
    inline std::vector<RayGrid::Strip> strips_from_bounds(const plane_bounds_t& b)
    {
        std::vector<RayGrid::Strip> strips;
        strips.reserve(5);
        strips.push_back(RayGrid::Strip{0, {0, 1}});
        strips.push_back(RayGrid::Strip{1, {0, 1}});
        for (int p = 0; p < 3; ++p) {
            strips.push_back(RayGrid::Strip{2 + p, {b[p].first, b[p].second}});
        }
        return strips;
    }

    /// The blob shape from its bounds.  nudge as ClusterArrays' loader default
    /// (aux/inc/WireCellAux/ClusterArrays.h:98-117); it only affects corners,
    /// which the stepped / charge_stepped samplers read only for their empty-blob
    /// center fallback.
    inline RayGrid::Blob shape_from_bounds(const RayGrid::Coordinates& coords,
                                           const plane_bounds_t& b, double nudge = 1e-3)
    {
        RayGrid::Blob shape;
        for (const auto& strip : strips_from_bounds(b)) {
            shape.add(coords, strip, nudge);
        }
        return shape;
    }

    /// One wire's reconstructed slice activity.
    struct WireValue {
        int wire{0};
        double charge{0};
        double error{0};
    };

    /// Compose the activity of wires [lo, hi] (INCLUSIVE) of one plane in one
    /// time slice.  Precedence:
    ///   1. live   -- live(wire, q, err) returns true: the ctpc entry, kept as is;
    ///   2. dead   -- otherwise dead(wire) returns true: (dead_charge, dead_error);
    ///   3. absent -- otherwise the wire is left out, exactly as imaging leaves a
    ///                live channel with no signal out of the slice.
    /// Live wins over dead on purpose: the dead-wind registry is an x-range
    /// UNION over slices (PointTreeBuilding add_dead_winds, +-0.1 cm) plus the
    /// hand-declared regions, so it can cover a (slice, wire) that imaging held
    /// live.  Never call Grouping::get_wire_charge() for "absent": it returns
    /// {0, 1e12} for a missing wire, which BlobSampler's is_plane_bad would read
    /// as a dead boundary wire.
    template <typename LiveFn, typename DeadFn>
    std::vector<WireValue> compose_activity(int lo, int hi, LiveFn&& live, DeadFn&& dead)
    {
        std::vector<WireValue> out;
        if (hi < lo) return out;
        out.reserve(hi - lo + 1);
        for (int w = lo; w <= hi; ++w) {
            double q = 0, err = 0;
            if (live(w, q, err)) {
                out.push_back(WireValue{w, q, err});
            }
            else if (dead(w)) {
                out.push_back(WireValue{w, dead_charge, dead_error});
            }
        }
        return out;
    }

    /// The four "scalar" arrays that depend on the "3d" points and are read by
    /// Facade::Blob (hash(), blob_less, sanity()).  Everything else in "scalar"
    /// (charge, wpid, slice and wire bounds, wire intervals) is a property of
    /// the blob, not of its points, and is kept.
    inline const std::vector<std::string>& point_dependent_scalar_names()
    {
        static const std::vector<std::string> names = {"center_x", "center_y", "center_z", "npoints"};
        return names;
    }

    /// Replace the point-dependent arrays of `scalar` by those of `fresh`
    /// (erase first: Dataset::add throws on a duplicate key).  Arrays absent
    /// from `fresh` are left untouched.
    inline void refresh_scalar(PointCloud::Dataset& scalar, const PointCloud::Dataset& fresh)
    {
        for (const auto& name : point_dependent_scalar_names()) {
            auto arr = fresh.get(name);
            if (!arr) continue;
            scalar.erase(name);
            scalar.add(name, *arr);
        }
    }

}  // namespace WireCell::Clus::ResampleLive

#endif
