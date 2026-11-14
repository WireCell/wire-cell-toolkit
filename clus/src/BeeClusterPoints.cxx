#include "WireCellClus/BeeClusterPoints.h"
#include "WireCellClus/Facade_Cluster.h"

#include <cmath>
#include <cstdint>
#include <unordered_map>

using namespace WireCell;
using namespace WireCell::Clus;

bool Clus::bee_cluster_in_scope(const Facade::Cluster& cluster,
                                const PointCloud::Tree::Scope& scope, int filter)
{
    const bool filter_scope = cluster.get_scope_filter(scope);
    if (filter == 1) return filter_scope;
    if (filter == -1) return !filter_scope;
    return true;
}

void Clus::bee_cluster_points(const Facade::Cluster& cluster,
                              const PointCloud::Tree::Scope& scope, bool with_charge,
                              const std::function<void(size_t, const Point&, double)>& fn)
{
    const auto& coords = scope.coords;
    const auto& sv = cluster.sv<double>(scope);
    const auto& spcs = sv.pcs();

    // Dead-channel threshold: uncertainty > 1e10 flags a dead wire
    // (matches PointTreeBuilding m_dead_threshold and Facade_Cluster::is_wire_dead).
    const double dead_threshold = 1e10;

    // For each scoped point cloud (one per blob), compute per-point wire charge.
    // Each point gets q = mean(Q_U, Q_V, Q_W) over the non-dead planes for the
    // specific wires that define this 3D point -- matching the prototype formula.
    // The global cluster point index is obtained via a KD-tree nearest-neighbour
    // lookup (exact match: every scoped-view point is also in the cluster PC).
    for (size_t spc_idx = 0; spc_idx < spcs.size(); ++spc_idx) {
        const auto& spc = spcs[spc_idx];
        auto x = spc.get().get(coords[0])->elements<double>();
        auto y = spc.get().get(coords[1])->elements<double>();
        auto z = spc.get().get(coords[2])->elements<double>();

        const size_t size = x.size();
        for (size_t ind = 0; ind < size; ++ind) {
            const Point pt(x[ind], y[ind], z[ind]);
            double point_charge = 0.0;
            if (with_charge) {
                // charge_value() caches all per-plane charge vectors on first call,
                // so subsequent lookups are O(1) vector reads.
                const size_t pt_idx = cluster.get_closest_point_index(pt);
                double sum = 0.0;
                int nplanes = 0;
                for (int plane = 0; plane < 3; ++plane) {
                    if (!cluster.is_wire_dead(pt_idx, plane, dead_threshold)) {
                        sum += cluster.charge_value(pt_idx, plane);
                        ++nplanes;
                    }
                }
                point_charge = (nplanes > 0) ? sum / nplanes : 0.0;
            }
            fn(spc_idx, pt, point_charge);
        }
    }
}

std::vector<size_t> Clus::bee_touching_clusters(const std::vector<std::vector<float>>& clusters_xyz_cm,
                                                const std::vector<float>& cloud_x,
                                                const std::vector<float>& cloud_y,
                                                const std::vector<float>& cloud_z,
                                                double touch_cm)
{
    std::vector<size_t> out;
    const size_t ncloud = std::min({cloud_x.size(), cloud_y.size(), cloud_z.size()});
    if (!ncloud || !(touch_cm > 0)) return out;

    // A hash grid of the cloud with cell = touch_cm: a point within touch_cm of a
    // cloud point lies in one of the 27 cells around its own.
    const double cell = touch_cm;
    const double t2 = touch_cm * touch_cm;
    auto key = [](int64_t i, int64_t j, int64_t k) -> uint64_t {
        const uint64_t m = (1u << 21) - 1;
        return ((uint64_t(i) & m) << 42) | ((uint64_t(j) & m) << 21) | (uint64_t(k) & m);
    };
    auto icell = [cell](double v) -> int64_t { return static_cast<int64_t>(std::floor(v / cell)); };
    std::unordered_map<uint64_t, std::vector<uint32_t>> grid;
    grid.reserve(ncloud);
    for (size_t i = 0; i < ncloud; ++i) {
        grid[key(icell(cloud_x[i]), icell(cloud_y[i]), icell(cloud_z[i]))].push_back(static_cast<uint32_t>(i));
    }

    for (size_t ic = 0; ic < clusters_xyz_cm.size(); ++ic) {
        const auto& p = clusters_xyz_cm[ic];
        bool touched = false;
        for (size_t ip = 0; ip + 2 < p.size() && !touched; ip += 3) {
            const double px = p[ip], py = p[ip + 1], pz = p[ip + 2];
            const int64_t ci = icell(px), cj = icell(py), ck = icell(pz);
            for (int64_t di = -1; di <= 1 && !touched; ++di) {
                for (int64_t dj = -1; dj <= 1 && !touched; ++dj) {
                    for (int64_t dk = -1; dk <= 1 && !touched; ++dk) {
                        auto it = grid.find(key(ci + di, cj + dj, ck + dk));
                        if (it == grid.end()) continue;
                        for (const auto j : it->second) {
                            const double dx = px - cloud_x[j], dy = py - cloud_y[j], dz = pz - cloud_z[j];
                            if (dx * dx + dy * dy + dz * dz <= t2) { touched = true; break; }
                        }
                    }
                }
            }
        }
        if (touched) out.push_back(ic);
    }
    return out;
}
