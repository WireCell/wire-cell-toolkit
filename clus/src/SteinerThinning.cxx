#include "WireCellClus/SteinerThinning.h"

#include <array>
#include <cmath>
#include <map>

using namespace WireCell;

std::vector<size_t> Clus::Steiner::thin_by_min_separation(
    const std::vector<std::pair<Point, size_t>>& ordered,
    double min_separation)
{
    std::vector<size_t> kept;
    kept.reserve(ordered.size());

    // <= 0 is the OFF path and must be an exact pass-through: the caller's
    // byte-identity guarantee rests on it.
    if (min_separation <= 0.0) {
        for (const auto& [pt, id] : ordered) {
            (void) pt;
            kept.push_back(id);
        }
        return kept;
    }

    const double r2 = min_separation * min_separation;

    // Uniform grid of cell edge min_separation over the KEPT points only, so a
    // candidate need only consult its own cell and the 26 around it.  std::floor
    // rather than a cast, which truncates toward zero: PDVD spans x in
    // [-339.9, +339.9] cm and y in [-336.4, +336.4], so half the detector has
    // negative coordinates, and a cast would fuse (-r, +r) into one
    // double-width cell.  The +-1 sweep still covers that cell's reach, so the
    // result would not change -- but the cell edge would no longer BE
    // min_separation, and the "27 cells is enough" argument above would stop
    // being true of the code as written.  Keep the invariant, not the accident.
    using cell_t = std::array<long, 3>;
    auto cell_of = [min_separation](const WireCell::Point& p) -> cell_t {
        return cell_t{ (long) std::floor(p.x() / min_separation),
                       (long) std::floor(p.y() / min_separation),
                       (long) std::floor(p.z() / min_separation) };
    };

    // Ordered map keyed by the cell triple.  Lookup only -- the container is
    // never iterated, so no ordering enters the result.
    std::map<cell_t, std::vector<WireCell::Point>> grid;

    for (const auto& [pt, id] : ordered) {
        const cell_t c = cell_of(pt);

        bool suppressed = false;
        for (long dx = -1; dx <= 1 && !suppressed; ++dx) {
            for (long dy = -1; dy <= 1 && !suppressed; ++dy) {
                for (long dz = -1; dz <= 1 && !suppressed; ++dz) {
                    auto it = grid.find(cell_t{c[0] + dx, c[1] + dy, c[2] + dz});
                    if (it == grid.end()) continue;
                    for (const auto& q : it->second) {
                        const double ddx = pt.x() - q.x();
                        const double ddy = pt.y() - q.y();
                        const double ddz = pt.z() - q.z();
                        // Strict: a terminal exactly min_separation away survives.
                        if (ddx*ddx + ddy*ddy + ddz*ddz < r2) {
                            suppressed = true;
                            break;
                        }
                    }
                }
            }
        }
        if (suppressed) continue;

        kept.push_back(id);
        grid[c].push_back(pt);
    }

    return kept;
}
