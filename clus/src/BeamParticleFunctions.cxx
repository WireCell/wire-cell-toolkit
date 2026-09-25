// doc pdvd/120: the selection rules of the beam-particle PR stage.
#include "WireCellClus/BeamParticleFunctions.h"

#include <algorithm>
#include <cmath>
#include <map>

using namespace WireCell;
using namespace WireCell::Clus::PR;

BeamBundlePick WireCell::Clus::PR::beam_particle_pick_bundle(const std::vector<BeamBundleRow>& mains,
                                                             const std::vector<BeamFlashInfo>& flashes,
                                                             double lo, double hi)
{
    BeamBundlePick pick;
    if (lo >= hi) { pick.why = "window-off"; return pick; }

    // gid -> longest in-window main of that bundle (int-keyed => sorted iteration)
    std::map<int, double> longest;
    for (const auto& m : mains) {
        if (m.gid < 0) continue;
        if (m.t0 < lo || m.t0 >= hi) continue;
        ++pick.n_in_window;
        auto it = longest.find(m.gid);
        if (it == longest.end()) longest.emplace(m.gid, m.length);
        else it->second = std::max(it->second, m.length);
    }
    pick.n_gids = static_cast<int>(longest.size());
    if (longest.empty()) { pick.why = "none"; return pick; }

    // the brightest valid flash among the in-window gids
    double best_pe = -1; int best_gid = -1; bool tie = false;
    for (const auto& f : flashes) {
        if (!f.valid || f.pe < 0) continue;
        if (!longest.count(f.gid)) continue;
        if (f.pe > best_pe) { best_pe = f.pe; best_gid = f.gid; tie = false; }
        else if (f.pe == best_pe && f.gid != best_gid) tie = true;
    }
    if (best_gid >= 0 && !tie) { pick.gid = best_gid; pick.why = "brightest"; return pick; }

    // fallback: the bundle holding the longest main; ties -> smallest gid (map order)
    double best_len = -1;
    for (const auto& [gid, len] : longest) {
        if (len > best_len) { best_len = len; pick.gid = gid; }
    }
    pick.why = "longest";
    return pick;
}

int WireCell::Clus::PR::beam_particle_pick_main(const std::vector<BeamMainCand>& cands, double min_length)
{
    int best = -1;
    for (size_t i = 0; i < cands.size(); ++i) {
        const auto& c = cands[i];
        if (c.length < min_length) continue;
        if (best < 0) { best = static_cast<int>(i); continue; }
        const auto& b = cands[best];
        if (c.dist < b.dist) best = static_cast<int>(i);
        else if (c.dist == b.dist) {
            if (c.length > b.length) best = static_cast<int>(i);
            else if (c.length == b.length && c.cluster_id < b.cluster_id) best = static_cast<int>(i);
        }
    }
    return best;
}

BeamEntryChoice WireCell::Clus::PR::beam_particle_choose_entry(const Point& a, const Point& b,
                                                               const Point& nominal,
                                                               const Vector& beam_dir, double tie_tol)
{
    BeamEntryChoice c;
    const double da = (a - nominal).magnitude();
    const double db = (b - nominal).magnitude();
    const double bmag = beam_dir.magnitude();
    auto cos_of = [&](const Point& e, const Point& x) {
        const Vector d = x - e;
        const double m = d.magnitude();
        if (m <= 0 || bmag <= 0) return 0.0;
        return d.dot(beam_dir) / (m * bmag);
    };
    bool a_first = da <= db;
    if (std::abs(da - db) < tie_tol) {
        const double ca = cos_of(a, b), cb = cos_of(b, a);
        if (ca != cb) { a_first = ca > cb; c.tie_by_dir = true; }
    }
    c.entry = a_first ? a : b;
    c.exit = a_first ? b : a;
    c.dist = a_first ? da : db;
    c.cos_beam = cos_of(c.entry, c.exit);
    return c;
}
