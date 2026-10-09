#include "WireCellClus/GapFill.h"

#include <set>

using namespace WireCell::Clus::GapFill;

bool WireCell::Clus::GapFill::covers(const Box& b, int slice, const std::array<int, 3>& wire, int wire_margin,
                                     int slice_margin)
{
    if (slice < b.smin - slice_margin || slice >= b.smax + slice_margin) return false;
    for (int p = 0; p < 3; ++p) {
        if (wire[p] < b.w[p].first - wire_margin || wire[p] >= b.w[p].second + wire_margin) return false;
    }
    return true;
}

std::vector<Gap> WireCell::Clus::GapFill::find_gaps(const std::vector<PathPt>& pts, const std::vector<bool>& covered,
                                                    int min_slices, double max_length)
{
    std::vector<Gap> gaps;
    const size_t n = pts.size();
    if (covered.size() != n) return gaps;
    bool seen_cover = false;   // a covered on-face sample before the run
    size_t run_first = 0;
    bool in_run = false;
    for (size_t i = 0; i < n; ++i) {
        if (!pts[i].on_face) {
            // leaving the face ends any run without closing it
            in_run = false;
            seen_cover = false;
            continue;
        }
        if (covered[i]) {
            if (in_run && seen_cover) {
                Gap g;
                g.first = run_first;
                g.last = i - 1;
                g.length = pts[g.last].s - pts[g.first].s;
                std::set<int> sl;
                for (size_t k = g.first; k <= g.last; ++k) sl.insert(pts[k].slice);
                g.nslices = (int) sl.size();
                if (g.nslices >= min_slices && (max_length <= 0 || g.length <= max_length)) gaps.push_back(g);
            }
            in_run = false;
            seen_cover = true;
            continue;
        }
        if (!in_run) {
            in_run = true;
            run_first = i;
        }
    }
    return gaps;
}

void OwnerMap::add(int plane, int slice, int lo, int hi, double charge)
{
    if (plane < 0 || plane > 2 || hi <= lo) return;
    const double share = (charge > 0 ? charge : 0.0) / (hi - lo);
    auto& m = m_q[plane];
    for (int w = lo; w < hi; ++w) m[{slice, w}] += share;
}

bool OwnerMap::is_other(int plane, int slice, int wire) const
{
    if (plane < 0 || plane > 2) return false;
    return m_q[plane].count({slice, wire}) > 0;
}

double OwnerMap::other_share(int plane, int slice, int wire) const
{
    if (plane < 0 || plane > 2) return 0.0;
    auto it = m_q[plane].find({slice, wire});
    return it == m_q[plane].end() ? 0.0 : it->second;
}

size_t OwnerMap::size() const { return m_q[0].size() + m_q[1].size() + m_q[2].size(); }
