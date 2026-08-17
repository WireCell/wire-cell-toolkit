// doc sbnd_xin/docs/pr/50: unit tests for path_scan_vertex_kink, the
// image-path corner detector behind the main-vertex kink-consistency snap
// (172230-class near-vertex robustness).  Synthetic point paths only -- the
// full pass needs an event; its validation record is the pr/50 A/B arms.

#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus::PR;

namespace {

    // Straight run from `start` along `dir` (unit), `n` points at `step`.
    std::vector<Point> straight(Point start, Vector dir, int n, double step)
    {
        std::vector<Point> pts;
        for (int i = 0; i < n; i++) pts.push_back(start + dir * (i * step));
        return pts;
    }

    // Elbow: straight along +z for arm1_len, then turn into `dir2`.
    std::vector<Point> elbow(double arm1_len, Vector dir2, double arm2_len, double step)
    {
        auto pts = straight(Point(0, 0, 0), Vector(0, 0, 1), static_cast<int>(arm1_len / step) + 1, step);
        const Point corner = pts.back();
        const int n2 = static_cast<int>(arm2_len / step);
        for (int i = 1; i <= n2; i++) pts.push_back(corner + dir2 * (i * step));
        return pts;
    }

    // The pr/50 C++ defaults (NeutrinoPatternBase.h m_vks_* members).
    const double D_MIN = 0.5 * units::cm;
    const double D_MAX = 5 * units::cm;
    const double SKIRT = 0.3 * units::cm;
    const double BASELINE = 2 * units::cm;
    const double ANGLE = 25;
    const double MIN_ARM = 1.5 * units::cm;

    std::vector<VertexKinkScanResult> scan(const std::vector<Point>& pts)
    {
        return path_scan_vertex_kink(pts, D_MIN, D_MAX, SKIRT, BASELINE, ANGLE, MIN_ARM);
    }

}  // namespace

TEST_CASE("pr50 path_scan_vertex_kink: straight path yields no candidates")
{
    auto pts = straight(Point(0, 0, 0), Vector(0, 0, 1), 21, 0.5 * units::cm);
    CHECK(scan(pts).empty());
}

TEST_CASE("pr50 path_scan_vertex_kink: 90-deg elbow at 2.5 cm is found at the corner")
{
    // The 172230 shape: vertex-side stub too short for the symmetric PCA
    // windows, so the chord fallback must carry the vertex-side direction.
    auto pts = elbow(2.5 * units::cm, Vector(1, 0, 0), 5 * units::cm, 0.5 * units::cm);
    auto cands = scan(pts);
    REQUIRE(!cands.empty());
    // The strongest candidate sits at the corner arclength with ~90 deg turn.
    const auto best = *std::max_element(cands.begin(), cands.end(),
        [](const auto& a, const auto& b) { return a.turn_deg < b.turn_deg; });
    CHECK(best.arc == doctest::Approx(2.5 * units::cm).epsilon(0.25));
    CHECK(best.turn_deg > 60);
}

TEST_CASE("pr50 path_scan_vertex_kink: corner outside the scan window is ignored")
{
    // Inside d_min: corner at 0.25 cm (indices there are below D_MIN).
    {
        auto pts = elbow(0.25 * units::cm, Vector(1, 0, 0), 6 * units::cm, 0.25 * units::cm);
        for (const auto& c : scan(pts)) CHECK(c.arc >= D_MIN);
        // No candidate may claim the sub-d_min corner itself.
        for (const auto& c : scan(pts)) CHECK(c.arc > 0.4 * units::cm);
    }
    // Beyond d_max: corner at 7 cm.
    {
        auto pts = elbow(7 * units::cm, Vector(1, 0, 0), 5 * units::cm, 0.5 * units::cm);
        for (const auto& c : scan(pts)) CHECK(c.arc <= D_MAX);
        CHECK(scan(pts).empty());
    }
}

TEST_CASE("pr50 path_scan_vertex_kink: corner without outward support is never claimed")
{
    // Only 1 cm of path beyond the corner (< MIN_ARM = 1.5 cm): the corner
    // index (and everything past it) is excluded.  A blurred PRE-corner
    // index with just enough support may still surface -- that is by
    // design; the snap pass's bendV / fit-miss guards do the filtering.
    auto pts = elbow(2.5 * units::cm, Vector(1, 0, 0), 1 * units::cm, 0.5 * units::cm);
    for (const auto& c : scan(pts)) {
        CHECK(c.arc < 2.5 * units::cm - 0.01 * units::cm);
    }
}

TEST_CASE("pr50 path_scan_vertex_kink: mm-scale zigzag noise does not fire")
{
    // Straight path with alternating +-1 mm transverse displacement -- the
    // wide-baseline PCA must average it away (the whole point of not using
    // point-to-point angles).
    auto pts = straight(Point(0, 0, 0), Vector(0, 0, 1), 21, 0.5 * units::cm);
    for (size_t i = 0; i < pts.size(); i++) {
        pts[i] = pts[i] + Vector(0, (i % 2) ? 0.1 * units::cm : -0.1 * units::cm, 0);
    }
    CHECK(scan(pts).empty());
}

TEST_CASE("pr50 path_scan_vertex_kink: degenerate inputs are safe")
{
    CHECK(path_scan_vertex_kink({}, D_MIN, D_MAX, SKIRT, BASELINE, ANGLE, MIN_ARM).empty());
    auto three = straight(Point(0, 0, 0), Vector(0, 0, 1), 3, 0.5 * units::cm);
    CHECK(path_scan_vertex_kink(three, D_MIN, D_MAX, SKIRT, BASELINE, ANGLE, MIN_ARM).empty());
}
