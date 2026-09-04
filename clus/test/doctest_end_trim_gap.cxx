/** doc pdvd/38 -- the detached-island rule of the trajectory end trim.
 *
 * `TrackFitting::examine_end_ps_vec` needs a Grouping, a DetectorVolumes and a
 * PCTransformSet, so it is not reachable from a unit test.  The rule it applies
 * is, though: `Clus::scan_tip_island` is the whole of the new logic and is
 * exercised here on synthetic paths with a synthetic support predicate, with
 * the caller's decision spelled out as `detached_island()` so the two clauses
 * are tested exactly as `examine_end_ps_vec` combines them.
 * The knob that gates it (`TrackFitting` parameter `end_trim_gap_len`) is
 * checked to default to 0 = off.
 */
#include "WireCellClus/EndTrimGap.h"
#include "WireCellClus/TrackFitting.h"
#include "WireCellClus/TrackFittingPresets.h"
#include "WireCellUtil/Point.h"
#include "WireCellUtil/Units.h"

#include "WireCellUtil/doctest.h"

#include <cmath>
#include <list>
#include <set>
#include <vector>

using namespace WireCell;
using WireCell::Clus::scan_tip_island;

static const double STEP = 0.6 * units::cm;      // the fitted-trajectory spacing

/// A straight path along z at STEP, so path length == index distance * STEP.
static std::list<Point> straight_path(size_t n)
{
    std::list<Point> out;
    for (size_t i = 0; i < n; ++i) out.push_back(Point(0, 0, i * STEP));
    return out;
}

/// index of a point on such a path
static size_t idx_of(const Point& p) { return static_cast<size_t>(std::lround(p.z() / STEP)); }

/// exactly the decision examine_end_ps_vec makes
static bool detached_island(const Clus::TipIsland& isl, double L)
{
    return isl.detached && isl.gap_len > L
        && isl.gap_len  > isl.supported_len
        && isl.body_len > isl.supported_len;
}

TEST_CASE("clus end trim gap: a short tip stretch across a long gap is an island")
{
    auto path = straight_path(200);
    const double L = 5 * units::cm;

    SUBCASE("the cluster-109 shape: a 1-point island holding a 54 cm chord") {
        // supported at the tip, then 90 unsupported points (~54 cm), then body
        auto good = [](const Point& p) { const size_t i = idx_of(p); return i == 0 || i > 90; };
        auto isl = scan_tip_island(path.begin(), path.end(), good);
        CHECK(isl.npoints == 1);
        CHECK(isl.supported_len == doctest::Approx(0.0));
        CHECK(isl.gap_len == doctest::Approx(91 * STEP));   // 54.6 cm
        CHECK(isl.detached);
        CHECK(detached_island(isl, L));
    }

    SUBCASE("the cluster-53 shape: a 2-point island -- the immediate neighbour PASSES") {
        // This is why the rule cannot be "is the next point supported?":
        // the island is longer than one point.
        auto good = [](const Point& p) { const size_t i = idx_of(p); return i <= 1 || i > 110; };
        auto isl = scan_tip_island(path.begin(), path.end(), good);
        CHECK(isl.npoints == 2);
        CHECK(isl.supported_len == doctest::Approx(STEP));
        CHECK(isl.gap_len == doctest::Approx(110 * STEP));
        CHECK(detached_island(isl, L));
    }

    SUBCASE("a real track end beyond a dead region: the stretch dwarfs the gap, KEPT") {
        // 100 supported points (60 cm), then a 20-point (12 cm) dead region,
        // then the body.  The gap clears the threshold but not the stretch.
        auto good = [](const Point& p) { const size_t i = idx_of(p); return i < 100 || i >= 120; };
        auto isl = scan_tip_island(path.begin(), path.end(), good);
        CHECK(isl.npoints == 100);
        CHECK(isl.supported_len == doctest::Approx(99 * STEP));
        CHECK(isl.gap_len == doctest::Approx(21 * STEP));
        CHECK(isl.detached);
        CHECK(isl.gap_len > L);                 // the threshold alone would fire
        CHECK_FALSE(detached_island(isl, L));   // the comparison clauses save it
    }

    SUBCASE("the chewing case: a long gap but no body bigger than the island") {
        // Measured on real trajectories (doc 38 sec 2): a path whose supported
        // runs are ALL short.  Without the body clause the rule pops the tip,
        // the pop loop eats the gap, the next short run becomes the tip, and
        // the whole trajectory is amputated one island at a time.
        // 2 supported, a 100-point gap, then only 2 supported again.
        auto good = [](const Point& p) {
            const size_t i = idx_of(p);
            return i <= 1 || (i >= 102 && i <= 103);
        };
        auto isl = scan_tip_island(path.begin(), path.end(), good);
        CHECK(isl.npoints == 2);
        CHECK(isl.supported_len == doctest::Approx(STEP));
        CHECK(isl.gap_len == doctest::Approx(101 * STEP));
        CHECK(isl.detached);
        CHECK(isl.gap_len > L);                     // threshold fires
        CHECK(isl.gap_len > isl.supported_len);     // gap clause fires
        CHECK(isl.body_len == doctest::Approx(STEP));
        CHECK_FALSE(isl.body_len > isl.supported_len);
        CHECK_FALSE(detached_island(isl, L));       // the body clause saves it
    }

    SUBCASE("body_len stops measuring as soon as it beats the island") {
        // cost bound: the body is 90 points but only the first step or two are
        // ever walked, because the island is one point long.
        auto good = [](const Point& p) { const size_t i = idx_of(p); return i == 0 || i > 90; };
        auto isl = scan_tip_island(path.begin(), path.end(), good);
        CHECK(isl.supported_len == doctest::Approx(0.0));
        CHECK(isl.body_len == doctest::Approx(STEP));   // one step, then it wins
        CHECK(isl.body_len > isl.supported_len);
        CHECK(detached_island(isl, L));
    }

    SUBCASE("a short gap inside a track is not an island at any stretch length") {
        auto good = [](const Point& p) { const size_t i = idx_of(p); return i != 1 && i != 2; };
        auto isl = scan_tip_island(path.begin(), path.end(), good);
        CHECK(isl.npoints == 1);
        CHECK(isl.gap_len == doctest::Approx(3 * STEP));    // 1.8 cm
        CHECK(isl.detached);
        CHECK_FALSE(detached_island(isl, L));
    }

    SUBCASE("no body beyond the gap: detached is false, so nothing is popped") {
        // (the same sentinel the chewing case relies on, in its extreme form)
        // The conservative sentinel.  Without it the caller would pop the tip,
        // then the next point, and amputate the whole path one point at a time.
        auto good = [](const Point& p) { return idx_of(p) == 0; };
        auto isl = scan_tip_island(path.begin(), path.end(), good);
        CHECK(isl.npoints == 1);
        CHECK_FALSE(isl.detached);
        CHECK_FALSE(detached_island(isl, L));
        CHECK_FALSE(detached_island(isl, 0.1 * units::cm));
    }

    SUBCASE("everything supported: no gap, nothing popped") {
        auto all_good = [](const Point&) { return true; };
        auto isl = scan_tip_island(path.begin(), path.end(), all_good);
        CHECK(isl.npoints == 200);
        CHECK(isl.gap_len == doctest::Approx(0.0));
        CHECK_FALSE(isl.detached);
        CHECK_FALSE(detached_island(isl, L));
    }

    SUBCASE("walking inward from the far end uses the reverse iterators") {
        // island at the FAR end: supported at 199, unsupported 109..198, body below
        auto good = [](const Point& p) { const size_t i = idx_of(p); return i == 199 || i < 109; };
        auto fwd = scan_tip_island(path.begin(), path.end(), good);
        CHECK(fwd.npoints == 109);                       // the body, from the near end
        CHECK_FALSE(detached_island(fwd, L));
        auto rev = scan_tip_island(path.rbegin(), path.rend(), good);
        CHECK(rev.npoints == 1);
        CHECK(rev.gap_len == doctest::Approx(91 * STEP));
        CHECK(detached_island(rev, L));
    }

    SUBCASE("the scan cap bounds the predicate calls and reads as not detached") {
        size_t calls = 0;
        auto good = [&calls](const Point& p) { ++calls; return idx_of(p) == 0 || idx_of(p) > 90; };
        auto isl = scan_tip_island(path.begin(), path.end(), good, 7);
        CHECK(calls <= 7);
        CHECK_FALSE(isl.detached);
        CHECK_FALSE(detached_island(isl, L));
    }

    SUBCASE("degenerate ranges") {
        std::list<Point> empty_path;
        auto all_good = [](const Point&) { return true; };
        auto e = scan_tip_island(empty_path.begin(), empty_path.end(), all_good);
        CHECK(e.npoints == 0);
        CHECK_FALSE(e.detached);
        auto one = straight_path(1);
        auto o = scan_tip_island(one.begin(), one.end(), all_good);
        CHECK(o.npoints == 1);
        CHECK_FALSE(o.detached);
    }

    SUBCASE("lengths follow the polyline, not the chord") {
        // an L-bend: the gap turns a corner, so the chord would understate it
        std::list<Point> bent;
        for (size_t i = 0; i < 6; ++i) bent.push_back(Point(0, 0, i * STEP));
        for (size_t i = 1; i <= 5; ++i) bent.push_back(Point(0, i * STEP, 5 * STEP));
        auto good = [](const Point& p) {
            return p.z() == doctest::Approx(0.0) || p.y() == doctest::Approx(5 * STEP);
        };
        auto isl = scan_tip_island(bent.begin(), bent.end(), good);
        CHECK(isl.npoints == 1);
        CHECK(isl.gap_len == doctest::Approx(10 * STEP));   // chord would be 7.07
        CHECK(isl.detached);
    }
}

TEST_CASE("clus knob defaults: TrackFitting end_trim_gap_len is off")
{
    // doc pdvd/38.  0 = the legacy tip-only trim.  The knob is a LENGTH in WCT
    // internal units; the PDVD operating point lives in pdvd_track_fitting.json.
    Clus::TrackFitting tf;
    CHECK(tf.get_parameter("end_trim_gap_len") == doctest::Approx(0.0));

    tf.set_parameter("end_trim_gap_len", 5.0 * units::cm);
    CHECK(tf.get_parameter("end_trim_gap_len") == doctest::Approx(5.0 * units::cm));

    auto preset = Clus::TrackFittingPresets::create_with_current_values();
    CHECK(preset.get_parameters().end_trim_gap_len == doctest::Approx(0.0));
}
