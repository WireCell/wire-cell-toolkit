// doc sbnd_xin/docs/pr/51 round 5 -- steiner gap penalty: unit tests for
// the pure scan core (gap_edge_bad_fraction, the per-edge penalty input).
// The classifier is injected, so no cluster/grouping is needed; the full
// flavor build (ensure_steiner_gap_graph) needs an event and its validation
// record is the round-5 A/B arms and on-census.

#include "WireCellClus/PRSegmentFunctions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>

using namespace WireCell;
using namespace WireCell::Clus::PR;

namespace {

    const double STEP = 0.3 * units::cm;    // m_sgp_sample_step default
    const double ALPHA = 0.25;              // m_sgp_dead_alpha default

    // A 3 cm edge along z.
    const Point A(0, 0, 0);
    const Point B(0, 0, 3 * units::cm);
    const double LEN = 3 * units::cm;

}  // namespace

TEST_CASE("steiner gap penalty: uniform classifications")
{
    // All live => no penalty input.
    CHECK(gap_edge_bad_fraction(A, B, LEN, STEP, ALPHA,
                                [](const Point&) { return 0; }) == doctest::Approx(0.0));
    // All unsupported => full penalty input.
    CHECK(gap_edge_bad_fraction(A, B, LEN, STEP, ALPHA,
                                [](const Point&) { return 2; }) == doctest::Approx(1.0));
    // All dead => exactly alpha (dead is traversable but not free -- the
    // round-4 H2 finding: is_good_point blesses dead regions, the penalty
    // must not).
    CHECK(gap_edge_bad_fraction(A, B, LEN, STEP, ALPHA,
                                [](const Point&) { return 1; }) == doctest::Approx(ALPHA));
}

TEST_CASE("steiner gap penalty: half-gap geometry")
{
    // Support vanishes past the midpoint: bad_fraction ~ the unsupported
    // half.  10 intervals => 11 samples, 5 of them (z > 1.5 cm) unsupported.
    auto half = [](const Point& p) { return p.z() > 1.5 * units::cm ? 2 : 0; };
    CHECK(gap_edge_bad_fraction(A, B, LEN, STEP, ALPHA, half) == doctest::Approx(5.0 / 11.0));
}

TEST_CASE("steiner gap penalty: nsteps rounding floor")
{
    // length ~ step => nsteps = max(1, round(1)) = 1 => exactly two samples
    // (the endpoints).  One unsupported endpoint => 0.5.
    auto far_end = [](const Point& p) { return p.z() > 0.2 * units::cm ? 2 : 0; };
    CHECK(gap_edge_bad_fraction(A, Point(0, 0, 0.3 * units::cm), 0.3 * units::cm, STEP, ALPHA,
                                far_end) == doctest::Approx(0.5));
    // Degenerate: zero length still samples both (coincident) endpoints and
    // divides by n, never by zero.
    CHECK(gap_edge_bad_fraction(A, A, 0.0, STEP, ALPHA,
                                [](const Point&) { return 2; }) == doctest::Approx(1.0));
}

TEST_CASE("steiner gap penalty: chord-vs-arc arithmetic at the shipped scales")
{
    // The H1 mechanism in numbers (doc pr/51 round 4 sec "Root cause"): a
    // turn with equal legs L and interior angle theta offers Dijkstra a
    // supported detour of 2L against a gap-spanning chord of
    // 2L*sin(theta/2); the unpenalized weight cannot beat the chord for any
    // theta > ~150 deg because its charge factor is capped at 1.5x.  With
    // the gap penalty at scale s, a FULLY-unsupported chord costs
    // w * (1 + s), so the supported detour wins whenever
    //     1 / sin(theta/2) < 1 + s.
    auto detour_over_chord = [](double theta_deg) {
        return 1.0 / std::sin(theta_deg / 2.0 * M_PI / 180.0);
    };
    const double chord_factor_s2 = 1.0 + 2.0 * 1.0;  // scale 2, bad_fraction 1

    // Round-4 table rows: 150 deg (1.04), 120 deg (1.16), 90 deg (1.41),
    // 60 deg (2.00) -- every one is structurally lost today (cap 1.5 only
    // holds above ~150 deg with charge exactly zero at both endpoints) and
    // every one is won at scale 2.
    for (double theta : {150.0, 120.0, 90.0, 60.0}) {
        CHECK(detour_over_chord(theta) < chord_factor_s2);
    }
    // Scale 1 already covers the measured target events (detour/chord well
    // under 2): the worst measured reroute, 506746's 65.5 cm main track,
    // pays +0.38 cm on 76 cm -- ratio 1.005.
    CHECK(detour_over_chord(120.0) < 1.0 + 1.0);
    // The geometry the penalty deliberately does NOT override: a detour
    // longer than (1 + s) x chord keeps the chord even when unsupported
    // (no infinite-weight edges; the graph stays connected).
    CHECK(detour_over_chord(30.0) > chord_factor_s2);
}

// doc sbnd_xin/docs/pr/51 round 6 -- weak-charge deficit term.

TEST_CASE("weak charge deficit: thresholded endpoint form")
{
    const double QREF = 2000.0;  // m_sgp_weak_qref default
    // Both endpoints at/above the reference => no deficit.
    CHECK(weak_charge_deficit(QREF, QREF, QREF) == doctest::Approx(0.0));
    CHECK(weak_charge_deficit(5000.0, 3000.0, QREF) == doctest::Approx(0.0));
    // Both chargeless => full deficit.
    CHECK(weak_charge_deficit(0.0, 0.0, QREF) == doctest::Approx(1.0));
    // Linear below the threshold, averaged over the two endpoints.
    CHECK(weak_charge_deficit(1000.0, 1000.0, QREF) == doctest::Approx(0.5));
    CHECK(weak_charge_deficit(0.0, QREF, QREF) == doctest::Approx(0.5));
    CHECK(weak_charge_deficit(500.0, 3000.0, QREF) == doctest::Approx(0.375));
    // Degenerate qref never divides by zero and never penalizes.
    CHECK(weak_charge_deficit(0.0, 0.0, 0.0) == doctest::Approx(0.0));
    CHECK(weak_charge_deficit(0.0, 0.0, -1.0) == doctest::Approx(0.0));
}

TEST_CASE("weak charge deficit: 131357 chord-vs-corner arithmetic")
{
    // The round-6 target in numbers (doc pr/51 round 6, work-pr51r5-s2probe
    // arms): 18259-131357's trunk takes a 3.69 cm fully-SUPPORTED chord
    // whose interior charge is ~674-1276, cutting a corner whose real route
    // is ~5.4 cm at charge ~1900-5400.  The round-5 term sees bad = 0 on
    // both (all points live) -- scale-invariant.  The weak term at qref
    // 2000 prices the chord (endpoint charges ~1000 => deficit ~0.5) but
    // not the corner route (>= 2000 => deficit 0), so at weak_scale s the
    // chord costs ~3.69*(1 + s*0.5) against the corner's ~5.4: the corner
    // wins for s > ~0.93, comfortably at the grid scales 2-5.
    const double QREF = 2000.0;
    const double chord_deficit = weak_charge_deficit(1000.0, 1000.0, QREF);
    const double corner_deficit = weak_charge_deficit(2400.0, 4000.0, QREF);
    CHECK(corner_deficit == doctest::Approx(0.0));
    for (double s : {2.0, 3.0, 5.0}) {
        CHECK(3.69 * (1.0 + s * chord_deficit) > 5.4 * (1.0 + s * corner_deficit));
    }
    // And the round-5 term alone cannot: bad = 0 on both routes.
    CHECK(3.69 * (1.0 + 5.0 * 0.0) < 5.4);
}
