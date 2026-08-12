#include "WireCellClus/PRSegment.h"
#include "WireCellClus/PRGraph.h"
#include "WireCellClus/PRSegmentFunctions.h"

#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus;
using namespace WireCell::Clus::PR;

TEST_CASE("clus pr segment basic") {
    PR::Segment seg;

    CHECK(! seg.descriptor_valid());
}

// F7 regression: find_vertices must not dereference wcpts().front() when wcpts
// is empty.  Before the fix, calling find_vertices on a segment built only from
// fits (no wcpts) caused undefined behaviour.  After the fix the code falls back
// to fits().front() if wcpts is empty.
TEST_CASE("clus pr segment find_vertices empty wcpts fallback") {
    Graph g;

    // Two vertices at known positions along X.
    auto vtx1 = make_vertex(g);
    vtx1->wcpt().point = Point(0, 0, 0);

    auto vtx2 = make_vertex(g);
    vtx2->wcpt().point = Point(10*units::cm, 0, 0);

    // Segment with fits but NO wcpts.
    auto seg = make_segment(g, vtx1, vtx2);
    REQUIRE(seg->wcpts().empty());

    // Add a fit closer to vtx1 so find_vertices can order the pair.
    Fit f;
    f.point = Point(1*units::cm, 0, 0);
    f.index = 0;
    seg->fits().push_back(f);

    // Must not crash and must return both vertices non-null.
    auto [va, vb] = find_vertices(g, seg);
    REQUIRE(va != nullptr);
    REQUIRE(vb != nullptr);

    // va should be closer to the fit point (i.e. vtx1 at x=0 cm).
    CHECK(va.get() == vtx1.get());
    CHECK(vb.get() == vtx2.get());
}

// F5 regression: clear_fit(nullptr) must not attempt to populate paf (no DV),
// so every fit's paf must remain at the default {-1,-1}.
// (Testing with a non-null DV requires a Facade::Cluster attached to the segment
// because create_segment_fit_point_cloud internally accesses the cluster's
// grouping.  That is covered by integration tests.)
TEST_CASE("clus pr segment clear_fit paf with null dv") {
    // Segment must be owned by a shared_ptr for shared_from_this().
    auto seg = std::make_shared<Segment>();

    // Populate wcpts at two known positions.
    WCPoint wp1; wp1.point = Point(0, 0, 0);
    WCPoint wp2; wp2.point = Point(5*units::cm, 0, 0);
    seg->wcpts({wp1, wp2});

    // clear_fit(nullptr) should reset fits to match wcpts but leave paf at default.
    seg->clear_fit(nullptr);

    REQUIRE(seg->fits().size() == 2);
    for (const auto& fit : seg->fits()) {
        CHECK(fit.paf.first  == -1);
        CHECK(fit.paf.second == -1);
    }

    // Points must be copied from wcpts.
    CHECK(seg->fits()[0].point == Point(0, 0, 0));
    CHECK(seg->fits()[1].point == Point(5*units::cm, 0, 0));

    // All other numeric fields must be reset to their defaults.
    for (const auto& fit : seg->fits()) {
        CHECK(fit.dQ           == doctest::Approx(-1.0));
        CHECK(fit.dx           == doctest::Approx( 0.0));
        CHECK(fit.reduced_chi2 == doctest::Approx(-1.0));
        CHECK(fit.index        == -1);
        CHECK(fit.flag_fix     == false);
    }
}

// F6 regression: both segment_track_direct_length and segment_track_max_deviation
// must interpret n2 = -1 as "last fit index" (fits.size()-1), not as 0.
// Before the fix, segment_track_direct_length clamped n2<0 to 0, so calling
// (seg, 0, -1) on a 5-fit segment measured from index 0 to 0 (length 0).
TEST_CASE("clus pr segment track functions n2 minus1 clamp") {
    // Build a segment with 5 collinear fits at x = 0..4 cm, y=z=0.
    auto seg = std::make_shared<Segment>();
    std::vector<Fit> fits(5);
    for (int i = 0; i < 5; ++i) {
        fits[i].point = Point(i * units::cm, 0, 0);
        fits[i].index = i;
    }
    seg->fits(fits);

    const double full_len = 4*units::cm;   // x=0 to x=4 cm

    // -- segment_track_direct_length --

    // Default (n1=-1, n2=-1): full segment.
    CHECK(segment_track_direct_length(seg) == doctest::Approx(full_len));

    // (0, -1) must equal (0, 4): length from first to last fit.
    CHECK(segment_track_direct_length(seg, 0, -1) == doctest::Approx(full_len));

    // (2, -1) must equal (2, 4): partial length from index 2 to end.
    CHECK(segment_track_direct_length(seg, 2, -1) == doctest::Approx(2*units::cm));

    // Explicit endpoints must agree.
    CHECK(segment_track_direct_length(seg, 0, 4) ==
          doctest::Approx(segment_track_direct_length(seg, 0, -1)));

    // -- segment_track_max_deviation --
    // Fits are collinear so max deviation is 0 in all sub-ranges.

    CHECK(segment_track_max_deviation(seg)        == doctest::Approx(0.0));
    CHECK(segment_track_max_deviation(seg, 0, -1) == doctest::Approx(0.0));
    CHECK(segment_track_max_deviation(seg, 2, -1) == doctest::Approx(0.0));

    // Explicit endpoints must agree.
    CHECK(segment_track_max_deviation(seg, 0, 4) ==
          doctest::Approx(segment_track_max_deviation(seg, 0, -1)));
}

// doc sbnd_xin/docs/pr/32 §11 F2 (was P3).
//
// The prototype's ProtoSegment::is_shower_trajectory() opens with
// `flag_shower_trajectory = false` (ProtoSegment.cxx:544) and only sets it true
// if the test fires (:608), so every call RE-CACHES the label and a segment that
// no longer qualifies is demoted.  The toolkit's port only ever calls set_flags,
// so kShowerTrajectory is monotone -- which is why the two improve_vertex gates
// recompute instead of reading it.
//
// PR::g_shower_traj_refresh_flag restores the prototype's semantics.  These two
// cases pin BOTH polarities, so reverting the refresh block in
// segment_is_shower_trajectory fails the second one.
TEST_CASE("clus pr32 F2 shower-trajectory flag is monotone by default") {
    Graph g;
    auto vtx1 = make_vertex(g);
    vtx1->wcpt().point = Point(0, 0, 0);
    auto vtx2 = make_vertex(g);
    vtx2->wcpt().point = Point(60*units::cm, 0, 0);
    auto seg = make_segment(g, vtx1, vtx2);

    // Longer than the function's 50 cm ceiling, so it returns false on the
    // early-out path -- the same path the prototype still clears on.
    for (int i = 0; i <= 6; ++i) {
        Fit f;
        f.point = Point(i * 10*units::cm, 0, 0);
        f.index = i;
        seg->fits().push_back(f);
    }
    REQUIRE(segment_track_length(seg, 0) > 50*units::cm);

    seg->set_flags(SegmentFlags::kShowerTrajectory);

    const bool saved = g_shower_traj_refresh_flag;
    g_shower_traj_refresh_flag = false;
    CHECK(segment_is_shower_trajectory(seg, 10*units::cm, 50000/units::cm) == false);
    // Today's behaviour: the negative answer does NOT demote the label.
    CHECK(seg->flags_any(SegmentFlags::kShowerTrajectory) == true);
    g_shower_traj_refresh_flag = saved;
}

TEST_CASE("clus pr32 F2 refresh clears the flag on a negative answer") {
    Graph g;
    auto vtx1 = make_vertex(g);
    vtx1->wcpt().point = Point(0, 0, 0);
    auto vtx2 = make_vertex(g);
    vtx2->wcpt().point = Point(60*units::cm, 0, 0);
    auto seg = make_segment(g, vtx1, vtx2);

    for (int i = 0; i <= 6; ++i) {
        Fit f;
        f.point = Point(i * 10*units::cm, 0, 0);
        f.index = i;
        seg->fits().push_back(f);
    }

    seg->set_flags(SegmentFlags::kShowerTrajectory);

    const bool saved = g_shower_traj_refresh_flag;
    g_shower_traj_refresh_flag = true;
    CHECK(segment_is_shower_trajectory(seg, 10*units::cm, 50000/units::cm) == false);
    // Prototype behaviour: the label is demoted, so a stored-flag read
    // downstream now agrees with the test.
    CHECK(seg->flags_any(SegmentFlags::kShowerTrajectory) == false);
    // Other flags must survive -- unset_flags is bit-masked, not clear_flags.
    seg->set_flags(SegmentFlags::kShowerTopology);
    CHECK(segment_is_shower_trajectory(seg, 10*units::cm, 50000/units::cm) == false);
    CHECK(seg->flags_any(SegmentFlags::kShowerTopology) == true);
    g_shower_traj_refresh_flag = saved;
}

// ---------------------------------------------------------------------------
// doc pr/47 sec 8 (O1): segment_cathode_wide_kink_accepts -- the wide-baseline
// skirt-excluded PCA turn angle across a cathode crossing.  Synthetic fit
// trajectories with the production-like 0.6 cm point spacing.
// ---------------------------------------------------------------------------
namespace {

// Polyline of Fits: arm A travels along dir_a and ends STEP/2 before x=0,
// arm B starts STEP/2 after x=0 along dir_b.  Both arms arm_len long.
// offset_b translates the whole B arm (models the transverse cathode
// mismatch, which shifts but does not rotate the downstream arm).
std::vector<Fit> two_arm_track(const Vector& dir_a, const Vector& dir_b,
                               double arm_len_a, double arm_len_b,
                               const Vector& offset_b = Vector(0, 0, 0))
{
    const double STEP = 0.6 * units::cm;
    std::vector<Fit> fits;
    auto ua = dir_a * (1.0 / dir_a.magnitude());
    auto ub = dir_b * (1.0 / dir_b.magnitude());
    // Arm A: from -arm_len_a (measured along ua) up to the crossing.
    const int na = static_cast<int>(arm_len_a / STEP);
    Point junction(0, 0, 0);
    for (int k = na; k >= 1; k--) {
        Fit f;
        f.point = junction - ua * (k * STEP - STEP / 2);
        fits.push_back(f);
    }
    // Arm B: from the crossing outward.
    const int nb = static_cast<int>(arm_len_b / STEP);
    for (int k = 1; k <= nb; k++) {
        Fit f;
        f.point = junction + ub * (k * STEP - STEP / 2) + offset_b;
        fits.push_back(f);
    }
    return fits;
}

} // namespace

TEST_CASE("clus pr47 wide kink: straight through-going track never fires") {
    auto fits = two_arm_track(Vector(1, 0.2, 0.1), Vector(1, 0.2, 0.1),
                              30 * units::cm, 30 * units::cm);
    auto acc = segment_cathode_wide_kink_accepts(fits, 0, 25.0,
                                                 3 * units::cm, 15 * units::cm);
    CHECK(acc.empty());
}

TEST_CASE("clus pr47 wide kink: 35 deg kink at the crossing fires at 25, not 40") {
    const double th = 35.0 * M_PI / 180.0;
    auto fits = two_arm_track(Vector(1, 0, 0),
                              Vector(std::cos(th), std::sin(th), 0),
                              30 * units::cm, 30 * units::cm);
    auto acc25 = segment_cathode_wide_kink_accepts(fits, 0, 25.0,
                                                   3 * units::cm, 15 * units::cm);
    REQUIRE(acc25.size() == 1);
    // Contract: interior index, stepped clear of the |x| < 0.5 cm slab.
    CHECK(acc25[0] > 0);
    CHECK(acc25[0] + 1 < fits.size());
    CHECK(std::abs(fits[acc25[0]].point.x()) >= 0.5 * units::cm);
    // Well below the crossing-adjacent band edge (5 cm veto band).
    CHECK(std::abs(fits[acc25[0]].point.x()) < 2.0 * units::cm);

    auto acc40 = segment_cathode_wide_kink_accepts(fits, 0, 40.0,
                                                   3 * units::cm, 15 * units::cm);
    CHECK(acc40.empty());
}

TEST_CASE("clus pr47 wide kink: pure transverse offset does not fake a kink") {
    // Straight track, downstream arm shifted 1.5 cm in y across the gap --
    // the data-scale cathode mismatch.  A translation does not rotate a PCA
    // axis, so the turn angle stays ~0 and nothing fires.
    auto fits = two_arm_track(Vector(1, 0.2, 0.1), Vector(1, 0.2, 0.1),
                              30 * units::cm, 30 * units::cm,
                              Vector(0, 1.5 * units::cm, 0));
    auto acc = segment_cathode_wide_kink_accepts(fits, 0, 25.0,
                                                 3 * units::cm, 15 * units::cm);
    CHECK(acc.empty());
}

TEST_CASE("clus pr47 wide kink: bends inside the skirt are excluded") {
    // Straight arms beyond 3 cm; corrupt every point within the 3 cm skirt
    // with a large transverse zig.  The windowed PCA never sees them.
    auto fits = two_arm_track(Vector(1, 0, 0), Vector(1, 0, 0),
                              30 * units::cm, 30 * units::cm);
    for (auto& f : fits) {
        if (std::abs(f.point.x()) < 3 * units::cm) {
            f.point += Vector(0, 1.0 * units::cm, 0);
        }
    }
    auto acc = segment_cathode_wide_kink_accepts(fits, 0, 25.0,
                                                 3 * units::cm, 15 * units::cm);
    CHECK(acc.empty());
}

// ---------------------------------------------------------------------------
// doc pr/48 sec 6: segment_wide_turn_angle -- the cathode-free wide-baseline
// PCA turn angle at one fit index (the R2 statistic of the two-end break).
// ---------------------------------------------------------------------------

TEST_CASE("clus pr48 wide turn: straight track reads ~0, bent track reads the bend") {
    auto straight = two_arm_track(Vector(1, 0.2, 0.1), Vector(1, 0.2, 0.1),
                                  30 * units::cm, 30 * units::cm);
    const size_t mid_s = straight.size() / 2;
    CHECK(segment_wide_turn_angle(straight, mid_s, 3 * units::cm, 15 * units::cm) < 3.0);

    const double th = 40.0 * M_PI / 180.0;
    auto bent = two_arm_track(Vector(1, 0, 0),
                              Vector(std::cos(th), std::sin(th), 0),
                              30 * units::cm, 30 * units::cm);
    // The junction sits at the index where arm A ends.
    const size_t junc = static_cast<size_t>((30 * units::cm) / (0.6 * units::cm)) - 1;
    const double ang = segment_wide_turn_angle(bent, junc, 3 * units::cm, 15 * units::cm);
    CHECK(ang > 35.0);
    CHECK(ang < 45.0);
}

TEST_CASE("clus pr48 wide turn: terminus-adjacent index and short arms read 0") {
    const double th = 40.0 * M_PI / 180.0;
    auto bent = two_arm_track(Vector(1, 0, 0),
                              Vector(std::cos(th), std::sin(th), 0),
                              30 * units::cm, 4 * units::cm);
    // Downstream arm only 4 cm: < 3 points beyond the 3 cm skirt.
    const size_t junc = static_cast<size_t>((30 * units::cm) / (0.6 * units::cm)) - 1;
    CHECK(segment_wide_turn_angle(bent, junc, 3 * units::cm, 15 * units::cm) == 0.0);
    // Index 0 / last index have no upstream/downstream window at all.
    CHECK(segment_wide_turn_angle(bent, 0, 3 * units::cm, 15 * units::cm) == 0.0);
    CHECK(segment_wide_turn_angle(bent, bent.size() - 1, 3 * units::cm, 15 * units::cm) == 0.0);
    // Out-of-range index declines.
    CHECK(segment_wide_turn_angle(bent, bent.size(), 3 * units::cm, 15 * units::cm) == 0.0);
}

TEST_CASE("clus pr47 wide kink: short arm cannot fire; angle 0 disables") {
    const double th = 35.0 * M_PI / 180.0;
    // Downstream arm only 4 cm: past the 3 cm skirt there are < 3 points.
    auto fits_short = two_arm_track(Vector(1, 0, 0),
                                    Vector(std::cos(th), std::sin(th), 0),
                                    30 * units::cm, 4 * units::cm);
    CHECK(segment_cathode_wide_kink_accepts(fits_short, 0, 25.0,
                                            3 * units::cm, 15 * units::cm).empty());
    // angle_cut <= 0 is the OFF switch, even for a hard kink.
    auto fits_kink = two_arm_track(Vector(1, 0, 0),
                                   Vector(std::cos(th), std::sin(th), 0),
                                   30 * units::cm, 30 * units::cm);
    CHECK(segment_cathode_wide_kink_accepts(fits_kink, 0, 0.0,
                                            3 * units::cm, 15 * units::cm).empty());
}

// doc sbnd_xin/docs/pr/72 round 2 -- es3_stub_suppress, the pure predicate
// examine_structure_3 (NeutrinoStructureExaminer.cxx) uses, when
// es3_stub_guard is on, to decline merging a genuine near-vertex track
// stub into a shower/track trunk.  Numbers below are the real 18255-196649
// measurements (worktree pinned at efef4535): a 6.28cm stub, terminal
// (deg_short=1) at the true neutrino vertex, meeting a 33.19cm trunk at
// ang10=16.450deg / ang3=20.993deg.

TEST_CASE("clus pr72 es3 stub guard: 196649 shape suppresses with default params") {
    Es3StubGuardParams p;  // shipped defaults (117-evt grid-scan pick): 7cm / 2.0 / 15deg / 1.0 / terminal-required
    CHECK(es3_stub_suppress(6.28 * units::cm, 33.19 * units::cm, 20.993, 16.450,
                             /*deg_short=*/1, /*nfit_short=*/11, /*nfit_long=*/56, p));
}

TEST_CASE("clus pr72 es3 stub guard: non-terminal short arm is NOT suppressed") {
    // Same shape as 196649, but the short arm's far end has degree 2 (it
    // chains further into the graph instead of ending at a free vertex
    // candidate) -- the real near-miss pattern seen at evt235435/423981 in
    // the 117-event census, which the terminal-degree requirement exists to
    // spare.
    Es3StubGuardParams p;
    CHECK_FALSE(es3_stub_suppress(6.28 * units::cm, 33.19 * units::cm, 20.993, 16.450,
                                   /*deg_short=*/2, 11, 56, p));
    // With require_terminal explicitly off, the same non-terminal junction
    // DOES suppress -- proves the toggle acts in both directions, not just
    // "off never fires".
    p.require_terminal = false;
    CHECK(es3_stub_suppress(6.28 * units::cm, 33.19 * units::cm, 20.993, 16.450,
                             2, 11, 56, p));
}

TEST_CASE("clus pr72 es3 stub guard: degenerate arm (nfit<2) never suppresses, either side") {
    // segment_track_length returns 0.0 for an arm with < 2 fit points
    // (PRSegmentFunctions.cxx), which would otherwise read as "the short
    // arm" and divide len_long/len_short by (near-)zero.  The guard must
    // catch degeneracy on EITHER arm, not just the short one.
    Es3StubGuardParams p;
    CHECK_FALSE(es3_stub_suppress(0.0, 33.19 * units::cm, 20.993, 16.450,
                                   1, /*nfit_short=*/1, 56, p));
    CHECK_FALSE(es3_stub_suppress(6.28 * units::cm, 33.19 * units::cm, 20.993, 16.450,
                                   1, 11, /*nfit_long=*/1, p));
    // Below the absolute length floor (0.5cm), even with nfit>=2 on both
    // sides: a near-zero-length "stub" is a fit artifact, not a real track.
    CHECK_FALSE(es3_stub_suppress(0.2 * units::cm, 33.19 * units::cm, 20.993, 16.450,
                                   1, 3, 56, p));
}

TEST_CASE("clus pr72 es3 stub guard: threshold boundaries") {
    Es3StubGuardParams p;  // stub_max=7cm len_ratio=2.0 ang3_min=15deg ang_ratio=1.0

    // stub_max: strictly less-than.
    CHECK(es3_stub_suppress(6.99 * units::cm, 33.19 * units::cm, 20.993, 16.450, 1, 11, 56, p));
    CHECK_FALSE(es3_stub_suppress(7.01 * units::cm, 33.19 * units::cm, 20.993, 16.450, 1, 11, 56, p));

    // len_ratio: strictly greater-than (long/short).
    CHECK(es3_stub_suppress(6.28 * units::cm, 6.28 * 2.01 * units::cm, 20.993, 16.450, 1, 11, 56, p));
    CHECK_FALSE(es3_stub_suppress(6.28 * units::cm, 6.28 * 1.99 * units::cm, 20.993, 16.450, 1, 11, 56, p));

    // ang3_min: strictly greater-than.
    CHECK(es3_stub_suppress(6.28 * units::cm, 33.19 * units::cm, 15.01, 10.0, 1, 11, 56, p));
    CHECK_FALSE(es3_stub_suppress(6.28 * units::cm, 33.19 * units::cm, 14.99, 10.0, 1, 11, 56, p));

    // ang_ratio: ang3 must strictly exceed ang_ratio * ang10.
    CHECK(es3_stub_suppress(6.28 * units::cm, 33.19 * units::cm, 20.0, 19.9, 1, 11, 56, p));
    CHECK_FALSE(es3_stub_suppress(6.28 * units::cm, 33.19 * units::cm, 20.0, 20.1, 1, 11, 56, p));
}
