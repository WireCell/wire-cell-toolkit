/** doc pr/82 sec 12.7 -- organize_ps_path must survive an empty path.
 *
 * REPRODUCER FIRST.  Before the guard, SBND data event 54629 (2000-event
 * MCP2025C batch, work-mcp2k-cb0816) died with a hard, deterministic SIGSEGV:
 *
 *   #6  D3Vector<double>::x ()                        util/inc/WireCellUtil/D3Vector.h:108
 *   #8  TrackFitting::organize_ps_path (...)          clus/src/TrackFitting.cxx:2071
 *   #9  TrackFitting::do_single_tracking (...)        clus/src/TrackFitting.cxx:8870
 *   #10 PatternAlgorithms::find_other_segments (...)  clus/src/NeutrinoOtherSegments.cxx:561
 *
 * The failure needs two things to line up, and neither is wrong alone:
 *
 *   1. examine_end_ps_vec DELIBERATELY returns an empty list when the whole
 *      path was drained as face-invalid -- TrackFitting.cxx:1985-1993 says so:
 *      "returning an empty list lets the caller (organize_ps_path) fall back to
 *      the original pts rather than handing back an out-of-detector point".
 *   2. organize_ps_path implements that fallback as `if (ps_vec.size() <= 1)
 *      ps_vec = pts;` -- which silently assumes `pts` is non-empty.  At the
 *      second call site (:8870) `pts` has just been rebuilt from `ptss`, and
 *      `ptss` can come back empty.
 *
 * Both empty => ps_vec empty => ps_vec.front() reads unmapped memory.
 *
 * These cases reach the guard before any member (m_dv, m_pcts, m_grouping) or
 * the segment is touched, which is why they need no fixture -- and it is also
 * why a null segment here is deliberate rather than lazy: if the guard were
 * removed, the empty case would dereference through examine_end_ps_vec and this
 * test would crash the runner instead of failing politely.  That is the
 * intended revert behaviour: this test cannot pass without the fix.
 *
 * ---------------------------------------------------------------------------
 * doc pdvd/30: a SECOND, sibling empty-list read in the same function is NOT
 * covered here, and deliberately so.  examine_end_ps_vec's flag_start block can
 * leave ps_list empty on a NON-empty input (it pops every point that fails
 * is_good_point, and its S1.7 else-branch re-inserts temp_start only when that
 * point has a valid face); the flag_end block below it then read ps_list.back()
 * on that empty list.  Worse than a bare UB read: the resulting sentinel
 * garbage becomes temp_end, and the symmetric S1.7 else-branch push_back()s it
 * whenever contained_by() reports a valid face -- inventing an out-of-detector
 * point.  Now guarded (`flag_end && !ps_list.empty()`).
 *
 * It gets no reproducer in this file because it CANNOT have a fixture-free one:
 * unlike the cases above, draining the list requires m_dv->contained_by(),
 * m_pcts->pc_transform() and m_grouping->is_good_point() to be live, so a null
 * segment and a default-constructed TrackFitting crash long before the drain.
 * A real reproducer needs geometry injected through TrackFittingTestHarness
 * (the doc pr/98 friend seam used by doctest_update_association.cxx) plus a
 * grouping whose is_good_point() returns false everywhere.  That fixture is
 * worth building; it was out of scope for the round that found the bug, and is
 * recorded here rather than papered over with a test that could not fail.
 */
#include "WireCellClus/TrackFitting.h"
#include "WireCellUtil/Units.h"

#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace WireCell::Clus;

TEST_CASE("clus trackfitting organize_ps_path empty path is a no-op")
{
    TrackFitting tf;
    std::shared_ptr<PR::Segment> no_segment;   // never dereferenced: see above

    std::vector<WireCell::Point> pts;          // the evt54629 condition
    REQUIRE(pts.empty());

    tf.organize_ps_path(no_segment, pts, 1.2 * units::cm, 0.6 * units::cm);

    // Survived, and left the (empty) path alone rather than inventing points.
    CHECK(pts.empty());

    // The end_point_limit=0 overload is the one that actually crashed
    // (TrackFitting.cxx:8870); cover it explicitly.
    tf.organize_ps_path(no_segment, pts, 1.2 * units::cm, 0);
    CHECK(pts.empty());
}

TEST_CASE("clus trackfitting examine_end_ps_vec empty path returns empty")
{
    TrackFitting tf;
    std::shared_ptr<PR::Segment> no_segment;

    const std::vector<WireCell::Point> pts;
    auto out = tf.examine_end_ps_vec(no_segment, pts, true, true);
    CHECK(out.empty());

    // flag combinations all take the same early exit; :1949 (front) and :1998
    // (back) are each unguarded, so both flags matter.
    CHECK(tf.examine_end_ps_vec(no_segment, pts, true, false).empty());
    CHECK(tf.examine_end_ps_vec(no_segment, pts, false, true).empty());
    CHECK(tf.examine_end_ps_vec(no_segment, pts, false, false).empty());
}
