// doc pr/62: the S7 "long-edge corridor connectivity" kill predicate of the
// "relaxed_strict_img_2d_rescue_long" / "..._long2" graph flavors.
//
// S7 covers exactly the distance band S6 (two_d_connectivity_bad, doc pr/56)
// declines to judge: candidates at or above s6_dis_cap (30cm). Unlike S6,
// this band carries NO hand-scan labels yet -- doc pr/62 measured only that
// on evt 142421 cluster 7, 14 of 210 candidate pairs survive S1-S3 at
// >=30cm, three of which (30.03, 36.87, 37.39cm) are the MST bridges that
// wrongly reconnect a 17-point island. The cases below pin the seam with S6,
// the fail-OPEN budget posture (the one property this predicate must never
// silently lose, since S6's own breaker fails CLOSED and it would be an easy
// "harmonization" mistake to copy that here), and the structural properties
// (no-vote sentinels, W-never-excused, monotonicity, strictness parameter)
// rather than fitted physics thresholds -- there is nothing to fit yet.

#include "WireCellClus/Graphs.h"

#include "WireCellUtil/doctest.h"

using WireCell::Clus::Graphs::S7CorridorInput;
using WireCell::Clus::Graphs::long_corridor_bad;

namespace {

struct Build {
    S7CorridorInput in;
    Build()
    {
        for (int p = 0; p < 3; ++p) in.has_plane[p] = true;
    }
    Build& dis(double v) { in.dis_cm = v; return *this; }
    Build& gap(bool u, bool v, bool w)
    {
        in.gap[0] = u; in.gap[1] = v; in.gap[2] = w; return *this;
    }
    Build& budget(bool u, bool v, bool w)
    {
        in.budget_hit[0] = u; in.budget_hit[1] = v; in.budget_hit[2] = w; return *this;
    }
    Build& has(bool u, bool v, bool w)
    {
        in.has_plane[0] = u; in.has_plane[1] = v; in.has_plane[2] = w; return *this;
    }
    Build& excuse(bool u, bool v) { in.excuse_u = u; in.excuse_v = v; return *this; }
    Build& gap_cm(double u, double v, double w)
    {
        in.gap_cm[0] = u; in.gap_cm[1] = v; in.gap_cm[2] = w; return *this;
    }
    operator const S7CorridorInput&() const { return in; }
};

}  // namespace

TEST_CASE("long_corridor_bad: seam with S6 -- identically false below 30cm")
{
    // All three planes gapped, nothing excused: would kill decisively at or
    // above the seam, so this isolates the floor itself.
    for (const double d : {0.5, 5.0, 15.0, 29.9, 29.99, 29.999}) {
        CHECK(!long_corridor_bad(Build().dis(d).gap(true, true, true)));
    }
    // Exactly at and just above the shipped seam (must equal
    // connect_graph_relaxed_strict.cxx's s6_dis_cap / s7_dis_floor_cm --
    // 30.0cm -- so S6 and S7 partition the distance axis with no gap and no
    // overlap).
    CHECK(long_corridor_bad(Build().dis(30.0).gap(true, false, false)));
    CHECK(long_corridor_bad(Build().dis(30.03).gap(true, false, false)));  // the measured 142421 bridge
}

TEST_CASE("long_corridor_bad: fails OPEN on an exhausted search budget")
{
    // A named regression guard: S6's own breaker (s6_planes_connected)
    // fails CLOSED (exhaustion == kill) because at S6's <30cm range that is
    // the safe default. S7 must NOT be "harmonized" to match -- at 30cm+ an
    // exhausted budget on a large sparse real object must abstain, not kill.
    CHECK(!long_corridor_bad(Build().dis(40.0).gap(false, false, true).budget(false, false, true)));
    // Clearing budget_hit with the same gap makes it a real, decisive kill.
    CHECK(long_corridor_bad(Build().dis(40.0).gap(false, false, true).budget(false, false, false)));
}

TEST_CASE("long_corridor_bad: a plane with no seeds on one side never votes")
{
    // No plane has data at all -> no votes -> survives regardless of gap[].
    CHECK(!long_corridor_bad(Build().dis(40.0).has(false, false, false).gap(true, true, true)));
    // Only U lacks data; V and W both have data but neither is gapped -> survives.
    CHECK(!long_corridor_bad(Build().dis(40.0).has(false, true, true).gap(true, false, false)));
    // Only U lacks data; W has data and IS gapped -> still kills via W.
    CHECK(long_corridor_bad(Build().dis(40.0).has(false, true, true).gap(true, false, true)));
}

TEST_CASE("long_corridor_bad: W is never excused; U/V are")
{
    // W gap alone, no excuse possible for W -- kills.
    CHECK(long_corridor_bad(Build().dis(40.0).gap(false, false, true)));
    // U gap, excused -> survives (nothing else gapped).
    CHECK(!long_corridor_bad(Build().dis(40.0).gap(true, false, false).excuse(true, false)));
    // V gap, excused -> survives.
    CHECK(!long_corridor_bad(Build().dis(40.0).gap(false, true, false).excuse(false, true)));
    // Both U and V excused, W gapped -> still killed.
    CHECK(long_corridor_bad(Build().dis(40.0).gap(true, true, true).excuse(true, true)));
}

TEST_CASE("long_corridor_bad: min_gapped_planes strictness is one argument away")
{
    // Single unexcused U gap: kills at the default (1), survives at 2.
    const S7CorridorInput single = Build().dis(40.0).gap(true, false, false);
    CHECK(long_corridor_bad(single, 1));
    CHECK(!long_corridor_bad(single, 2));
    // Two unexcused gaps (U+V): kills at both operating points.
    const S7CorridorInput both = Build().dis(40.0).gap(true, true, false);
    CHECK(long_corridor_bad(both, 1));
    CHECK(long_corridor_bad(both, 2));
}

TEST_CASE("long_corridor_bad: monotone -- more gaps, fewer excuses, or less abstention never un-kill")
{
    for (int g = 0; g < 8; ++g) {
        const bool gu = g & 1, gv = g & 2, gw = g & 4;
        for (int e = 0; e < 4; ++e) {
            const bool eu = e & 1, ev = e & 2;
            const auto base = Build().dis(40.0).gap(gu, gv, gw).excuse(eu, ev);
            const bool verdict = long_corridor_bad(base);
            // Adding a gap on a currently-clean plane can only add a kill.
            if (!gu) CHECK(verdict <= long_corridor_bad(Build().dis(40.0).gap(true, gv, gw).excuse(eu, ev)));
            if (!gv) CHECK(verdict <= long_corridor_bad(Build().dis(40.0).gap(gu, true, gw).excuse(eu, ev)));
            if (!gw) CHECK(verdict <= long_corridor_bad(Build().dis(40.0).gap(gu, gv, true).excuse(eu, ev)));
            // Withdrawing an excuse on a currently-excused plane can only add a kill.
            if (eu) CHECK(verdict <= long_corridor_bad(Build().dis(40.0).gap(gu, gv, gw).excuse(false, ev)));
            if (ev) CHECK(verdict <= long_corridor_bad(Build().dis(40.0).gap(gu, gv, gw).excuse(eu, false)));
        }
    }
    // Setting budget_hit on a gapped, unexcused plane can only REMOVE a kill
    // (it converts a vote into an abstention), never add one.
    const S7CorridorInput wgap = Build().dis(40.0).gap(false, false, true);
    CHECK(long_corridor_bad(wgap));
    CHECK(!long_corridor_bad(Build().dis(40.0).gap(false, false, true).budget(false, false, true)));
    // Likewise clearing has_plane on the only gapped plane can only remove a kill.
    CHECK(!long_corridor_bad(Build().dis(40.0).gap(false, false, true).has(true, true, false)));
}

TEST_CASE("long_corridor_bad: gap_floor_cm is inactive at its shipped default")
{
    // gap_cm all zero (as if unmeasured) still kills at the default floor
    // (0.0) -- the floor is censused but deliberately unfitted (doc pr/62:
    // no scan yet justifies turning it on).
    const S7CorridorInput z = Build().dis(40.0).gap(false, false, true).gap_cm(0, 0, 0);
    CHECK(long_corridor_bad(z));
    CHECK(long_corridor_bad(z, /*min_gapped_planes=*/1, /*dis_floor_cm=*/30.0, /*gap_floor_cm=*/0.0));
    // Raising the floor above the measured gap length removes the kill.
    CHECK(!long_corridor_bad(z, 1, 30.0, /*gap_floor_cm=*/3.0));
    // A gap length at or above a raised floor still kills.
    const S7CorridorInput big = Build().dis(40.0).gap(false, false, true).gap_cm(0, 0, 5.0);
    CHECK(long_corridor_bad(big, 1, 30.0, 3.0));
}

TEST_CASE("long_corridor_bad: long_check=false is a byte-identical no-op")
{
    // connect_graph_relaxed_strict's long_check parameter defaults to false;
    // when false the S7 block is never entered (long_corridor_bad is never
    // called), so every existing flavor's behavior is untouched by this
    // predicate's mere existence. This is a documentation check, not a call
    // into production code -- the actual off-state proof is the C++
    // if(long_check) guard itself plus the doc pr/62 off-gate.
    CHECK(!long_corridor_bad(Build().dis(0.0)));    // no gaps at all, any distance
    CHECK(!long_corridor_bad(Build().dis(100.0)));
}
