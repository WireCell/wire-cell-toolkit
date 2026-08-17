// doc pr/64 round 4: the S6 W-plane long-track exception predicate of the
// "relaxed_strict_img_2d_rescue_long_wtrack" graph flavor.
//
// Every named case below is a REAL candidate from the owner's 899-label
// hand scan (docs/pr/pr57r6-truth.tsv) or from the pr/64 1000-event PR-data
// census, with feature values measured by
// sbnd_xin/scripts/analysis/pr64/wgate_sweep.py.  The tightened operating
// point (Tmax < 1.7 cm OR ab_global < 6 deg on top of the base R2w cuts) is
// the owner's choice: it protects evt122660's good pair at the cost of the
// marginal 64959/172656 recoveries.  If a constant in
// Graphs::two_d_w_track_ok changes, these cases say which owner verdicts
// flip.

#include "WireCellClus/Graphs.h"

#include "WireCellUtil/doctest.h"

using WireCell::Clus::Graphs::S6WTrackInput;
using WireCell::Clus::Graphs::two_d_w_track_ok;

namespace {

struct Build {
    S6WTrackInput in;
    Build()
    {
        // A killed candidate whose sole voting plane is W, with both
        // features measured -- each case then sets what it pins.
        in.w_gap = true;
        in.w_sole_vote = true;
    }
    Build& dis(double v) { in.dis_cm = v; return *this; }
    Build& wgap(bool v) { in.w_gap = v; return *this; }
    Build& sole(bool v) { in.w_sole_vote = v; return *this; }
    Build& np(int v) { in.npmin = v; return *this; }
    Build& lmin(double v) { in.lmin_cm = v; return *this; }
    Build& tmax(double v) { in.tmax_cm = v; return *this; }
    Build& ab(double v) { in.ab_global_deg = v; return *this; }
    Build& dead_w(int v) { in.dead_w = v; return *this; }
    operator const S6WTrackInput&() const { return in; }
};

}  // namespace

TEST_CASE("s6 w-track: guards and sentinels")
{
    // Default-constructed input (all sentinels): never revives.  This is
    // also the documentation no-op case -- with w_track_excuse false the
    // predicate is never even consulted.
    CHECK(!two_d_w_track_ok(S6WTrackInput{}));
    // Not the sole voting plane: a perfect track pair stays killed -- an
    // induction plane's unexcused verdict is never overridden here.
    CHECK(!two_d_w_track_ok(Build().sole(false).dis(1.0).np(1000).lmin(50).tmax(0.5).ab(1)));
    // Unmeasured transverse RMS (tmax sentinel 1e9) fails the thinness cut.
    CHECK(!two_d_w_track_ok(Build().dis(1.0).np(1000).lmin(50).ab(1)));
    // Unmeasured global axis (ab sentinel 90 deg) fails the angle cut.
    CHECK(!two_d_w_track_ok(Build().dis(1.0).np(1000).lmin(50).tmax(0.5)));
}

TEST_CASE("s6 w-track: the three pr/64 target events revive")
{
    // evt174224 1-2 (owner-labelled bad): the plain 6-wire W-only charge
    // hole in an 84 cm track -- the motivating case the local-axis rescue
    // missed (local 20.7 deg vs global 4.2 deg).
    CHECK(two_d_w_track_ok(Build().dis(1.65).np(1295).lmin(84.4).tmax(1.38).ab(4.2)));
    // evt276836 1-3 (census): drift-parallel track, all planes prolonged,
    // only W convicts.  Beyond the 5 cm rescue population on purpose.
    CHECK(two_d_w_track_ok(Build().dis(5.97).np(191).lmin(21.2).tmax(0.38).ab(10.1)));
    // evt276836 6-10 (census): shorter piece of the same track.
    CHECK(two_d_w_track_ok(Build().dis(1.0).np(59).lmin(7.0).tmax(0.49).ab(3.8)));
    // evt314507 4-5 (census): the original pr/64 W-gap exhibit.
    CHECK(two_d_w_track_ok(Build().dis(2.0).np(77).lmin(14.0).tmax(0.42).ab(0.9)));
    // Two more owner-labelled bad pairs the gate recovers.
    CHECK(two_d_w_track_ok(Build().dis(1.0).np(130).lmin(15.9).tmax(0.62).ab(0.8)));   // evt60017 0-1
    CHECK(two_d_w_track_ok(Build().dis(1.0).np(234).lmin(26.9).tmax(1.17).ab(6.5)));   // evt407280 2-3
    // evt287517 0-1: fat-ish (tmax 1.85 >= 1.7) but ultra-collinear -- the
    // A_tight branch of the tightening is what keeps it.
    CHECK(two_d_w_track_ok(Build().dis(1.0).np(588).lmin(32.1).tmax(1.85).ab(5.9)));
}

TEST_CASE("s6 w-track: owner-labelled good separations stay killed")
{
    // evt122660 13-16 (good): passes every BASE R2w cut (lmin 10.5, np 170,
    // tmax 1.73 < 2, ab 7.6 < 25) -- the tightening exists exactly for it:
    // tmax 1.73 >= 1.7 AND ab 7.6 >= 6.
    CHECK(!two_d_w_track_ok(Build().dis(1.70).np(170).lmin(10.5).tmax(1.73).ab(7.6)));
}

TEST_CASE("s6 w-track: marginal bad pairs deliberately forgone")
{
    // evt64959 0-1 (bad, forgone): tmax 1.97 and ab 21.9 -- inside the base
    // cuts, outside the tightening.  The owner chose zero good regressions
    // over recovering these two.
    CHECK(!two_d_w_track_ok(Build().dis(1.0).np(79).lmin(6.6).tmax(1.97).ab(21.9)));
    // evt172656 0-1 (bad, forgone): tmax exactly 1.70 (not < 1.7), ab 18.4.
    CHECK(!two_d_w_track_ok(Build().dis(1.0).np(487).lmin(26.9).tmax(1.70).ab(18.4)));
    // evt137238 0-1 (bad, out of scope since pr/57): fat and kinked.
    CHECK(!two_d_w_track_ok(Build().dis(1.0).np(1431).lmin(45.4).tmax(4.78).ab(48.3)));
    // evt288287 2-3 (bad, out of scope): huge but tmax 14.1 -- a shower.
    CHECK(!two_d_w_track_ok(Build().dis(1.0).np(1554).lmin(114.6).tmax(14.13).ab(9.0)));
}

TEST_CASE("s6 w-track: R2w threshold boundaries")
{
    const auto base = [](){ return Build().dis(1.0).np(500).lmin(20.0).tmax(1.0).ab(3.0); };
    // lmin strictly > 6.0
    CHECK(!two_d_w_track_ok(base().lmin(6.0)));
    CHECK(two_d_w_track_ok(base().lmin(6.01)));
    // npmin >= 50
    CHECK(!two_d_w_track_ok(base().np(49)));
    CHECK(two_d_w_track_ok(base().np(50)));
    // tmax strictly < 2.0 (and 1.9 needs the A_tight branch)
    CHECK(!two_d_w_track_ok(base().tmax(2.0)));
    CHECK(two_d_w_track_ok(base().tmax(1.9)));       // ab 3 < 6 carries it
    CHECK(!two_d_w_track_ok(base().tmax(1.9).ab(10)));  // neither tightening branch
    CHECK(two_d_w_track_ok(base().tmax(1.6).ab(10)));   // T_tight branch
    // ab strictly < 25 (with tmax thin enough for T_tight)
    CHECK(!two_d_w_track_ok(base().ab(25.0)));
    CHECK(two_d_w_track_ok(base().ab(24.9)));
}

TEST_CASE("s6 w-track: R2d dead-W band")
{
    // evt60669 0-1 (bad): W gap through a dead band, np=46 -- revives even
    // though an induction plane also voted (NOT sole-vote-gated; the
    // sole-vote sole-W variant is already handled by the shipped rescue).
    CHECK(two_d_w_track_ok(Build().sole(false).dis(0.68).np(46).lmin(3.5).dead_w(8)));
    // Boundaries: dead_w >= 3, np >= 20, dis < 3.0, and the W gap itself.
    CHECK(two_d_w_track_ok(Build().sole(false).dis(2.9).np(20).dead_w(3)));
    CHECK(!two_d_w_track_ok(Build().sole(false).dis(2.9).np(20).dead_w(2)));
    CHECK(!two_d_w_track_ok(Build().sole(false).dis(2.9).np(19).dead_w(3)));
    CHECK(!two_d_w_track_ok(Build().sole(false).dis(3.0).np(20).dead_w(3)));
    CHECK(!two_d_w_track_ok(Build().sole(false).wgap(false).dis(2.9).np(20).dead_w(3)));
    // dead_w sentinel -1 (unmeasured, dis > 5 cm population) never revives.
    CHECK(!two_d_w_track_ok(Build().sole(false).dis(2.9).np(20).dead_w(-1)));
}
