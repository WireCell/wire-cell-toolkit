// doc pr/57 round 6: the S6 rescue predicate of the
// "relaxed_strict_img_2d_rescue" graph flavor.
//
// Every case below is a REAL candidate from the owner's 2026-08-10 hand
// scan of the S6 separation displays (doc pr/57 sec 14), with the feature
// values measured by sbnd_xin/scripts/analysis/pr57/oc56_fit.py from the
// work-pr57r4-scan* dumps.  Each rescue branch is pinned by a bad-labelled
// pair it must rescue AND by its nearest good-labelled neighbour it must
// leave killed -- the pairs that actually bracketed the fitted thresholds.
// If a constant in Graphs::two_d_rescue_ok changes, these cases say which
// owner verdicts flip.

#include "WireCellClus/Graphs.h"

#include "WireCellUtil/doctest.h"

using WireCell::Clus::Graphs::S6RescueInput;
using WireCell::Clus::Graphs::two_d_rescue_ok;

namespace {

struct Build {
    S6RescueInput in;
    Build()
    {
        for (int p = 0; p < 3; ++p) in.has_plane[p] = true;
    }
    Build& dis(double v) { in.dis_cm = v; return *this; }
    Build& gap(bool u, bool v, bool w)
    {
        in.gap[0] = u; in.gap[1] = v; in.gap[2] = w; return *this;
    }
    Build& excuse(bool u, bool v) { in.excuse_u = u; in.excuse_v = v; return *this; }
    Build& ab(double v) { in.ab_local_deg = v; return *this; }
    Build& np(int v) { in.npmin = v; return *this; }
    Build& lmin(double v) { in.lmin_cm = v; return *this; }
    Build& close(int u, int v, int w)
    {
        in.close_mx[0] = u; in.close_mx[1] = v; in.close_mx[2] = w; return *this;
    }
    Build& cov(double u, double v, double w)
    {
        in.cov[0] = u; in.cov[1] = v; in.cov[2] = w; return *this;
    }
    Build& slope(double u, double v) { in.slope[0] = u; in.slope[1] = v; return *this; }
    Build& ext(double u, double v) { in.ext_med[0] = u; in.ext_med[1] = v; return *this; }
    Build& ov(double u, double v, double w)
    {
        in.ov[0] = u; in.ov[1] = v; in.ov[2] = w; return *this;
    }
    Build& dead_w(int v) { in.dead_w = v; return *this; }
    operator const S6RescueInput&() const { return in; }
};

}  // namespace

TEST_CASE("s6 rescue: guards")
{
    // No gap at all: nothing to rescue.
    CHECK(!two_d_rescue_ok(Build().dis(1.0).gap(false, false, false).np(1000).lmin(50)));
    // Beyond the 5 cm display population: abstain even on a perfect case.
    CHECK(!two_d_rescue_ok(Build().dis(5.1).gap(false, true, false).ab(2).np(1000)
                               .lmin(50).close(1, 3, 1).cov(1, 1, 1)));
}

TEST_CASE("s6 rescue: dead-W band (owner case a)")
{
    // evt59003 1-2 (bad): U+V gapped, W clean, dead-W band of 6 wires,
    // np=481 -- the SBND mid-TPC dead column distorting both inductions.
    CHECK(two_d_rescue_ok(Build().dis(1.91).gap(true, true, false).ab(5).np(481)
                              .lmin(41.1).close(2, 3, 1).cov(0.91, 0.90, 0.79)
                              .slope(1.6, 1.8).ext(12, 12).dead_w(6)));
    // evt60669 0-1 (bad): W-only gap THROUGH the dead band (cov_w 0.47),
    // np=46 -- the dead branch has no length floor and a 20-point floor.
    CHECK(two_d_rescue_ok(Build().dis(0.68).gap(false, false, true).ab(24).np(46)
                              .lmin(3.5).close(1, 1, 2).cov(1, 1, 0.47).dead_w(8)));
    // evt316553 0-2 (good): also dead_w=6, but U coverage 0.58 -- a REAL
    // hole below the 0.65 dead-branch floor.  Must stay killed.
    CHECK(!two_d_rescue_ok(Build().dis(3.31).gap(true, true, false).ab(3).np(47)
                               .lmin(11.6).close(5, 4, 1).cov(0.58, 0.89, 0.76)
                               .slope(1.6, 0.8).ext(32, 12).dead_w(6)));
    // evt282217-shape (good): dead band but only 15 points -- below the floor.
    CHECK(!two_d_rescue_ok(Build().dis(4.71).gap(false, true, false).ab(46).np(15)
                               .lmin(2.4).close(1, 3, 1).cov(1, 0.83, 1)
                               .slope(4.0, 2.8).ext(16, 16).dead_w(6)));
}

TEST_CASE("s6 rescue: W-robustness -- tiny W closure needs collinearity")
{
    // evt286241 0-1 (bad): W-only gap closing at (2,2), ab=14, covered.
    CHECK(two_d_rescue_ok(Build().dis(0.76).gap(false, false, true).ab(14).np(66)
                              .lmin(6.2).close(1, 1, 2).cov(1, 1, 0.91)));
    // evt61579 3-4 (good): nearly identical -- close_w=2, np=59 -- but
    // ab=20, just outside the 15-degree W collinearity.  Stays killed.
    CHECK(!two_d_rescue_ok(Build().dis(0.69).gap(false, false, true).ab(20).np(59)
                               .lmin(6.6).close(1, 1, 2).cov(1, 0.94, 0.86)));
    // evt60017 0-1 (bad): V+W gap, W closes at (2,2), ab=2, and the V gap
    // is prolonged-signal (slope 6.5 on a >5.5cm/130-point pair).
    CHECK(two_d_rescue_ok(Build().dis(1.21).gap(false, true, true).excuse(true, false)
                              .ab(2).np(130).lmin(15.9).close(1, 3, 2)
                              .cov(0.78, 0.83, 0.84).slope(9.0, 6.5).ext(44, 20)
                              .dead_w(1)));
    // evt280972 11-16 (good): same shape but V slope 2.1 (not prolonged)
    // and V coverage 0.88 (below the 0.90 collinear floor).  Stays killed.
    CHECK(!two_d_rescue_ok(Build().dis(1.56).gap(false, true, true).ab(5).np(863)
                               .lmin(24.5).close(1, 4, 2).cov(1.0, 0.88, 1.0)
                               .slope(2.0, 2.1).ext(20, 32)));
}

TEST_CASE("s6 rescue: direction-consistent substantial pairs (owner case b)")
{
    // evt348471 0-2 (bad): V-only gap, ab=13, full coverage, 472 points.
    CHECK(two_d_rescue_ok(Build().dis(1.66).gap(false, true, false).ab(13).np(472)
                              .lmin(51.4).close(1, 4, 1).cov(1, 1, 1)
                              .slope(2.8, 3.4).ext(12, 46)));
    // evt287621 1-2 (bad): U coverage only 0.59 (contaminated span) but
    // ab=2 and np=68 -- the W-clean track-like branch reaches down to 0.55.
    CHECK(two_d_rescue_ok(Build().dis(1.95).gap(true, false, false).ab(2).np(68)
                              .lmin(21.1).close(5, 1, 1).cov(0.59, 1, 0.96)));
    // evt169598-shape (good): 9-point fragment, ab=82 -- no branch.
    CHECK(!two_d_rescue_ok(Build().dis(1.35).gap(true, false, false).ab(82).np(9)
                               .lmin(1.0).close(3, 1, 1).cov(1, 1, 1)
                               .slope(3.2, 4.7).ext(12, 34)));
}

TEST_CASE("s6 rescue: prolonged-signal tiers")
{
    // evt315849 0-1 (bad): V slope 10.5 at np=85 -- the np>=15 tier.
    CHECK(two_d_rescue_ok(Build().dis(1.39).gap(false, true, false).ab(29).np(85)
                              .lmin(9.3).close(1, 4, 1).cov(1, 0.89, 1)
                              .slope(4.9, 10.5).ext(20, 26)));
    // evt282385 0-2 (bad): 6-point stub, ab=9, V slope 20, full coverage --
    // the micro tier (np>=5, collinear, slope>=15, cov>=0.95).
    CHECK(two_d_rescue_ok(Build().dis(1.43).gap(false, true, false).ab(9).np(6)
                              .lmin(1.9).close(1, 4, 1).cov(1, 1, 1)
                              .slope(2.7, 20.0).ext(12, 52)));
    // evt53361 0-1 (good): slopes 6.6/6.1 on a 4.9cm pair -- below both the
    // np>=15 tier (needs 8) and outside the substantial tier (needs
    // Lmin > 5.5).  Stays killed.
    CHECK(!two_d_rescue_ok(Build().dis(1.79).gap(true, true, false).ab(80).np(61)
                               .lmin(4.9).close(4, 2, 1).cov(1, 1, 1)
                               .slope(6.6, 6.1).ext(16, 20)));
}

TEST_CASE("s6 rescue: big pair, tight gap")
{
    // evt234638 2-6 (bad): np=1247, dis=1.38, V closes at 3, full coverage
    // -- non-collinear (ab=68) but the big-tight branch carries it.
    CHECK(two_d_rescue_ok(Build().dis(1.38).gap(false, true, false).ab(68).np(1247)
                              .lmin(39.1).close(1, 3, 1).cov(0.94, 1, 1)
                              .slope(2.1, 1.3).ext(48, 32).ov(0.3, 0.0, 0.5)));
    // evt348691 0-1 (good): np=1371 and U closes at 3, but dis=2.28 -- the
    // 2 cm tightness bound is exactly what separates this owner-good pair
    // from the five owner-bad pairs the branch was fitted on.
    CHECK(!two_d_rescue_ok(Build().dis(2.28).gap(true, false, false).ab(65).np(1371)
                               .lmin(53.9).close(3, 1, 1).cov(0.91, 1, 1)
                               .slope(2.6, 2.4).ext(36, 34).ov(0.0, 0.5, 0.5)));
}

TEST_CASE("s6 rescue: two-view co-location")
{
    // evt290729 0-1 (bad): stub at a track head -- footprints overlap fully
    // in the connected V and W views; U-side evidence is hopeless (ab=39,
    // cov 0.76, close 5) and only the overlap branch reaches it.
    CHECK(two_d_rescue_ok(Build().dis(1.51).gap(true, false, false).ab(39).np(112)
                              .lmin(9.4).close(5, 1, 1).cov(0.76, 1, 1)
                              .slope(2.2, 2.1).ext(12, 12).ov(0.0, 1.0, 1.0)));
    // Same shape with a W-plane gap: the branch is W-anchored and must NOT
    // fire (U+V overlap alone is common at shower junctions labelled good).
    CHECK(!two_d_rescue_ok(Build().dis(1.51).gap(true, false, true).ab(39).np(112)
                               .lmin(9.4).close(5, 1, 4).cov(0.76, 1, 0.7)
                               .slope(2.2, 2.1).ext(12, 12).ov(0.0, 1.0, 1.0)));
}
