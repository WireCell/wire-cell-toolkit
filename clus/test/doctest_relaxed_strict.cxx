// doc pr/53 round 6: the "relaxed_strict" graph flavor's kill predicate
// (S1 interior-step denominator + S2 ">= 2" floor, caps unchanged).
//
// The reference cases are the owner-flagged SBND over-clustering pairs
// replayed per-plane in doc pr/53 (sbnd_xin/docs/pr/53_*.md sec 6, sec 16):
// both live pairs (18255-422851, 18255-521075) sit at num_bad1[0]=2 on a
// 3-step neck -- exactly one step short of the legacy "> 2" floor, which is
// arithmetically unreachable below ~3cm because the last 1cm sample is the
// far endpoint itself and good by construction.

#include "WireCellClus/Graphs.h"

#include "WireCellUtil/doctest.h"

using WireCell::Clus::Graphs::relaxed_strict_bad;
using WireCell::Clus::Graphs::relaxed_img_bad;

// The legacy relaxed test, restated locally as the comparison baseline
// (connect_graph_relaxed.cxx:268 shape).
static bool legacy_bad(int nbad, int num_steps, int cap = 7)
{
    return nbad > cap || (nbad > 2 && nbad >= 0.75 * num_steps);
}

TEST_CASE("relaxed_strict kill predicate: sub-3cm blind spot closed")
{
    // 422851 / 521075: 3-step neck (~2.3-2.8cm), 2 of the 2 interior steps
    // bad.  Legacy passes (blind spot); strict kills.
    CHECK(!legacy_bad(2, 3));
    CHECK(relaxed_strict_bad(2, 3));

    // A clean sub-3cm neck must survive both.
    CHECK(!legacy_bad(0, 3));
    CHECK(!relaxed_strict_bad(0, 3));
    CHECK(!relaxed_strict_bad(1, 3));  // 1 of 2 interior bad: below 0.75 ratio? 1 >= 1.5 false
}

TEST_CASE("relaxed_strict kill predicate: 1-interior-step necks survive")
{
    // num_steps=2 (~1-2cm): only one interior sample exists, so a single bad
    // step never reaches the ">= 2" floor -- no new kill possible below the
    // floor's reach, same as legacy.
    CHECK(!legacy_bad(1, 2));
    CHECK(!relaxed_strict_bad(1, 2));
}

TEST_CASE("relaxed_strict kill predicate: long-path behavior unchanged in kind")
{
    // Caps unchanged: > 7 (single counter) and > 9 (sum tests) kill outright.
    CHECK(legacy_bad(8, 100));
    CHECK(relaxed_strict_bad(8, 100));
    CHECK(legacy_bad(10, 100, 9));
    CHECK(relaxed_strict_bad(10, 100, 9));

    // Below cap on a long path the 0.75 ratio governs and S1's one-step
    // denominator shift is negligible: 7 bad of 11 steps survives both
    // (7 < 0.75*11 = 8.25 legacy; 7 < 0.75*10 = 7.5 strict).
    CHECK(!legacy_bad(7, 11));
    CHECK(!relaxed_strict_bad(7, 11));

    // Ratio kill on a mid-length path: 6 bad of 8 steps.
    // legacy: 6 > 2 && 6 >= 0.75*8 = 6 -> kill; strict: 6 >= 0.75*7 = 5.25 -> kill.
    CHECK(legacy_bad(6, 8));
    CHECK(relaxed_strict_bad(6, 8));
}

TEST_CASE("relaxed_strict kill predicate: strict is a superset of legacy")
{
    // Every path legacy kills, strict also kills (no kill is ever relaxed).
    for (int steps = 1; steps <= 20; ++steps) {
        for (int nbad = 0; nbad <= steps; ++nbad) {
            if (legacy_bad(nbad, steps)) {
                CHECK(relaxed_strict_bad(nbad, steps));
            }
            if (legacy_bad(nbad, steps, 9)) {
                CHECK(relaxed_strict_bad(nbad, steps, 9));
            }
        }
    }
}

// doc pr/53 round 7 sec 18: relaxed_img_bad (S5, "3D-image support").  The
// reference cases below are the operating point's own justification --
// quoted from the round-7 offline threshold scan
// (sbnd_xin/scripts/analysis/pr53/threshold_scan.py) run against every
// OC53CENSUS-S closest-pair edge on the 27 round-6 mover events plus the
// owner's round-7 hand-scan events, NOT guessed:
//   - every owner-flagged edge measured max-contiguous-ghost-run 5-10 on an
//     edge under 12cm (269774 j=1,k=3: run=10 dis=11.87cm; 269774
//     j=13,k=14: run=7 dis=8.40cm; 71372 j=1,k=13: run=5 dis=6.06cm
//     (the tightest target); 71372 j=2,k=11: run=7 dis=7.54cm; 463565
//     j=2,k=7: run=7 dis=8.21cm);
//   - a raw ghost-step COUNT or RATIO does not separate these from the
//     surviving-edge background at all (background reaches comparable
//     counts/ratios); the longest CONTIGUOUS run does (background median 0,
//     p90 2 vs owner-target min 5);
//   - the one confirmed false-positive class at long run length is a
//     closest-pair edge whose "nearest points between components" skims the
//     low-density corona of one large blob rather than crossing a genuine
//     gap (evt 52672 j=0,k=2: run=19 but dis=45.89cm, spot-checked visually
//     -- doc pr/53 round 7 sec 18.2) -- excluded by the edge-length cap,
//     which costs nothing on the targets (longest is 11.87cm).
TEST_CASE("relaxed_img kill predicate: owner-flagged target edges all kill")
{
    CHECK(relaxed_img_bad(10, 11.87));  // 269774 j=1 k=3
    CHECK(relaxed_img_bad(7, 8.40));    // 269774 j=13 k=14
    CHECK(relaxed_img_bad(5, 6.06));    // 71372 j=1 k=13 (tightest target)
    CHECK(relaxed_img_bad(7, 7.54));    // 71372 j=2 k=11
    CHECK(relaxed_img_bad(7, 8.21));    // 463565 j=2 k=7
}

TEST_CASE("relaxed_img kill predicate: corona false-positive excluded by the length cap")
{
    // evt 52672 j=0 k=2: 19-step contiguous ghost run would fire the
    // run-floor alone, but the 45.89cm edge is nearly 4x the longest owner
    // target and the spot-check figure shows the "closest pair" running
    // along one large blob's boundary, not a gap between two objects.
    CHECK(!relaxed_img_bad(19, 45.89));
    // The run alone, ignoring length, would have fired -- confirm the cap
    // is doing the work, not the floor.
    CHECK(relaxed_img_bad(19, 14.9));  // same run, under the cap: fires
}

TEST_CASE("relaxed_img kill predicate: short runs and short edges survive")
{
    // Background median run is 0, p90 is 2 -- well under the floor of 4.
    CHECK(!relaxed_img_bad(0, 5.0));
    CHECK(!relaxed_img_bad(2, 5.0));
    CHECK(!relaxed_img_bad(3, 5.0));   // one short of the floor
    CHECK(relaxed_img_bad(4, 5.0));    // exactly at the floor
    // At the length cap boundary (dis_cap_cm is an exclusive "<").
    CHECK(relaxed_img_bad(4, 14.99));
    CHECK(!relaxed_img_bad(4, 15.0));
    CHECK(!relaxed_img_bad(4, 15.01));
}

TEST_CASE("relaxed_img kill predicate: image_check=false is a byte-identical no-op")
{
    // connect_graph_relaxed_strict's image_check parameter defaults to
    // false; when false the S5 block is never entered (max_ghost_run stays
    // 0, relaxed_img_bad is never called), so "relaxed_strict_img" with
    // image_check=false and "relaxed_strict" are the same code path. This
    // is a documentation check, not a call -- relaxed_img_bad(0, dis) never
    // kills regardless of dis (0 < the floor for any floor >= 1), which is
    // consistent with "no ghost steps observed" rather than "test skipped";
    // the actual off-state proof is the C++ if(image_check) guard itself
    // plus the round-7 off-gate (doc pr/53 round 7 sec 18.4).
    CHECK(!relaxed_img_bad(0, 1.0));
    CHECK(!relaxed_img_bad(0, 100.0));
}
