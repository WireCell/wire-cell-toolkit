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
