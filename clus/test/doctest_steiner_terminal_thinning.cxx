// Pins the contract of the cross-blob Steiner terminal thinning, doc pdvd/37 R1.
//
// WHAT THE KNOB IS FOR.  find_peak_point_indices runs once per blob
// (find_steiner_terminals), and its map_index_charge holds only that blob's
// points, so the local-maximum test skips every out-of-blob neighbour and can
// never remove a blob's LAST candidate.  Terminal density is therefore floored
// at the candidate-bearing blob density -- 1.02 terminals per such blob,
// measured in doc pdvd/31 round 6 -- which makes the spacing a function of the
// time-slice pitch and the track's drift alignment, not of any physical scale.
// thin_by_min_separation does not touch selection; it filters the selected set
// at a real length.
//
// WHY A FREE FUNCTION.  The Grapher needs detector volumes, a PC transform set
// and a retiler to exist at all, so a Grapher-level test would exercise the
// fixture rather than the rule.  The rule is pure geometry over an ordering the
// caller owns, so it is tested as pure geometry.  What the C++ default is (0,
// i.e. OFF) is pinned separately in doctest_clus_knob_defaults.cxx.
//
// WHAT THIS FILE DOES NOT SAY.  It says nothing about what PDVD runs.  The
// operating point lives in pdvd/wct-pr-perevt.jsonnet (wcp-porting-img) and is
// gated by the compiled-config proof, not here.

#include "WireCellUtil/doctest.h"
#include "WireCellClus/SteinerThinning.h"
#include "WireCellUtil/Units.h"

#include <algorithm>

using namespace WireCell;
using namespace WireCell::Clus::Steiner;

namespace {
    // Terminals laid along x, spaced `pitch`, in admission order 0,1,2,...
    // -- i.e. as if charge decreased monotonically down the track.  This is the
    // shape the knob actually meets: PDVD's production spacing is ~0.5 cm along
    // a straight muon (doc pdvd/37 sec.3.1).
    std::vector<std::pair<Point, size_t>> line(size_t n, double pitch)
    {
        std::vector<std::pair<Point, size_t>> v;
        for (size_t i = 0; i < n; ++i) {
            v.emplace_back(Point{i * pitch, 0.0, 0.0}, i);
        }
        return v;
    }
}

TEST_CASE("steiner thinning: a non-positive separation is an exact pass-through")
{
    // The OFF path is what every byte-identity claim about this branch rests
    // on, so it is pinned as identity of the SEQUENCE, not just of the set.
    const auto in = line(6, 0.5 * units::cm);
    const std::vector<size_t> all{0, 1, 2, 3, 4, 5};

    CHECK(thin_by_min_separation(in, 0.0) == all);
    CHECK(thin_by_min_separation(in, -1.0 * units::cm) == all);

    // ... and an empty input is not a special case.
    CHECK(thin_by_min_separation({}, 0.5 * units::cm).empty());
}

TEST_CASE("steiner thinning: the exclusion is STRICT at the radius")
{
    // Two terminals exactly min_separation apart are BOTH kept.  The guarantee
    // the survivors satisfy is "no two closer than R", not "no two at R".
    // Written so the comparison is exact in IEEE: the same product forms both
    // sides of the `d2 < r2` test.
    const double R = 1.0 * units::cm;

    std::vector<std::pair<Point, size_t>> at_R{
        {Point{0.0, 0.0, 0.0}, 0},
        {Point{R, 0.0, 0.0}, 1},
    };
    CHECK(thin_by_min_separation(at_R, R) == std::vector<size_t>{0, 1});

    // A hair inside and the second one goes.
    std::vector<std::pair<Point, size_t>> inside{
        {Point{0.0, 0.0, 0.0}, 0},
        {Point{0.999 * R, 0.0, 0.0}, 1},
    };
    CHECK(thin_by_min_separation(inside, R) == std::vector<size_t>{0});
}

TEST_CASE("steiner thinning: the first terminal offered can never be suppressed")
{
    // This is the invariant the CALL SITE depends on, and the reason phase 3b
    // sits before phase 4 rather than after it.  create_steiner_tree adds
    // get_extreme_points_for_reference unconditionally AFTER this pass: the
    // extreme points pin the ends of the tree and must not be suppressible by
    // anything.  Were the thinning ever moved below phase 4, an extreme point
    // arriving mid-order would be dropped by exactly the rule below.
    //
    // Offered first == kept, whatever else is on top of it.
    const double R = 1.0 * units::cm;
    std::vector<std::pair<Point, size_t>> pile;
    for (size_t i = 0; i < 20; ++i) {
        pile.emplace_back(Point{0.01 * i * units::cm, 0.0, 0.0}, i);
    }
    const auto kept = thin_by_min_separation(pile, R);
    REQUIRE(kept.size() == 1);
    CHECK(kept[0] == 0);
}

TEST_CASE("steiner thinning: the caller's order decides which member of a crowd survives")
{
    // Same four positions, reversed priority -> the other end survives.  The
    // function carries no tie-break of its own; at the production call site the
    // order is decreasing calc_charge_wcp, so the survivor is the point phase 1
    // itself ranks highest.
    const double R = 1.0 * units::cm;
    const double p = 0.3 * units::cm;

    std::vector<std::pair<Point, size_t>> fwd{
        {Point{0 * p, 0, 0}, 0}, {Point{1 * p, 0, 0}, 1},
        {Point{2 * p, 0, 0}, 2}, {Point{3 * p, 0, 0}, 3},
    };
    auto rev = fwd;
    std::reverse(rev.begin(), rev.end());

    CHECK(thin_by_min_separation(fwd, R) == std::vector<size_t>{0});
    CHECK(thin_by_min_separation(rev, R) == std::vector<size_t>{3});
}

TEST_CASE("steiner thinning: a thinned line keeps its span and its spacing bound")
{
    // The two properties quoted to the owner: thinning can only REMOVE, and it
    // cannot open a hole -- every dropped terminal was within R of a kept one,
    // so a gap of G can grow to at most G + 2R.  On a uniform line of pitch
    // 0.5 cm at R = 0.5 cm nothing is dropped at all (0.5 is not < 0.5); at
    // R = 1.0 cm every other one goes and the survivors are 1.0 cm apart.
    const double pitch = 0.5 * units::cm;
    const auto in = line(11, pitch);   // 0 .. 5 cm

    const auto keep_half = thin_by_min_separation(in, 0.5 * units::cm);
    CHECK(keep_half.size() == 11);

    const auto keep_cm = thin_by_min_separation(in, 1.0 * units::cm);
    CHECK(keep_cm == std::vector<size_t>{0, 2, 4, 6, 8, 10});
    // Span preserved: both ends survive, so the cluster's extent is unchanged.
    CHECK(keep_cm.front() == 0);
    CHECK(keep_cm.back() == 10);
}

TEST_CASE("steiner thinning: negative and large coordinates behave like any other")
{
    // Half of PDVD sits at negative x and y (the detector spans +-339.91 cm in
    // x, +-336.4 in y), and the grid is keyed by floor(coord / R).  A pair
    // straddling the origin, and one far out in the negative corner, must be
    // suppressed exactly as a pair near it would be.
    const double R = 1.0 * units::cm;

    std::vector<std::pair<Point, size_t>> straddle{
        {Point{-0.4 * units::cm, 0, 0}, 0},
        {Point{+0.4 * units::cm, 0, 0}, 1},   // 0.8 cm away -> dropped
        {Point{+1.6 * units::cm, 0, 0}, 2},   // 2.0 cm from #0  -> kept
    };
    CHECK(thin_by_min_separation(straddle, R) == std::vector<size_t>{0, 2});

    std::vector<std::pair<Point, size_t>> corner{
        {Point{-339.91 * units::cm, -336.4 * units::cm, 0.05 * units::cm}, 0},
        {Point{-339.41 * units::cm, -336.4 * units::cm, 0.05 * units::cm}, 1},
        {Point{-337.91 * units::cm, -336.4 * units::cm, 0.05 * units::cm}, 2},
    };
    CHECK(thin_by_min_separation(corner, R) == std::vector<size_t>{0, 2});
}

TEST_CASE("steiner thinning: the guarantee holds pairwise over a dense 3-D cloud")
{
    // The grid is an accelerator, never an arbiter.  Brute-force the invariant
    // over a cloud dense enough that most candidates are suppressed and every
    // one of the 27 neighbour cells gets used: NO two survivors are closer than
    // R, and every dropped point had a survivor within R that was admitted
    // BEFORE it (which is what bounds the hole a removal can open).
    const double R = 0.5 * units::cm;
    std::vector<std::pair<Point, size_t>> cloud;
    size_t id = 0;
    for (int i = -6; i <= 6; ++i) {
        for (int j = -6; j <= 6; ++j) {
            for (int k = -6; k <= 6; ++k) {
                cloud.emplace_back(Point{0.21 * i * units::cm,
                                               0.21 * j * units::cm,
                                               0.21 * k * units::cm}, id++);
            }
        }
    }
    const auto kept = thin_by_min_separation(cloud, R);
    REQUIRE(kept.size() > 1);
    REQUIRE(kept.size() < cloud.size());

    auto pos = [&cloud](size_t i) { return cloud[i].first; };

    for (size_t a = 0; a < kept.size(); ++a) {
        for (size_t b = a + 1; b < kept.size(); ++b) {
            CHECK((pos(kept[a]) - pos(kept[b])).magnitude() >= R);
        }
    }

    std::vector<size_t> rank(cloud.size(), cloud.size());
    for (size_t r = 0; r < kept.size(); ++r) rank[kept[r]] = r;

    for (const auto& [pt, i] : cloud) {
        if (rank[i] != cloud.size()) continue;   // survived
        bool covered = false;
        for (size_t r = 0; r < kept.size() && !covered; ++r) {
            if (kept[r] > i) break;              // admitted after this point
            if ((pt - pos(kept[r])).magnitude() < R) covered = true;
        }
        CHECK(covered);
    }
}
