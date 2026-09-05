// doc pdvd/40 round 3 -- unit tests for the pure decision core of
// ImproveCluster_1::remove_bad_blobs (BadBlobRuns.h).  The caller's adjacency
// and support tests are exercised on events (arms d41base / d41fixNN); the
// vote and the run bound are exercised here on synthetic index graphs.

#include "WireCellClus/BadBlobRuns.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <algorithm>

using namespace WireCell;
using namespace WireCell::Clus::BadBlobRuns;

namespace {
    // A chain of n blobs along z, 1 cm apart, linked i-(i+1).
    struct Chain {
        int n;
        std::vector<std::pair<int, int>> edges;
        std::vector<Point> centers;
        std::vector<bool> supported;
        std::vector<int> slice;
        explicit Chain(int n_) : n(n_), supported(n_, true), slice(n_)
        {
            for (int i = 0; i < n; ++i) {
                centers.emplace_back(0, 0, i * units::cm);
                slice[i] = i;
                if (i) edges.emplace_back(i - 1, i);
            }
        }
    };
}

TEST_CASE("bad blob runs: (i) a 30 cm unsupported middle in a supported chain")
{
    Chain c(41);                          // 0..40 cm
    for (int i = 5; i <= 35; ++i) c.supported[i] = false;   // 31 blobs, span 30 cm
    auto r20 = analyze(c.n, c.edges, c.supported, c.centers, 20 * units::cm, c.slice);
    CHECK(r20.ncomp == 1);
    CHECK(r20.removed_by_vote.empty());
    REQUIRE(r20.runs.size() == 1);
    CHECK(r20.runs[0].span == doctest::Approx(30 * units::cm));
    CHECK(r20.runs[0].nslices == 31);
    CHECK(r20.removed_by_run.size() == 31);
    CHECK(r20.removed_by_run.front() == 5);
    CHECK(r20.removed_by_run.back() == 35);
    auto r40 = analyze(c.n, c.edges, c.supported, c.centers, 40 * units::cm, c.slice);
    CHECK(r40.removed_by_run.empty());
    // The historical vote sees one component and removes nothing (hole a).
    CHECK(legacy_component_vote(c.n, c.edges, c.supported).empty());
}

TEST_CASE("bad blob runs: (ii) a detached unsupported component goes regardless of bound")
{
    Chain c(10);
    c.edges.erase(std::remove(c.edges.begin(), c.edges.end(), std::make_pair(4, 5)), c.edges.end());
    for (int i = 5; i < 10; ++i) c.supported[i] = false;   // component 5..9, span 4 cm
    auto r = analyze(c.n, c.edges, c.supported, c.centers, 0.0, c.slice);
    CHECK(r.ncomp == 2);
    CHECK(r.removed_by_vote == std::vector<int>{5, 6, 7, 8, 9});
    CHECK(r.removed_by_run.empty());
    CHECK(r.runs.empty());                 // removed components contribute no runs
    CHECK(legacy_component_vote(c.n, c.edges, c.supported) == r.removed_by_vote);
}

TEST_CASE("bad blob runs: (iii) and (v) bound 0 equals the legacy vote")
{
    Chain c(20);
    for (int i = 3; i < 17; ++i) c.supported[i] = false;    // 13 cm unsupported, one component
    auto r = analyze(c.n, c.edges, c.supported, c.centers, 0.0, c.slice);
    CHECK(r.removed_by_vote.empty());
    CHECK(r.removed_by_run.empty());
    REQUIRE(r.runs.size() == 1);           // still reported, for the census
    CHECK(r.runs[0].span == doctest::Approx(13 * units::cm));
    CHECK(legacy_component_vote(c.n, c.edges, c.supported).empty());
    auto r10 = analyze(c.n, c.edges, c.supported, c.centers, 10 * units::cm, c.slice);
    CHECK(r10.removed_by_run.size() == 14);
}

TEST_CASE("bad blob runs: (iv) the vote is the historical any-blob component vote, unchanged")
{
    Chain c(10);
    c.edges.erase(std::remove(c.edges.begin(), c.edges.end(), std::make_pair(4, 5)), c.edges.end());
    // component A = 0..4 with blob 0 unsupported; component B = 5..9 all unsupported
    c.supported[0] = false;
    for (int i = 5; i < 10; ++i) c.supported[i] = false;
    // The historical loop tests A's blobs in order until one is supported
    // (blob 1), so A is kept and only B goes -- a first-blob reading of that
    // loop is wrong, see BadBlobRuns.h.
    CHECK(legacy_component_vote(c.n, c.edges, c.supported) == std::vector<int>{5, 6, 7, 8, 9});
    auto r = analyze(c.n, c.edges, c.supported, c.centers, 20 * units::cm, c.slice);
    CHECK(r.removed_by_vote == legacy_component_vote(c.n, c.edges, c.supported));
    CHECK(r.removed_by_run.empty());
    REQUIRE(r.runs.size() == 1);           // blob 0: a 1-blob run of span 0 inside kept A
    CHECK(r.runs[0].blobs == std::vector<int>{0});
    CHECK(r.runs[0].span == doctest::Approx(0.0));
    // A component whose ONLY supported blob is its last one is still kept.
    for (int i = 0; i < 5; ++i) c.supported[i] = (i == 4);
    CHECK(legacy_component_vote(c.n, c.edges, c.supported) == std::vector<int>{5, 6, 7, 8, 9});
}

TEST_CASE("bad blob runs: same-slice adjacency joins a column into one run")
{
    // 6 blobs in ONE slice spread 25 cm along y; without same-slice edges each
    // is its own run of span 0 and no bound fires.
    const int n = 6;
    std::vector<Point> centers;
    std::vector<bool> sup(n, false);
    std::vector<int> slice(n, 7);
    for (int i = 0; i < n; ++i) centers.emplace_back(0, i * 5 * units::cm, 0);
    std::vector<std::pair<int, int>> none;
    auto r0 = analyze(n, none, sup, centers, 20 * units::cm, slice);
    CHECK(r0.ncomp == 6);
    CHECK(r0.removed_by_vote.size() == 6);   // every lone blob unsupported => voted out
    std::vector<std::pair<int, int>> chain;
    for (int i = 1; i < n; ++i) chain.emplace_back(i - 1, i);
    // anchor with one supported blob so the component survives the vote
    sup[0] = true;
    auto r1 = analyze(n, chain, sup, centers, 20 * units::cm, slice);
    CHECK(r1.ncomp == 1);
    REQUIRE(r1.runs.size() == 1);
    CHECK(r1.runs[0].span == doctest::Approx(20 * units::cm));
    CHECK(r1.runs[0].nslices == 1);
    CHECK(r1.removed_by_run.empty());        // 20 is not > 20
    auto r2 = analyze(n, chain, sup, centers, 19 * units::cm, slice);
    CHECK(r2.removed_by_run.size() == 5);
}

TEST_CASE("bad blob runs: runs are reported longest first")
{
    Chain c(30);
    for (int i = 2; i <= 4; ++i) c.supported[i] = false;      // 2 cm
    for (int i = 10; i <= 20; ++i) c.supported[i] = false;    // 10 cm
    auto r = analyze(c.n, c.edges, c.supported, c.centers, 0.0, c.slice);
    REQUIRE(r.runs.size() == 2);
    CHECK(r.runs[0].span == doctest::Approx(10 * units::cm));
    CHECK(r.runs[1].span == doctest::Approx(2 * units::cm));
    CHECK(r.runs[0].center.z() == doctest::Approx(15 * units::cm));
}
