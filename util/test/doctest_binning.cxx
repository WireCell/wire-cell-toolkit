#include "WireCellUtil/Binning.h"
#include "WireCellUtil/doctest.h"

#include <vector>
#include <cmath>

using namespace WireCell;

TEST_SUITE("binning") {

TEST_CASE("subset basic") {
    Binning a(10, 0, 10);

    SUBCASE("xmin interior, xmax at boundary") {
        auto s = subset(a, 0.5, 10);
        CHECK(s.nbins()   == a.nbins());
        CHECK(s.binsize() == a.binsize());
        CHECK(s.min()     == a.min());
        CHECK(s.max()     == a.max());
    }

    SUBCASE("xmin on edge, xmax at boundary") {
        // xmin=1 lands exactly on edge(1); the returned subset begins at
        // that edge, dropping bin 0.
        auto s = subset(a, 1, 10);
        CHECK(s.nbins()   == a.nbins() - 1);
        CHECK(s.binsize() == a.binsize());
        CHECK(s.min()     == a.min() + 1);
        CHECK(s.max()     == a.max());
    }

    SUBCASE("xmin below range is clamped") {
        auto s = subset(a, -1, 10);
        CHECK(s.nbins()   == a.nbins());
        CHECK(s.binsize() == a.binsize());
        CHECK(s.min()     == a.min());
        CHECK(s.max()     == a.max());
    }

    SUBCASE("both args above range returns empty") {
        Binning b(122, 1.528e+06, 1.772e+06);
        auto s = subset(b, 1791637.8683138976, 1795821.8405835067);
        CHECK(s.nbins() == 0);
    }
}

TEST_CASE("subset contains both bounds") {
    // Bug #473: when xmax falls strictly inside a bin (not on an edge), the old
    // code set hi = floor(bin(xmax)), making the subset end at edge(hi) which
    // is to the *left* of xmax.  xmax was therefore excluded.  The fix must
    // ensure the returned Binning::inside(xmax) is true.
    Binning a(10, 0, 10);
    auto s = subset(a, 0.5, 5.3);
    CHECK(s.inside(0.5));  // xmin is covered
    CHECK(s.inside(5.3));  // xmax must also be covered (fails before the fix)
}

TEST_CASE("gcumulative") {
    CHECK(std::abs(gcumulative( 0,  0,  1) - 0.5) < 1e-6);
    CHECK(std::abs(gcumulative(10,  0,  1) - 1.0) < 1e-6);
    CHECK(std::abs(gcumulative( 5,  5, 30) - 0.5) < 1e-6);
    CHECK(std::abs(gcumulative(300, 5, 30) - 1.0) < 1e-6);
}

TEST_CASE("gbounds") {
    CHECK(std::abs(gbounds(-1.0,  1.0, 0, 1) - 0.682689) < 1e-6);
    CHECK(std::abs(gbounds( 1.0, -1.0, 0, 1) - 0.682689) < 1e-6);  // args swapped: same result
    CHECK(std::abs(gbounds( 0.0, 10.0, 0, 1) - 0.5)      < 1e-6);
    CHECK(std::abs(gbounds(-10.0, 0.0, 0, 1) - 0.5)      < 1e-6);
}

TEST_CASE("gaussian bin integration") {
    Binning bins(20, -10, 10);
    std::vector<double> g(bins.nbins(), 0);
    const double total = gaussian(g.begin(), bins, 0, 1);
    double sum = 0;
    for (auto v : g) sum += v;
    CHECK(std::abs(total - sum) < 1e-6);   // total returned == sum of bin values
    CHECK(std::abs(total - 1.0) < 1e-6);   // near-complete normalisation
}

}  // TEST_SUITE("binning")
