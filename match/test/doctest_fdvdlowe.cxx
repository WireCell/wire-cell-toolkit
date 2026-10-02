// Tests of the FD-VD low-energy matcher's numeric core (FdvdLowE.h): numpy-order
// helpers, the flash grouping (equal-time flashes share a group), the response
// curve and the strict one-to-one decision with its float32 cut semantics.

#include "WireCellUtil/doctest.h"

#include "WireCellMatch/FdvdLowE.h"

#include <cmath>
#include <limits>
#include <vector>

using namespace WireCell::Match;
namespace L = WireCell::Match::FdvdLowE;

TEST_CASE("fdvd lowe np_interp matches numpy.interp")
{
    // expected values: numpy 2.1.1 np.interp(x, xp, fp, left=-7, right=99)
    const std::vector<double> xp{0.0, 1.0, 2.0, 4.0}, fp{0.0, 10.0, 15.0, 16.0};
    CHECK(L::np_interp(-1.0, xp, fp, -7, 99) == -7.0);
    CHECK(L::np_interp(0.0, xp, fp, -7, 99) == 0.0);
    CHECK(L::np_interp(0.5, xp, fp, -7, 99) == 5.0);
    CHECK(L::np_interp(1.0, xp, fp, -7, 99) == 10.0);
    CHECK(L::np_interp(3.0, xp, fp, -7, 99) == 15.5);
    CHECK(L::np_interp(4.0, xp, fp, -7, 99) == 16.0);
    CHECK(L::np_interp(5.0, xp, fp, -7, 99) == 99.0);
}

TEST_CASE("fdvd lowe np_sum is pairwise")
{
    // a[i] = sin(i) 10^(i % 7): numpy 2.1.1 a.sum() (also checked in doctest_fdvdophitfinder)
    std::vector<double> a(300);
    for (size_t i = 0; i < a.size(); ++i) a[i] = std::sin((double) i) * std::pow(10.0, (double) (i % 7));
    CHECK(L::np_sum(a.data(), 300) == 0x1.ab4159e844d66p+20);
}

TEST_CASE("fdvd lowe response")
{
    L::Resp R;
    R.xc = {1.0, 2.0};
    R.y = {0.5, 1.5};
    R.f = {0.2, 0.9};
    R.tail = 0.75;
    CHECK(R(0.5) == 0.0);      // left of the first bin: 0
    CHECK(R(1.5) == 1.0);
    CHECK(R(4.0) == 3.0);      // above the last bin: p * tail
    CHECK(R.fire(0.5) == 0.0);
    CHECK(R.fire(9.0) == 0.9);
}

static std::vector<double> flash_row(int nlit, double pe)
{
    std::vector<double> r(L::NCH, 0.0);
    for (int k = 0; k < nlit; ++k) r[k] = pe;
    return r;
}

TEST_CASE("fdvd lowe flash groups")
{
    L::Constants C;
    // three flashes: two at the same time merge (anchor = brighter), one with 2 lit OpDets is not stored
    std::vector<double> t_ns{100e3, 100e3, 500e3};
    std::vector<double> pe;
    for (auto [n, v] : std::vector<std::pair<int, double>>{{4, 2.0}, {6, 3.0}, {2, 50.0}}) {
        auto r = flash_row(n, v);
        pe.insert(pe.end(), r.begin(), r.end());
    }
    auto G = L::build_groups(t_ns, pe, C);
    REQUIRE(G.t.size() == 1);
    CHECK(G.t[0] == 100.0);
    CHECK(G.nflash[0] == 2);
    CHECK(G.pe[0] == 5.0f);    // 3 + 2 on the shared OpDets
    CHECK(G.pe[5] == 3.0f);
}

TEST_CASE("fdvd lowe strict decision")
{
    L::Calibration cal;
    cal.rwin_lo = 0.5;
    cal.rwin_hi = 2.0;
    const auto arm = L::default_arms()[0];   // P5 E100, Qc 100 ke, KS 0.3, centroid 200 cm
    std::vector<L::Cluster> cl(3);
    cl[0].Q = 150e3;
    cl[1].Q = 150e3;
    cl[2].Q = 90e3;   // below Qc: its pairs never count
    L::Features Fe;
    Fe.gt = {10.0, 20.0};
    Fe.gtot = {500.0, 500.0};
    Fe.gnpd = {10, 10};
    auto pair = [&](int k, int g, float ks) {
        Fe.k.push_back(k);
        Fe.g.push_back(g);
        Fe.r.push_back(1.0f);
        Fe.ks.push_back(ks);
        Fe.dc.push_back(50.0f);
        Fe.npp.push_back(6.0f);
    };
    pair(0, 0, 0.1f);
    pair(1, 0, 0.5f);   // fails KS: does not make group 0 ambiguous
    pair(1, 1, 0.1f);
    pair(2, 1, 0.1f);   // below Qc: does not make group 1 ambiguous
    std::vector<bool> all(Fe.k.size(), true);
    auto d = L::decide(Fe, all, cl, arm, cal);
    REQUIRE(d.size() == 2);
    CHECK(d.at(0) == 0);
    CHECK(d.at(1) == 1);
    // a second passing cluster on group 0 makes both unmatched
    Fe.ks[1] = 0.2f;
    d = L::decide(Fe, std::vector<bool>(Fe.k.size(), true), cl, arm, cal);
    CHECK(d.size() == 0);   // cluster 1 now has two passing groups, group 0 two passing clusters
    // NumPy 2 compares float32 features with float32(cut): ks = float32(0.3) passes KS <= 0.3
    Fe.ks = {0.3f, 0.9f, 0.9f, 0.9f};
    d = L::decide(Fe, std::vector<bool>(Fe.k.size(), true), cl, arm, cal);
    CHECK(d.size() == 1);
    CHECK(d.at(0) == 0);
}
