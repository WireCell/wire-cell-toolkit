/** The array-level pieces of the cascade deghosting core (moved from img/test/doctest_cascade_deghosting.cxx with
    the code, pdvd doc 130 phase 2): in-slice pairs, the coverage guards, the iso fallback and the Steiner repair.
    No blobs, no components, no model.
*/
#include "WireCellAux/CascadeGraph.h"
#include "WireCellAux/CellSteiner.h"

#include "WireCellUtil/doctest.h"

#include <cmath>
#include <vector>

using namespace WireCell;
using namespace WireCell::Aux;

TEST_CASE("in-slice pairs: overlap or abut on all planes, within a group")
{
    std::vector<int64_t> g{0, 0, 0, 1};
    std::vector<std::array<int, 6>> r{{0, 4, 0, 4, 0, 4},     // a
                                      {4, 8, 2, 6, 3, 9},     // abuts a on U, overlaps V and W
                                      {9, 12, 0, 4, 0, 4},    // U gap to a (and to b: 8 < 9)
                                      {0, 4, 0, 4, 0, 4}};    // same ranges as a, other group
    auto p = Cascade::inslice_pairs(g, r);
    REQUIRE(p.size() == 1);
    CHECK(p[0][0] == 0);
    CHECK(p[0][1] == 1);
}

TEST_CASE("coverage guard keeps a wire's last explanation")
{
    Cascade::Level lev;
    lev.blobs.resize(3);
    lev.wq = {1, 1};
    lev.wplane = {0, 0};
    // node 0 -> wire 0; node 1 -> wires 0, 1; node 2 -> wire 1
    lev.bw_src = {0, 1, 1, 2};
    lev.bw_dst = {0, 0, 1, 1};
    size_t ng = 0;
    auto pr = Cascade::guard_prune(lev, {1, 0, 2}, ng);
    CHECK(pr[1]);          // first: both wires still have another node
    CHECK(!pr[0]);         // wire 0 would lose its last node
    CHECK(!pr[2]);
    CHECK(ng == 2);
}

TEST_CASE("final guard: every charged wire keeps its best node; uncharged wires are not guarded")
{
    Cascade::Level lev;
    lev.blobs.resize(4);
    lev.wq = {1, 1, 0};
    lev.wplane = {0, 0, 0};
    // node 0 -> wire 0; node 1 -> wires 0, 1; node 2 -> wire 1; node 3 -> wire 2 (no charge)
    lev.bw_src = {0, 1, 1, 2, 3};
    lev.bw_dst = {0, 0, 1, 1, 2};
    std::vector<float> logit{0.5f, 2.0f, 1.0f, -5.0f};
    {
        std::vector<bool> keep(4, false);
        CHECK(Cascade::final_guard(lev, logit, keep) == 1);     // wire 0 -> node 1 (best), which also covers wire 1
        CHECK(keep == std::vector<bool>{false, true, false, false});
    }
    {
        std::vector<bool> keep{false, false, true, false};         // wire 1 covered, wire 0 not
        CHECK(Cascade::final_guard(lev, logit, keep) == 1);
        CHECK(keep == std::vector<bool>{false, true, true, false});
    }
    {
        std::vector<bool> keep{true, false, true, false};          // already covered: nothing added
        CHECK(Cascade::final_guard(lev, logit, keep) == 0);
        CHECK(keep == std::vector<bool>{true, false, true, false});
    }
    {
        std::vector<float> tie{1.0f, 1.0f, 1.0f, 1.0f};              // ties: the lowest node index
        std::vector<bool> keep(4, false);
        CHECK(Cascade::final_guard(lev, tie, keep) == 2);          // wire 0 -> node 0, then wire 1 -> node 1
        CHECK(keep == std::vector<bool>{true, true, false, false});
    }
}

TEST_CASE("iso fallback: dense ambiguous slices keep down to t_keep, others untouched")
{
    // slice 0 face 0: 12 nodes over 2 charged wires (6 per wire), 11 ambiguous -> triggers
    // slice 0 face 1: 12 nodes but only 6 ambiguous -> no trigger
    // slice 1 face 0: 12 ambiguous nodes over 12 charged wires (1 per wire) -> no trigger
    std::vector<int> face, sidx;
    std::vector<float> logit;
    for (int i = 0; i < 12; ++i) { face.push_back(0); sidx.push_back(0); logit.push_back(i == 0 ? 2.0f : -1.0f - 0.1f * i); }
    for (int i = 0; i < 12; ++i) { face.push_back(1); sidx.push_back(0); logit.push_back(i < 6 ? -1.0f : -5.0f); }
    for (int i = 0; i < 12; ++i) { face.push_back(0); sidx.push_back(1); logit.push_back(-1.0f); }
    std::vector<float> wq{10, 20, 0, 0};            // slice 0: 2 charged wires (+2 empty)
    std::vector<int> wsidx{0, 0, 0, 0};
    for (int w = 0; w < 12; ++w) { wq.push_back(5); wsidx.push_back(1); }
    const size_t N = logit.size();
    std::vector<bool> keep0(N, false);
    for (size_t i = 0; i < N; ++i) keep0[i] = logit[i] >= -0.1877f;

    Cascade::IsoParams p;
    p.nmin = 10; p.mmin = 4; p.amin = 0.8; p.t_keep = -1.5;
    auto keep = keep0;
    auto r = Cascade::iso_fallback(face, sidx, wq, wsidx, logit, keep, p);
    CHECK(r.nslices == 1);
    // slice 0 face 0: logits -1.1 .. -2.1; kept those >= -1.5 (i = 1..5), node 0 was kept already
    for (int i = 0; i < 12; ++i) CHECK(keep[i] == (i <= 5));
    CHECK(r.nadded == 5);
    for (size_t i = 12; i < N; ++i) CHECK(keep[i] == keep0[i]);

    // nmin above the group size, or mmin above the density: nothing changes
    for (auto q : {10000.0, 7.0}) {
        auto p2 = p;
        if (q > 100) p2.nmin = (int) q; else p2.mmin = q;
        auto k2 = keep0;
        auto r2 = Cascade::iso_fallback(face, sidx, wq, wsidx, logit, k2, p2);
        CHECK(r2.nslices == 0);
        CHECK(r2.nadded == 0);
        CHECK(k2 == keep0);
    }

    // node order does not matter: reverse every node array
    std::vector<int> rf(face.rbegin(), face.rend()), rs(sidx.rbegin(), sidx.rend());
    std::vector<float> rl(logit.rbegin(), logit.rend());
    std::vector<bool> rk(keep0.rbegin(), keep0.rend());
    auto rr = Cascade::iso_fallback(rf, rs, wq, wsidx, rl, rk, p);
    CHECK(rr.nslices == r.nslices);
    CHECK(rr.nadded == r.nadded);
    for (size_t i = 0; i < N; ++i) CHECK(rk[N - 1 - i] == keep[i]);

    // size mismatch throws
    std::vector<bool> bad(N - 1, false);
    CHECK_THROWS(Cascade::iso_fallback(face, sidx, wq, wsidx, logit, bad, p));
}

TEST_CASE("steiner repair: bridge within budget, not beyond, weak island dropped")
{
    auto chain = [](size_t n) {
        std::vector<std::array<int64_t, 2>> e;
        for (size_t i = 0; i + 1 < n; ++i) e.push_back({(int64_t) i, (int64_t) i + 1});
        return e;
    };
    const float hi = 4.6f, lo = (float) std::log(0.2 / 0.8);   // P = 0.99, P = 0.2
    Cascade::SteinerParams par;   // p_term 0.8, q_floor 1e4, budget 5.5
    {
        // 0 1 [2 3 4] 5 6 : three P = 0.2 cells between two confident fragments; node 7 an isolated weak island
        std::vector<float> logit{hi, hi, lo, lo, lo, hi, hi, 0.4f};
        std::vector<float> qhat(8, 0.1f);
        std::vector<bool> keep{true, true, false, false, false, true, true, true};
        auto e = chain(7);   // node 7 unconnected
        auto res = Cascade::steiner_repair(8, e, logit, qhat, keep, par);
        CHECK(res.nbridges == 1);
        CHECK(res.nadded == 3);
        for (int i = 0; i < 7; ++i) CHECK(res.keep[i]);
        CHECK(!res.keep[7]);
        CHECK(res.nweak_cells == 1);
    }
    {
        // a gap of five P = 0.2 cells costs more than the budget
        std::vector<float> logit{hi, hi, lo, lo, lo, lo, lo, hi, hi};
        std::vector<float> qhat(9, 1.0f);
        std::vector<bool> keep{true, true, false, false, false, false, false, true, true};
        auto res = Cascade::steiner_repair(9, chain(9), logit, qhat, keep, par);
        CHECK(res.nbridges == 0);
        CHECK(res.nadded == 0);
    }
}
