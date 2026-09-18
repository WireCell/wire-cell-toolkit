/** doc pdvd/114 -- the blank-plane admission policy for Steiner terminal
    candidates (WireCellClus/SteinerBlankPlane.h) and its knob defaults.

    The policy is a pure function of a blob's candidate set, a per-point
    zero-plane count and a per-point "three-plane candidate nearby" predicate,
    so it is pinned here without a cluster.  What is pinned:
      * the C++ defaults ("wcp", radius 0) and the mode parser (a typo is
        refused, not mapped to the legacy path);
      * "wcp" returns its input unchanged, whatever the predicates say;
      * "prefer3" drops the candidates with a zero plane only when the blob
        holds a three-plane candidate, and leaves an all-two-plane blob whole
        (the gap-jump guarantee);
      * "nearby" drops a two-plane candidate only when the predicate fires,
        and never a three-plane one;
      * "prefer3+nearby" is the conjunction;
      * the survivors keep the caller's ordering (a std::set with
        std::greater, as at the production site).
 */

#include "WireCellClus/SteinerBlankPlane.h"

#include "WireCellUtil/doctest.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"

#include <functional>
#include <map>
#include <set>
#include <vector>

using namespace WireCell::Clus::Steiner;

namespace {
    using cand_t = std::set<std::pair<double, size_t>, std::greater<std::pair<double, size_t>>>;

    cand_t make(std::initializer_list<std::pair<double, size_t>> il)
    {
        cand_t c;
        for (const auto& p : il) c.insert(p);
        return c;
    }

    std::set<size_t> ids(const cand_t& c)
    {
        std::set<size_t> out;
        for (const auto& p : c) out.insert(p.second);
        return out;
    }
}

TEST_CASE("steiner blank plane: defaults and parser")
{
    // the component's defaults: "wcp" / 0 => a config without the keys is byte-identical
    WireCell::PluginManager::instance().add("WireCellClus");
    auto icfg = WireCell::Factory::lookup<WireCell::IConfigurable>("CreateSteinerGraph", "doc114_blank_plane_probe");
    REQUIRE(icfg);
    auto cfg = icfg->default_configuration();
    REQUIRE(cfg.isMember("terminal_blank_plane_mode"));
    CHECK(cfg["terminal_blank_plane_mode"].asString() == "wcp");
    REQUIRE(cfg.isMember("terminal_blank_plane_radius"));
    CHECK(cfg["terminal_blank_plane_radius"].asDouble() == 0.0);

    BlankPlaneMode m = BlankPlaneMode::prefer3;
    CHECK(parse_blank_plane_mode("wcp", m));
    CHECK(m == BlankPlaneMode::wcp);
    CHECK(parse_blank_plane_mode("prefer3", m));
    CHECK(m == BlankPlaneMode::prefer3);
    CHECK(parse_blank_plane_mode("nearby", m));
    CHECK(m == BlankPlaneMode::nearby);
    CHECK(parse_blank_plane_mode("prefer3+nearby", m));
    CHECK(m == BlankPlaneMode::prefer3_nearby);
    // a typo is refused and leaves the mode untouched
    CHECK_FALSE(parse_blank_plane_mode("prefer", m));
    CHECK(m == BlankPlaneMode::prefer3_nearby);
    CHECK_FALSE(parse_blank_plane_mode("", m));
}

TEST_CASE("steiner blank plane: the policy on one blob")
{
    // point ids 1..6: 1,2 three-plane; 3,4,5 one zero plane; 6 two zero planes.
    std::map<size_t, int> nz = {{1, 0}, {2, 0}, {3, 1}, {4, 1}, {5, 1}, {6, 2}};
    auto nzero = [&](size_t i) { return nz.at(i); };
    // the nearby predicate fires for 4 and 6 only
    std::set<size_t> near = {4, 6};
    auto near3 = [&](size_t i) { return near.count(i) > 0; };
    // the one-blank point 3 is the brightest, as at doc 112's h1 vertex 850
    const cand_t blob = make({{44343.0, 3}, {30000.0, 1}, {20000.0, 4}, {12000.0, 2}, {9000.0, 5}, {8000.0, 6}});

    SUBCASE("wcp returns the input unchanged")
    {
        auto out = apply_blank_plane_policy(blob, BlankPlaneMode::wcp, nzero, near3);
        CHECK(out == blob);
    }
    SUBCASE("prefer3 keeps only the three-plane candidates when the blob has any")
    {
        auto out = apply_blank_plane_policy(blob, BlankPlaneMode::prefer3, nzero, near3);
        CHECK(ids(out) == std::set<size_t>{1, 2});
        // ordering preserved: the brightest survivor first
        CHECK(out.begin()->second == 1);
    }
    SUBCASE("prefer3 leaves an all-two-plane blob whole (gap jumps kept)")
    {
        const cand_t gap = make({{5000.0, 3}, {4000.0, 4}, {3000.0, 6}});
        auto out = apply_blank_plane_policy(gap, BlankPlaneMode::prefer3, nzero, near3);
        CHECK(out == gap);
    }
    SUBCASE("nearby drops a zero-plane candidate only where the predicate fires")
    {
        auto out = apply_blank_plane_policy(blob, BlankPlaneMode::nearby, nzero, near3);
        CHECK(ids(out) == std::set<size_t>{1, 2, 3, 5});
    }
    SUBCASE("nearby never drops a three-plane candidate even if the predicate would fire")
    {
        std::set<size_t> all = {1, 2, 3, 4, 5, 6};
        auto fire = [&](size_t i) { return all.count(i) > 0; };
        auto out = apply_blank_plane_policy(blob, BlankPlaneMode::nearby, nzero, fire);
        CHECK(ids(out) == std::set<size_t>{1, 2});
    }
    SUBCASE("prefer3+nearby is the conjunction")
    {
        auto out = apply_blank_plane_policy(blob, BlankPlaneMode::prefer3_nearby, nzero, near3);
        CHECK(ids(out) == std::set<size_t>{1, 2});
        const cand_t gap = make({{5000.0, 3}, {4000.0, 4}, {3000.0, 6}});
        auto out2 = apply_blank_plane_policy(gap, BlankPlaneMode::prefer3_nearby, nzero, near3);
        CHECK(ids(out2) == std::set<size_t>{3});
    }
    SUBCASE("an empty candidate set stays empty under every mode")
    {
        const cand_t none;
        for (auto m : {BlankPlaneMode::wcp, BlankPlaneMode::prefer3, BlankPlaneMode::nearby, BlankPlaneMode::prefer3_nearby}) {
            CHECK(apply_blank_plane_policy(none, m, nzero, near3).empty());
        }
    }
}
