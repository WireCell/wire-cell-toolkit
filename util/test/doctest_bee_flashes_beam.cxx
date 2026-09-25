// doc pdvd/119: Bee::Flashes op_beam label.  The array is written only when
// set_beam() is called, so every existing op JSON stays byte-identical.

#include "WireCellUtil/Bee.h"
#include "WireCellUtil/doctest.h"

#include <vector>

using namespace WireCell;

namespace {
    Bee::Flashes two_flashes()
    {
        Bee::Flashes f("protodunevd", "op", 39305, 0, 157312);
        f.append(84.3, {1.0, 2.0}, 3.0, {40}, {1.5, 1.5}, 0);
        f.append(430.6, {10.0, 20.0}, 30.0, {}, {}, 0);
        return f;
    }
}

TEST_CASE("bee flashes op_beam absent unless set")
{
    auto f = two_flashes();
    const auto j = f.asJson();
    CHECK(f.size() == 2);
    CHECK_FALSE(j.isMember("op_beam"));
    // Pre-existing keys are unchanged by the new API's existence.
    CHECK(j.isMember("op_t"));
    CHECK(j.isMember("op_cluster_ids"));
    CHECK_FALSE(j.isMember("op_t1"));
    CHECK_FALSE(j.isMember("op_flash_group"));
}

TEST_CASE("bee flashes op_beam written row-aligned when set")
{
    auto f = two_flashes();
    f.set_beam(std::vector<int>{0, 1});
    const auto j = f.asJson();
    REQUIRE(j.isMember("op_beam"));
    REQUIRE(j["op_beam"].isArray());
    REQUIRE(j["op_beam"].size() == j["op_t"].size());
    CHECK(j["op_beam"][0].asInt() == 0);
    CHECK(j["op_beam"][1].asInt() == 1);
}

TEST_CASE("bee flashes op_beam all-zero means labelled, no beam flash")
{
    auto f = two_flashes();
    f.set_beam(std::vector<int>{0, 0});
    const auto j = f.asJson();
    REQUIRE(j.isMember("op_beam"));
    CHECK(j["op_beam"].size() == 2);
    CHECK(j["op_beam"][0].asInt() + j["op_beam"][1].asInt() == 0);
}
