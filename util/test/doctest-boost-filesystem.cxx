#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"
#include <boost/filesystem.hpp>

using namespace boost::filesystem;

TEST_SUITE("boost filesystem") {

TEST_CASE("executable path operations") {
    // Use /proc/self/exe as a known-to-exist regular file (Linux)
    path me("/proc/self/exe");
    REQUIRE(exists(me));
    REQUIRE(is_regular_file(me));
    auto p = me.parent_path();
    auto f = me.filename();
    auto me2 = p / f;
    REQUIRE(me == me2);
    spdlog::debug("{}", me.string());
    spdlog::debug("{}", absolute(me).string());
    spdlog::debug("{}", canonical(me).string());
}

TEST_CASE("nonexistent path operations") {
    path you = "path/to/you/";
    spdlog::debug("{}", you.string());
    spdlog::debug("{}", absolute(you).string());

    bool caught = false;
    try {
        spdlog::debug("{}", canonical(you).string());
    }
    catch (const filesystem_error& err) {
        caught = true;
        spdlog::debug("no canonical version of {}", you.string());
    }
    CHECK(caught);

    auto combined = you / "another";
    CHECK(combined.string() == "path/to/you/another");
    spdlog::debug("{}", combined.string());
}

}  // TEST_SUITE("boost filesystem")
