#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include "WireCellUtil/Testing.h"
#include <boost/filesystem.hpp>
#include <iostream>
using namespace boost::filesystem;

TEST_CASE("util boost filesystem with format")
{
    SUBCASE("absolute") {
        path base("/frame/{ident}");
        path rel("plane/{plane}");
        path tot = base / rel;
        path resolved(fmt::format(tot.string(), fmt::arg("ident", 42), fmt::arg("plane","U")));
        CHECK(resolved.string() == "/frame/42/plane/U");
    }
    SUBCASE("relative") {
        path base("frame/{ident}");
        path rel("plane/{plane}");
        path tot = base / rel;
        path resolved(fmt::format(tot.string(), fmt::arg("ident", 42), fmt::arg("plane","U")));
        CHECK(resolved.string() == "frame/42/plane/U");
    }
    SUBCASE("extra") {
        path base("frame/{ident}");
        path rel("/plane/{plane}");
        path tot = base / rel;
        path resolved(fmt::format(tot.string(), fmt::arg("ident", 42), fmt::arg("plane","U")));
        CHECK(resolved.string() == "frame/42/plane/U");
        std::cerr << resolved << "\n";
    }
    SUBCASE("ending") {
        path base("frame/{ident}/");
        path rel("plane/{plane}");
        path tot = base / rel;
        path resolved(fmt::format(tot.string(), fmt::arg("ident", 42), fmt::arg("plane","U")));
        CHECK(resolved.string() == "frame/42/plane/U");
        std::cerr << resolved << "\n";
    }
}

TEST_CASE("util boost filesystem string")
{
    SUBCASE("absolute") {
        path base("/frame/{ident}");
        std::string rel = "plane/{plane}";
        path tot = base / rel;
        path resolved(fmt::format(tot.string(), fmt::arg("ident", 42), fmt::arg("plane","U")));
        CHECK(resolved.string() == "/frame/42/plane/U");
    }
    SUBCASE("relative") {
        path base("frame/{ident}");
        std::string rel = "plane/{plane}";
        path tot = base / rel;
        path resolved(fmt::format(tot.string(), fmt::arg("ident", 42), fmt::arg("plane","U")));
        CHECK(resolved.string() == "frame/42/plane/U");
    }
    SUBCASE("extra") {
        path base("frame/{ident}");
        std::string rel = "/plane/{plane}";
        path tot = base / rel;
        path resolved(fmt::format(tot.string(), fmt::arg("ident", 42), fmt::arg("plane","U")));
        CHECK(resolved.string() == "frame/42/plane/U");
        std::cerr << resolved << "\n";
    }
    SUBCASE("ending") {
        path base("frame/{ident}/");
        std::string rel = "plane/{plane}";
        path tot = base / rel;
        path resolved(fmt::format(tot.string(), fmt::arg("ident", 42), fmt::arg("plane","U")));
        CHECK(resolved.string() == "frame/42/plane/U");
        std::cerr << resolved << "\n";
    }
}

TEST_CASE("util boost filesystem playground")
{
    path pw("/etc/passwd");
    Assert(exists(pw));
    Assert(is_regular_file(pw));
    auto p = pw.parent_path();
    auto f = pw.filename();
    auto pw2 = p / f;
    Assert(pw == pw2);
    std::cerr << pw << std::endl;
    std::cerr << absolute(pw).string() << std::endl;
    std::cerr << canonical(pw).string() << std::endl;
    
    path you = "path/to/you/";
    std::cerr << you << std::endl;
    std::cerr << absolute(you).string() << std::endl;
    try {
        std::cerr << canonical(you).string() << std::endl;
    }
    catch (const filesystem_error& err) {
        std::cerr << "no canonical version of " << you << std::endl;
    }

    std::cerr << you / "another" << std::endl;
}
