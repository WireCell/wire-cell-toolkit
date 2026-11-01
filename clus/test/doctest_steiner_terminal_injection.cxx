// pr_ml doc 36 sec 2.2: the terminal file of CreateSteinerGraph's terminals_inject_file and the matched-count selection.

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellClus/SteinerTerminalInjection.h"

#include <set>
#include <sstream>

using namespace WireCell;
using WireCell::Clus::Steiner::InjectedTerminals;
using WireCell::Clus::Steiner::keep_best_terminals;
using WireCell::Clus::Steiner::parse_terminal_injection;

namespace {
    auto parse(const std::string& text)
    {
        std::istringstream ss(text);
        return parse_terminal_injection(ss);
    }
}

TEST_CASE("steiner terminal injection: parse")
{
    SUBCASE("records, comments and blank lines")
    {
        const auto inj = parse("# header\n\nT 7 3 0.5\nT 7 1 2 # trailing\nN 7 1\nT -2 0 1\n");
        REQUIRE(inj.size() == 2);
        CHECK(inj.at(7).score.size() == 2);
        CHECK(inj.at(7).score.at(3) == 0.5);
        CHECK(inj.at(7).keep == 1);
        CHECK(inj.at(-2).keep == -1);   // no N record: keep all
    }
    SUBCASE("an empty file holds no cluster, so every cluster keeps the rule") { CHECK(parse("").empty()); }
    SUBCASE("malformed records are refused, not skipped")
    {
        CHECK_THROWS_AS(parse("T 7 3\n"), ValueError);          // no score
        CHECK_THROWS_AS(parse("T 7 -1 1\n"), ValueError);       // negative index
        CHECK_THROWS_AS(parse("T 7 3 1 9\n"), ValueError);      // trailing field
        CHECK_THROWS_AS(parse("T 7 3 1\nT 7 3 2\n"), ValueError);   // listed twice
        CHECK_THROWS_AS(parse("N 7 -1\n"), ValueError);
        CHECK_THROWS_AS(parse("X 7 3 1\n"), ValueError);
    }
}

TEST_CASE("steiner terminal injection: keep the n best")
{
    InjectedTerminals inj;
    inj.score = {{1, 0.1}, {2, 0.9}, {3, 0.9}, {4, 0.5}};
    const std::set<size_t> all = {1, 2, 3, 4};

    SUBCASE("no N record keeps all") { CHECK(keep_best_terminals(all, inj) == all); }
    SUBCASE("n at or above the size keeps all")
    {
        inj.keep = 4;
        CHECK(keep_best_terminals(all, inj) == all);
    }
    SUBCASE("best scores win")
    {
        inj.keep = 3;
        CHECK(keep_best_terminals(all, inj) == std::set<size_t>{2, 3, 4});
    }
    SUBCASE("equal scores break by ascending index")
    {
        inj.keep = 1;
        CHECK(keep_best_terminals(all, inj) == std::set<size_t>{2});
    }
    SUBCASE("zero keeps none")
    {
        inj.keep = 0;
        CHECK(keep_best_terminals(all, inj).empty());
    }
    SUBCASE("a terminal with no score is an error when ranking is needed")
    {
        inj.keep = 1;
        const std::set<size_t> extra = {1, 2, 9};
        CHECK_THROWS_AS(keep_best_terminals(extra, inj), ValueError);
    }
}
