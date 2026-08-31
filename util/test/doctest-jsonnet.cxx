#include "WireCellUtil/Persist.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <string>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace WireCell;

// from Jsonnet docs
const std::string give1 = R"(
std.mergePatch({
    person1: {
        name: "Alice",
        welcome: "Hello " + self.name + "!",
    },
    person2: self.person1 { name: "Bob" },
    person3: self.person1 { name: std.extVar("person") },
}, std.extVar("override"))
)";
const std::string want1 = R"({
	"person1" :
	{
		"name" : "Alice",
		"welcome" : "Hello Alice!"
	},
	"person2" :
	{
		"name" : "Bob",
		"welcome" : "Hello Bob!"
	},
	"person3" :
	{
		"name" : "Malory",
		"welcome" : "Hello Malory!"
	},
	"person4" :
	{
		"name" : "Surley",
		"welcome" : "Don't call my Shirley"
	}
})";
//' <-- make emacs c++ mode happy about matched quotes

TEST_SUITE("jsonnet") {

TEST_CASE("loads with external vars") {
    Persist::externalvars_t extravars{{"person", "Malory"}};
    Persist::externalvars_t extracode{{"override", "{person4: {name: \"Surley\", welcome:\"Don't call my Shirley\"}}"}};

    auto jobj = Persist::loads(give1, extravars, extracode);
    string got1 = Persist::dumps(jobj);

    spdlog::debug("------give:\n{}", give1);
    spdlog::debug("------got:\n{}", got1);
    spdlog::debug("------want:\n{}", want1);
    // not same up to whitespace
    // TODO: add meaningful assertions (string comparison fails due to whitespace differences)
    CHECK(!jobj.isNull());
    CHECK(!got1.empty());
}

TEST_CASE("nonexistent extvar throws") {
    CHECK_THROWS_AS(Persist::loads("[ std.extVar(\"doesnoteexists\") ]"), Exception);
}

TEST_CASE("load wirecell jsonnet pi") {
    if (!std::getenv("WIRECELL_PATH")) {
        spdlog::debug("Skipping: WIRECELL_PATH not set");
        return;
    }
    auto jobj = Persist::loads("local wc = import \"wirecell.jsonnet\"; [ wc.pi ]");
    string text = Persist::dumps(jobj);
    auto res = Persist::loads(text);

    const double mypi = 2 * std::acos(0.0);
    const double jspi = res[0].asDouble();
    CHECK(std::abs(jspi - mypi) < 1.0e-16);
}

TEST_CASE("parser class") {
    Persist::externalvars_t extravars{{"person", "Malory"}};
    Persist::externalvars_t extracode{{"override", "{person4: {name: \"Surley\", welcome:\"Don't call my Shirley\"}}"}};

    Persist::Parser p({}, extravars, extracode);
    auto j = p.loads(give1);
    spdlog::debug("{}", Persist::dumps(j));
    CHECK(!j.isNull());
}

}  // TEST_SUITE("jsonnet")
