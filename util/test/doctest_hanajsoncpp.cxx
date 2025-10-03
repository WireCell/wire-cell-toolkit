#include "WireCellUtil/doctest.h"
#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellUtil/Type.h"

#include <string>
#include <iostream>

using WireCell::type;
using namespace WireCell::HanaJsonCPP;

struct MyConfig {

    int number = 42;
    std::string title = "default title";
    bool isActive = true;
    double weight = 75.5;

};

BOOST_HANA_ADAPT_STRUCT(MyConfig, number, title, isActive, weight);

TEST_CASE("util hana jsoncpp") {

    MyConfig mc;
    Json::Value j = to_json(mc);
    // std::cerr << j << "\n";
    CHECK(j["number"].asInt() == 42);

    auto j2 = j;
    j2["number"] = 69;
    MyConfig mc2;
    from_json(mc2, j2);
    // std::cerr << to_json(mc2) << "\n";
    CHECK(mc2.number == 69);

}


