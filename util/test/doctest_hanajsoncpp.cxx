#include "WireCellUtil/doctest.h"
#include "WireCellUtil/HanaJsonCPP.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/Type.h"

#include <string>
#include <variant>
#include <iostream>

using WireCell::type;
using WireCell::Configuration;
using namespace WireCell::HanaJsonCPP;

struct MyConfig {

    int number = 42;
    std::string title = "default title";
    bool isActive = true;
    double weight = 75.5;

};
BOOST_HANA_ADAPT_STRUCT(MyConfig, number, title, isActive, weight);

struct SubConfig {
    int sub_id = 1;
    std::string sub_name = "default_sub";
};
BOOST_HANA_ADAPT_STRUCT(SubConfig, sub_id, sub_name);

struct MainConfig {
    int root_id = 1000;
    std::string root_title = "Master Config";
    SubConfig primary_sub = SubConfig{900, "Primary Default"};
    std::vector<SubConfig> sub_list;
    std::vector<int> simple_list = {1, 2};
};
BOOST_HANA_ADAPT_STRUCT(MainConfig, root_id, root_title, primary_sub, sub_list, simple_list);

// Poor-man's variant.
struct IOS {
    int number = 0;
    std::string name = "";
};
struct ConfigWithVariants {
    
    std::vector<IOS> ios{0};
};
BOOST_HANA_ADAPT_STRUCT(IOS, number, name);
BOOST_HANA_ADAPT_STRUCT(ConfigWithVariants, ios);

struct ConfigWithConfig {
    int number=0;
    Configuration json={};
};
BOOST_HANA_ADAPT_STRUCT(ConfigWithConfig, number, json);


TEST_SUITE("util hana jsoncpp") {

    TEST_CASE("hana jsoncpp vector of fake variant") {

        ConfigWithVariants vov;
        auto jvov = to_json(vov);
        CHECK(jvov["ios"].empty());
        jvov["ios"][0]["number"] = 42;
        jvov["ios"][1]["name"] = "the answer";
        from_json(vov, jvov);
        REQUIRE(vov.ios.size() == 2);
        CHECK(vov.ios[0].number == 42);
        CHECK(vov.ios[1].name == "the answer");
    }

    TEST_CASE("hana jsoncpp basics")
    {
        MyConfig mc;
        Json::Value j = to_json(mc);
        // std::cerr << j << "\n";
        CHECK(j["number"].asInt() == 42);
        CHECK(j["isActive"].asBool() == true);

        auto j2 = j;
        j2["number"] = 69;
        j2["isActive"] = false;
        MyConfig mc2;
        from_json(mc2, j2);
        // std::cerr << to_json(mc2) << "\n";
        CHECK(mc2.number == 69);
        CHECK(mc2.isActive == false);
    }

    TEST_CASE("hana jsoncpp with residual")
    {
        MyConfig mc;
        auto jmc = to_json(mc);
        jmc["extra"] = "this is not in the struct";
        auto jres = from_json(mc, jmc);
        CHECK(jres.isMember("extra"));
        CHECK(jres["extra"].asString() == "this is not in the struct");
        CHECK(! jres.isMember("number"));
        CHECK(! jres.isMember("title"));
        CHECK(! jres.isMember("isActive"));
        CHECK(! jres.isMember("weight"));
    }

    TEST_CASE("hana jsoncpp bigger test") {

        MainConfig cfg;

        cfg.sub_list.push_back(SubConfig{10, "Component A"});
        cfg.sub_list.push_back(SubConfig{20, "Component B"});

        std::cout << "--- Initial Configuration ---\n";
        std::cout << "Primary Sub ID: " << cfg.primary_sub.sub_id << std::endl;
        std::cout << "Sub List size: " << cfg.sub_list.size() << std::endl;

        // --- Serialization ---
        Json::Value json_val = to_json(cfg);
        // ... (rest of main function remains the same)
        Json::StreamWriterBuilder writerBuilder;
        writerBuilder["commentStyle"] = "None";
        writerBuilder["indentation"] = "  ";
        std::string document = Json::writeString(writerBuilder, json_val);
        std::cout << "\n--- Serialized JSON ---\n" << document << std::endl;

        // --- Deserialization Test ---
        Json::CharReaderBuilder readerBuilder;
        Json::Value new_json_val;
        std::string errs;

        std::istringstream stream(R"({
        "root_id": 999,
        "primary_sub": {
            "sub_id": 500,
            "sub_name": "Override Sub"
        },
        "sub_list": [
            {"sub_id": 101, "sub_name": "New Item 1"},
            {"sub_id": 102, "sub_name": "New Item 2"}
        ],
        "simple_list": [10, 20, 30],
        "root_title": "Fully Overridden"
    })");

        bool parse_okay = Json::parseFromStream(readerBuilder, stream, &new_json_val, &errs);
        CHECK(parse_okay);

        MainConfig new_cfg;
        from_json(new_cfg, new_json_val);

        std::cout << "\n--- Deserialized Configuration ---\n";
        std::cout << "Root ID: " << new_cfg.root_id << std::endl;
        std::cout << "Primary Sub ID: " << new_cfg.primary_sub.sub_id << std::endl;
        std::cout << "Sub List size: " << new_cfg.sub_list.size() << std::endl;
        std::cout << "Sub List[0] Name: " << new_cfg.sub_list[0].sub_name << std::endl;
        std::cout << "Simple List[2]: " << new_cfg.simple_list[2] << std::endl;

    }

    TEST_CASE("hana jsoncpp config with config") {
        ConfigWithConfig cwc{9};
        cwc.json["key"] = "value";

        Configuration cfg = to_json(cwc);
        REQUIRE(cfg["number"].isInt());
        CHECK(cfg["number"].asInt() == 9);
        REQUIRE(cfg["json"].isObject());
        CHECK(cfg["json"]["key"].isString());
        CHECK(cfg["json"]["key"].asString() == "value");

        ConfigWithConfig cwc2;
        from_json(cwc2, cfg);
        CHECK(cwc2.number == 9);
        REQUIRE(cwc2.json.isObject());
        REQUIRE(cwc2.json["key"].isString());
        CHECK(cwc2.json["key"].asString() == "value");


    }

}



