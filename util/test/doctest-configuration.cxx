#include "WireCellUtil/Logging.h"
#include "WireCellUtil/ConfigurationTesting.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/TimeKeeper.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Point.h"

using spdlog::debug;
using namespace WireCell;
using namespace WireCell::ConfigurationTesting;

TEST_SUITE("util configuration") {

    TEST_CASE("configuration omnibus")
    {
        std::string json = R"(
{
"my_int": 1,
"my_float": 6.9,
"my_string": "hello",
"my_struct": { "x": 1, "y": 2, "z": 3 },
"my_array" : [ "one", "two", "three" ],
"my_array_dict" : [ {"a":1, "b":2}, {"a":10, "b":20}],
"data1" : { "data2": { "a":1, "b":2, "c":3 } }
}
)";
        std::string extra_json = R"(
{
"my_int" : 2,
"data1" : { "data2": { "d":4 } },
"data3": { "data4" : 4 }
}
)";

        Configuration cfg = Persist::loads(json);
        Configuration extra_cfg = Persist::loads(extra_json);

        CHECK(get(cfg, "my_int", 0) == 1);
        CHECK(get(cfg, "my_float", 0.0) == 6.9);

        CHECK(convert<std::string>(cfg["my_string"]) == "hello");

        //cerr << "my_string=" << get(cfg, "my_string", std::string("")) << endl;
        CHECK(get(cfg, "my_string", std::string("")) == "hello");

        put(cfg, "a.b.c", 42);
        // cerr << cfg << endl;
        const int n42get = get(cfg, "a.b.c", 0);
        CHECK(42 == n42get);
        const int n42ind = cfg["a"]["b"]["c"].asInt();
        CHECK(42 == n42ind);

        CHECK(get<int>(cfg, "my_struct.x") == 1);

        auto nums = get<std::vector<std::string> >(cfg, "my_array");
        // for (auto anum : nums) {
        //     cerr << anum << endl;
        // }

        Configuration a = branch(cfg, "data1.data2.a");
        CHECK(convert<int>(a) == 1);

        Configuration other;
        update(other, cfg);
        Configuration last = update(other, extra_cfg);
        // cerr << "other:\n" << other << endl;
        // cerr << "last:\n" << last << endl;
        CHECK(last["a"]["b"]["c"] == 42);
        CHECK(last["data3"]["data4"] == 4);
    }


    TEST_CASE("configuration loads") {
        Configuration cfg = Persist::loads(tracks_json());

        for (auto comp : cfg) {
            CHECK(get<std::string>(comp, "type") == "TrackDepos");
            Configuration data = comp["data"];
            CHECK(get<double>(data, "step_size") == 1.0);
            Configuration tracks = data["tracks"];
            for (auto track : tracks) {
                CHECK(get<double>(track, "charge") < 0);
                /*Ray ray =*/ get<Ray>(track, "ray");
            }
        }
    }

    TEST_CASE("configuration sizes")
    {
        {
            Configuration cfg;
            std::vector<size_t> sizes = {2,4,8};
            assign(cfg, sizes);
            CHECK(cfg.size() == 3);
            for (size_t ind=0; ind<sizes.size(); ++ind) {
                auto c = cfg[(int)ind];
                // std::cerr << ind << " " << c << "\n";
                CHECK(c.asUInt64() == sizes[ind]);
            }
        }
        {
            Configuration cfg;
            std::map<std::string, size_t> named = {{"two",2},{"four",4},{"eight",8}};
            assign(cfg, named);
            CHECK(cfg.size() == 3);
            for (const auto& [nam,siz] : named) {
                auto c = cfg[nam];
                // std::cerr << nam << " " << c << "\n";
                CHECK(c.asUInt64() == siz);
            }
        }
        {
            Configuration cfg;
            std::vector<size_t> sizes = {2,4,8};
            std::map<std::string, std::vector<size_t>> nameds = {{"foo", sizes}, {"bar", {1,2,3}}};
            assign(cfg, nameds);
            CHECK(cfg.size() == 2);
            for (const auto& [nam,sizs] : nameds) {
                auto cs = cfg[nam];
                for (size_t ind=0; ind<sizes.size(); ++ind) {
                    auto c = cs[(int)ind];
                    // std::cerr << nam << " " << ind << " " << c << "\n";
                    CHECK(c.asUInt64() == sizs[ind]);
                }
            }
        }    
    }


    TEST_CASE("configuration hash speed")
    {
        Configuration cfg = Persist::loads(tracks_json());
        TimeKeeper kp("configuration hash speed");

        const int ntimes = 10000;

        std::hash<Configuration> chash;
        std::hash<std::string> shash;

        for (int count=0; count<ntimes; ++count) {
            chash(cfg);
        }
        kp("... std::hash'ed object");
        for (int count=0; count<ntimes; ++count) {
            shash(Persist::dumps(cfg));
        }
        kp("... std::hash'ed json text");

        debug("repeated {} times\n{}", ntimes, kp.summary());

    }

    TEST_CASE("configuration various")
    {

        Configuration cfg;
        cfg["chirp"] = "bad";
        CHECK("bad" == cfg["chirp"].asString());
        cfg["bad"] = "bad";

        Configuration top;
        top["maskmap"] = cfg;

        auto jmm = top["maskmap"];

        std::unordered_map<std::string, std::string> mm;
        for (auto name : jmm.getMemberNames()) {
            mm[name] = jmm[name].asString();
            debug("{} {}", name, mm[name]);
        }

        auto newthing = top["maskmap"]["adc"];
        CHECK(newthing.isNull() == true);

    }


    TEST_CASE("configuration override") {
        Configuration a, b, c;

        a["tag"] = Json::nullValue;
        b["tag"] = "youreit";
        c = a;
        b = update(a,b);
        std::cerr << "a: " << a << "\n";
        std::cerr << "b: " << b << "\n";
        std::cerr << "c: " << c << "\n";
        CHECK(c["tag"] == Json::nullValue);
        CHECK(a["tag"] == b["tag"]);
        CHECK(a["tag"] == "youreit");
    }

}
