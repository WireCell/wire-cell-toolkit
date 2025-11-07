#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/doctest.h"

#include <spdlog/spdlog.h>
using spdlog::debug;

TEST_CASE("jsoncpp zero is zero") {

    Json::Value arr = Json::arrayValue;
    const double q = 0;
    CHECK(arr.size() == 0);
    CHECK(arr.empty() == true);
    arr.append(q);
    CHECK(arr.size() == 1);
    CHECK(arr.empty() == false);

    CHECK(! arr[0].isNull());
    CHECK(arr[0].isDouble());

    Json::Value arr2 = Json::arrayValue;
    arr2.append(arr[0]);

    CHECK(arr2.size() == 1);

    CHECK(! arr2[0].isNull());
    CHECK(arr2[0].isDouble());


}

TEST_CASE("jsoncpp object order") {
    Json::Value arr = Json::objectValue;
    arr["z"] = 1;
    arr["y"] = 2;
    arr["x"] = 3;
    arr["a"] = 12;
    arr["b"] = 11;
    arr["c"] = 10;
    for (const auto& name : arr.getMemberNames()) {
        debug("{} {}", name, arr[name].asInt());
    }
}
