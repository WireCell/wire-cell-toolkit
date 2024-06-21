#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Configuration.h"
#include "WireCellUtil/doctest.h"

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
