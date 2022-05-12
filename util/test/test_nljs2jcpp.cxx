// provided at top level by WireCellUtil
#include "WireCellUtil/nljs2jcpp.hpp"

#include <cassert>

int main()
{
    nlohmann::json nljs{{"a",42}};
    auto jval = nljs.get<Json::Value>();
    assert (jval["a"].asInt() == 42);
    nlohmann::json nljs2 = jval;
    assert (nljs2["a"] == 42);
    
    return 0;
}
