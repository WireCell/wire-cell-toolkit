#include "WireCellGen/Cfg/TrackDepos/Nljs.hpp"
#include "WireCellUtil/nljs2jcpp.hpp"

//using namespace WireCell;
using TDConfig = WireCellGen::Cfg::TrackDepos::Config;

using nljs_t = nlohmann::json;

int main()
{
    TDConfig td{};
    nljs_t nljs = td;
    Json::Value jval2 = nljs.get<Json::Value>();
    nljs_t nljs2 = jval2;
    auto td2 = nljs2.get<TDConfig>();

    return 0;
}
