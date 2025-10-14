#include "WireCellUtil/Configuration.h"

#include <boost/container_hash/hash.hpp>

using namespace WireCell;
using namespace std;

WireCell::Configuration WireCell::branch(WireCell::Configuration cfg, const std::string& dotpath)
{
    std::vector<std::string> path;
    boost::algorithm::split(path, dotpath, boost::algorithm::is_any_of("."));
    for (auto name : path) {
        cfg = cfg[name];
    }
    return cfg;
}

// http://stackoverflow.com/a/23860017
WireCell::Configuration WireCell::update(WireCell::Configuration& a, WireCell::Configuration& b)
{
    if (a.isNull()) {
        a = b;
        return b;
    }
    if (!a.isObject() || !b.isObject()) {
        return a;
    }

    for (const auto& key : b.getMemberNames()) {
        if (a[key].isObject()) {
            update(a[key], b[key]);
        }
        else {
            a[key] = b[key];
        }
    }
    return a;
}
WireCell::Configuration WireCell::update(WireCell::Configuration& a, const WireCell::Configuration& b)
{
    auto c = b;
    return update(a, c);
}

/// Append array b onto end of a and return a.
WireCell::Configuration WireCell::append(Configuration& a, Configuration& b)
{
    Configuration ret(Json::arrayValue);
    for (auto x : a) {
        ret.append(x);
    }
    for (auto x : b) {
        ret.append(x);
    }
    return ret;
}

size_t WireCell::hash(const Configuration& cfg)
{
    std::stringstream ss;
    ss << cfg;
    return std::hash<std::string>{}(ss.str());
}
