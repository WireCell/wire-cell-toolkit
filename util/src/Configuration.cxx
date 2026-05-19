#include "WireCellUtil/Configuration.h"

#include <boost/container_hash/hash.hpp>

using namespace WireCell;
using namespace std;

WireCell::Configuration WireCell::branch(const WireCell::Configuration& cfg, const std::string& dotpath)
{
    std::vector<std::string> path;
    boost::algorithm::split(path, dotpath, boost::algorithm::is_any_of("."));
    const WireCell::Configuration* cur = &cfg;
    for (const auto& name : path) {
        cur = &((*cur)[name]);
    }
    return *cur;
}

// ABI-compatibility shim for binaries (e.g. libWireCellLarsoft.so on cvmfs)
// compiled against the older by-value branch() signature. Defined inside a
// namespace block so the qualified definition is a declaration too; not
// re-exposed in the header, so new code keeps using the const-ref overload.
namespace WireCell {
    Configuration branch(Configuration cfg, const std::string& dotpath)
    {
        std::vector<std::string> path;
        boost::algorithm::split(path, dotpath, boost::algorithm::is_any_of("."));
        const Configuration* cur = &cfg;
        for (const auto& name : path) {
            cur = &((*cur)[name]);
        }
        return *cur;
    }
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
