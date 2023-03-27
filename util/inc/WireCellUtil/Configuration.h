#ifndef WIRECELL_CONFIGURATION
#define WIRECELL_CONFIGURATION

#include <json/json.h>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

#include <sstream>
#include <vector>

namespace WireCell {

    /** The Wire Cell Toolkit configuration layer uses Json::Value
        objects for its transient data model.

        The Configuration type is a recursive in that one
        Configuration object may contain others.

        The WCT assumes an Object Configuration Protocol is
        implemented by the "client" code that uses the toolkit.  See
        for example the reference implementation in the wire-cell
        command line program.  This protocol consists of:

        1) Creation of the list of Configuration objects,
        typically via WireCell::Persist::load().

        2) Iterating this list and for each suitable Configuration
        object instantiating the appropriately matching IConfigurable
        via the WireCell::NamedFactory facility.

        3) Obtaining the default Configuration object from the
        IConfigurable and merging the user-supplied on top of it.

        4) Passing the resulting Configuration object to the
        IConfigurable::configure() method.

        A suitable Configuration object must have two top level keys
        and may have a third optional key:

       - type :: must match the "component class name" (as given as
         first argument to the WIRECELL_FACTORY() macro in the
         implementation file)

       - data :: an object which follows a schema which is specific to
         each IConfigurable implementation.

       - name :: an optional Instance Name.  If not given, the default instance of the type will be used.
     */

    typedef Json::Value Configuration;

    /// For persistence use WireCell::Persist::load() and
    /// WireCell::Persist::dump().

    /// The following functions provide some access methods which add
    /// some value beyond what Json::Value provides including some
    /// support for basic WCT types.

    /// Convert a configuration value to a particular type.
    template <typename T>
    T convert(const Configuration& cfg, const T& def = T())
    {
        return def;
    }
    template <>
    inline bool convert<bool>(const Configuration& cfg, const bool& def)
    {
        if (cfg.isNull()) return def;
        return cfg.asBool();
    }
    template <>
    inline int convert<int>(const Configuration& cfg, const int& def)
    {
        if (cfg.isNull()) return def;
        return cfg.asInt();
    }
    template <>
    inline float convert<float>(const Configuration& cfg, const float& def)
    {
        if (cfg.isNull()) return def;
        return cfg.asDouble();
    }
    template <>
    inline double convert<double>(const Configuration& cfg, const double& def)
    {
        if (cfg.isNull()) return def;
        return cfg.asDouble();
    }
    template <>
    inline std::string convert<std::string>(const Configuration& cfg, const std::string& def)
    {
        if (cfg.isNull()) return def;
        return cfg.asString();
    }
    template <>
    inline  // fixme: ignores default
        std::vector<std::string>
        convert<std::vector<std::string> >(const Configuration& cfg, const std::vector<std::string>& def)
    {
        std::vector<std::string> ret;
        for (auto v : cfg) {
            ret.push_back(convert<std::string>(v));
        }
        return ret;
    }
    template <>
    inline  // fixme: ignores default
        std::vector<int>
        convert<std::vector<int> >(const Configuration& cfg, const std::vector<int>& def)
    {
        std::vector<int> ret(cfg.size());
        std::transform(cfg.begin(), cfg.end(), ret.begin(),
                       [](const auto& e) { return e.asInt(); });
        return ret;
    }
    template <>
    inline  // fixme: ignores default
        std::vector<float>
        convert<std::vector<float> >(const Configuration& cfg, const std::vector<float>& def)
    {
        std::vector<float> ret(cfg.size());
        std::transform(cfg.begin(), cfg.end(), ret.begin(),
                       [](const auto& e) { return e.asFloat(); });
        return ret;
    }
    template <>
    inline  // fixme: ignores default
        std::vector<double>
        convert<std::vector<double> >(const Configuration& cfg, const std::vector<double>& def)
    {
        std::vector<double> ret(cfg.size());
        std::transform(cfg.begin(), cfg.end(), ret.begin(),
                       [](const auto& e) { return e.asDouble(); });
        return ret;
    }
    // For Point and Ray converters, see Point.h 


    /// Follow a dot.separated.path and return the branch there.
    Configuration branch(Configuration cfg, const std::string& dotpath);

    /// Merge dictionary b into a, return a
    Configuration update(Configuration& a, Configuration& b);

    /// Return an array which is composed of the array b appended to the array a.
    // fixme: this should be called "extend".
    Configuration append(Configuration& a, Configuration& b);

    /// Return dictionary in given list if it value at dotpath matches
    template <typename T>
    Configuration find(Configuration& lst, const std::string& dotpath, const T& val)
    {
        for (auto ent : lst) {
            auto maybe = branch(ent, dotpath);
            if (maybe.isNull()) {
                continue;
            }
            if (convert<T>(maybe) == val) {
                return maybe;
            }
        }
        return Configuration();
    }

    /// Get value in configuration at the dotted path from or return default.
    template <typename T>
    T get(Configuration cfg, const std::string& dotpath, const T& def = T())
    {
        return convert(branch(cfg, dotpath), def);
    }

    /// Set a configuration value with a typed value
    template <typename T>
    void assign(Configuration& cfg, const T& val)
    {
        cfg = val;
    }

    template <typename T>
    void assign(Configuration& cfg, const std::vector<T>& val)
    {
        const size_t size = val.size();
        cfg = Json::arrayValue;
        cfg.resize(size);
        for (size_t ind=0; ind<size; ++ind) {
            assign(cfg[(int)ind], val[ind]);
        }
    }

    template <typename KeyType, typename ValueType>
    void assign(Configuration& cfg, const std::map<KeyType, ValueType>& val)
    {
        cfg = Json::objectValue;
        for (const auto& [k,v] : val) {
            assign(cfg[k], v);
        }
    }

    /// Put value in configuration at the dotted path.
    template <typename T>
    void put(Configuration& cfg, const std::string& dotpath, const T& val)
    {
        Configuration* ptr = &cfg;
        std::vector<std::string> path;
        boost::algorithm::split(path, dotpath, boost::algorithm::is_any_of("."));
        for (auto name : path) {
            ptr = &(*ptr)[name];
        }
        assign(*ptr, val);
    }

}  // namespace WireCell

#endif
