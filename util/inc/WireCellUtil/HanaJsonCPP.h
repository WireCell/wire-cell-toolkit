#ifndef WIRECELL_UTIL_HANAJSONCPP
#define WIRECELL_UTIL_HANAJSONCPP

/** Serialize between json and c++ struct

    It is recommended that configurable components define an explicit
    configuration struct and then use the to_json() and from_json() functions
    defined here to serialize between the struct and the JSON object in the
    IConfigurable::configure() and IConfigurable::default_configuration()
    methods.  Defining the C++ config struct helps you and others reason about
    configuration information.

    See util/tests/doctest_hanajsoncpp.cxx for one example.
*/
#include "boost/hana.hpp"
#include "json/json.h"          // JsonCPP::Value aka WireCell::Configuration

#include <iostream>
#include <string>

namespace WireCell::HanaJsonCPP {

    // Generic JSON serialization function for any Hana-reflectable struct
    template <typename T>
    Json::Value to_json(const T& obj) {
        Json::Value root;
        boost::hana::for_each(boost::hana::accessors<T>(), [&](auto const& accessor) {
            auto name_hana = boost::hana::first(accessor); // boost::hana::string compile-time literal
            const char* name_cstr = boost::hana::to<char const*>(name_hana); // to C-string for Json::Value key
            auto const& value = boost::hana::second(accessor)(obj); // Get value by reference

            if constexpr (boost::hana::is_a<int, decltype(value)>) {
                root[name_cstr] = value;
            } else if constexpr (boost::hana::is_a<std::string, decltype(value)>) {
                root[name_cstr] = value;
            } else if constexpr (boost::hana::is_a<bool, decltype(value)>) {
                root[name_cstr] = value;
            } else if constexpr (boost::hana::is_a<double, decltype(value)>) {
                root[name_cstr] = value;
            }
            // Add more `else if constexpr` for other types as needed
        });
        return root;
    }

    // Generic JSON deserialization function for any Hana-reflectable struct
    template <typename T>
    void from_json(T& obj, const Json::Value& root) {
        boost::hana::for_each(boost::hana::accessors<T>(), [&](auto const& accessor) {
            auto name_hana = boost::hana::first(accessor);
            const char* name_cstr = boost::hana::to<char const*>(name_hana);
            auto& value = boost::hana::second(accessor)(obj); // Get mutable reference to value

            if (root.isMember(name_cstr)) {
                if constexpr (boost::hana::is_a<int, decltype(value)>) {
                    value = root[name_cstr].asInt();
                } else if constexpr (boost::hana::is_a<std::string, decltype(value)>) {
                    value = root[name_cstr].asString();
                } else if constexpr (boost::hana::is_a<bool, decltype(value)>) {
                    value = root[name_cstr].asBool();
                } else if constexpr (boost::hana::is_a<double, decltype(value)>) {
                    value = root[name_cstr].asDouble();
                }
                // Add more `else if constexpr` for other types as needed
            }
        });
    }
}


#endif
