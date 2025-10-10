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

#define IS_HANA_STRUCT(T) boost::hana::is_a<boost::hana::Struct, T>()

namespace WireCell::HanaJsonCPP {

    namespace hana = boost::hana;

// --- Helper Traits (SFINAE) ---

// 1. Trait to detect if T is a specialization of std::vector
    template <typename T>
    struct is_vector : std::false_type {};

    template <typename T, typename A>
    struct is_vector<std::vector<T, A>> : std::true_type {};

    template <typename T>
    constexpr bool is_vector_v = is_vector<T>::value;

// 2. Trait to detect if T is a Hana-reflectable struct.
// CORRECTED: Use hana::Struct<T>::value syntax.
    template <typename T>
    constexpr bool is_hana_struct_v = hana::Struct<T>::value;

// --- Forward declaration of the JSON handlers for recursion ---
    template <typename T>
    Json::Value to_json(const T& obj);

    template <typename T>
    void from_json(T& obj, const Json::Value& root);


// --- 2. Generic JSON Serialization ---
    template <typename T>
    Json::Value to_json(const T& obj) {
        Json::Value root;
        hana::for_each(hana::accessors<T>(), [&](auto const& accessor) {
            auto name_hana = hana::first(accessor);
            const char* name_cstr = hana::to<char const*>(name_hana);
            auto const& value = hana::second(accessor)(obj);
            using ValueType = std::decay_t<decltype(value)>;

            if constexpr (is_hana_struct_v<ValueType>) {
                // Case 1: Nested Hana Struct -> Recursive call
                root[name_cstr] = to_json(value);
            } else if constexpr (is_vector_v<ValueType>) {
                // Case 2: Vector
                if (value.empty()) return;

                using InnerType = typename ValueType::value_type;
            
                if constexpr (is_hana_struct_v<InnerType>) {
                    // Case 2a: Vector of Hana Structs -> Iterate and recurse
                    for (const auto& item : value) {
                        root[name_cstr].append(to_json(item));
                    }
                } else {
                    // Case 2b: Vector of Primitives
                    for (const auto& item : value) {
                        if constexpr (std::is_arithmetic_v<InnerType> || std::is_same_v<InnerType, std::string>) {
                            root[name_cstr].append(item);
                        }
                    }
                }
            } else if constexpr (std::is_arithmetic_v<ValueType> || std::is_same_v<ValueType, std::string>) {
                // Handle primitives
                root[name_cstr] = value;
            }
        });
        return root;
    }

// --- 3. Generic JSON Deserialization ---
    template <typename T>
    void from_json(T& obj, const Json::Value& root) {
        hana::for_each(hana::accessors<T>(), [&](auto const& accessor) {
            auto name_hana = hana::first(accessor);
            const char* name_cstr = hana::to<char const*>(name_hana);
            auto& value = hana::second(accessor)(obj);
            using ValueType = std::decay_t<decltype(value)>;

            if (!root.isMember(name_cstr)) {
                return;
            }

            auto const& json_field = root[name_cstr];

            if constexpr (is_hana_struct_v<ValueType>) {
                // Case 1: Nested Hana Struct -> Recursive call
                from_json(value, json_field);

            } else if constexpr (is_vector_v<ValueType>) {
                // Case 2: Vector
                if (!json_field.isArray()) return;
            
                value.clear();
                using InnerType = typename ValueType::value_type;
            
                if constexpr (is_hana_struct_v<InnerType>) {
                    // Case 2a: Vector of Hana Structs -> Iterate and recurse
                    for (Json::Value::ArrayIndex i = 0; i < json_field.size(); ++i) {
                        InnerType temp_item{};
                        from_json(temp_item, json_field[i]);
                        value.push_back(std::move(temp_item));
                    }
                } else {
                    // Case 2b: Vector of Primitives
                    for (Json::Value::ArrayIndex i = 0; i < json_field.size(); ++i) {
                        if constexpr (std::is_same_v<InnerType, int>) {
                            value.push_back(json_field[i].asInt());
                        } else if constexpr (std::is_same_v<InnerType, std::string>) {
                            value.push_back(json_field[i].asString());
                        } else if constexpr (std::is_same_v<InnerType, bool>) {
                            value.push_back(json_field[i].asBool());
                        } else if constexpr (std::is_same_v<InnerType, double>) {
                            value.push_back(json_field[i].asDouble());
                        }
                    }
                }
            } else if constexpr (std::is_same_v<ValueType, int>) {
                value = json_field.asInt();
            } else if constexpr (std::is_same_v<ValueType, std::string>) {
                value = json_field.asString();
            } else if constexpr (std::is_same_v<ValueType, bool>) {
                value = json_field.asBool();
            } else if constexpr (std::is_same_v<ValueType, double>) {
                value = json_field.asDouble();
            }
        });
    }

}


#endif
