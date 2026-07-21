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

    // Trait to detect if T is a specialization of std::vector
    template <typename T>
    struct is_vector : std::false_type {};

    template <typename T, typename A>
    struct is_vector<std::vector<T, A>> : std::true_type {};

    template <typename T>
    constexpr bool is_vector_v = is_vector<T>::value;

    // Trait to detect if T is a Hana-reflectable struct.
    template <typename T>
    constexpr bool is_hana_struct_v = hana::Struct<T>::value;

    // --- Forward declaration of the JSON handlers for recursion ---

    /// Return JSON built from type.
    template <typename T>
    Json::Value to_json(const T& obj);

    /// Set type from JSON, return portion of root that was unused.
    template <typename T>
    Json::Value from_json(T& obj, const Json::Value& root);

    /// Convert a single element value (scalar, nested vector, or Hana struct) to JSON.
    template <typename T>
    Json::Value elem_to_json(const T& val);

    /// Deserialize a single element value from JSON. Returns true on success.
    template <typename T>
    bool elem_from_json(T& val, const Json::Value& j);

    // --- Element-level helpers (handle scalars, nested vectors, and Hana structs) ---

    template <typename T>
    Json::Value elem_to_json(const T& val) {
        if constexpr (is_hana_struct_v<T>) {
            return to_json(val);
        } else if constexpr (is_vector_v<T>) {
            Json::Value arr(Json::arrayValue);
            for (const auto& item : val) {
                arr.append(elem_to_json(item));
            }
            return arr;
        } else if constexpr (std::is_arithmetic_v<T> || std::is_same_v<T, std::string>) {
            return Json::Value(val);
        }
        return Json::Value();
    }

    template <typename T>
    bool elem_from_json(T& val, const Json::Value& j) {
        if constexpr (is_hana_struct_v<T>) {
            from_json(val, j);
            return true;
        } else if constexpr (is_vector_v<T>) {
            if (!j.isArray()) return false;
            val.clear();
            using Item = typename T::value_type;
            bool all_ok = true;
            for (Json::Value::ArrayIndex i = 0; i < j.size(); ++i) {
                Item item{};
                if (!elem_from_json(item, j[i])) all_ok = false;
                val.push_back(std::move(item));
            }
            return all_ok;
        } else if constexpr (std::is_same_v<T, int>) {
            if (j.isInt()) { val = j.asInt(); return true; }
        } else if constexpr (std::is_same_v<T, std::string>) {
            if (j.isString()) { val = j.asString(); return true; }
        } else if constexpr (std::is_same_v<T, bool>) {
            if (j.isBool()) { val = j.asBool(); return true; }
        } else if constexpr (std::is_same_v<T, double>) {
            if (j.isNumeric()) { val = j.asDouble(); return true; }
        }
        return false;
    }

    // --- Generic JSON Serialization ---
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
                // Case 2: Vector (scalars, nested vectors, or Hana structs)
                if (value.empty()) return;

                for (const auto& item : value) {
                    root[name_cstr].append(elem_to_json(item));
                }
            } else if constexpr (std::is_same_v<ValueType, Json::Value>) {
                // Json::Value attribute -> copy the entire JSON tree
                root[name_cstr] = value;
            } else if constexpr (std::is_arithmetic_v<ValueType> || std::is_same_v<ValueType, std::string>) {
                // Handle primitives
                root[name_cstr] = value;
            }
        });
        return root;
    }

    /**
     * @param obj The C++ struct to deserialize into.
     * @param root The source JSON object.
     * @return A Json::Value containing attributes from 'root' that were not
     *         consumed by the deserialization process.
     */
    template <typename T>
    Json::Value from_json(T& obj, const Json::Value& root) {
        // Start with a copy of the root. We remove members as they are consumed.
        Json::Value unused = root;

        hana::for_each(hana::accessors<T>(), [&](auto const& accessor) {
            auto name_hana = hana::first(accessor);
            const char* name_cstr = hana::to<char const*>(name_hana);
            auto& value = hana::second(accessor)(obj);
            using ValueType = std::decay_t<decltype(value)>;

            if (!unused.isMember(name_cstr)) {
                return; // Attribute not present in JSON input
            }

            auto const& json_field = root[name_cstr];

            if constexpr (is_hana_struct_v<ValueType>) {
                // Case 1: Nested Hana Struct -> Recursive call
                Json::Value nested_unused = from_json(value, json_field);

                if (nested_unused.empty() || nested_unused.isNull()) {
                    // If the nested struct consumed everything, remove the top-level key entirely.
                    unused.removeMember(name_cstr);
                } else {
                    // If the nested struct left fragments, update the top-level unused tracker with those fragments.
                    unused[name_cstr] = std::move(nested_unused);
                }

            } else if constexpr (is_vector_v<ValueType>) {
                // Case 2: Vector (scalars, nested vectors, or Hana structs)
                if (!json_field.isArray()) {
                    // Field exists but is not an array; leave it in 'unused'.
                    return;
                }

                value.clear();
                using InnerType = typename ValueType::value_type;
                bool all_consumed = true;
                for (Json::Value::ArrayIndex i = 0; i < json_field.size(); ++i) {
                    InnerType item{};
                    if (!elem_from_json(item, json_field[i])) all_consumed = false;
                    value.push_back(std::move(item));
                }
                if (all_consumed) {
                    unused.removeMember(name_cstr);
                }
            } else if constexpr (std::is_same_v<ValueType, Json::Value>) {
                // Json::Value attribute -> copy the parsed JSON data. Fully consumed.
                value = json_field;
                unused.removeMember(name_cstr);
            } else {
                // Handle primitives
                bool consumed = false;
                if constexpr (std::is_same_v<ValueType, int>) {
                    if (json_field.isInt()) { value = json_field.asInt(); consumed = true; }
                } else if constexpr (std::is_same_v<ValueType, std::string>) {
                    if (json_field.isString()) { value = json_field.asString(); consumed = true; }
                } else if constexpr (std::is_same_v<ValueType, bool>) {
                    if (json_field.isBool()) { value = json_field.asBool(); consumed = true; }
                } else if constexpr (std::is_same_v<ValueType, double>) {
                    if (json_field.isNumeric()) { value = json_field.asDouble(); consumed = true; }
                }

                if (consumed) {
                    unused.removeMember(name_cstr);
                }
            }
        });
        
        return unused;
    }

    // // --- Generic JSON Deserialization ---
    // template <typename T>
    // Json::Value from_json(T& obj, const Json::Value& root) {
    //     hana::for_each(hana::accessors<T>(), [&](auto const& accessor) {
    //         auto name_hana = hana::first(accessor);
    //         const char* name_cstr = hana::to<char const*>(name_hana);
    //         auto& value = hana::second(accessor)(obj);
    //         using ValueType = std::decay_t<decltype(value)>;

    //         if (!root.isMember(name_cstr)) {
    //             return root;
    //         }

    //         auto const& json_field = root[name_cstr];

    //         if constexpr (is_hana_struct_v<ValueType>) {
    //             // Case 1: Nested Hana Struct -> Recursive call
    //             from_json(value, json_field);

    //         } else if constexpr (is_vector_v<ValueType>) {
    //             // Case 2: Vector
    //             if (!json_field.isArray()) return;
            
    //             value.clear();
    //             using InnerType = typename ValueType::value_type;
            
    //             if constexpr (is_hana_struct_v<InnerType>) {
    //                 // Case 2a: Vector of Hana Structs -> Iterate and recurse
    //                 for (Json::Value::ArrayIndex i = 0; i < json_field.size(); ++i) {
    //                     InnerType temp_item{};
    //                     from_json(temp_item, json_field[i]);
    //                     value.push_back(std::move(temp_item));
    //                 }
    //             } else {
    //                 // Case 2b: Vector of Primitives
    //                 for (Json::Value::ArrayIndex i = 0; i < json_field.size(); ++i) {
    //                     if constexpr (std::is_same_v<InnerType, int>) {
    //                         value.push_back(json_field[i].asInt());
    //                     } else if constexpr (std::is_same_v<InnerType, std::string>) {
    //                         value.push_back(json_field[i].asString());
    //                     } else if constexpr (std::is_same_v<InnerType, bool>) {
    //                         value.push_back(json_field[i].asBool());
    //                     } else if constexpr (std::is_same_v<InnerType, double>) {
    //                         value.push_back(json_field[i].asDouble());
    //                     }
    //                 }
    //             }
    //         } else if constexpr (std::is_same_v<ValueType, Json::Value>) {
    //             // Json::Value attribute -> copy the parsed JSON data
    //             value = json_field;
    //         } else if constexpr (std::is_same_v<ValueType, int>) {
    //             value = json_field.asInt();
    //         } else if constexpr (std::is_same_v<ValueType, std::string>) {
    //             value = json_field.asString();
    //         } else if constexpr (std::is_same_v<ValueType, bool>) {
    //             value = json_field.asBool();
    //         } else if constexpr (std::is_same_v<ValueType, double>) {
    //             value = json_field.asDouble();
    //         }
    //     });
    // }

}


#endif
