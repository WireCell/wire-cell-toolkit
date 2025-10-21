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
                // Case 2: Vector
                if (!json_field.isArray()) {
                    // If the field exists but is not an array, we skip conversion,
                    // leaving it in 'unused'.
                    return;
                }
            
                value.clear();
                using InnerType = typename ValueType::value_type;
                Json::Value array_unused(Json::arrayValue);
            
                if constexpr (is_hana_struct_v<InnerType>) {
                    // Case 2a: Vector of Hana Structs -> Iterate and recurse
                    bool all_consumed = true;
                    for (Json::Value::ArrayIndex i = 0; i < json_field.size(); ++i) {
                        InnerType temp_item{};
                        // Recurse on the item
                        Json::Value item_unused = from_json(temp_item, json_field[i]);
                        value.push_back(std::move(temp_item));

                        if (!item_unused.empty()) {
                            array_unused.append(std::move(item_unused));
                            all_consumed = false;
                        } else {
                            // Use null placeholder to maintain index sync if needed,
                            // although JsonCPP arrays skip nulls on output unless asked otherwise.
                            // For simplicity, we track only real unused parts.
                            // However, we track if *anything* was left unused.
                            array_unused.append(Json::nullValue); 
                        }
                    }

                    if (all_consumed) {
                        unused.removeMember(name_cstr); // Entire array was fully consumed
                    } else {
                        // Clean up null placeholders before assigning array_unused back
                        Json::Value cleaned_array_unused(Json::arrayValue);
                        for (const auto& item : array_unused) {
                            if (!item.isNull()) {
                                cleaned_array_unused.append(item);
                            }
                        }
                        
                        // Assign the cleaned list of unused fragments
                        unused[name_cstr] = std::move(cleaned_array_unused);
                    }

                } else {
                    // Case 2b: Vector of Primitives
                    bool all_consumed = true;
                    for (Json::Value::ArrayIndex i = 0; i < json_field.size(); ++i) {
                        // Assuming primitives are always consumed if they match the type check
                        if constexpr (std::is_same_v<InnerType, int>) {
                            if (json_field[i].isInt()) {
                                value.push_back(json_field[i].asInt());
                            } else { all_consumed = false; }
                        } else if constexpr (std::is_same_v<InnerType, std::string>) {
                            if (json_field[i].isString()) {
                                value.push_back(json_field[i].asString());
                            } else { all_consumed = false; }
                        } else if constexpr (std::is_same_v<InnerType, bool>) {
                            if (json_field[i].isBool()) {
                                value.push_back(json_field[i].asBool());
                            } else { all_consumed = false; }
                        } else if constexpr (std::is_same_v<InnerType, double>) {
                            if (json_field[i].isNumeric()) {
                                value.push_back(json_field[i].asDouble());
                            } else { all_consumed = false; }
                        }
                    }
                    if (all_consumed) {
                        // Primitive vectors are fully consumed if the type was an array
                        unused.removeMember(name_cstr);
                    }
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
