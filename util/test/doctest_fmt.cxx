// test libfmt that we get via spdlog

#include "WireCellUtil/Fmt.h"
#include "WireCellUtil/doctest.h"
#include <iostream>

using namespace WireCell;

TEST_CASE("fmtlib named parameters") {
    // --- Case 1: Scalar parameter ---
    Configuration scalar_str = "World";
    CHECK("Hello, World!" == Fmt::format("Hello, {}!", scalar_str));

    Configuration scalar_int = 42;
    CHECK("The answer is 42." == Fmt::format("The answer is {}.", scalar_int));


    Configuration scalar_float = 3.14159;
    CHECK("Value of Pi: 3.14" == Fmt::format("Value of Pi: {:.2f}", scalar_float));

    // --- Case 2: Array parameter ---
    Json::Value array_params;
    array_params.append("one");
    array_params.append(2);
    array_params.append(true);
    CHECK("Got one, then 2, then true" == Fmt::format("Got {}, then {}, then {}", array_params));

    // --- Case 3: Object parameter (named arguments) ---
    Json::Value object_params;
    object_params["name"] = "Leia";
    object_params["planet"] = "Alderaan";
    object_params["age"] = 19;
    CHECK("Character: Leia, Planet: Alderaan, Age: 19" == Fmt::format("Character: {name}, Planet: {planet}, Age: {age}", object_params));

    // --- Error Case Example ---
    try {
        Json::Value bad_params;
        bad_params["data"].append(1); // Nested array
        Fmt::format("This will fail: {data}", bad_params);
    } catch (const std::invalid_argument& e) {
        
    }

}
