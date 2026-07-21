#include <chrono>
#include <thread>
#include <iostream>
#include "WireCellSpng/TorchnvTools.h"

// Name collision for "CHECK" between torch and doctest.
#undef CHECK
#include "WireCellUtil/doctest.h"

TEST_CASE("spng nvtx basic functionality") {
    std::cerr << "Running NVTX basic functionality test..." << std::endl;
    
    // Test basic NVTX scoped range - should not crash or throw
    auto test_function = []() {
        NVTX_SCOPED_RANGE("test_range");
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    };
    
    REQUIRE_NOTHROW(test_function());
    std::cerr << "Basic NVTX test completed successfully!" << std::endl;
}


