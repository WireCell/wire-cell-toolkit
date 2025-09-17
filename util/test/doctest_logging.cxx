// Some basic tests of logging.
// 
// Note, wcb generates a main() with some WCT related logging init code.
// see, eg the generated file:  build/util/wcdoctest-util.cxx

#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <string>
#include <chrono>

using namespace WireCell;

TEST_CASE("logging params")
{
    auto b = Log::logger("params");
    // Three styles
    b->debug("x={} y={}", fmt::arg("x", 42), fmt::arg("y", 6.9));

    // can not mix "automatic" and "manual"
    // b->debug("x={} y={1}", fmt::arg("x", 42), fmt::arg("y", 6.9));

    // Positional
    b->debug("x={0} y={1}", fmt::arg("x", 42), fmt::arg("y", 6.9));

    // Keyword
    b->debug("x={x} y={y}", fmt::arg("x", 42), fmt::arg("y", 6.9));

    // Can mix positional and keyword.
    b->debug("x={0} y={y}", fmt::arg("x", 42), fmt::arg("y", 6.9));
    b->debug("x={x} y={1}", fmt::arg("x", 42), fmt::arg("y", 6.9));

}

TEST_CASE("logging various")
{
    spdlog::level::level_enum active_level = (spdlog::level::level_enum)SPDLOG_ACTIVE_LEVEL;
    spdlog::info("compiled active level {} ({})",
                 spdlog::level::to_short_c_str(active_level),
                 SPDLOG_ACTIVE_LEVEL);

    auto b = Log::logger("before");

    Log::set_level("info");     // overrides SPDLOG_LEVEL from env

    auto a = Log::logger("after");

    Log::set_level("debug", "special");

    auto l = Log::logger("notshared", false);
    REQUIRE(l != spdlog::default_logger());
    Log::set_pattern("special pattern: %v", "notshared");

    auto s = Log::logger("special");

    l->error("test error l logger");
    b->error("test error b logger");
    a->error("test error a logger");
    s->error("test error s logger");
    spdlog::error("error default logger");

    l->info("info l logger");
    b->info("info b logger");
    a->info("info a logger");
    s->info("info s logger");
    spdlog::info("info default logger");

    l->debug("debug l logger");
    b->debug("debug b logger");
    a->debug("debug a logger");
    s->debug("debug s logger");
    spdlog::debug("debug default logger");

    SPDLOG_DEBUG("log from default debug CPP macro, compile --with-spdlog-active-level=debug to see");
    SPDLOG_LOGGER_DEBUG(s, "log from debug CPP macro, compile --with-spdlog-active-level=debug to see");
    SPDLOG_TRACE("log from default trace CPP macro, compile --with-spdlog-active-level=trace to see");
    SPDLOG_LOGGER_TRACE(s, "log from trace CPP macro, compile --with-spdlog-active-level=trace to see");

    auto t0 = std::chrono::high_resolution_clock::now();
    const int nlookups = 100000;
    for (int count = 0; count < nlookups; ++count) {
        auto l = Log::logger("lookup");
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
    spdlog::info("{} in {} us, {:.3f} MHz", nlookups, us.count(), double(nlookups) / us.count());

}
