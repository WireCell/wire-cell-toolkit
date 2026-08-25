// Some basic tests of logging.
// 
// Note, wcb generates a main() with some WCT related logging init code.
// see, eg the generated file:  build/util/wcdoctest-util.cxx

#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"

#include <string>
#include <chrono>
#include <vector>
#include <fstream>
#include <filesystem>

#include <unistd.h>

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

// A log FILE must not tear.
//
// Log::logger(name, share_sinks=false) makes a fresh sink per logger so each
// can carry its own pattern -- Aux::Logger::set_name() does exactly that for
// every configured component.  With spdlog's basic_file_sink that also means a
// fresh FILE* per logger, so N components leave N independent 4096-byte stdio
// buffers appending to one file.  They drain at different times, so one
// component's buffer lands in the MIDDLE of another's line.
//
// Seen in production: an SBND imaging job with 98 component loggers tore 70 of
// 5240 lines, every splice on an exact 4096-byte file offset, with the spliced
// fragment up to 26 seconds older than the line it cut.  Downstream log parsing
// was patched four separate times to tolerate it before the cause was found.
//
// The test writes enough through enough loggers to cross many buffer
// boundaries and then requires that every line carry exactly one record
// header.  Two headers on one line is the splice.
TEST_CASE("logging file is not torn")
{
    const std::string fname =
        (std::filesystem::temp_directory_path() /
         ("wct-doctest-logging-tear-" + std::to_string(getpid()) + ".log")).string();
    std::filesystem::remove(fname);

    Log::add_file(fname, "debug");

    // The sentinel opens each record and cannot occur in a payload, so
    // counting it per line counts records per line.
    const std::string sentinel = "|@|";
    const size_t nloggers = 32, nmsgs = 200;

    std::vector<Log::logptr_t> logs;
    for (size_t i = 0; i < nloggers; ++i) {
        const std::string name = "tear/" + std::to_string(i);
        auto l = Log::logger(name, false);   // unique sinks, as Aux::Logger does
        l->set_level(spdlog::level::debug);
        // A per-logger pattern -- the reason unique sinks exist at all.
        l->set_pattern(sentinel + " " + name + " %v");
        logs.push_back(l);
    }

    // Round-robin so the per-logger buffers fill at different rates, which is
    // what a real job does.  Messages are long enough that 32 x 200 of them
    // cross buffer boundaries many times over.
    const std::string filler(180, 'x');
    for (size_t m = 0; m < nmsgs; ++m) {
        for (size_t i = 0; i < nloggers; ++i) {
            logs[i]->debug("msg {} of logger {} {}", m, i, filler);
        }
    }
    for (auto& l : logs) {
        l->flush();
    }

    // Other loggers in this test binary share the file, so count RECORDS by
    // sentinel rather than counting lines: a line carrying two sentinels is a
    // splice, and the sentinel total must account for every message written.
    std::ifstream fp(fname);
    REQUIRE(fp.good());
    std::string line;
    size_t nrecords = 0, ntorn = 0;
    std::string first_torn;
    while (std::getline(fp, line)) {
        size_t nsent = 0;
        for (size_t p = line.find(sentinel); p != std::string::npos;
             p = line.find(sentinel, p + sentinel.size())) {
            ++nsent;
        }
        nrecords += nsent;
        if (nsent > 1) {
            ++ntorn;
            if (first_torn.empty()) {
                first_torn = line.substr(0, 160);
            }
        }
    }
    fp.close();

    spdlog::info("tear check: {} records, {} spliced lines, first: {}",
                 nrecords, ntorn, first_torn);

    // No line may carry two records, and nothing may be lost.
    CHECK(ntorn == 0);
    CHECK(nrecords == nloggers * nmsgs);

    if (ntorn == 0) {
        std::filesystem::remove(fname);
    }
}
