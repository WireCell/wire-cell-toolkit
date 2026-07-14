#include "WireCellUtil/doctest.h"
#include "WireCellUtil/TimeKeeper.h"
#include "WireCellUtil/Logging.h"

#include <unistd.h>

using namespace WireCell;

TEST_SUITE("timekeeper") {

TEST_CASE("timekeeper") {
    TimeKeeper tk("test_testing");
    spdlog::debug("{}", tk("sleeping"));
    sleep(1);
    spdlog::debug("{}", tk("awake"));

    TimeKeeper::deltat dt = tk.last_duration();
    REQUIRE_MESSAGE(dt.seconds() == 1, "Bad sleep.");

    spdlog::debug("TimeKeeper summary: {}", tk.summary());
}

}  // TEST_SUITE("timekeeper")
