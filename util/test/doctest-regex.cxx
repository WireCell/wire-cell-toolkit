#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <regex>

void retest(const std::string& res, const std::string& str, bool should_match = true)
{
    std::regex re(res);
    std::smatch smatch;
    bool matched = std::regex_match(str, smatch, re);
    CHECK(matched == should_match);
    spdlog::debug("regex_match(\"{}\", \"{}\") {}", str, res,
                  (should_match ? "matches" : "differs"));
    for (size_t ind = 0; ind < smatch.size(); ++ind) {
        spdlog::debug("\t{}: {}", ind, smatch[ind].str());
    }
}

void test_frames(const std::string& res, bool named = false)
{
    retest(res, "frames", false);
    retest(res, "frames/", false);
    retest(res, "frames/0", true);
    retest(res, "frames/NAME", named);
    retest(res, "frames/0/traces", false);
    retest(res, "frames/0/traces/", false);
    retest(res, "frames/0/traces/0", false);
    retest(res, "frame/0", false);
    retest(res, "frame", false);
}

TEST_SUITE("regex") {

TEST_CASE("frames any non-slash segment") {
    test_frames("^frames/[^/]+", true);
}

TEST_CASE("frames digits only segment") {
    test_frames("^frames/[[:digit:]]+$");
}

TEST_CASE("frames alphanumeric segment") {
    test_frames("^frames/[[:alnum:]]+$", true);
}

}  // TEST_SUITE("regex")
