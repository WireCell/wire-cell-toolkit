#include "WireCellUtil/Testing.h"
#include <regex>
#include <iostream>

void retest(const std::string& res, const std::string& str, bool should_match = true)
{
    std::regex re(res);
    std::smatch smatch;
    bool matched = std::regex_match(str, smatch, re);
    Assert(matched == should_match);
    std::cerr << "regex_match(\"" << str << "\", \"" << res << "\") "
              << (should_match ? "matches\n" : "differs\n");
    for (size_t ind=0; ind<smatch.size(); ++ind) {
        std::cerr << "\t" << ind << ": " << smatch[ind] << "\n";
    }
}

void test_frames(const std::string& res, bool named=false)
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

int main()
{
    test_frames("^frames/[^/]+", true);
    test_frames("^frames/[[:digit:]]+$");
    test_frames("^frames/[[:alnum:]]+$", true);
}
