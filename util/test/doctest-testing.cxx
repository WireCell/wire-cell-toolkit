/** Test the WireCellUtil/Testing helpers.
 */

#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/Testing.h"

using namespace WireCell;

TEST_SUITE("testing") {

TEST_CASE("assert passes when condition is true") {
    Assert(true);
    AssertMsg(true, "this assert should not throw");
    CHECK(true);
}

TEST_CASE("AssertionError thrown and catchable") {
    bool caught = false;
    try {
        AssertMsg(false, "this assert should be caught");
    }
    catch (AssertionError& e) {
        spdlog::debug("Caught:\n{}", errstr(e));
        caught = true;
    }
    CHECK(caught);
}

}  // TEST_SUITE("testing")
