#include "WireCellUtil/Logging.h"
#include "WireCellUtil/String.h"

#include "WireCellUtil/doctest.h"

using namespace WireCell::String;

TEST_CASE("string strip") {
    
    CHECK(strip("") == "");
    CHECK(strip(" ") == "");
    CHECK(strip("blah") == "blah");
    CHECK(strip("blah ") == "blah");
    CHECK(strip(" blah") == "blah");
    CHECK(strip(" blah ") == "blah");
    CHECK(strip(" \t\nblah\n\t \t\n") == "blah");
    CHECK(strip(" \t\nblah blah\n\t \t\n") == "blah blah");
}

