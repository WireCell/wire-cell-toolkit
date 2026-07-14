#include "WireCellUtil/RandomIter.h"
#include "WireCellUtil/doctest.h"
#include "WireCellUtil/Logging.h"

#include <vector>

using namespace WireCell;

TEST_SUITE("randomiter") {

TEST_CASE("arithmetic and iteration") {
    typedef int value;
    typedef std::vector<value> container;

    container array(4);
    array[0] = 1;
    array[1] = 2;
    array[2] = 3;
    array[3] = 4;

    typedef RandomIter<container, value> iterator;

    iterator beg(array), end(true, array);
    int count = 0;
    for (iterator it = beg; it != end; ++it) {
        spdlog::debug("randomiter value: {}", *it);
        ++count;
    }
    CHECK(count == 4);

    beg += 2;
    spdlog::debug("after +=2: {}", *beg);
    CHECK(3 == *beg);

    beg -= 2;
    spdlog::debug("after -=2: {}", *beg);
    CHECK(1 == *beg);
}

}  // TEST_SUITE("randomiter")
