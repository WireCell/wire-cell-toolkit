#include "WireCellUtil/Testing.h"
#include "WireCellUtil/doctest.h"

#include <list>

TEST_CASE("list iterator stability")
{
    using List = std::list<int>;
    using Iter = List::iterator;
    List li = {3,2,1,0};
    Iter i3 = li.begin();
    Iter i2 = std::next(i3);
    Iter i1 = std::next(i2);
    Iter i0 = std::next(i1);
    li.sort();

    Iter i = li.begin();
    REQUIRE(*i == 0);
    ++i;
    REQUIRE(*i == 1);
    ++i;
    REQUIRE(*i == 2);
    ++i;
    REQUIRE(*i == 3);
    ++i;

    REQUIRE(*i0 == 0);
    REQUIRE(*i1 == 1);
    REQUIRE(*i2 == 2);
    REQUIRE(*i3 == 3);


    REQUIRE(i0 == li.begin());
    REQUIRE(i1 == std::next(i0));
    REQUIRE(i3 == std::prev(li.end()));
    REQUIRE(std::next(i3) == li.end());

    Iter i4 = li.end();
    li.push_back(4);
    REQUIRE(i4 != std::prev(li.end()));
    REQUIRE(4 == *std::prev(li.end()));

}
