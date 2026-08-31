#include "WireCellUtil/Quantity.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;
using namespace std;

TEST_SUITE("quantity") {

TEST_CASE("arithmetic operations") {
    Quantity a1(5, 1), b1(2, 3);

    CHECK(Quantity() == 0.0);
    CHECK((a1 * b1) == 10.0);
    CHECK((a1 / b1) == 2.5);
    CHECK((a1 + b1) == 7.0);
    CHECK((a1 - b1) == 3.0);
    CHECK((-a1) == -5.0);
}

TEST_CASE("comparison operators") {
    Quantity a1(5, 1), b1(2, 3);
    Quantity a2(5, 2);

    CHECK(a1 < 10.0);
    CHECK(a1 > 4.0);
    CHECK(a1 == 5.0);
    CHECK(a1 == a2);
    CHECK(a1 != b1);
    CHECK(a1 != 3.0);

    CHECK(b1 < a1);
    CHECK(a1 > b1);
}

}  // TEST_SUITE("quantity")
