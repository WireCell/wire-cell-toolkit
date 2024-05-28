#include "WireCellUtil/Dtype.h"
#include "WireCellUtil/doctest.h"

using namespace WireCell;

template<typename T>
void check(size_t size)
{
    CHECK(dtype_size(dtype<T>()) == size);
}

TEST_CASE("dtype sizes")
{
    check<int>(4);
    check<float>(4);
    check<double>(8);
    check<int8_t>(1);
    check<int16_t>(2);
    check<int32_t>(4);
    check<int64_t>(8);
    check<char>(1);
    check<std::byte>(1);
}

TEST_CASE("dtype various") 
{
    CHECK(dtype<std::string>().empty());
    using myint = int64_t;
    using myuint = uint64_t;

    CHECK(dtype<myint>() == "i8");
    CHECK(dtype<myuint>() == "u8");

    CHECK(dtype_size("c16") == 16);

    CHECK(dtype(typeid(myint)) == "i8");
}
