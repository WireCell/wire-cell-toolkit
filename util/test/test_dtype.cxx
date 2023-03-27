#include "WireCellUtil/Dtype.h"
#include "WireCellUtil/Testing.h"

using namespace WireCell;

int main()
{
    Assert(dtype<std::string>().empty());
    using myint = int64_t;
    using myuint = uint64_t;

    Assert(dtype<myint>() == "i8");
    Assert(dtype<myuint>() == "u8");

    Assert(dtype_size("c16") == 16);

    Assert(dtype(typeid(myint)) == "i8");

    return 0;
}
