#include "WireCellUtil/Type.h"
#include "WireCellUtil/Testing.h"

#include <iostream>
#include <vector>
using namespace std;
using namespace WireCell;

int main()
{
    int i;
    vector<int> vi;

    Assert(typeid(int) == typeid(i));


    cerr << "int: " << type(i) << endl;
    cerr << "vector<int>: \"" << type(vi) << "\"\n";

    std::vector<std::string> type_strings = {
        "std::vector<int, std::allocator<int> >", // GCC, newer clang
        "std::__1::vector<int, std::__1::allocator<int> >" // older clang
    };
    for (const auto& tstr : type_strings) {
        if (tstr == type(vi)) {
            return 0;
        }
    }
    cerr << "Failed to match any expected type string\n";
    return 1;
}
