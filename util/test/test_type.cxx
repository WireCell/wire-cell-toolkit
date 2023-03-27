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

// the two don't give exactly the same pattern....
#if defined(__clang__)
    AssertMsg("std::__1::vector<int, std::__1::allocator<int> >" == type(vi), "Clang demangling fails");
#else
    AssertMsg("std::vector<int, std::allocator<int> >" == type(vi), "GCC demangling fails");
#endif
}
