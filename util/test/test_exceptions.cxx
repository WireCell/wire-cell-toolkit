#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/String.h"

#include <iostream>

using namespace WireCell;
using namespace std;

using WireCell::String::format;

int main()
{
    int value = 42;
    std::string omg = "omg";
    const double number = 6.9;

    cerr << format("something") << endl;
    // cerr << format("no percent", value) << endl;
    cerr << format("value=%d", value) << endl;
    cerr << format("value=%d, omg=\"%s\", number=%f", value, omg, number) << endl;

    try {
        THROW(ValueError() << errmsg{format("some error with value=%d msg=\"%s\"", value, omg)});
    }
    catch (ValueError& e) {
        cerr << "Caught:\n";
        cerr << e.what() << "\n";
    }
}
