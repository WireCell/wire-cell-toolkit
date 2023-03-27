#include "WireCellUtil/WireSchema.h"
#include <iostream>

using namespace WireCell;
using namespace WireCell::WireSchema;

int main(int argc, char* argv[])
{
    StoreDB storedb;
    auto& plane = get_append(storedb, 0, 0, 0, 0);

    Ray pitch(Point(0,0,0), Point(0,0,1));
    Ray bounds(Point(0,0,0), Point(0,10,10));

    std::cerr << "pitch=" << pitch << "\n"
              << "bounds="<< bounds << "\n";

    int nwires = generate(storedb, plane, pitch, bounds);

    std::cerr << "nwires=" << nwires << "\n";
    for (const auto& wire : storedb.wires) {
        std::cerr << wire.ident << "," << wire.channel << "," << wire.segment << " " << Ray(wire.tail,wire.head) << "\n";
    }

    Store store(std::make_shared<StoreDB>(storedb));
    try {
        validate(store);
    }
    catch (ValueError& err) {
        std::cerr << "invalid\n";
        return -1;
    }
    std::cerr << "valid\n";

    return 0;
}
