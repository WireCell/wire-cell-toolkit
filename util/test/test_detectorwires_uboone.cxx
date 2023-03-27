#include <WireCellUtil/DetectorWires.h>
#include <WireCellUtil/Testing.h>
#include <vector>
#include <iostream>

using namespace WireCell;
using namespace WireCell::WireSchema;

int main(int argc, char* argv[])
{
    std::vector<size_t> expected_nwires = {2400, 2400, 3456};
    StoreDB storedb;
    for (size_t ind=0; ind<3; ++ind) {
        auto& plane = get_append(storedb, ind, 0, 0, 0);
        size_t nwires = DetectorWires::plane(storedb, plane, DetectorWires::uboone[ind]);
        Assert(nwires == expected_nwires[ind]);
    }
    DetectorWires::flat_channels(storedb);
    Store store(std::make_shared<StoreDB>(storedb));
    try {
        validate(store);
    }
    catch (ValueError& err) {
        std::cerr << "invalid\n";
        return -1;
    }
    std::cerr << "valid\n";

    std::string fname = argv[0];
    fname += ".json.bz2";
    dump(fname.c_str(), store);
    
}
