#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Units.h"

#include <iostream>
using namespace std;
using namespace WireCell;

void test_read(const char* filename)
{
    /// BIG FAT NOTE: don't directly load() from user code.  Use
    /// IWireSchema component named WireSchemaFile instead.
    cerr << "Test loading file twice, should read but once\n";
    auto store1 = WireSchema::load(filename);
    auto store2 = WireSchema::load(filename);

    Assert(store1.db().get() == store2.db().get());

    cerr << "Total store has:\n";
    cerr << "\t" << store1.detectors().size() << " detectors\n";
    cerr << "\t" << store1.anodes().size() << " anodes\n";
    cerr << "\t" << store1.faces().size() << " faces\n";
    cerr << "\t" << store1.planes().size() << " planes\n";
    cerr << "\t" << store1.wires().size() << " wires\n";

    Assert(store1.wires().size() > 0);

    auto& a0 = store1.anode(0);
    cerr << "Got anode 0 with " << store1.faces(a0).size() << " faces\n";

    const int bogusid = 0xdeadbeaf;
    try {
        store1.anode(bogusid);
        cerr << "Unexpectedly got bogus anode ID  " << bogusid << "\n";
        Assert(false);
    }
    catch (WireCell::Exception& e) {
        cerr << "Correctly caught exception with unlikely anode id " << bogusid << ":\n"
             << "----\n"
             << e.what() << "----\n";
    }

    for (const auto& anode : store1.anodes()) {
        for (const auto& face : store1.faces(anode)) {
            for (const auto& plane : store1.planes(face)) {
                auto bb = store1.bounding_box(plane);
                auto wp = store1.wire_pitch(plane);
                auto ex = wp.first.cross(wp.second);

                auto wires = store1.wires(plane);
                std::set<int> channels;
                for (const auto& w : wires) {
                    channels.insert(w.channel);
                }
                auto mm = std::minmax_element(channels.begin(), channels.end());

                cerr << "anode:" << anode.ident << " face:" << face.ident << " plane:" << plane.ident << ":\n"
                     << "\tdrift: " << ex << "\n"
                     << "\twire:  " << wp.first << "\n"
                     << "\tpitch: " << wp.second << "\n"
                     << "\tbb: " << bb.bounds() / units::cm << "cm\n"
                     << "\tchids: "<<channels.size()<<" in [" << *mm.first << "..." << *mm.second << "] reading "<<wires.size()<<" wires\n";
            }
            auto wplane = store1.planes(face).back();
            auto wires = store1.wires(wplane);

            const double zmax = wires.back().head.z();
            const double zmin = wires.front().head.z();
            const double dd = zmax - zmin;
            const int n = wires.size();
            cerr << "\tW plane: " << n << " wires, "
                 << "[" << zmin / units::m << "," << zmax / units::m << "]m, "
                 << "extent=" << dd / units::cm << "cm, pitch=" << dd / (n - 1) / units::mm << "mm\n";

        }
    }
}

int main(int argc, char* argv[])
{
    const char* filename = "microboone-celltree-wires-v2.1.json.bz2";
    if (argc > 1) {
        filename = argv[1];
    }
    test_read(filename);

    return 0;
}
