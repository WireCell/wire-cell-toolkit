#include "WireCellUtil/WireSchema.h"

#include <iostream>

using namespace WireCell;


int main(int argc, char* argv[])
{
    std::string fname = "microboone-celltree-wires-v2.1.json.bz2";
    if (argc > 1) {
        fname = argv[1];
    }
    auto store = WireSchema::load(fname);
    const auto anode = store.anodes()[0];
    for (const auto& face : store.faces(anode)) {

        for (const auto& plane : store.planes(face)) {

            auto [W,P] = store.wire_pitch(plane);
            auto pitch = store.mean_pitch(plane);
            
            std::cerr << "f" << face.ident << " p" << plane.ident
                      << " pitch="<<pitch
                      << " phat=" << P
                      << "\n";
        }
    }
    return 0;
}
