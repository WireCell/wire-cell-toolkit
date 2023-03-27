
#include "WireCellIface/WirePlaneId.h"
#include "WireCellUtil/Testing.h"

#include <iostream>

using namespace WireCell;
using namespace std;

void test_roundtrip(WirePlaneLayer_t layer, int face, int apa)
{
    WirePlaneId wpid(layer, face, apa);
    Assert(wpid.layer() == layer);
    Assert(wpid.face() == face);
    Assert(wpid.apa() == apa);
}

int main()
{
    WirePlaneLayer_t layers[] = {kUnknownLayer, kUlayer, kVlayer, kWlayer};

    WirePlaneId u(kUlayer), v(kVlayer), w(kWlayer);

    cerr << "u.ident=" << u.ident() << " v.ident=" << v.ident() << " w.ident=" << w.ident() << endl;
    cerr << "u=" << u << " v=" << v << " w=" << w << endl;

    Assert(u.ident() == 1);
    Assert(v.ident() == 2);
    Assert(w.ident() == 4);

    Assert(u.ilayer() == 1);
    Assert(v.ilayer() == 2);
    Assert(w.ilayer() == 4);

    Assert(u.layer() == kUlayer);
    Assert(v.layer() == kVlayer);
    Assert(w.layer() == kWlayer);

    Assert(u.index() == 0);
    Assert(v.index() == 1);
    Assert(w.index() == 2);

    for (int ilayer = 0; ilayer < 4; ++ilayer) {
        WirePlaneLayer_t layer = layers[ilayer];
        for (int face = 0; face < 2; ++face) {
            for (int apa = 0; apa < 3; ++apa) {
                cerr << "Raw: " << ilayer << " " << layer << " " << face << " " << apa << endl;

                test_roundtrip(layer, face, apa);

                WirePlaneId wpid(layer, face, apa);

                cerr << "\twpid=" << wpid << endl;

                cerr << "\tident=" << wpid.ident() << " ilayer=" << wpid.ilayer() << " layer=" << wpid.layer()
                     << " index=" << wpid.index() << endl;

                if (ilayer) {
                    AssertMsg(wpid.valid(), "known layer should give true wpid");
                }
                else {
                    AssertMsg(!wpid.valid(), "unknown layer should give false wpid");
                }

                Assert(ilayer - 1 == wpid.index());
                Assert(face == wpid.face());
                Assert(apa == wpid.apa());
            }
        }
    }

    {
        std::vector<int> packed = {0, 8, 32, 56, 80, 88};
        for (auto p : packed) {
            WirePlaneId wpid(p);
            cerr << "wpid=" << wpid.ident() << " packed=" << p << " layer=" << wpid.layer() << " face=" << wpid.face() << " apa=" << wpid.apa() << "\n";
            Assert(p == wpid.ident());
        }
    }
    return 0;
}
