#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/Stats.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Math.h"

using namespace WireCell;
using namespace WireCell::WireSchema;
using namespace WireCell::String;

void dump(const Stats& s, const std::string& lab, int nbins=10)
{
    Binning bins(nbins,0,0);
    auto vec = s.histogram(bins);
    std::cerr << "\t" << lab << "\t: " << s << "\n\t\t: " << bins << ":";
    for (auto v : vec) std::cerr << " " << v;
    std::cerr << "\n";
}

void dump(const Store& store)
{
    const auto& detectors = store.detectors();
    // std::cerr << detectors.size() << " detectors\n";
    for (const auto& detector : detectors) {

        const auto anodes = store.anodes(detector);
        // std::cerr << "detector " << detector.ident << ", " << anodes.size() << " anodes:\n";

        for (const auto& anode : anodes) {

            const auto faces = store.faces(anode);
            // std::cerr << "anode " << anode.ident << ", " << faces.size() << " faces:\n";

            for (const auto& face : faces) {

                const auto planes = store.planes(face);
                // std::cerr << "face " << face.ident << ",  " << planes.size() << " planes:\n";
            
                for (const auto& plane : planes) {

                    Stats X{true}, P{true}, DX{true}, DY{true}, DZ{true};
                    Ray rprev;
                    Wire wprev;
                    Ray rprev2;
                    Wire wprev2;
                    int pnum=0;

                    const auto wires = store.wires(plane);
                    // std::cerr << "plane " << plane.ident << ", " << wires.size() << " wires:\n";

                    for (const auto& wire : wires) {
                        const Ray rnext(wire.tail, wire.head);

                        X(wire.tail.x()/units::mm);
                        X(wire.head.x()/units::mm);

                        if (pnum) {
                            const auto pray = ray_pitch(rprev, rnext);
                            const auto pvec = ray_vector(pray);
                            const double pmag = pvec.magnitude()/units::mm;

                            P(pmag);
                            const auto pdir = pvec.norm();
                            const double dx = 180/pi*acos(pdir.x());
                            const double dy = 180/pi*acos(pdir.y());
                            const double dz = 180/pi*acos(pdir.z());
                            DX(dx);
                            DY(dy);
                            DZ(dz);
                            bool do_dump=false;
                            if (pmag < 1e-4) {
                                std::cerr << "TINY PITCH:\n";
                                do_dump=true;
                            }
                            if (std::abs(dy-52.542) < 2 or std::abs(dy-1.739) < 2) {
                                std::cerr << "WEIRD ANGLE:\n";
                                do_dump=true;
                            }
                            if (do_dump) {
                                std::cerr << " aid=" << anode.ident
                                          << " fid=" << face.ident
                                          << " pid=" << plane.ident 
                                          << " wid=" << wire.ident
                                          << " chn=" << wire.channel
                                          << " <-- "
                                          << " wid=" << wprev.ident
                                          << " chn=" << wprev.channel
                                          << " <-- "
                                          << " wi2=" << wprev2.ident
                                          << " ch2=" << wprev2.channel
                                          << "\n";
                                std::cerr << "\tnext: " << rnext << "\n";
                                std::cerr << "\tprev: " << rprev << "\n";
                                std::cerr << "\tpre2: " << rprev2 << "\n";
                            }
                        }
                        rprev2 = rprev;
                        wprev2 = wprev;
                        rprev = rnext;
                        wprev = wire;
                        ++pnum;
                    }
                    std::cerr << "d:" << detector.ident
                              << " a:" << anode.ident
                              << " f:" << face.ident
                              << " p:" << plane.ident
                              << " nwires:" << wires.size() << ":\n";
                    dump(X, "X(mm)");
                    dump(DX, "DX(deg)");
                    dump(DY, "DY(deg)");
                    dump(DZ, "DZ(deg)");
                    dump(P, "P(mm)");
                }
            }
        }
    }
}

// assume same order
void compare(const Store& a, const Store& b)
{
    const auto& detectors = a.detectors();
    for (const auto& detector : detectors) {
        const auto anodes = a.anodes(detector);
        for (const auto& anode : anodes) {
            const auto faces = a.faces(anode);
            for (const auto& face : faces) {
                for (int iplane : face.planes) {
                    const auto aplane = a.planes()[iplane];
                    const auto bplane = b.planes()[iplane];
                    Assert(aplane.ident == bplane.ident);

                    const auto was = a.wires(aplane);
                    const auto wbs = b.wires(bplane);
                    Assert(was.size() == wbs.size());

                    Stats A{true}, P{true}, H{true}, T{true};

                    for (size_t wind=0; wind<was.size(); ++wind) {
                        const auto& wa = was[wind];
                        const auto& wb = wbs[wind];

                        const Ray ra(wa.tail, wa.head);
                        const Ray rb(wb.tail, wb.head);

                        const auto pray = ray_pitch(ra, rb);
                        const auto pvec = ray_vector(pray);
                        const double pmag = pvec.magnitude()/units::mm;
                        const double ang = 180/pi*acos(ray_unit(ra).dot(ray_unit(rb)));
                        P(pmag/units::mm);
                        A(ang);
                        H((wa.head - wb.head).magnitude()/units::mm);
                        T((wa.tail - wb.tail).magnitude()/units::mm);
                    }
                    std::cerr << "d:" << detector.ident
                              << " a:" << anode.ident
                              << " f:" << face.ident
                              << " p:" << aplane.ident
                              << " nwires:" << was.size() << ":\n";
                    dump(T, "T(mm)");
                    dump(H, "H(mm)");
                    dump(P, "P(mm)");
                    dump(A, "A(deg)");
                }
            }
        }
    }
                    
}


int main(int argc, char* argv[])
{
    std::string fname = "microboone-celltree-wires-v2.1.json.bz2";
    if (argc > 1) {
        fname = argv[1];
    }

    std::vector<Store> stores;
    for (Correction cor=Correction::load; cor<Correction::ncorrections; ++cor) {
        auto store = WireSchema::load(fname.c_str(), (Correction)cor);
        std::cerr << fname << " correction: " << (int)cor << ":\n";
        dump(store);
        stores.push_back(store);
        std::cerr << "\n";
    }

    std::cerr << "size of correction:\n";
    // We use the order-corrected as the baseline instead of "raw" so
    // that we gain wire correspondence.
    compare(stores[1], stores[3]);

    return 0;
}
