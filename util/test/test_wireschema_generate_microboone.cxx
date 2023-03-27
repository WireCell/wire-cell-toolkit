/*
  Generate microboone wires using numbers gleended from microboone
  larsoft wires dump.

  Take numbers from:

  wirecell-util wire-summary ~/opt/wire-cell-data/microboone-celltree-wires-v2.1.json.bz2 | jq '.[].anodes[].faces[].plan
es[].pray'
  wirecell-util wire-summary ~/opt/wire-cell-data/microboone-celltree-wires-v2.1.json.bz2 | jq '.[].anodes[].faces[].plan
es[].bb'

 */

#include <WireCellUtil/WireSchema.h>
#include <WireCellUtil/Testing.h>
#include <vector>
#include <iostream>

using namespace WireCell;
using namespace WireCell::WireSchema;



int main(int argc, char* argv[])
{
    std::vector<Ray> boundes = {
        Ray(Point( 0, -1155.1, 0.352608),
            Point( 0,  1174.5, 10369.6)),
        Ray(Point(-3, -1155.1, 0.352608),
            Point(-3,  1174.5, 10369.6)),
        Ray(Point(-6, -1155.3, 2.5),
            Point(-6,  1174.7, 10367.5))
    };
    std::vector<Ray> file_pitch_rays = {
        Ray(Point(0, 1173.0155, 2.919594),
            Point(0, 1170.4148, 4.420849440518711)),
        Ray(Point(-3, -1153.615, 2.919594),
            Point(-3, -1151.0147600126693, 4.420849440518712)),
        Ray(Point(-6, 9.7, 2.5),
            Point(-6, 9.7, 5.5))
    };
    std::vector<double> angles = {-60,60,0};
    std::vector<Ray> pitches(3);

    const double pi = 3.14159265;
    const double mbpitch = 3.0;

    for (size_t ind=0; ind<3; ++ind) {
        const auto& fpitch = file_pitch_rays[ind];

        const auto fpdir = ray_unit(fpitch);
        const auto fang = atan2(fpdir[1], fpdir[2]) * 180 / pi;
        const auto wang = angles[ind];

        const double rad = wang*pi/180;
        const Vector pdir(0, sin(rad), cos(rad));
        const Vector pvec = mbpitch*pdir;

        const auto& cen = fpitch.first;
        const Ray pray(cen, cen+pvec);
        pitches[ind] = pray;
        const auto gang = atan2(pdir[1], pdir[2]) * 180/pi;

        std::cerr << ind << ": angle:"
                  << " file=" << fang
                  << " want=" << wang
                  << " got=" << gang
                  << "\n";
    }

    int errors = 0;
    for (const auto& p : pitches) {
        const auto len = ray_length(p);
        const auto err = std::abs(len-3);
        std::cerr << "len: " << len << " err:" << err << std::endl;
        if (err > 1e-6) {
            ++errors;
        }
    }
    Assert(!errors);

    std::vector<size_t> expected_nwires = {2400, 2400, 3456};

    StoreDB storedb;
    errors = 0;
    for (int iplane=0; iplane<3; ++iplane) {
        auto& plane = get_append(storedb, iplane, 0, 0, 0);
        
        const Ray& bounds = boundes[iplane];
        const Ray& pitch = pitches[iplane];

        size_t iwire0 = plane.wires.size();
        size_t nwires = generate(storedb, plane, pitch, bounds);

        std::cerr << iplane << ": N: " << nwires << " =?= " << expected_nwires[iplane] << "\n";
        const auto& wf = storedb.wires[plane.wires[iwire0]];
        const auto& wl = storedb.wires[plane.wires[iwire0+nwires-1]];
        std::cerr << "wire["<<iwire0<<"]:\t" << Ray(wf.tail, wf.head) << "\n"
                  << "wire["<<iwire0+nwires-1<<"]:\t" << Ray(wl.tail, wl.head) << "\n";
        if (nwires != expected_nwires[iplane]) ++errors;
    }
    Assert(errors == 0);

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
