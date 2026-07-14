/*
  Generate microboone wires using numbers gleended from microboone
  larsoft wires dump.
 */

#include "WireCellUtil/WireSchema.h"
#include "WireCellUtil/Logging.h"
#include "WireCellUtil/doctest.h"
#include <vector>

using namespace WireCell;
using namespace WireCell::WireSchema;

TEST_SUITE("wireschema generate microboone") {

TEST_CASE("pitch angles and lengths") {
    std::vector<Ray> file_pitch_rays = {
        Ray(Point(0, 1173.0155, 2.919594),
            Point(0, 1170.4148, 4.420849440518711)),
        Ray(Point(-3, -1153.615, 2.919594),
            Point(-3, -1151.0147600126693, 4.420849440518712)),
        Ray(Point(-6, 9.7, 2.5),
            Point(-6, 9.7, 5.5))
    };
    std::vector<double> angles = {-60, 60, 0};
    std::vector<Ray> pitches(3);

    const double pi = 3.14159265;
    const double mbpitch = 3.0;

    for (size_t ind = 0; ind < 3; ++ind) {
        const auto& fpitch = file_pitch_rays[ind];
        const auto fpdir = ray_unit(fpitch);
        const auto fang = atan2(fpdir[1], fpdir[2]) * 180 / pi;
        const auto wang = angles[ind];

        const double rad = wang * pi / 180;
        const Vector pdir(0, sin(rad), cos(rad));
        const Vector pvec = mbpitch * pdir;

        const auto& cen = fpitch.first;
        const Ray pray(cen, cen + pvec);
        pitches[ind] = pray;
        const auto gang = atan2(pdir[1], pdir[2]) * 180 / pi;

        spdlog::debug("{}: angle: file={} want={} got={}", ind, fang, wang, gang);
    }

    for (size_t ind = 0; ind < 3; ++ind) {
        const auto len = ray_length(pitches[ind]);
        const auto err = std::abs(len - 3.0);
        spdlog::debug("plane {}: len={} err={}", ind, len, err);
        CHECK(err <= 1e-6);
    }
}

TEST_CASE("generate wires and validate") {
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
    std::vector<double> angles = {-60, 60, 0};
    std::vector<Ray> pitches(3);

    const double pi = 3.14159265;
    const double mbpitch = 3.0;

    for (size_t ind = 0; ind < 3; ++ind) {
        const auto& fpitch = file_pitch_rays[ind];
        const double rad = angles[ind] * pi / 180;
        const Vector pdir(0, sin(rad), cos(rad));
        const Vector pvec = mbpitch * pdir;
        const auto& cen = fpitch.first;
        pitches[ind] = Ray(cen, cen + pvec);
    }

    std::vector<size_t> expected_nwires = {2400, 2400, 3456};

    StoreDB storedb;
    for (int iplane = 0; iplane < 3; ++iplane) {
        auto& plane = get_append(storedb, iplane, 0, 0, 0);

        const Ray& bounds = boundes[iplane];
        const Ray& pitch = pitches[iplane];

        size_t iwire0 = plane.wires.size();
        size_t nwires = generate(storedb, plane, pitch, bounds);

        const auto& wf = storedb.wires[plane.wires[iwire0]];
        const auto& wl = storedb.wires[plane.wires[iwire0 + nwires - 1]];
        spdlog::debug("{}: N={} expected={}", iplane, nwires, expected_nwires[iplane]);
        spdlog::debug("first wire tail=({},{},{}) head=({},{},{})",
                      wf.tail[0], wf.tail[1], wf.tail[2],
                      wf.head[0], wf.head[1], wf.head[2]);
        spdlog::debug("last wire  tail=({},{},{}) head=({},{},{})",
                      wl.tail[0], wl.tail[1], wl.tail[2],
                      wl.head[0], wl.head[1], wl.head[2]);
        CHECK(nwires == expected_nwires[iplane]);
    }

    Store store(std::make_shared<StoreDB>(storedb));
    bool valid = true;
    try {
        validate(store);
    }
    catch (ValueError&) {
        valid = false;
    }
    CHECK(valid);
}

}  // TEST_SUITE("wireschema generate microboone")
