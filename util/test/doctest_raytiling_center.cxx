// Blob::add nudges each candidate corner "toward the blob center" before testing strip inclusion (62f4f1e9b).  The center
// is the mean over the old corners of Coordinates::pitch_location(corner, layer), and pitch_location() is 0 when the layer
// is one of the corner's own two layers.  In a U x V blob every corner is a U x V crossing, so the U and V centers come
// out as 0 and the nudge points toward ray 0: a real corner just inside a low edge is rejected, and prune() then cuts the
// blob short.  FD-VD's 30 degree wires (U/V pitch 7.65 mm, W 5.10 mm) put about half of all U x V crossings within 0.02
// pitch of a W boundary, so its tiling and the cascade cut lose area (pr_ml doc 34 sec 8.9).
//
// The geometry is FD-VD 1x8x14 v7 CRM 94 (dunevd10kt_3view_30deg_v7_refactored_1x8x14.json.bz2, default WireSchemaFile
// pitch correction): the box and wire 0 / 1 rays as gen/src/AnodeFace.cxx get_raypairs() builds them.

#include "WireCellUtil/RayTiling.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <vector>

using namespace WireCell;
using namespace WireCell::RayGrid;

static Coordinates crm94_coords()
{
    const Point ll(0, 3393.302110930631, 16518.749999999996);
    const Point ur(0, 5071.6978916064818, 18007.903670416323);
    const Point lr(0, ll.y(), ur.z()), ul(0, ur.y(), ll.z());
    ray_pair_vector_t raypairs;
    raypairs.push_back(ray_pair_t(Ray(ll, ul), Ray(lr, ur)));
    raypairs.push_back(ray_pair_t(Ray(ll, lr), Ray(ul, ur)));
    const std::vector<std::vector<Point>> w01 = {
        {{0, 3393.341192068684, 16526.31226062006}, {0, 3397.7370723271279, 16518.698284211274},
         {0, 3393.3412471421666, 16541.612273703642}, {0, 3406.5696285799727, 16518.699778760136}},
        {{0, 5067.2629286431666, 16518.698283649217}, {0, 5071.6588077363658, 16526.312260730745},
         {0, 5058.4303731800665, 16518.699777741858}, {0, 5071.6587511113412, 16541.612274709834}},
        {{0, 3393, 16521.299999999996}, {0, 5072, 16521.299999999996},
         {0, 3393, 16526.399999999998}, {0, 5072, 16526.399999999998}}};
    for (const auto& w : w01) {
        const Ray wray0(w[0], w[1]), wray1(w[2], w[3]);
        const auto pitvec = ray_vector(ray_pitch(wray0, wray1));
        raypairs.push_back(ray_pair_t(Ray(wray0.first - 0.5 * pitvec, wray0.second - 0.5 * pitvec),
                                      Ray(wray0.first + 0.5 * pitvec, wray0.second + 0.5 * pitvec)));
    }
    return Coordinates(raypairs);
}

// The upper half of the first bisection of an FD-VD L0 blob (evt1012 CRM 94): U [39,44) V [193,199) W [63,71).  Its
// polygon has 5 vertices and spans W 63..71; U ray 39 and V ray 193 cross at w = 62.9886, and the two vertices on W ray 63
// lie 0.0080 / 0.0076 pitch inside the V / U strips.
static Blob crm94_half(const Coordinates& coords, double nudge, bool exact_center)
{
    Blob blob;
    blob.add(coords, Strip{0, {0, 1}}, nudge, exact_center);
    blob.add(coords, Strip{1, {0, 1}}, nudge, exact_center);
    blob.add(coords, Strip{2, {39, 44}}, nudge, exact_center);
    blob.add(coords, Strip{3, {193, 199}}, nudge, exact_center);
    blob.add(coords, Strip{4, {63, 71}}, nudge, exact_center);
    return blob;
}

TEST_CASE("raytiling blob corners near a triple crossing")
{
    const auto coords = crm94_coords();
    // the near-triple crossing: U ray 39 x V ray 193 sits 0.0114 pitch below W ray 63
    const double w = coords.pitch_relative(coords.pitch_location({2, 39}, {3, 193}, 4), 4);
    CHECK(w == doctest::Approx(62.9886).epsilon(1e-5));

    SUBCASE("legacy add: the defect, kept as is")
    {
        // the 4-argument add with exact_center=false is the 3-argument add
        const Blob b3 = [&] {
            Blob b;
            for (const auto& s : {Strip{0, {0, 1}}, Strip{1, {0, 1}}, Strip{2, {39, 44}}, Strip{3, {193, 199}},
                                  Strip{4, {63, 71}}}) {
                b.add(coords, s, 0.01);
            }
            return b;
        }();
        const Blob b4 = crm94_half(coords, 0.01, false);
        REQUIRE(b3.corners().size() == b4.corners().size());
        for (size_t i = 0; i < b3.corners().size(); ++i) {
            const auto &c3 = b3.corners()[i], &c4 = b4.corners()[i];
            CHECK(c3.first.layer == c4.first.layer);
            CHECK(c3.first.grid == c4.first.grid);
            CHECK(c3.second.layer == c4.second.layer);
            CHECK(c3.second.grid == c4.second.grid);
        }
        // the two vertices on W ray 63 are pushed out, prune cuts W to [70, 71)
        blobs_t blobs{b4};
        CHECK(blobs[0].corners().size() == 3);
        prune(coords, blobs, 0.01);
        REQUIRE(blobs.size() == 1);
        CHECK(blobs[0].strips()[4].bounds.first == 70);
        CHECK(blobs[0].strips()[4].bounds.second == 71);
    }

    SUBCASE("exact_center keeps the polygon")
    {
        for (double nudge : {0.01, 1e-3}) {
            blobs_t blobs{crm94_half(coords, nudge, true)};
            CHECK(blobs[0].corners().size() == 5);
            prune(coords, blobs, nudge);
            REQUIRE(blobs.size() == 1);
            CHECK(blobs[0].strips()[2].bounds.first == 39);
            CHECK(blobs[0].strips()[2].bounds.second == 44);
            CHECK(blobs[0].strips()[3].bounds.first == 193);
            CHECK(blobs[0].strips()[3].bounds.second == 199);
            CHECK(blobs[0].strips()[4].bounds.first == 63);
            CHECK(blobs[0].strips()[4].bounds.second == 71);
        }
    }

    SUBCASE("make_blobs: exact_center reaches the tiling")
    {
        // one active strip per layer, the blob above: the legacy tiling prunes W to [70, 71), exact_center keeps it
        auto activities = [&] {
            activities_t acts;
            const int lo[5] = {0, 0, 39, 193, 63}, hi[5] = {1, 1, 44, 199, 71};
            for (int l = 0; l < 5; ++l) {
                acts.push_back(Activity(l, hi[l] - lo[l], 1.0, lo[l]));
            }
            return acts;
        }();
        const auto legacy = make_blobs(coords, activities, 0.01);
        const auto exact = make_blobs(coords, activities, 0.01, true);
        REQUIRE(legacy.size() == 1);
        REQUIRE(exact.size() == 1);
        CHECK(legacy[0].strips()[4].bounds.first == 70);
        CHECK(exact[0].strips()[4].bounds.first == 63);
        CHECK(exact[0].strips()[4].bounds.second == 71);
    }
}
