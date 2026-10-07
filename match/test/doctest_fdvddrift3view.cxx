// Tests of the crop arithmetic of FdvdDriftRegressor3View (fdvd_sim doc 37): the footprint-masked crop of one
// view around its charge centroid, and the tick shift that moves an induction crop to the collection origin.

#include "WireCellUtil/doctest.h"

#include "WireCellMatch/FdvdDriftRegressor3View.h"

#include <vector>

using WireCell::Match::FdvdDriftRegressor3View;
using box_t = FdvdDriftRegressor3View::box_t;

TEST_CASE("fdvd drift 3view crop: footprint mask, centroid, window")
{
    const int nrow = 40;
    const size_t ntick = 200;
    const int nch = 16, ntk = 32;
    std::vector<float> img((size_t) nrow * ntick, 0.0f);
    auto at = [&](int r, int t) -> float& { return img[(size_t) r * ntick + t]; };
    at(20, 100) = 10.0f;    // inside the footprint
    at(21, 102) = 30.0f;    // inside
    at(30, 100) = 1000.0f;  // outside: must not enter the centroid or the crop
    std::vector<float> out((size_t) nch * ntk, -1.0f);
    long c0 = -99, t0 = -99;
    // box wires [20, 22), ticks [100, 103); dilation 1 wire, 2 ticks
    const bool ok = FdvdDriftRegressor3View::view_crop(img.data(), nrow, ntick, {box_t{20, 22, 100, 103}}, 1, 2, nch, ntk,
                                                       out.data(), c0, t0);
    REQUIRE(ok);
    // centroid: wire (20*10 + 21*30)/40 = 20.75 -> 21; tick (100*10 + 102*30)/40 = 101.5 -> 102 (half to even)
    CHECK(c0 == 21 - nch / 2);
    CHECK(t0 == 102 - ntk / 2);
    CHECK(out[(20 - c0) * ntk + (100 - t0)] == 10.0f);
    CHECK(out[(21 - c0) * ntk + (102 - t0)] == 30.0f);
    float sum = 0;
    for (float v : out) sum += v;
    CHECK(sum == 40.0f);
}

TEST_CASE("fdvd drift 3view crop: empty footprint, window clipped at the frame edge")
{
    const int nrow = 10;
    const size_t ntick = 50;
    const int nch = 16, ntk = 32;
    std::vector<float> img((size_t) nrow * ntick, 0.0f);
    std::vector<float> out((size_t) nch * ntk, 7.0f);
    long c0 = 5, t0 = 5;
    CHECK(!FdvdDriftRegressor3View::view_crop(img.data(), nrow, ntick, {box_t{2, 4, 10, 12}}, 2, 12, nch, ntk, out.data(), c0, t0));
    CHECK(c0 == 0);
    CHECK(t0 == 0);
    for (float v : out) REQUIRE(v == 0.0f);

    img[1 * ntick + 3] = 5.0f;   // near the corner: the window starts before the frame
    REQUIRE(FdvdDriftRegressor3View::view_crop(img.data(), nrow, ntick, {box_t{1, 2, 3, 4}}, 2, 12, nch, ntk, out.data(), c0, t0));
    CHECK(c0 == 1 - nch / 2);
    CHECK(t0 == 3 - ntk / 2);
    CHECK(out[(1 - c0) * ntk + (3 - t0)] == 5.0f);
}

TEST_CASE("fdvd drift 3view tick shift")
{
    const int nch = 2, ntk = 8;
    std::vector<float> src((size_t) nch * ntk), dst((size_t) nch * ntk, -1.0f);
    for (int i = 0; i < nch * ntk; ++i) src[i] = (float) (i + 1);
    FdvdDriftRegressor3View::shift_ticks(src.data(), 0, nch, ntk, dst.data());
    CHECK(dst == src);
    FdvdDriftRegressor3View::shift_ticks(src.data(), 3, nch, ntk, dst.data());   // later: content moves right
    for (int r = 0; r < nch; ++r) {
        for (int t = 0; t < ntk; ++t) CHECK(dst[r * ntk + t] == (t < 3 ? 0.0f : src[r * ntk + t - 3]));
    }
    FdvdDriftRegressor3View::shift_ticks(src.data(), -2, nch, ntk, dst.data());  // earlier: content moves left
    for (int r = 0; r < nch; ++r) {
        for (int t = 0; t < ntk; ++t) CHECK(dst[r * ntk + t] == (t >= ntk - 2 ? 0.0f : src[r * ntk + t + 2]));
    }
    FdvdDriftRegressor3View::shift_ticks(src.data(), ntk + 1, nch, ntk, dst.data());
    for (float v : dst) CHECK(v == 0.0f);
}
