#include "WireCellMatch/PhotonLibraryModel.h"

#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/cnpy.h"
#include "WireCellUtil/doctest.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>

using WireCell::Match::PhotonLibraryModel;

TEST_CASE("photonlibrarymodel trilinear")
{
    // synthetic 2x3x4 grid, 2 channels: vis(ch0) = x + 10y + 100z in index
    // space (trilinear-exact), vis(ch1) = 7.
    const int nx = 2, ny = 3, nz = 4, nch = 2;
    std::vector<float> vis(nx * ny * nz * nch);
    for (int ix = 0; ix < nx; ++ix)
        for (int iy = 0; iy < ny; ++iy)
            for (int iz = 0; iz < nz; ++iz) {
                const size_t base = ((ix * ny + iy) * nz + iz) * nch;
                vis[base + 0] = ix + 10.f * iy + 100.f * iz;
                vis[base + 1] = 7.f;
            }

    const std::string dir = "/tmp";
    const std::string npy = dir + "/doctest_photlib_vis.npy";
    const std::string meta = dir + "/doctest_photlib_meta.json";
    cnpy::npy_save(npy, vis.data(),
                   {(size_t) nx, (size_t) ny, (size_t) nz, (size_t) nch});
    {
        std::ofstream f(meta);
        // origin (-5, 0, 10), steps (2, 4, 8) cm
        f << "{\"origin_cm\": [-5, 0, 10], \"step_cm\": [2, 4, 8], "
          << "\"n\": [" << nx << ", " << ny << ", " << nz << "], "
          << "\"nchan\": " << nch << ", "
          << "\"vis_npy\": \"doctest_photlib_vis.npy\"}";
    }

    PhotonLibraryModel plm(meta);
    CHECK(plm.nchan() == (size_t) nch);

    std::vector<double> v;

    // exact at grid nodes
    plm.visibilities(v, {-5., 0., 10.});
    CHECK(v[0] == doctest::Approx(0.0));
    CHECK(v[1] == doctest::Approx(7.0));
    plm.visibilities(v, {-3., 8., 34.});   // node (1, 2, 3)
    CHECK(v[0] == doctest::Approx(1. + 20. + 300.));

    // trilinear at a generic interior point: index (0.5, 0.25, 0.75)
    plm.visibilities(v, {-4., 1., 16.});
    CHECK(v[0] == doctest::Approx(0.5 + 10. * 0.25 + 100. * 0.75));
    CHECK(v[1] == doctest::Approx(7.0));

    // outside points clamp to the boundary
    plm.visibilities(v, {-100., -100., -100.});
    CHECK(v[0] == doctest::Approx(0.0));
    plm.visibilities(v, {100., 100., 100.});
    CHECK(v[0] == doctest::Approx(1. + 20. + 300.));

    // mask skips channels
    std::vector<unsigned int> mask{0, 1};
    plm.visibilities(v, {-4., 1., 16.}, &mask);
    CHECK(v[0] == doctest::Approx(0.0));
    CHECK(v[1] == doctest::Approx(7.0));

    std::remove(npy.c_str());
    std::remove(meta.c_str());
}

// Exercises the optional chan_pos_cm field (QLMatching's channel-order
// cross-check against m_opdets relies on has_positions()/position() -- see
// QLMatching.cxx's library-mode configure() block).
TEST_CASE("photonlibrarymodel chan_pos_cm")
{
    const int nx = 2, ny = 2, nz = 2, nch = 2;
    std::vector<float> vis(nx * ny * nz * nch, 1.f);

    const std::string dir = "/home/xqian/tmp";
    const std::string npy = dir + "/doctest_photlib_pos_vis.npy";
    const std::string meta_ok = dir + "/doctest_photlib_pos_meta_ok.json";
    const std::string meta_bad = dir + "/doctest_photlib_pos_meta_bad.json";
    cnpy::npy_save(npy, vis.data(),
                   {(size_t) nx, (size_t) ny, (size_t) nz, (size_t) nch});

    const std::string header =
        "{\"origin_cm\": [0, 0, 0], \"step_cm\": [1, 1, 1], "
        "\"n\": [2, 2, 2], \"nchan\": 2, "
        "\"vis_npy\": \"doctest_photlib_pos_vis.npy\", ";

    {
        std::ofstream f(meta_ok);
        f << header << "\"chan_pos_cm\": [[1.0, 2.0, 3.0], [-1.0, 0.0, 5.0]]}";
    }
    {
        // length mismatch (1 entry, nchan=2) -- must raise, not silently truncate.
        std::ofstream f(meta_bad);
        f << header << "\"chan_pos_cm\": [[1.0, 2.0, 3.0]]}";
    }

    PhotonLibraryModel plm_ok(meta_ok);
    CHECK(plm_ok.has_positions());
    CHECK(plm_ok.position(0).x() == doctest::Approx(1.0));
    CHECK(plm_ok.position(1).z() == doctest::Approx(5.0));

    CHECK_THROWS_AS(PhotonLibraryModel{meta_bad}, WireCell::ValueError);

    // absent chan_pos_cm (the format all library files shipped before this
    // check was added use) -- has_positions() stays false, no error.
    const std::string meta_nopos = dir + "/doctest_photlib_pos_meta_none.json";
    {
        std::ofstream f(meta_nopos);
        f << header.substr(0, header.size() - 2) << "}";
    }
    PhotonLibraryModel plm_none(meta_nopos);
    CHECK_FALSE(plm_none.has_positions());

    std::remove(npy.c_str());
    std::remove(meta_ok.c_str());
    std::remove(meta_bad.c_str());
    std::remove(meta_nopos.c_str());
}

// Optional check against a production library file whose meta JSON carries a
// "selfcheck" block of python-computed trilinear values (see
// pdvd/photlib/export_wct_photlib.py).  Enable with e.g.
//   PHOTLIB_META=/path/pdvd-photlib-vis-v5-128nm.json ./wcdoctest-match
TEST_CASE("photonlibrarymodel selfcheck file")
{
    const char* env = std::getenv("PHOTLIB_META");
    if (!env) return;

    PhotonLibraryModel plm(env);
    const auto meta = WireCell::Persist::load(env);
    REQUIRE(meta["selfcheck"].isArray());
    std::vector<double> v;
    for (const auto& sc : meta["selfcheck"]) {
        const WireCell::Point p(sc["p_cm"][0].asDouble(), sc["p_cm"][1].asDouble(),
                                sc["p_cm"][2].asDouble());
        plm.visibilities(v, p);
        CHECK(v.at(sc["ch"].asInt())
              == doctest::Approx(sc["vis"].asDouble()).epsilon(1e-6));
    }
}
