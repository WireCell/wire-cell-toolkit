// Tests of FdvdAdjOpFlashFinder, the AdjOpHits (duneopdet AdjOpHitsUtils::CalcAdjOpHits)
// port: the plane assignment, the scans and their kept quirks, and the component's
// tensor interface (fdvd_sim doc 21).

#include "WireCellUtil/doctest.h"

#include "WireCellFlash/FdvdAdjOpFlashFinder.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include <array>
#include <filesystem>
#include <fstream>

using namespace WireCell;
using Flash::FdvdAdjParams;

namespace {
    // 4 OpDets: 0, 1, 2 on the cathode plane 100 cm apart in z, 3 on membrane +y.
    std::vector<std::vector<double>> dist4()
    {
        std::vector<std::vector<double>> d(4, std::vector<double>(4, 0.0));
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) d[a][b] = 100.0 * std::abs(a - b);
        for (int a = 0; a < 3; ++a) d[a][3] = d[3][a] = 50.0;
        return d;
    }
    const std::vector<int> plane4{0, 0, 0, 1};

    FdvdAdjParams wide()
    {
        FdvdAdjParams p;
        p.min_time_us = 1.0;
        p.max_time_us = 1.6;
        p.radius_cm = 150.0;
        p.nhit = 0;
        p.pe = 0.0;
        p.trigger_pe = 0.0;
        return p;
    }
}  // namespace

TEST_CASE("fdvd adj: plane assignment (GetOpHitPlane order, open +-0.1 cm band)")
{
    FdvdAdjParams p;
    CHECK(Flash::fdvd_adj_plane(-327.5, 0, 0, p) == 0);
    CHECK(Flash::fdvd_adj_plane(-327.41, 0, 0, p) == 0);
    CHECK(Flash::fdvd_adj_plane(-327.4, 0, 0, p) == -1);   // the band is open
    CHECK(Flash::fdvd_adj_plane(0, 743.302, 0, p) == 1);
    CHECK(Flash::fdvd_adj_plane(0, -743.302, 0, p) == 2);
    CHECK(Flash::fdvd_adj_plane(0, 0, 2188.38, p) == 3);
    CHECK(Flash::fdvd_adj_plane(0, 0, -96.5, p) == 4);
    CHECK(Flash::fdvd_adj_plane(-327.5, 743.302, 0, p) == 0);   // cathode tested first
}

TEST_CASE("fdvd adj: brightest seed collects same-plane hits inside radius and window")
{
    // t [ns], pe, opdet; hit 0 is the time-first hit (position 0 of the sorted list)
    std::vector<double> t{0.0, 1000.0, 1500.0, 1600.0, 2700.0, 1200.0};
    std::vector<double> pe{1.0, 10.0, 3.0, 2.0, 4.0, 5.0};
    std::vector<int> od{0, 0, 1, 2, 0, 3};
    auto cl = Flash::fdvd_calc_adj_ophits(t, pe, od, dist4(), plane4, wide());
    // sorted: 0(0) 1(1000) 5(1200) 2(1500) 3(1600) 4(2700)
    // hit 0 seeds first: forward hit 1 is brighter -> demoted
    // hit 1 seeds: forward 5 (other plane, skip), 2 (100 cm, in), 3 (200 cm > 150, out), 4 (|dt| 1700 > 1600, stop);
    //   backward: nothing, sorted position 0 (hit 0, 1000 ns <= 1 us) is never reached
    // hit 5 (membrane) alone; hit 3: backward finds the claimed, brighter hit 2 (duplicates on) -> demoted;
    // hit 4: backward hit 3 at 1100 ns > 1 us -> alone
    REQUIRE(cl.size() == 3);
    CHECK(cl[0] == Flash::FdvdAdjCluster{1, 2});
    CHECK(cl[1] == Flash::FdvdAdjCluster{5});
    CHECK(cl[2] == Flash::FdvdAdjCluster{4});
}

TEST_CASE("fdvd adj: without duplicates a claimed hit stays with its first flash")
{
    std::vector<double> t{1000.0, 1500.0, 1600.0};
    std::vector<double> pe{10.0, 3.0, 2.0};
    std::vector<int> od{0, 1, 2};
    auto p = wide();
    p.hit_duplicates = false;
    auto cl = Flash::fdvd_calc_adj_ophits(t, pe, od, dist4(), plane4, p);
    REQUIRE(cl.size() == 2);
    CHECK(cl[0] == Flash::FdvdAdjCluster{0, 1});
    CHECK(cl[1] == Flash::FdvdAdjCluster{2});
    p.hit_duplicates = true;      // the claimed, brighter hit 1 now demotes the seed hit 2
    cl = Flash::fdvd_calc_adj_ophits(t, pe, od, dist4(), plane4, p);
    REQUIRE(cl.size() == 1);
    CHECK(cl[0] == Flash::FdvdAdjCluster{0, 1});
}

TEST_CASE("fdvd adj: window edges are inclusive (break on |dt| > window), backward demotion un-claims")
{
    std::vector<double> t{0.0, 500.0, 2100.0};
    std::vector<double> pe{1.0, 2.0, 5.0};
    std::vector<int> od{0, 0, 0};
    auto cl = Flash::fdvd_calc_adj_ophits(t, pe, od, dist4(), plane4, wide());
    // hit 0: forward 1 brighter -> demoted.  hit 1: forward 2 at exactly 1600 ns, brighter -> demoted.
    // hit 2: backward window 1000 ns: 1 is 1600 away -> stop.  -> {2}; then 0 and 1 stay unclaimed but were visited.
    REQUIRE(cl.size() == 1);
    CHECK(cl[0] == Flash::FdvdAdjCluster{2});
    auto p = wide();
    p.max_time_us = 1.5999;
    cl = Flash::fdvd_calc_adj_ophits(t, pe, od, dist4(), plane4, p);
    REQUIRE(cl.size() == 2);
    CHECK(cl[0] == Flash::FdvdAdjCluster{1});
    CHECK(cl[1] == Flash::FdvdAdjCluster{2});
}

TEST_CASE("fdvd adj: no-plane hits (-1) count as one plane; nhit and trigger_pe cuts")
{
    std::vector<std::vector<double>> d(2, std::vector<double>(2, 10.0));
    std::vector<int> pl{-1, -1};
    std::vector<double> t{0.0, 100.0, 200.0};
    std::vector<double> pe{3.0, 1.0, 0.5};
    std::vector<int> od{0, 1, 0};
    auto p = wide();
    auto cl = Flash::fdvd_calc_adj_ophits(t, pe, od, d, pl, p);
    REQUIRE(cl.size() == 1);
    CHECK(cl[0] == Flash::FdvdAdjCluster{0, 1, 2});
    p.pe = 0.8;           // member threshold drops hit 2
    p.nhit = 3;
    cl = Flash::fdvd_calc_adj_ophits(t, pe, od, d, pl, p);
    CHECK(cl.empty());
    p.nhit = 0;
    p.trigger_pe = 5.0;   // no seed
    cl = Flash::fdvd_calc_adj_ophits(t, pe, od, d, pl, p);
    CHECK(cl.empty());
}

namespace {
    std::string write_json(const std::string& name, const std::string& body)
    {
        auto p = std::filesystem::temp_directory_path() / name;
        std::ofstream f(p);
        f << body;
        return p.string();
    }
}  // namespace

TEST_CASE("fdvd adj: component tensors")
{
    // OpDets 0, 1 on the cathode (x = -3275 mm), 1000 mm apart; OpChannel 10*od + {0,1}
    const auto geom = write_json("doctest_fdvdadj_geom.json",
                                 "{\"opdets\":[{\"opdet\":0,\"x\":-3275.0,\"y\":0.0,\"z\":0.0},"
                                 "{\"opdet\":1,\"x\":-3275.0,\"y\":0.0,\"z\":1000.0}]}");
    const auto chmap = write_json("doctest_fdvdadj_chmap.json",
                                  "{\"channels\":[{\"opch\":0,\"opdet\":0},{\"opch\":1,\"opdet\":0},"
                                  "{\"opch\":10,\"opdet\":1},{\"opch\":11,\"opdet\":1}]}");
    Flash::FdvdAdjOpFlashFinder ff;
    auto cfg = ff.default_configuration();
    cfg["nchan"] = 2;
    cfg["geom_file"] = geom;
    cfg["channel_map_file"] = chmap;
    cfg["min_time_us"] = 1.0;
    cfg["max_time_us"] = 1.6;
    cfg["radius_cm"] = 600.0;
    cfg["nhit"] = 0;
    cfg["pe"] = 0.0;
    cfg["trigger_pe"] = 0.0;
    ff.configure(cfg);
    // rows {opch, peak, width, area, amplitude, PE, start, flash_id, f2t}
    std::vector<std::array<double, 9>> rows{{1, 1100, 50, 0, 0, 4.0, 1000, -1, 0},
                                           {10, 1300, 50, 0, 0, 9.0, 1200, -1, 0},
                                           {0, 1400, 50, 0, 0, 2.0, 1300, -1, 0}};
    std::vector<double> flat;
    for (const auto& r : rows) flat.insert(flat.end(), r.begin(), r.end());
    Configuration md;
    md["name"] = "ophits";
    auto tv = std::make_shared<ITensor::vector>();
    tv->push_back(std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{rows.size(), (size_t) 9}, flat.data(), md));
    auto in = std::make_shared<Aux::SimpleTensorSet>(7, Configuration{}, ITensor::shared_vector(tv));
    ITensorSet::pointer out;
    REQUIRE(ff(in, out));
    REQUIRE(out);
    CHECK(out->ident() == 7);
    ITensor::pointer opflash, memb;
    for (const auto& t : *out->tensors()) {
        if (t->metadata()["name"].asString() == "opflash") opflash = t;
        if (t->metadata()["name"].asString() == "adj_membership") memb = t;
    }
    REQUIRE(opflash);
    REQUIRE(opflash->shape() == ITensor::shape_t{1, 3});
    const double* f = (const double*) opflash->data();
    // hit 0 is demoted by hit 1 and, at sorted position 0, never reached by hit 1's backward scan
    CHECK(f[0] == 1200.0);          // start time of the max-PE member (hit 1)
    CHECK(f[1] == 2.0);             // OpDet 0: hit 2 only
    CHECK(f[2] == 9.0);             // OpDet 1
    REQUIRE(memb);
    CHECK(memb->shape() == ITensor::shape_t{2, 2});
    ITensorSet::pointer eos;
    REQUIRE(ff(nullptr, eos));
    CHECK(eos == nullptr);
}
