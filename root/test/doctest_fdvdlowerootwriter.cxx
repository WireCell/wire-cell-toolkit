// FdvdLowERootWriter: a hand-made FdvdLowEQLMatching tensor set (2 clusters, 2 flash groups, 3 pairs,
// 1 decision, 1 drift row) goes through the writer, which must pass it on unchanged and write
// trees whose rows read back as the tables (match vectors per arm x spec, t0 from the group time,
// x corrected with vcm).

#include "WireCellUtil/doctest.h"

#include "WireCellRoot/FdvdLowERootWriter.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "TFile.h"
#include "TTree.h"

#include <filesystem>
#include <vector>

using namespace WireCell;

static ITensor::pointer tab(const std::string& name, std::vector<double>& v, size_t nc)
{
    Configuration md;
    md["name"] = name;
    return std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{v.size() / nc, nc}, v.data(), md);
}

TEST_CASE("fdvd lowe root writer")
{
    const auto path = (std::filesystem::temp_directory_path() / "doctest_fdvdlowerootwriter.root").string();
    Root::FdvdLowERootWriter w;
    auto cfg = w.default_configuration();
    cfg["output_filename"] = path;
    cfg["event"] = 42;
    w.configure(cfg);

    // clusters: tree index, rep, nblob, Q, x_app, t_lo, t_hi, dt, mu, sigma
    std::vector<double> C{7, 100, 3, 150e3, 10.0, -100, 2000, 31.1, 120.0, 30.0,
                          9, 200, 5, 60e3, -50.0, -400, 1500, 31.1, NAN, NAN};
    // groups: window index, group index, t, tot, npd, pe[184]
    std::vector<double> G;
    for (int g = 0; g < 2; ++g) {
        G.insert(G.end(), {(double) g, (double) (g + 3), 100.0 * (g + 1), 500.0, 9.0});
        for (int ch = 0; ch < 184; ++ch) G.push_back(ch == g ? 500.0 : 0.0);
    }
    // pairs: k, g, r, ks, dc, npp, pass(none), pass(veto3)
    std::vector<double> P{0, 0, 1.0, 0.1, 20, 6, 1, 1, 0, 1, 2.0, 0.4, 80, 6, 1, 0, 1, 1, 0.9, 0.2, 30, 7, 1, 1};
    // decisions: arm, spec, k, g
    std::vector<double> D{1, 0, 0, 1};
    std::vector<double> R{100, 120.0, 30.0};
    auto tv = std::make_shared<ITensor::vector>();
    tv->push_back(tab("clusters", C, 10));
    tv->push_back(tab("groups", G, 189));
    tv->push_back(tab("pairs", P, 8));
    tv->push_back(tab("decisions", D, 4));
    tv->push_back(tab("drift", R, 3));
    Configuration md;
    md["producer"] = "FdvdLowEQLMatching";
    md["arms"][0] = "P5E100";
    md["arms"][1] = "P8E0";
    md["specs"][0] = "none";
    md["specs"][1] = "veto3";
    md["n_tree_clusters"] = 2;
    md["n_flash"] = 4;
    auto in = std::make_shared<Aux::SimpleTensorSet>(100, md, ITensor::shared_vector(tv));

    ITensorSet::pointer out;
    REQUIRE(w(in, out));
    CHECK(out == in);   // pass-through
    ITensorSet::pointer eos;
    REQUIRE(w(nullptr, eos));
    CHECK(!eos);

    TFile f(path.c_str());
    REQUIRE(!f.IsZombie());
    auto* tc = (TTree*) f.Get("T_cluster");
    auto* tm = (TTree*) f.Get("T_match");
    auto* te = (TTree*) f.Get("T_event");
    auto* tp = (TTree*) f.Get("T_pair");
    REQUIRE(tc);
    REQUIRE(tm);
    REQUIRE(te);
    REQUIRE(tp);
    CHECK(tc->GetEntries() == 2);
    CHECK(tm->GetEntries() == 1);
    CHECK(tp->GetEntries() == 2);   // pairs "qc": the pair of the 60 ke cluster is not written
    int event = -1;
    te->SetBranchAddress("event", &event);
    te->GetEntry(0);
    CHECK(event == 42);

    std::vector<int>* match_g = nullptr;
    std::vector<double>* match_x = nullptr;
    double Q = 0;
    tc->SetBranchAddress("match_g", &match_g);
    tc->SetBranchAddress("match_x", &match_x);
    tc->SetBranchAddress("Q", &Q);
    tc->GetEntry(0);
    CHECK(Q == 150e3);
    REQUIRE(match_g->size() == 4);   // 2 arms x 2 specs
    CHECK((*match_g)[2] == 1);       // arm 1, spec 0
    CHECK((*match_g)[0] == -1);
    CHECK((*match_x)[2] == doctest::Approx(10.0 + 0.160563 * 200.0));

    double t0 = 0;
    int tree_index = -1;
    tm->SetBranchAddress("t0", &t0);
    tm->SetBranchAddress("tree_index", &tree_index);
    tm->GetEntry(0);
    CHECK(t0 == 200.0);
    CHECK(tree_index == 7);
    f.Close();
    std::filesystem::remove(path);
}
