// ai-helper issue 39: RootFlatTreeVisitor::flatten_file is a lossless re-layout
// of every row tree of a per-event file into one flat, CAF-style entry.

#include "WireCellUtil/doctest.h"

#include "WireCellRoot/RootFlatTreeVisitor.h"

#include "TFile.h"
#include "TLeaf.h"
#include "TTree.h"

#include <cstdio>
#include <string>
#include <vector>

using namespace WireCell;

namespace {
    std::vector<double> leafvals(TTree* t, const char* name)
    {
        TLeaf* l = t->GetLeaf(name);
        REQUIRE_MESSAGE(l, "missing leaf " << name);
        std::vector<double> v;
        for (int i = 0; i < l->GetLen(); ++i) v.push_back(l->GetValue(i));
        return v;
    }
    std::string leaftype(TTree* t, const char* name)
    {
        TLeaf* l = t->GetLeaf(name);
        REQUIRE_MESSAGE(l, "missing leaf " << name);
        return l->GetTypeName();
    }
}

TEST_CASE("RootFlatTreeVisitor flattens row trees losslessly")
{
    const std::string fname = "doctest_root_flat_tree.root";
    {
        TFile f(fname.c_str(), "RECREATE");
        // Trun: one row, ints, a double, a string
        TTree trun("Trun", "Trun");
        Int_t runNo = 301, subRunNo = 55, eventNo = 17;
        Double_t scale = 0.5;
        std::string infile = "/path/reco1.root";
        trun.Branch("runNo", &runNo, "runNo/I");
        trun.Branch("subRunNo", &subRunNo, "subRunNo/I");
        trun.Branch("eventNo", &eventNo, "eventNo/I");
        trun.Branch("scale", &scale, "scale/D");
        trun.Branch("input_file", &infile);
        trun.Fill();
        // T_rows: three rows, a float scalar and a vector<float> (one empty)
        TTree trows("T_rows", "T_rows");
        Float_t x = 0;
        Int_t id = 0;
        std::vector<float> payload;
        trows.Branch("x", &x, "x/F");
        trows.Branch("id", &id, "id/I");
        trows.Branch("payload", &payload);
        const std::vector<std::vector<float>> pls = {{1.5f, 2.5f}, {}, {7.f, 8.f, 9.f}};
        for (int i = 0; i < 3; ++i) {
            x = 10.f + i;
            id = 100 + i;
            payload = pls[i];
            trows.Fill();
        }
        // T_nested: one row with a vector<vector<int>>
        TTree tnest("T_nested", "T_nested");
        std::vector<std::vector<int>> vv = {{1, 2}, {}, {3}};
        tnest.Branch("vv", &vv);
        tnest.Fill();
        // T_empty: zero rows
        TTree tempty("T_empty", "T_empty");
        Double_t e = 0;
        tempty.Branch("e", &e, "e/D");
        f.Write();
        f.Close();
    }

    std::string err;
    REQUIRE_MESSAGE(Root::RootFlatTreeVisitor::flatten_file(fname, "recTreeWireCell", "wc", false, err), err);

    TFile f(fname.c_str());
    auto* t = dynamic_cast<TTree*>(f.Get("recTreeWireCell"));
    REQUIRE(t);
    CHECK(t->GetEntries() == 1);
    REQUIRE(f.Get("T_rows"));  // row trees kept
    t->GetEntry(0);

    CHECK(leafvals(t, "wc.run") == std::vector<double>{301});
    CHECK(leafvals(t, "wc.subrun") == std::vector<double>{55});
    CHECK(leafvals(t, "wc.evt") == std::vector<double>{17});
    CHECK(leafvals(t, "wc.valid") == std::vector<double>{1});

    CHECK(leafvals(t, "wc.Trun..length") == std::vector<double>{1});
    CHECK(leafvals(t, "wc.Trun.scale") == std::vector<double>{0.5});
    CHECK(leaftype(t, "wc.Trun.scale") == "Double_t");
    {
        auto c = leafvals(t, "wc.Trun.input_file");
        std::string s;
        for (double ch : c) s.push_back((char) ch);
        CHECK(s == "/path/reco1.root");
        CHECK(leaftype(t, "wc.Trun.input_file") == "Char_t");
    }

    CHECK(leafvals(t, "wc.T_rows..length") == std::vector<double>{3});
    CHECK(leafvals(t, "wc.T_rows.x") == std::vector<double>{10, 11, 12});
    CHECK(leaftype(t, "wc.T_rows.x") == "Float_t");
    CHECK(leafvals(t, "wc.T_rows.id") == std::vector<double>{100, 101, 102});
    CHECK(leafvals(t, "wc.T_rows.payload..length") == std::vector<double>{2, 0, 3});
    CHECK(leafvals(t, "wc.T_rows.payload..idx") == std::vector<double>{0, 2, 2});
    CHECK(leafvals(t, "wc.T_rows.payload..totarraysize") == std::vector<double>{5});
    CHECK(leafvals(t, "wc.T_rows.payload") == std::vector<double>{1.5, 2.5, 7, 8, 9});

    CHECK(leafvals(t, "wc.T_nested..length") == std::vector<double>{1});
    CHECK(leafvals(t, "wc.T_nested.vv..length") == std::vector<double>{3});
    CHECK(leafvals(t, "wc.T_nested.vv..idx") == std::vector<double>{0});
    CHECK(leafvals(t, "wc.T_nested.vv.v..length") == std::vector<double>{2, 0, 1});
    CHECK(leafvals(t, "wc.T_nested.vv.v..idx") == std::vector<double>{0, 2, 2});
    CHECK(leafvals(t, "wc.T_nested.vv.v") == std::vector<double>{1, 2, 3});

    CHECK(leafvals(t, "wc.T_empty..length") == std::vector<double>{0});
    CHECK(leafvals(t, "wc.T_empty.e").empty());

    f.Close();
    std::remove(fname.c_str());
}

// PR 544 review, item 2: a missing file is an error, not a created file with a
// "valid" entry; an entry without Trun (no header RSE) is valid = 0.
TEST_CASE("RootFlatTreeVisitor refuses a missing file")
{
    const std::string fname = "doctest_root_flat_tree_missing.root";
    std::remove(fname.c_str());
    std::string err;
    CHECK_FALSE(Root::RootFlatTreeVisitor::flatten_file(fname, "recTreeWireCell", "wc", false, err));
    CHECK(err.find("no such file") != std::string::npos);
    FILE* fp = std::fopen(fname.c_str(), "r");
    CHECK_MESSAGE(fp == nullptr, "flatten_file created " << fname);
    if (fp) { std::fclose(fp); std::remove(fname.c_str()); }
}

TEST_CASE("RootFlatTreeVisitor without Trun is not valid")
{
    const std::string fname = "doctest_root_flat_tree_notrun.root";
    {
        TFile f(fname.c_str(), "RECREATE");
        TTree t("T_rows", "T_rows");
        Int_t id = 7;
        t.Branch("id", &id, "id/I");
        t.Fill();
        t.Write();
        f.Close();
    }
    std::string err;
    REQUIRE_MESSAGE(Root::RootFlatTreeVisitor::flatten_file(fname, "recTreeWireCell", "wc", false, err), err);
    TFile f(fname.c_str());
    auto* t = dynamic_cast<TTree*>(f.Get("recTreeWireCell"));
    REQUIRE(t);
    t->GetEntry(0);
    CHECK(leafvals(t, "wc.valid") == std::vector<double>{0});
    CHECK(leafvals(t, "wc.run") == std::vector<double>{-1});
    CHECK(leafvals(t, "wc.T_rows.id") == std::vector<double>{7});
    f.Close();
    std::remove(fname.c_str());
}
