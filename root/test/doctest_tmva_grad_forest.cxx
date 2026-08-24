/** TmvaGradForest must reproduce TMVA::Reader bit-for-bit on the class of
 *  weight file it accepts (sbnd_xin/docs/76 round 2).  A small forest is
 *  written in the xgboost2TMVA layout (one line, the attribute spellings
 *  and number formats of the production files), booked by both readers and
 *  evaluated on inputs that include exact cut ties, cType=0 nodes, NaN and
 *  extreme values.  Also checks the refusals for files outside that class.
 */
#include "WireCellRoot/TmvaGradForest.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/doctest.h"

#include "TMVA/Reader.h"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

using namespace WireCell;
using WireCell::Root::TmvaGradForest;

namespace {

    struct XNode {
        int ivar{-1};
        float cut{0};
        int ctype{1};
        float res{0};
        int ntype{0};
        int left{-1}, right{-1};
    };

    struct XTree {
        std::vector<XNode> nodes;
        int root{0};
    };

    // xgboost2TMVA prints cuts/responses with %.9g-ish formats; mix in the
    // "0.0e+00" spelling it uses for zero.
    std::string fmt(float v)
    {
        if (v == 0.f) return "0.0e+00";
        char buf[64];
        std::snprintf(buf, sizeof buf, "%.9g", v);
        return buf;
    }

    void emit(std::string& out, const XTree& t, int idx, int depth, char pos)
    {
        const XNode& n = t.nodes[idx];
        out += "<Node Cut=\"" + fmt(n.cut) + "\" IVar=\"" + std::to_string(n.ivar) + "\" NCoef=\"0\" cType=\"" +
               std::to_string(n.ctype) + "\" depth=\"" + std::to_string(depth) + "\" nType=\"" +
               std::to_string(n.ntype) + "\" pos=\"" + pos + "\" purity=\"0.0e+00\" res=\"" + fmt(n.res) +
               "\" rms=\"0.0e+00\"";
        if (n.ntype != 0) {
            out += " />";
            return;
        }
        out += ">";
        emit(out, t, n.left, depth + 1, 'l');
        emit(out, t, n.right, depth + 1, 'r');
        out += "</Node>";
    }

    int grow(XTree& t, std::mt19937& rng, const std::vector<float>& cutpool, int nvar, int depth, int maxdepth)
    {
        XNode n;
        std::uniform_real_distribution<float> resd(-0.2f, 0.2f);
        std::uniform_int_distribution<int> coin(0, 1);
        const int idx = static_cast<int>(t.nodes.size());
        t.nodes.push_back(n);
        if (depth >= maxdepth || (depth > 1 && coin(rng))) {
            t.nodes[idx].ntype = -99;
            t.nodes[idx].res = resd(rng);
            return idx;
        }
        t.nodes[idx].ivar = std::uniform_int_distribution<int>(0, nvar - 1)(rng);
        t.nodes[idx].cut = cutpool[std::uniform_int_distribution<size_t>(0, cutpool.size() - 1)(rng)];
        t.nodes[idx].ctype = coin(rng);
        const int l = grow(t, rng, cutpool, nvar, depth + 1, maxdepth);
        const int r = grow(t, rng, cutpool, nvar, depth + 1, maxdepth);
        t.nodes[idx].left = l;
        t.nodes[idx].right = r;
        return idx;
    }

    std::string make_xml(const std::vector<std::string>& names, const std::vector<XTree>& trees,
                         const std::string& boost = "Grad", const std::string& anal = "1")
    {
        std::string out = "<?xml version=\"1.0\"?><MethodSetup Method=\"BDT::BDT\"><Variables NVar=\"" +
                          std::to_string(names.size()) + "\">";
        for (size_t i = 0; i < names.size(); ++i) {
            out += "<Variable Expression=\"" + names[i] + "\" Internal=\"" + names[i] + "\" Label=\"" + names[i] +
                   "\" Max=\"0.0e+00\" Min=\"0.0e+00\" Title=\"" + names[i] + "\" Type=\"F\" Unit=\"\" VarIndex=\"" +
                   std::to_string(i) + "\" />";
        }
        out += "</Variables><GeneralInfo><Info name=\"Creator\" value=\"xgboost2TMVA\" /><Info name=\"AnalysisType\" "
               "value=\"Classification\" /></GeneralInfo><Options><Option modified=\"No\" "
               "name=\"NodePurityLimit\">5.00e-01</Option><Option modified=\"Yes\" name=\"BoostType\">" +
               boost + "</Option></Options><Weights AnalysisType=\"" + anal + "\" NTrees=\"" +
               std::to_string(trees.size()) + "\">";
        for (size_t i = 0; i < trees.size(); ++i) {
            out += "<BinaryTree boostWeight=\"1.0e+00\" itree=\"" + std::to_string(i) + "\" type=\"DecisionTree\">";
            emit(out, trees[i], trees[i].root, 0, 's');
            out += "</BinaryTree>";
        }
        out += "</Weights></MethodSetup>";
        return out;
    }

    std::string write_tmp(const std::string& xml, const char* tag)
    {
        auto p = std::filesystem::temp_directory_path() / (std::string("doctest_tmva_grad_forest_") + tag + ".xml");
        std::ofstream(p) << xml;
        return p.string();
    }

}  // namespace

TEST_CASE("tmva grad forest matches TMVA::Reader bit for bit")
{
    const int nvar = 6;
    std::vector<std::string> names;
    for (int i = 0; i < nvar; ++i) names.push_back("v" + std::to_string(i));

    std::mt19937 rng(20260823);
    // Cut pool: values that the inputs below will hit EXACTLY (ties) plus
    // fractional ones printed with 9 significant digits.
    std::vector<float> cutpool = {0.5f, 1.5f, 2.5f, -1.f, 28.5f, 42.1889267f, 383.108154f, 0.00953744911f, 1e-6f};
    std::vector<XTree> trees(40);
    for (auto& t : trees) t.root = grow(t, rng, cutpool, nvar, 0, 6);

    const std::string xml = make_xml(names, trees);
    const std::string path = write_tmp(xml, "ok");

    TmvaGradForest fast;
    fast.load(path, names);
    CHECK(fast.nvars() == static_cast<size_t>(nvar));
    CHECK(fast.ntrees() == trees.size());
    size_t nn = 0;
    for (auto& t : trees) nn += t.nodes.size();
    CHECK(fast.nnodes() == nn);

    std::vector<float> vals(nvar, 0.f);
    TMVA::Reader reader("!V:Silent");
    for (int i = 0; i < nvar; ++i) reader.AddVariable(names[i], &vals[i]);
    reader.BookMVA("MyBDT", path);

    // Inputs: random, exact cut values (ties on >=), signed extremes, NaN.
    std::uniform_real_distribution<float> ud(-5.f, 50.f);
    std::uniform_int_distribution<int> mode(0, 3);
    int nties = 0;
    for (int trial = 0; trial < 4000; ++trial) {
        for (int i = 0; i < nvar; ++i) {
            switch (mode(rng)) {
            case 0:
                vals[i] = cutpool[std::uniform_int_distribution<size_t>(0, cutpool.size() - 1)(rng)];
                ++nties;
                break;
            case 1: vals[i] = (trial % 7 == 0) ? 1e30f : -1e30f; break;
            default: vals[i] = ud(rng);
            }
        }
        const double want = reader.EvaluateMVA("MyBDT");
        const double got = fast.evaluate(vals.data());
        REQUIRE(got == want);
    }
    CHECK(nties > 0);

    // NaN in any variable: TMVA::Reader returns -999 (and logs an error).
    vals[2] = std::nanf("");
    CHECK(fast.evaluate(vals.data()) == -999);
    CHECK(reader.EvaluateMVA("MyBDT") == -999);
    std::filesystem::remove(path);
}

TEST_CASE("tmva grad forest refuses files outside the mirrored class")
{
    std::vector<std::string> names = {"a", "b"};
    std::mt19937 rng(1);
    std::vector<float> cutpool = {0.5f, 1.5f};
    std::vector<XTree> trees(3);
    for (auto& t : trees) t.root = grow(t, rng, cutpool, 2, 0, 3);

    TmvaGradForest f;
    // the good file loads
    CHECK_NOTHROW(f.load_string(make_xml(names, trees), names));
    // variable name / count mismatch
    CHECK_THROWS_AS(f.load_string(make_xml(names, trees), {"a", "c"}), ValueError);
    CHECK_THROWS_AS(f.load_string(make_xml(names, trees), {"a"}), ValueError);
    // AdaBoost forests take TMVA's weighted-vote branch, not GradBoost
    CHECK_THROWS_AS(f.load_string(make_xml(names, trees, "AdaBoost"), names), ValueError);
    // classification trees return nType/purity, not the response
    CHECK_THROWS_AS(f.load_string(make_xml(names, trees, "Grad", "0"), names), ValueError);
    // input-variable transformations
    std::string x = make_xml(names, trees);
    x.insert(x.find("<Weights"), "<Transformations NTransformations=\"1\"><Transform Name=\"Gauss\"/></Transformations>");
    CHECK_THROWS_AS(f.load_string(x, names), ValueError);
    // a Fisher-cut node
    x = make_xml(names, trees);
    x.replace(x.find("NCoef=\"0\""), 9, "NCoef=\"2\"");
    CHECK_THROWS_AS(f.load_string(x, names), ValueError);
}
