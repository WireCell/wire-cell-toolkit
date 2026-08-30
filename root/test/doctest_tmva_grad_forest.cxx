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
#include <limits>
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

// ===========================================================================
// doc sbnd_xin/docs/85 sec 7.2 -- the log-odds transform, unclamped.
//
// The bug these cases pin: the two BDT scorers used to clamp the forest
// output to +-0.9999 before log10((1+v)/(1-v)), capping |score| at 4.30103.
// MicroBooNE's nue working point is nue_score > 7.0, so that cut selected
// zero events by construction, and every strongly-signal-like event piled up
// on one value.  Written against the clamped code these cases FAIL: the
// reachability case returns 4.30103 instead of ~8.3, and the two ceilings
// come back as 4.3009362 (float) / 4.3010083 (double) instead of 16.25562.
// ===========================================================================
TEST_CASE("bdt log-odds: the MicroBooNE nue working point is reachable")
{
    using WireCell::Root::bdt_log_odds_score;
    auto quiet = Log::logger("test.bdt_log_odds");

    // v = 0.99999998 is what a uB nue_score of 7.0 corresponds to; a forest
    // that separates this well must be able to say so.
    const double v7 = (std::pow(10.0, 7.0) - 1.0) / (std::pow(10.0, 7.0) + 1.0);
    CHECK(bdt_log_odds_score(v7, quiet, "test") == doctest::Approx(7.0).epsilon(1e-9));
    CHECK(bdt_log_odds_score(v7, quiet, "test") > 7.0 - 1e-6);

    // The removed clamp's ceiling is now an ordinary interior value, not a
    // wall: scores above it exist and stay ordered.
    const double s_clamp = std::log10((1.0 + 0.9999) / (1.0 - 0.9999));
    CHECK(s_clamp == doctest::Approx(4.30103).epsilon(1e-5));
    CHECK(bdt_log_odds_score(0.99999, quiet, "test") > s_clamp);
    CHECK(bdt_log_odds_score(0.999999, quiet, "test") >
          bdt_log_odds_score(0.99999, quiet, "test"));

    // Sign symmetry and the neutral point, so a sign slip cannot pass.
    CHECK(bdt_log_odds_score(0.0, quiet, "test") == doctest::Approx(0.0));
    CHECK(bdt_log_odds_score(-0.5, quiet, "test") ==
          doctest::Approx(-bdt_log_odds_score(0.5, quiet, "test")));
}

TEST_CASE("bdt log-odds: degenerate forest outputs stay finite")
{
    using WireCell::Root::bdt_log_odds_score;
    auto quiet = Log::logger("test.bdt_log_odds");

    // tanh(sum) rounds to exactly 1.0 in double once sum > ~18.4.  Unguarded
    // that is log10(2/0) = +inf, which would reach a ROOT float branch.
    const double hi = bdt_log_odds_score(1.0, quiet, "test");
    CHECK(std::isfinite(hi));
    CHECK(hi == doctest::Approx(16.25562).epsilon(1e-6));
    CHECK(static_cast<float>(hi) < std::numeric_limits<float>::max());

    // TmvaGradForest::evaluate returns -999 when any input variable is NaN
    // (mirroring TMVA::Reader::EvaluateMVA).  log10 of a negative ratio is
    // NaN; the guard must send it to the floor, the same direction the old
    // clamp sent it (-4.30103), not to NaN.
    const double lo = bdt_log_odds_score(-999.0, quiet, "test");
    CHECK(std::isfinite(lo));
    CHECK(lo == doctest::Approx(-16.25562).epsilon(1e-6));
    CHECK(lo < -15.0);   // below the "not filled" sentinel, as background-like

    CHECK(std::isfinite(bdt_log_odds_score(-1.0, quiet, "test")));
    CHECK(std::isfinite(bdt_log_odds_score(
        std::numeric_limits<double>::quiet_NaN(), quiet, "test")));
    CHECK(std::isfinite(bdt_log_odds_score(
        std::numeric_limits<double>::infinity(), quiet, "test")));

    // The guard is degenerate-only: it must not move an ordinary value.
    CHECK(bdt_log_odds_score(0.9999, quiet, "test") ==
          doctest::Approx(std::log10((1.0 + 0.9999) / (1.0 - 0.9999))));
}
