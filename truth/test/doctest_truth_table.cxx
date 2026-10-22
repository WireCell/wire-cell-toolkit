/** TruthTable round trip: a truth file written with the python layout
 *  (tables + uint8 string matrices, stored zip) loads with the schema
 *  columns, serves rows by name, remaps foreign process codes, and
 *  TruthTableSource emits the two tensors the SBND reader expects.
 *  A file written by truth/python/truth_table.py --selftest is loaded
 *  too when WCT_TRUTH_TEST_FILE names it.
 */

#include "WireCellTruth/TruthTable.h"
#include "WireCellTruth/TruthTableSource.h"
#include "WireCellTruth/Schema.h"

#include "WireCellAux/SimpleDepo.h"
#include "WireCellUtil/cnpy.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <boost/filesystem.hpp>
#include <cstdlib>
#include <cstring>

using namespace WireCell;

namespace {
    std::string tmpname(const std::string& base)
    {
        auto dir = boost::filesystem::temp_directory_path();
        const char* home = std::getenv("HOME");
        if (home) {
            boost::filesystem::path p(home);
            p /= "tmp";
            if (boost::filesystem::is_directory(p)) dir = p;
        }
        return (dir / (base + "-" + std::to_string(::getpid()) + ".npz")).string();
    }

    void save_strings(const std::string& fname, const std::string& name, const std::vector<std::string>& strs, const std::string& mode)
    {
        size_t L = 1;
        for (const auto& s : strs) L = std::max(L, s.size());
        std::vector<uint8_t> buf(strs.size() * L, 0);
        for (size_t i = 0; i < strs.size(); ++i) std::memcpy(&buf[i * L], strs[i].data(), strs[i].size());
        cnpy::npz_save<uint8_t>(fname, name, buf.data(), {strs.size(), L}, mode);
    }

    // Write a file with the python layout; process codes from a
    // *different* process list to exercise the remap.
    std::string write_file()
    {
        const auto fname = tmpname("doctest_truth_table");
        const auto& nuc = Truth::truth_nu_columns();
        const auto& pfc = Truth::truth_pf_columns();
        std::vector<double> nu(nuc.size(), 0.0);
        nu[Truth::column_index(nuc, "nu_idx")] = 0;
        nu[Truth::column_index(nuc, "pdg")] = 14;
        nu[Truth::column_index(nuc, "E")] = 1.5;
        cnpy::npz_save<double>(fname, Truth::truth_file_nu, nu.data(), {1, nuc.size()}, "w");

        std::vector<std::string> foreign = {"unknown", "Decay", "primary"};
        std::vector<double> pf(2 * pfc.size(), 0.0);
        auto set = [&](size_t r, const char* c, double v) { pf[r * pfc.size() + Truth::column_index(pfc, c)] = v; };
        set(0, "trackid", 1); set(0, "pdg", 13); set(0, "process", 2); set(0, "KE", 0.5); set(0, "nvis_q", 1e5);
        set(1, "trackid", 2); set(1, "pdg", 11); set(1, "process", 1); set(1, "mother_trackid", 1); set(1, "parent_trackid", 1);
        cnpy::npz_save<double>(fname, Truth::truth_file_pf, pf.data(), {2, pfc.size()}, "a");

        std::vector<int32_t> sed = {1, 1, 2, 7, 2, 2};
        cnpy::npz_save<int32_t>(fname, Truth::truth_file_sed_ids, sed.data(), {3, 2}, "a");
        std::vector<double> tq = {1, 1e5, 2, 3e3};
        cnpy::npz_save<double>(fname, Truth::truth_file_sed_trackq, tq.data(), {2, 2}, "a");

        save_strings(fname, "columns_nu", nuc, "a");
        save_strings(fname, "columns_pf", pfc, "a");
        save_strings(fname, "process_names", foreign, "a");
        save_strings(fname, "gen_labels", {"generator", "corsika"}, "a");
        save_strings(fname, "meta", {"{\"detector\": \"doctest\", \"id_mode\": \"index\"}"}, "a");
        return fname;
    }

    void check_table(const Truth::TruthTable& tt)
    {
        REQUIRE(tt.nu.size() == 1);
        REQUIRE(tt.pf.size() == 2);
        CHECK(tt.nu_value(0, "pdg") == 14);
        CHECK(tt.nu_value(0, "E") == doctest::Approx(1.5));
        CHECK(tt.pf_row(1) == 0);
        CHECK(tt.pf_row(2) == 1);
        CHECK(tt.pf_row(3) == -1);
        CHECK(tt.pf_value(1, "mother_trackid") == 1);
        CHECK(tt.pf_value(0, "pdg") == 13);
        CHECK(tt.sed_ids.size() == 3);
        CHECK(tt.sed_id(1) == std::make_pair(2, 7));
        CHECK(tt.sed_id(99) == std::make_pair(99, 99));
        CHECK(tt.sed_trackq.at(2) == doctest::Approx(3e3));
        CHECK(tt.source == "file");
        CHECK_THROWS_AS(tt.pf_value(0, "no_such_column"), KeyError);
    }
}

TEST_CASE("truth table loads the python layout and remaps process codes")
{
    const auto fname = write_file();
    auto tt = Truth::TruthTable::load(fname);
    check_table(tt);
    // the file coded "primary" as 2 and "Decay" as 1 in its own list
    CHECK(tt.pf_value(0, "process") == Truth::process_code("primary"));
    CHECK(tt.pf_value(1, "process") == Truth::process_code("Decay"));
    CHECK(tt.gen_labels == std::vector<std::string>{"generator", "corsika"});
    CHECK(tt.meta["detector"].asString() == "doctest");

    auto tens = tt.as_tensors(5, "truth/5");
    REQUIRE(tens.size() == 2);
    CHECK(tens[0]->metadata()["datatype"].asString() == "truth_nu");
    CHECK(tens[0]->metadata()["datapath"].asString() == "truth/5/nu");
    CHECK(tens[1]->metadata()["datatype"].asString() == "truth_pf");
    CHECK(tens[1]->shape() == ITensor::shape_t{2, Truth::truth_pf_columns().size()});
    CHECK(tens[1]->element_type() == typeid(double));
    CHECK(tens[1]->metadata()["columns"].size() == Truth::truth_pf_columns().size());
    const double* d = reinterpret_cast<const double*>(tens[1]->data());
    CHECK(d[1 * Truth::truth_pf_columns().size() + Truth::column_index(Truth::truth_pf_columns(), "trackid")] == 2);

    CHECK_THROWS_AS(Truth::TruthTable::load(fname + ".missing"), IOError);
    boost::filesystem::remove(fname);
}

TEST_CASE("truth table source serves one tensor set per ident then EOS")
{
    const auto fname = write_file();
    Truth::TruthTableSource src;
    auto cfg = src.default_configuration();
    cfg["inname"] = fname;      // no %d: the same file for every ident
    cfg["idents"] = Json::arrayValue;
    cfg["idents"].append(3);
    cfg["idents"].append(4);
    src.configure(cfg);

    ITensorSet::pointer ts;
    REQUIRE(src(ts));
    REQUIRE(ts);
    CHECK(ts->ident() == 3);
    CHECK(ts->tensors()->size() == 2);
    CHECK(ts->tensors()->at(0)->metadata()["datapath"].asString() == "truth/3/nu");
    REQUIRE(src(ts));
    REQUIRE(ts);
    CHECK(ts->ident() == 4);
    REQUIRE(src(ts));
    CHECK(ts == nullptr);       // EOS
    CHECK_FALSE(src(ts));       // nothing after EOS

    Truth::TruthTableSource bad;
    auto bcfg = bad.default_configuration();
    bcfg["inname"] = fname + ".missing_%d";
    bad.configure(bcfg);
    CHECK_THROWS_AS(bad(ts), IOError);

    Truth::TruthTableSource stop;
    auto scfg = stop.default_configuration();
    scfg["inname"] = fname + ".missing_%d";
    scfg["stop_on_missing"] = true;
    stop.configure(scfg);
    REQUIRE(stop(ts));
    CHECK(ts == nullptr);
    boost::filesystem::remove(fname);
}

TEST_CASE("truth table from depos when no file exists")
{
    IDepo::vector depos;
    auto prior1 = std::make_shared<Aux::SimpleDepo>(10 * units::us, Point(1, 2, 3), -100, nullptr, 0, 0, 7, 13, 2.0);
    auto prior1b = std::make_shared<Aux::SimpleDepo>(5 * units::us, Point(0, 2, 3), -50, nullptr, 0, 0, 7, 13, 1.0);
    auto prior2 = std::make_shared<Aux::SimpleDepo>(20 * units::us, Point(4, 5, 6), -30, nullptr, 0, 0, 9, 2212, 3.0);
    depos.push_back(std::make_shared<Aux::SimpleDepo>(1 * units::ms, Point(0, 2, 3), -100, prior1, 0, 0, 7, 0, 1.0));
    depos.push_back(std::make_shared<Aux::SimpleDepo>(1 * units::ms, Point(0, 2, 3), -50, prior1b, 0, 0, 7, 0, 1.0));
    depos.push_back(std::make_shared<Aux::SimpleDepo>(1 * units::ms, Point(0, 5, 6), -30, prior2, 0, 0, 9, 0, 1.0));
    auto tt = Truth::TruthTable::from_depos(depos);
    CHECK(tt.source == "depos");
    REQUIRE(tt.nu.size() == 1);
    CHECK(tt.nu_value(0, "nu_idx") == -1);
    REQUIRE(tt.pf.size() == 2);
    CHECK(tt.pf_value(0, "trackid") == 7);
    CHECK(tt.pf_value(0, "pdg") == 13);
    CHECK(tt.pf_value(0, "start_t") == doctest::Approx(5 * units::us));   // the earliest prior
    CHECK(tt.pf_value(0, "start_x") == doctest::Approx(0.0));
    CHECK(tt.pf_value(0, "nvis_q") == doctest::Approx(150));
    CHECK(tt.pf_value(1, "trackid") == 9);
    CHECK(tt.pf_value(1, "pdg") == 2212);
    CHECK(tt.pf_value(1, "process") == Truth::process_unknown);
    CHECK(tt.pf_row(9) == 1);
}

TEST_CASE("truth table loads a file written by truth_table.py when given")
{
    const char* fn = std::getenv("WCT_TRUTH_TEST_FILE");
    if (!fn) {
        MESSAGE("WCT_TRUTH_TEST_FILE not set, skipping the python file load");
        return;
    }
    auto tt = Truth::TruthTable::load(fn);
    check_table(tt);
    CHECK(tt.pf_value(0, "process") == Truth::process_code("primary"));
    CHECK(tt.pf_value(1, "process") == Truth::process_code("Decay"));
    CHECK(tt.gen_labels == std::vector<std::string>{"generator"});
    CHECK(tt.meta["detector"].asString() == "selftest");
}
