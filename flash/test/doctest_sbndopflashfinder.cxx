// End-to-end tests of SBNDOpFlashFinder on synthetic SBND-like hits: narrow single-photon
// tail hits after a few bright prompt hits, 8 us accumulator bins.  The component is made
// through the WCT factory by type name (as wire-cell does), so the test needs no link-time
// symbols of the class itself.

#include "WireCellUtil/doctest.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/ITensorSetFilter.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>

using namespace WireCell;

namespace {
    std::string geom_file()
    {
        static std::string path;
        if (path.empty()) {
            auto p = std::filesystem::temp_directory_path() / "doctest_sbndopflashfinder_geom.json";
            std::ofstream f(p);
            f << "{\"opdets\":[";
            for (int od = 0; od < 8; ++od) {
                if (od) f << ",";
                f << "{\"opdet\":" << od << ",\"x\":-2130.0,\"y\":" << od * 100.0 << ",\"z\":" << od * 500.0 << "}";
            }
            f << "]}";
            path = p.string();
        }
        return path;
    }

    // One ophit row: channel, time, width, area, amplitude, PE, start, flash id, fast/total.
    using Row = std::array<double, 9>;
    Row hit(int ch, double t, double pe) { return {double(ch), t, 30.0, pe, pe, pe, t - 10.0, -1.0, 0.0}; }

    // Prompt hits 60/50/40/30 PE on ch0-3 (or chans) at t, t+3, t+6, t+8.
    void add_prompt(std::vector<Row>& rows, double t, std::array<int, 4> ch = {0, 1, 2, 3})
    {
        const double dt[4] = {0, 3, 6, 8}, pe[4] = {60, 50, 40, 30};
        for (int i = 0; i < 4; ++i) rows.push_back(hit(ch[i], t + dt[i], pe[i]));
    }
    // n tail hits of 2 PE, one every step ns from t0, on the given channels in turn.
    void add_tail(std::vector<Row>& rows, double t0, int n, double step, std::array<int, 4> ch = {0, 1, 2, 3})
    {
        for (int i = 0; i < n; ++i) rows.push_back(hit(ch[i % 4], t0 + step * i, 2.0));
    }
    // Argon slow light: n hits of 2 PE from t0 on, spaced as an exponential with a 1.6 us
    // time constant (deterministic quantiles), on the given channels in turn.
    void add_exp_tail(std::vector<Row>& rows, double t0, int n, std::array<int, 4> ch = {0, 1, 2, 3})
    {
        for (int i = 0; i < n; ++i)
            rows.push_back(hit(ch[i % 4], t0 - 1600.0 * std::log(1.0 - (i + 0.5) / n), 2.0));
    }

    struct Out {
        std::vector<double> time, pe;
        std::vector<int> flash_of;  // per input row
    };
    // A fresh instance per call (factory instances are cached by type:name).
    std::string fresh_name()
    {
        PluginManager::instance().add("WireCellFlash");  // registers the flash factories
        static int n = 0;
        return "SBNDOpFlashFinder:doctest" + std::to_string(n++);
    }

    Out run(const Configuration& over, const std::vector<Row>& rows)
    {
        const auto tn = fresh_name();
        auto ff = Factory::lookup_tn<ITensorSetFilter>(tn);
        auto cf = Factory::find_tn<IConfigurable>(tn);
        auto cfg = cf->default_configuration();
        cfg["nchan"] = 8;
        cfg["geom_file"] = geom_file();
        for (const auto& key : over.getMemberNames()) cfg[key] = over[key];
        cf->configure(cfg);
        std::vector<double> flat;
        for (const auto& r : rows) flat.insert(flat.end(), r.begin(), r.end());
        Configuration md;
        md["name"] = "ophits";
        auto* tv = new ITensor::vector;
        tv->push_back(std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{rows.size(), (size_t) 9}, flat.data(), md));
        auto in = std::make_shared<Aux::SimpleTensorSet>(0, Configuration{}, ITensor::shared_vector(tv));
        ITensorSet::pointer out;
        REQUIRE((*ff)(in, out));
        REQUIRE(out);
        Out o;
        for (const auto& ten : *out->tensors()) {
            const auto name = ten->metadata()["name"].asString();
            const double* d = (const double*) ten->data();
            if (name == "opflash") {
                const size_t mcol = ten->shape()[1];
                for (size_t f = 0; f < ten->shape()[0]; ++f) {
                    o.time.push_back(d[f * mcol]);
                    double s = 0;
                    for (size_t c = 1; c < mcol; ++c) s += d[f * mcol + c];
                    o.pe.push_back(s);
                }
            }
            if (name == "ophits")
                for (size_t r = 0; r < ten->shape()[0]; ++r) o.flash_of.push_back(int(d[r * 9 + 7]));
        }
        return o;
    }
    Configuration only(std::initializer_list<std::pair<const char*, Json::Value>> kv)
    {
        Configuration c;
        for (const auto& [k, v] : kv) c[k] = v;
        return c;
    }
}

TEST_CASE("sbndopflashfinder prompt time and full PE of one flash")
{
    std::vector<Row> rows;
    add_prompt(rows, 1000.0);
    add_tail(rows, 1100.0, 100, 40.0);  // 200 PE of 2 PE hits, 1100-5060
    auto o = run(Configuration{}, rows);
    REQUIRE(o.time.size() == 1);
    CHECK(o.pe[0] == doctest::Approx(380.0));
    CHECK(o.time[0] == doctest::Approx(1001.5));  // mean(1000, 1003): 60+50 = 61 % of 180 PE
}

TEST_CASE("sbndopflashfinder prompt time falls back to the brightest bin")
{
    std::vector<Row> rows = {hit(0, 1000, 1), hit(1, 1001, 1), hit(2, 1002, 1), hit(3, 1025, 1)};
    for (int i = 0; i < 20; ++i) rows.push_back(hit(i % 4, 1030 + 5 * i, 1));
    auto o = run(Configuration{}, rows);
    REQUIRE(o.time.size() == 1);
    CHECK(o.time[0] == doctest::Approx(1000.0));
}

TEST_CASE("sbndopflashfinder pulse split separates a second pulse in the same bin")
{
    std::vector<Row> rows;
    add_prompt(rows, 1000.0);
    add_tail(rows, 1100.0, 100, 40.0);
    add_prompt(rows, 5000.0);
    auto o = run(Configuration{}, rows);
    REQUIRE(o.time.size() == 2);
    CHECK(o.time[0] == doctest::Approx(1001.5));
    CHECK(o.pe[0] == doctest::Approx(376.0));  // 180 + tail hits before 5000
    CHECK(o.time[1] == doctest::Approx(5001.5));
    CHECK(o.pe[1] == doctest::Approx(184.0));  // 180 + tail hits at 5020, 5060
    auto off = run(only({{"pulse_split", false}}), rows);
    REQUIRE(off.time.size() == 1);
}

TEST_CASE("sbndopflashfinder pulse split handles three pulses")
{
    std::vector<Row> rows;
    add_prompt(rows, 1000.0);
    add_tail(rows, 1100.0, 100, 40.0);
    add_prompt(rows, 3000.0);
    add_prompt(rows, 6000.0);
    auto o = run(Configuration{}, rows);
    REQUIRE(o.time.size() == 3);
    CHECK(o.time[1] == doctest::Approx(3001.5));
    CHECK(o.time[2] == doctest::Approx(6001.5));
    CHECK(o.pe[0] + o.pe[1] + o.pe[2] == doctest::Approx(740.0));
}

TEST_CASE("sbndopflashfinder pulse split ignores a one-PMT burst and a rising tail")
{
    std::vector<Row> rows;
    add_prompt(rows, 1000.0);
    add_tail(rows, 1100.0, 100, 40.0);
    auto one_pmt = rows;
    one_pmt.push_back(hit(0, 5600.0, 1000.0));  // the reco1 late artefact hit on one PMT
    CHECK(run(Configuration{}, one_pmt).time.size() == 1);
    // with the check off it splits off (min_fired_pds 0: the part lights one PMT only)
    CHECK(run(only({{"split_drop_brightest", false}, {"min_fired_pds", 0}}), one_pmt).time.size() == 2);
    auto rising = rows;  // 180 PE in the 500 ns before the burst: no dip
    for (int i = 0; i < 60; ++i) rising.push_back(hit(i % 4, 4500.0 + 8.0 * i, 3.0));
    add_prompt(rising, 5000.0);
    CHECK(run(Configuration{}, rising).time.size() == 1);
}

TEST_CASE("sbndopflashfinder join puts a tail piece back and late light no longer drops it")
{
    // A 1 PE hit at 0 anchors the bins; a flash at 7600 with a 1.6 us slow tail (600 PE): the
    // shifted bin [4000, 12000) claims it, and its tail after 12 us (~40 PE) is claimed by
    // [8000, 16000) as a separate piece, which must end up in the same flash.
    std::vector<Row> rows = {hit(4, 0.0, 1.0)};
    add_prompt(rows, 7600.0);
    add_exp_tail(rows, 7700.0, 300);
    auto o = run(Configuration{}, rows);
    REQUIRE(o.time.size() == 1);
    CHECK(o.time[0] == doctest::Approx(7601.5));
    CHECK(o.pe[0] > 780.0 - 10.0);  // all but a few stragglers after 16 us (a bin under 20 PE)
    CHECK(run(only({{"join_glow", "only"}}), rows).time.size() == 1);  // a glow piece: joined too
    auto nojoin = run(only({{"join", false}}), rows);  // the piece is then deleted as late light
    REQUIRE(nojoin.time.size() == 1);
    CHECK(nojoin.pe[0] < 780.0 - 20.0);
    auto neither = run(only({{"join", false}, {"remove_late_light", false}}), rows);
    CHECK(neither.time.size() == 2);
    // "never": the glow piece is not joined; step 5 then deletes it
    auto never = run(only({{"join_glow", "never"}}), rows);
    REQUIRE(never.time.size() == 1);
    CHECK(never.pe[0] < 780.0 - 20.0);
}

TEST_CASE("sbndopflashfinder join_glow decides a touching piece the glow cannot explain")
{
    // 60 PE of new light (no sharp burst) 10.5 us after a flash, touching its bin: the flash's
    // slow tail explains ~0 PE there.  "any" (default) joins it; "only" keeps it separate.
    std::vector<Row> rows = {hit(4, 0.0, 1.0)};
    add_prompt(rows, 7600.0);
    add_exp_tail(rows, 7700.0, 300);
    for (int i = 0; i < 30; ++i) rows.push_back(hit(i % 4, 18100.0 + 60.0 * i, 2.0));
    CHECK(run(Configuration{}, rows).time.size() == 1);
    CHECK(run(only({{"join_glow", "only"}}), rows).time.size() == 2);
}

TEST_CASE("sbndopflashfinder join does not re-join a split-off pulse")
{
    // Second pulse on a subset of the first's PMTs with less than half its PE: the join
    // conditions hold, but the two are parts of one split.
    std::vector<Row> rows;
    add_prompt(rows, 100.0);
    add_tail(rows, 1100.0, 150, 40.0);
    add_prompt(rows, 4100.0, {1, 2, 3, 0});
    auto o = run(Configuration{}, rows);
    CHECK(o.time.size() == 2);
}

TEST_CASE("sbndopflashfinder second flash near a bin edge comes out whole")
{
    // Flash A at 100, flash B at 7600 (0.5 us before the bin edge at 8100) on other PMTs, both
    // with a 1.6 us slow tail: every one of B's hits must be in B's flash.
    std::vector<Row> rows;
    add_prompt(rows, 100.0);
    add_exp_tail(rows, 200.0, 450);
    const size_t b0 = rows.size();
    add_prompt(rows, 7600.0, {4, 5, 6, 7});
    add_exp_tail(rows, 7700.0, 450, {4, 5, 6, 7});
    auto o = run(Configuration{}, rows);
    REQUIRE(o.time.size() == 2);
    CHECK(o.time[0] == doctest::Approx(101.5));
    CHECK(o.time[1] == doctest::Approx(7601.5));
    size_t in_b = 0;
    for (size_t r = b0; r < rows.size(); ++r) in_b += (o.flash_of[r] == 1);
    CHECK(in_b + 5 >= rows.size() - b0);  // all but a few stragglers in a last bin under 20 PE
}

TEST_CASE("sbndopflashfinder join keeps a new flash that starts after a bin edge")
{
    // A small real flash (burst) just after A's bin ends, on A's PMTs, dimmer: not a tail piece.
    std::vector<Row> rows;
    add_prompt(rows, 100.0);
    add_tail(rows, 1100.0, 175, 40.0);  // to 8060
    add_prompt(rows, 8300.0, {1, 2, 3, 0});
    auto o = run(Configuration{}, rows);
    REQUIRE(o.time.size() == 2);
    CHECK(o.time[1] == doctest::Approx(8301.5));
}

TEST_CASE("sbndopflashfinder join does not attach a bright flash to a dim earlier one")
{
    // A dim flash (180 PE) whose tail runs up to a much brighter flash 12 us later on the same
    // PMTs (as in data e37): the bright one has its own burst and must stay separate.
    std::vector<Row> rows;
    add_prompt(rows, 100.0);
    add_tail(rows, 1100.0, 280, 40.0);  // to 12260
    for (int i = 0; i < 4; ++i) rows.push_back(hit(i, 12300.0 + 2.0 * i, 3000.0));
    add_exp_tail(rows, 12400.0, 3000);
    auto o = run(Configuration{}, rows);
    REQUIRE(o.time.size() == 2);
    CHECK(o.time[0] == doctest::Approx(101.5));
    CHECK(o.time[1] == doctest::Approx(12303.0).epsilon(0.001));
    CHECK(o.pe[0] < 1000.0);
}

TEST_CASE("sbndopflashfinder late light keeps a bright flash right after a split-off one")
{
    // A dim pulse, then 0.9 us later a 3x brighter one with a long tail on the same PMTs: the
    // split cuts them; the second must not be deleted as the "late light" of the short first part.
    std::vector<Row> rows;
    add_prompt(rows, 100.0);
    for (int i = 0; i < 4; ++i) rows.push_back(hit(i, 1000.0 + 2.0 * i, 150.0));
    add_exp_tail(rows, 1100.0, 600);
    auto o = run(Configuration{}, rows);
    REQUIRE(o.time.size() == 2);
    CHECK(o.time[1] == doctest::Approx(1002.0).epsilon(0.01));
    CHECK(o.pe[1] > 1500.0);
}

TEST_CASE("sbndopflashfinder quality cut and bad config")
{
    std::vector<Row> rows = {hit(0, 100, 50)};  // one PMT only
    CHECK(run(Configuration{}, rows).time.empty());
    CHECK(run(only({{"min_fired_pds", 0}}), rows).time.size() == 1);
    auto cf = Factory::lookup_tn<IConfigurable>(fresh_name());
    auto cfg = cf->default_configuration();
    cfg["bin_width"] = 0;
    CHECK_THROWS(cf->configure(cfg));
    cfg = cf->default_configuration();
    cfg["join_glow"] = "sometimes";
    CHECK_THROWS(cf->configure(cfg));
}
