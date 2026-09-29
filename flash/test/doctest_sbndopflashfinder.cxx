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

    // Light travel-time inputs for the 8 test OpDets: box 0 = coated 0-2 + uncoated 3, box 1 =
    // coated 4-6 + uncoated 7; curve ratio 0.1 -> 190 cm ... 1.0 -> 10 cm (linear); SBND D, v.
    std::string lt_file()
    {
        static std::string path;
        if (path.empty()) {
            auto p = std::filesystem::temp_directory_path() / "doctest_sbndopflashfinder_lt.json";
            std::ofstream f(p);
            f << "{\"drift_cm\":201.3,\"v_vuv_cm_per_ns\":13.5,\"v_vis_cm_per_ns\":23.99,\"opdets\":[";
            for (int od = 0; od < 8; ++od) {
                if (od) f << ",";
                f << "{\"opdet\":" << od << ",\"type\":" << (od % 4 == 3 ? 2 : 1) << ",\"box\":" << od / 4 << "}";
            }
            f << "],\"curves\":{\"test\":{\"ratio\":[0.1,1.0],\"x_cm\":[190.0,10.0]}}}";
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
        std::vector<std::array<double, 4>> drift;  // "flash_drift" rows, if written
        bool have_drift{false};
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
            if (name == "flash_drift") {
                o.have_drift = true;
                for (size_t r = 0; r < ten->shape()[0]; ++r)
                    o.drift.push_back({d[r * 4], d[r * 4 + 1], d[r * 4 + 2], d[r * 4 + 3]});
            }
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

TEST_CASE("sbndopflashfinder a single-PMT fake hit cannot set the prompt time")
{
    // as the r713_s52_e6 flash: a giant one-hit fake on one PMT well after the real prompt light
    std::vector<Row> rows;
    add_prompt(rows, 1000.0);
    add_tail(rows, 1100.0, 100, 40.0);
    rows.push_back(hit(0, 4000.0, 1000.0));
    auto o = run(Configuration{}, rows);
    REQUIRE(o.time.size() == 1);
    CHECK(o.pe[0] == doctest::Approx(1380.0));  // the fake hit is still in the flash PE
    CHECK(o.time[0] == doctest::Approx(1001.5));
    // without the coincidence rule the brightest bin (the fake) sets the time
    o = run(only({{"prompt_min_hits", 1}, {"prompt_min_pe", 0.0}}), rows);
    REQUIRE(o.time.size() == 1);
    CHECK(o.time[0] == doctest::Approx(4000.0));
    // >= 5 OpDets: no bin qualifies -> brightest bin of all
    o = run(only({{"prompt_min_opdets", 5}}), rows);
    REQUIRE(o.time.size() == 1);
    CHECK(o.time[0] == doctest::Approx(4000.0));
    // >= 4 OpDets: the 4-PMT prompt bin still qualifies
    o = run(only({{"prompt_min_opdets", 4}}), rows);
    REQUIRE(o.time.size() == 1);
    CHECK(o.time[0] == doctest::Approx(1001.5));
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

namespace {
    // One flash: prompt hits on coated 0, 1, 2 (60/50/40 PE) and uncoated 3 (upe PE, none if 0);
    // prompt time 1001.5 ns as in the first test.
    std::vector<Row> lt_rows(double upe)
    {
        std::vector<Row> rows = {hit(0, 1000, 60), hit(1, 1003, 50), hit(2, 1006, 40)};
        rows.push_back(upe > 0 ? hit(3, 1008, upe) : hit(4, 1008, 30));
        return rows;
    }
    Configuration lt_cfg(int tpc, const char* fail = "wires")
    {
        return only({{"light_travel", true}, {"light_travel_file", lt_file()}, {"light_travel_curve", "test"},
                     {"light_travel_fail", fail}, {"tpc", tpc}});
    }
}

TEST_CASE("sbndopflashfinder light travel off: time and tensors unchanged")
{
    auto o = run(Configuration{}, lt_rows(30));
    REQUIRE(o.time.size() == 1);
    CHECK(o.time[0] == doctest::Approx(1001.5));
    CHECK_FALSE(o.have_drift);
}

TEST_CASE("sbndopflashfinder light travel: X from the PMT ratio, VUV branch")
{
    // ratio = 30 / 150 x 3 lit coated = 0.6 -> |X| = 190 - 180 x 0.5 / 0.9 = 90 cm > kink 44.0 cm
    // -> travel (201.3 - 90) / 13.5 = 8.2444 ns
    auto o = run(lt_cfg(0), lt_rows(30));
    REQUIRE(o.time.size() == 1);
    REQUIRE(o.drift.size() == 1);
    CHECK(o.drift[0][0] == doctest::Approx(-90.0));  // TPC 0: negative
    CHECK(o.drift[0][1] == doctest::Approx(111.3 / 13.5));
    CHECK(o.drift[0][2] == doctest::Approx(1001.5));
    CHECK(o.drift[0][3] == 1.0);
    CHECK(o.time[0] == doctest::Approx(1001.5 - 111.3 / 13.5));
    CHECK(run(lt_cfg(1), lt_rows(30)).drift[0][0] == doctest::Approx(90.0));
}

TEST_CASE("sbndopflashfinder light travel: near the cathode, visible-light branch")
{
    // ratio = 50 / 150 x 3 = 1.0 -> |X| = 10 cm < kink -> 10 / 13.5 + 201.3 / 23.99 ns
    auto o = run(lt_cfg(1), lt_rows(50));
    REQUIRE(o.drift.size() == 1);
    CHECK(o.drift[0][0] == doctest::Approx(10.0));
    CHECK(o.drift[0][1] == doctest::Approx(10.0 / 13.5 + 201.3 / 23.99));
    CHECK(o.drift[0][2] == doctest::Approx((1000.0 + 1003.0 + 1008.0) / 3));  // 60+50+50 PE > 60 %
    CHECK(o.time[0] == doctest::Approx(o.drift[0][2] - 10.0 / 13.5 - 201.3 / 23.99));
}

TEST_CASE("sbndopflashfinder light travel: X fails without uncoated light")
{
    // wires: the ratio -> 0 end of the curve, 190 cm -> 11.3 / 13.5 ns; X written 0, ok 0
    auto o = run(lt_cfg(0), lt_rows(0));
    REQUIRE(o.drift.size() == 1);
    CHECK(o.drift[0][0] == 0.0);
    CHECK(o.drift[0][3] == 0.0);
    CHECK(o.drift[0][1] == doctest::Approx(11.3 / 13.5));
    CHECK(o.time[0] == doctest::Approx(1001.5 - 11.3 / 13.5));
    // none: no correction
    auto n = run(lt_cfg(0, "none"), lt_rows(0));
    CHECK(n.drift[0][1] == 0.0);
    CHECK(n.time[0] == doctest::Approx(1001.5));
}

TEST_CASE("sbndopflashfinder light travel bad config")
{
    auto cf = Factory::lookup_tn<IConfigurable>(fresh_name());
    auto cfg = cf->default_configuration();
    cfg["light_travel"] = true;
    CHECK_THROWS(cf->configure(cfg));  // no file
    cfg["light_travel_file"] = lt_file();
    cfg["light_travel_curve"] = "test";
    cfg["light_travel_fail"] = "sometimes";
    CHECK_THROWS(cf->configure(cfg));
    cfg["light_travel_fail"] = "wires";
    cfg["light_travel_curve"] = "nope";
    CHECK_THROWS(cf->configure(cfg));
}

TEST_CASE("sbndopflashfinder fake_veto drops single-PMT fake hits after a bright pulse")
{
    // bright real pulse on ch0 (8000 PE, one 4.4 us wide hit) with the other PMTs, slow tail,
    // then on ch0 alone a 3000 PE fake at 4000 ns (the tail gives the other PMTs 8 PE around it)
    std::vector<Row> rows;
    add_prompt(rows, 1000.0);
    Row wide = hit(0, 1001.0, 8000.0);
    wide[2] = 4400.0;
    const size_t iwide = rows.size();
    rows.push_back(wide);
    add_tail(rows, 1100.0, 100, 40.0);
    rows.push_back(hit(0, 4000.0, 3000.0));
    const size_t ifake = rows.size() - 1;
    auto narrow = rows;  // the same, but the bright hit only 2.5 us wide: a split saturated pulse
    narrow[iwide][2] = 2500.0;
    auto o = run(Configuration{}, rows);  // off by default
    REQUIRE(o.time.size() == 1);
    CHECK(o.pe[0] == doctest::Approx(11380.0));
    CHECK(o.flash_of[ifake] == 0);

    o = run(only({{"fake_veto", true}}), rows);
    REQUIRE(o.time.size() == 1);
    CHECK(o.pe[0] == doctest::Approx(8380.0));
    CHECK(o.flash_of[ifake] == -1);
    CHECK(o.flash_of[0] == 0);        // the bright pulse itself is kept
    CHECK(o.time[0] == doctest::Approx(1001.0));

    // the same late hit with light on the other PMTs (3 x 150 PE) is real and kept
    auto seen = rows;
    for (int ch = 1; ch < 4; ++ch) seen.push_back(hit(ch, 4010.0, 150.0));
    o = run(only({{"fake_veto", true}}), seen);
    REQUIRE(o.time.size() >= 1);
    CHECK(o.flash_of[ifake] >= 0);
    // after a 2.5 us bright hit the late hit is the second piece of a saturated pulse: kept...
    o = run(only({{"fake_veto", true}}), narrow);
    CHECK(o.flash_of[ifake] == 0);
    // ...unless it is too bright for one PMT alone
    o = run(only({{"fake_veto", true}, {"fake_big_pe", 2000.0}}), narrow);
    CHECK(o.flash_of[ifake] == -1);
    // nor is a lone bright hit with no bright hit before it on its PMT
    std::vector<Row> lone;
    add_prompt(lone, 1000.0);
    add_tail(lone, 1100.0, 100, 40.0);
    lone.push_back(hit(0, 4000.0, 3000.0));
    o = run(only({{"fake_veto", true}}), lone);
    CHECK(o.flash_of.back() >= 0);
}
