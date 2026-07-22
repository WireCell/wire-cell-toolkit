// End-to-end tests of the OpFlashFinder slow-tail merge (flash_tail_merge,
// default OFF) and its independence from the legacy flash_refine satellite
// merge.  Synthetic ophits tensors run through the component; the scenario
// mirrors the PDVD split-flash defect (pdvd doc 23 §7d): a bright fast flash
// plus, ~1.25 us later, a dim flash made of very wide LAr slow-tail hits on
// the SAME OpDets.

#include "WireCellUtil/doctest.h"

#include "WireCellFlash/OpFlashFinder.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"

#include <array>
#include <filesystem>
#include <fstream>

using namespace WireCell;

// 8 OpDets in a z-line (single side, x>0): legacy adjacency = |Δod| <= 1.
static std::string geom_file()
{
    static std::string path;
    if (path.empty()) {
        auto p = std::filesystem::temp_directory_path() / "doctest_opflashfinder_geom.json";
        std::ofstream f(p);
        f << "{\"opdets\":[";
        for (int od = 0; od < 8; ++od) {
            if (od) f << ",";
            f << "{\"opdet\":" << od << ",\"x\":1000.0,\"y\":0.0,\"z\":" << od * 500.0 << "}";
        }
        f << "]}";
        path = p.string();
    }
    return path;
}

// One ophit row: {channel, peak_time, width, area, amplitude, pe, start_time,
// flash_id(-1), fast_to_total}.  Times/widths in WCT ns.
using HitRow = std::array<double, 9>;
static HitRow hit(int ch, double t, double width, double pe)
{
    return {double(ch), t, width, pe, pe, pe, t - 0.5 * width, -1.0, 0.0};
}

// The fast seed: 100 PE of narrow hits on ch0-3 at t=100.
static std::vector<HitRow> seed_hits()
{
    std::vector<HitRow> rows;
    for (int ch = 0; ch < 4; ++ch) rows.push_back(hit(ch, 100.0, 30.0, 25.0));
    return rows;
}

// A slow-tail candidate: wide hits at time t on ch0.. (pe each).
static void add_tail(std::vector<HitRow>& rows, double t, double width,
                     double pe_each, int ch_first = 0, int nch = 4)
{
    for (int ch = ch_first; ch < ch_first + nch; ++ch)
        rows.push_back(hit(ch, t, width, pe_each));
}

static ITensorSet::pointer run_finder(const Configuration& over,
                                      const std::vector<HitRow>& rows)
{
    Flash::OpFlashFinder ff;
    auto cfg = ff.default_configuration();
    cfg["nchan"] = 8;
    cfg["geom_file"] = geom_file();
    cfg["remove_late_light"] = false;  // isolate the merge logic
    for (const auto& key : over.getMemberNames()) cfg[key] = over[key];
    ff.configure(cfg);

    std::vector<double> flat;
    for (const auto& r : rows) flat.insert(flat.end(), r.begin(), r.end());
    Configuration md;
    md["name"] = "ophits";
    ITensor::vector* tensors = new ITensor::vector;
    tensors->push_back(std::make_shared<Aux::SimpleTensor>(
        ITensor::shape_t{rows.size(), (size_t) 9}, flat.data(), md));
    auto in = std::make_shared<Aux::SimpleTensorSet>(0, Configuration{},
                                                     ITensor::shared_vector(tensors));
    ITensorSet::pointer out = nullptr;
    REQUIRE(ff(in, out));
    REQUIRE(out != nullptr);
    return out;
}

static const ITensor::pointer named_tensor(const ITensorSet::pointer& ts, const std::string& name)
{
    for (const auto& ten : *ts->tensors()) {
        if (ten->metadata()["name"].asString() == name) return ten;
    }
    return nullptr;
}

TEST_CASE("opflashfinder tail merge off keeps the split pair")
{
    auto rows = seed_hits();
    add_tail(rows, 1350.0, 1600.0, 15.0);  // 60 PE, dt 1250 ns, on the seed's PDs

    auto out = run_finder(Configuration{}, rows);
    auto op = named_tensor(out, "opflash");
    REQUIRE(op->shape()[0] == 2);
    const double* M = (const double*) op->data();
    const size_t mcol = op->shape()[1];
    CHECK(M[0 * mcol] == doctest::Approx(100.0));   // fast member
    CHECK(M[1 * mcol] == doctest::Approx(1350.0));  // split-off tail
}

TEST_CASE("opflashfinder tail merge absorbs the slow tail, keeping the seed time")
{
    auto rows = seed_hits();
    add_tail(rows, 1350.0, 1600.0, 15.0);

    Configuration over;
    over["flash_tail_merge"] = true;
    auto out = run_finder(over, rows);

    auto op = named_tensor(out, "opflash");
    REQUIRE(op->shape()[0] == 1);
    const double* M = (const double*) op->data();
    const size_t mcol = op->shape()[1];
    CHECK(M[0] == doctest::Approx(100.0));  // seed (fast-peak) time kept, NOT PE-weighted
    double tot = 0;
    for (size_t c = 1; c < mcol; ++c) tot += M[c];
    CHECK(tot == doctest::Approx(160.0));

    const double* S = (const double*) named_tensor(out, "flash_summary")->data();
    CHECK(S[1] == doctest::Approx(160.0));  // total_pe
    CHECK(S[7] == doctest::Approx(8.0));    // nhits

    const double* H = (const double*) named_tensor(out, "ophits")->data();
    for (size_t r = 0; r < 8; ++r) CHECK(H[r * 9 + 7] == doctest::Approx(0.0));
}

TEST_CASE("opflashfinder tail merge leaves genuine pile-up alone")
{
    Configuration over;
    over["flash_tail_merge"] = true;

    SUBCASE("narrow second onset on the same PDs") {
        auto rows = seed_hits();
        add_tail(rows, 1350.0, 100.0, 15.0);  // narrow: a real second flash
        auto out = run_finder(over, rows);
        CHECK(named_tensor(out, "opflash")->shape()[0] == 2);
    }
    SUBCASE("wide hits but on PDs the seed never lit") {
        auto rows = seed_hits();
        add_tail(rows, 1350.0, 1600.0, 15.0, 4);  // ch4-7: unlit in seed
        auto out = run_finder(over, rows);
        CHECK(named_tensor(out, "opflash")->shape()[0] == 2);
    }
    SUBCASE("later flash brighter than the seed") {
        auto rows = seed_hits();
        add_tail(rows, 1350.0, 1600.0, 40.0);  // 160 PE > seed's 100
        auto out = run_finder(over, rows);
        CHECK(named_tensor(out, "opflash")->shape()[0] == 2);
    }
}

TEST_CASE("opflashfinder tail merge tolerates narrow passenger hits")
{
    // The late self-trigger wall/PMT records attach narrow, tiny-PE hits on
    // unlit PDs to the tail member; the PE-dominance fraction rides over them.
    auto rows = seed_hits();
    add_tail(rows, 1350.0, 1600.0, 15.0);
    rows.push_back(hit(5, 1350.0, 50.0, 1.5));  // 2.4% of the tail's PE

    Configuration over;
    over["flash_tail_merge"] = true;
    auto out = run_finder(over, rows);

    auto op = named_tensor(out, "opflash");
    REQUIRE(op->shape()[0] == 1);
    const double* M = (const double*) op->data();
    CHECK(M[0] == doctest::Approx(100.0));
    CHECK(M[1 + 5] == doctest::Approx(1.5));  // passenger PE lands on its PD
}

TEST_CASE("opflashfinder tail merge reaches past an intervening flash")
{
    // The cascade holds i while j walks: the tail merges into its seed even
    // with an unrelated flash in between (which itself stays unmerged).
    auto rows = seed_hits();
    add_tail(rows, 1200.0, 30.0, 20.0, 4);   // genuine flash, ch4-7, 80 PE
    add_tail(rows, 2300.0, 1600.0, 15.0);    // the seed's slow tail, 60 PE

    Configuration over;
    over["flash_tail_merge"] = true;
    auto out = run_finder(over, rows);

    auto op = named_tensor(out, "opflash");
    REQUIRE(op->shape()[0] == 2);
    const double* M = (const double*) op->data();
    const size_t mcol = op->shape()[1];
    CHECK(M[0 * mcol] == doctest::Approx(100.0));   // merged seed, time kept
    double tot0 = 0;
    for (size_t c = 1; c < mcol; ++c) tot0 += M[0 * mcol + c];
    CHECK(tot0 == doctest::Approx(160.0));
    CHECK(M[1 * mcol] == doctest::Approx(1200.0));  // intervening flash intact
}

TEST_CASE("opflashfinder legacy refine is unchanged and blind to the tail")
{
    Configuration over;
    over["flash_refine"] = true;
    over["refine_subset_merge"] = true;

    SUBCASE("PDHD-style dim subset satellite still merges, PE-weighted time") {
        auto rows = seed_hits();
        add_tail(rows, 1300.0, 30.0, 5.0, 0, 2);  // 10 PE on ch0-1: dim subset
        auto out = run_finder(over, rows);
        auto op = named_tensor(out, "opflash");
        REQUIRE(op->shape()[0] == 1);
        const double* M = (const double*) op->data();
        // Legacy semantics: construct_flash PE-weighted mean, not the seed time.
        CHECK(M[0] == doctest::Approx((100.0 * 100.0 + 10.0 * 1300.0) / 110.0));
    }
    SUBCASE("the tail scenario fails the legacy gate (PE ratio 0.6 > 0.5)") {
        auto rows = seed_hits();
        add_tail(rows, 1350.0, 1600.0, 15.0);
        auto out = run_finder(over, rows);
        CHECK(named_tensor(out, "opflash")->shape()[0] == 2);
    }
    SUBCASE("tail knob catches it on top of legacy refine") {
        auto rows = seed_hits();
        add_tail(rows, 1350.0, 1600.0, 15.0);
        over["flash_tail_merge"] = true;
        auto out = run_finder(over, rows);
        auto op = named_tensor(out, "opflash");
        REQUIRE(op->shape()[0] == 1);
        CHECK(((const double*) op->data())[0] == doctest::Approx(100.0));
    }
}
