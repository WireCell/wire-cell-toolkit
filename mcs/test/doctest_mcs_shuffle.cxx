// doc 80 sec 6.3: the RIGHT determinism test for this code is input-order
// invariance, not N-run repetition -- upstream is order-deterministic but
// order-SENSITIVE (ties in three unstable sorts propagate into segment
// membership).  The port adds id tie-breaks to every sort; this test permutes
// each fixture cloud by 20 fixed mt19937 permutations and requires
// byte-identical outputs.  Run it under `setarch x86_64 -R` as well (M4);
// nothing here reads addresses, so both must pass.

#include "WireCellMcs/MuonMCS.h"

#include "WireCellUtil/Persist.h"
#include "WireCellUtil/doctest.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <random>
#include <string>
#include <vector>

using namespace WireCell::Mcs;
using detail::VD;
using detail::VVD;

namespace {

    std::string data_dir()
    {
        const char* env = std::getenv("MCS_TEST_DATA");
        std::vector<std::string> cands;
        if (env) cands.push_back(env);
        cands.push_back("mcs/test/data");
        cands.push_back("../mcs/test/data");
        for (const auto& c : cands) {
            std::ifstream probe(c + "/mcs_golden_ref.json");
            if (probe.good()) return c;
        }
        FAIL("mcs test fixtures not found; run from the repo root or set MCS_TEST_DATA");
        return "";
    }

    double jd(const Json::Value& v)
    {
        if (v.isString()) return std::numeric_limits<double>::quiet_NaN();
        return v.asDouble();
    }

}  // namespace

TEST_CASE("mcs shuffle invariance")
{
    const std::vector<std::string> fixtures = {
        "mcs_golden_ref.json",
        "mcs_sbnd_evt166738.json",
        "mcs_sbnd_evt168526.json",
        "mcs_sbnd_evt168614.json",
        "mcs_sbnd_evt169356.json",
    };
    MuonMCS mcs;  // default = FIXED mode, the production configuration
    for (const auto& name : fixtures) {
        INFO("fixture ", name);
        auto fx = WireCell::Persist::load(data_dir() + "/" + name);
        VD vstart, vend;
        for (const auto& x : fx["input"]["start"]) vstart.push_back(jd(x));
        for (const auto& x : fx["input"]["end"]) vend.push_back(jd(x));
        VVD points;
        for (const auto& p : fx["input"]["points"]) {
            points.push_back({ jd(p[0]), jd(p[1]), jd(p[2]) });
        }

        McsResult ref = mcs.run(vstart, vend, points);
        for (unsigned seed = 1; seed <= 20; seed++) {
            INFO("seed ", seed);
            VVD shuffled = points;
            std::mt19937 rng(seed);
            std::shuffle(shuffled.begin(), shuffled.end(), rng);
            McsResult r = mcs.run(vstart, vend, shuffled);
            // byte-identical outputs required
            CHECK(r.mu_tracklen == ref.mu_tracklen);
            CHECK(r.emu_tracklen == ref.emu_tracklen);
            CHECK(r.emu_MCS == ref.emu_MCS);
            CHECK(r.ambiguity_MCS == ref.ambiguity_MCS);
            CHECK(r.nsegs == ref.nsegs);
            CHECK(r.bad_path == ref.bad_path);
        }
    }
}
