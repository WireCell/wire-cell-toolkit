// vertex_method (TaggerCheckNeutrino): the names the configuration accepts,
// the placement constant of dv2_split_move, how a DV2-MAX payload is read (the
// no-prediction fallback), and the configure-time failures of the dv2_*
// methods.  The placements on a real PR graph (split, long-muon sets, move)
// are in doctest_pattern_recognition.cxx ("pattern_recognition dv2 ...").

#include "WireCellClus/NeutrinoPatternBase.h"
#include "WireCellClus/TaggerCheckNeutrino.h"
#include "WireCellUtil/BuildConfig.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace WireCell;
using namespace WireCell::Clus::PR;
using VM = PatternAlgorithms::VertexMethod;
using Payload = PatternAlgorithms::Dv2Payload;

TEST_CASE("vertex_method: every accepted name parses and round-trips")
{
    for (const std::string name : {"original", "dv2_snap", "dv2_split", "dv2_split_move"}) {
        VM m = VM::original;
        REQUIRE(PatternAlgorithms::parse_vertex_method(name, m));
        CHECK(std::string(PatternAlgorithms::vertex_method_name(m)) == name);
    }
}

TEST_CASE("vertex_method: unknown and retired names are rejected and leave the method unchanged")
{
    // "scn" / "dv2max" were the values of the retired dl_vtx_engine knob.
    for (const std::string name : {"", "scn", "dv2max", "dv2", "snap", "split", "split_move", "Original", "dv2_split+move"}) {
        VM m = VM::dv2_split;
        CHECK_FALSE(PatternAlgorithms::parse_vertex_method(name, m));
        CHECK(m == VM::dv2_split);
    }
}

TEST_CASE("vertex_method: the default is the original chain, and only dv2_* methods are DV2-MAX")
{
    PatternAlgorithms pa;
    CHECK(pa.m_vertex_method == VM::original);
    CHECK_FALSE(pa.is_dv2());
    CHECK_FALSE(pa.m_dv2_target_valid);
    for (VM m : {VM::dv2_snap, VM::dv2_split, VM::dv2_split_move}) {
        pa.m_vertex_method = m;
        CHECK(pa.is_dv2());
    }
}

TEST_CASE("vertex_method: dv2_split_move reuses a PR-graph vertex within 1.5 cm (the MyFCN move limit)")
{
    CHECK(PatternAlgorithms::kDv2MoveReuseDistCm == doctest::Approx(1.5));
}

// ---- the DV2-MAX payload (determine_overall_main_vertex_DL keeps the
// traditional vertex for every status but ok) -------------------------------

TEST_CASE("vertex_method: an empty DV2-MAX payload is no prediction")
{
    const auto p = PatternAlgorithms::dv2_parse_payload({});
    CHECK(p.status == Payload::Status::no_prediction);
}

TEST_CASE("vertex_method: a DV2-MAX payload that is not [x, y, z, score] is not a prediction")
{
    for (size_t n : {1, 2, 3, 5, 8}) {
        CAPTURE(n);
        const auto p = PatternAlgorithms::dv2_parse_payload(std::vector<float>(n, 1.0f));
        CHECK(p.status == Payload::Status::bad_size);
    }
}

TEST_CASE("vertex_method: a DV2-MAX point with a NaN or infinite coordinate is rejected")
{
    const float bad[] = {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::infinity(),
                         -std::numeric_limits<float>::infinity()};
    for (int k = 0; k < 3; ++k) {
        for (float b : bad) {
            std::vector<float> v{10.0f, 20.0f, 30.0f, 0.9f};
            v[k] = b;
            CAPTURE(k);
            CAPTURE(b);
            CHECK(PatternAlgorithms::dv2_parse_payload(v).status == Payload::Status::nonfinite);
        }
    }
    // only the position is checked: a non-finite score does not reject the point
    const auto p = PatternAlgorithms::dv2_parse_payload({10.0f, 20.0f, 30.0f, std::numeric_limits<float>::quiet_NaN()});
    CHECK(p.status == Payload::Status::ok);
}

TEST_CASE("vertex_method: a finite DV2-MAX payload is the point (cm -> WCT units) and its score")
{
    const std::vector<float> v{-64.0228f, 135.6823f, 246.2073f, 0.848f};
    const auto p = PatternAlgorithms::dv2_parse_payload(v);
    REQUIRE(p.status == Payload::Status::ok);
    CHECK(p.P.x() == v[0] * units::cm);
    CHECK(p.P.y() == v[1] * units::cm);
    CHECK(p.P.z() == v[2] * units::cm);
    CHECK(p.score == double(v[3]));
}

// ---- configure-time failures (TaggerCheckNeutrino) --------------------------
// The vertex_method / dv2_bundle / model-load checks run before configure looks
// up any other component, so a bare TaggerCheckNeutrino reaches them.

static std::string configure_error(const Configuration& cfg)
{
    TaggerCheckNeutrino tcn;
    try {
        tcn.configure(cfg);
    }
    catch (const std::exception& e) {
        return e.what();
    }
    return "";
}

static Configuration dv2_config(const std::string& method, const std::string& bundle)
{
    TaggerCheckNeutrino tcn;
    Configuration cfg = tcn.default_configuration();
    cfg["vertex_method"] = method;
    cfg["dv2_bundle"] = bundle;
    return cfg;
}

TEST_CASE("vertex_method: configure rejects an unknown vertex_method")
{
    const auto msg = configure_error(dv2_config("dv2", "dv2max/sbnd/r7_cloud/bundle.json"));
    CHECK(msg.find("vertex_method must be original, dv2_snap, dv2_split or dv2_split_move") != std::string::npos);
}

TEST_CASE("vertex_method: configure fails for a dv2_* method without a bundle")
{
    for (const std::string method : {"dv2_snap", "dv2_split", "dv2_split_move"}) {
        for (const std::string bundle : {"", "/wct-doctest-no-such-dir/bundle.json", "wct-doctest-no-such-dir/bundle.json"}) {
            CAPTURE(method);
            CAPTURE(bundle);
            const auto msg = configure_error(dv2_config(method, bundle));
#ifdef HAVE_PYTHON_INC
            CHECK(msg.find("needs a DV2-MAX bundle") != std::string::npos);
            CHECK(msg.find("not found on WIRECELL_PATH") != std::string::npos);
#else
            CHECK(msg.find("needs the embedded Python") != std::string::npos);
#endif
        }
    }
}

#ifdef HAVE_PYTHON_INC
// $TMPDIR may name a directory that does not exist (e.g. a batch job's TMPDIR
// inside a container): fall back to /tmp, then the working directory.
static std::filesystem::path scratch_dir()
{
    std::error_code ec;
    auto p = std::filesystem::temp_directory_path(ec);
    if (!ec && std::filesystem::is_directory(p, ec)) return p;
    if (std::filesystem::is_directory("/tmp", ec)) return "/tmp";
    return std::filesystem::current_path();
}

TEST_CASE("vertex_method: configure fails when the DV2-MAX module cannot be imported")
{
    // An existing bundle file, so the import is the first thing to fail.
    const auto bundle = scratch_dir() / "wct_doctest_vertex_method_bundle.json";
    { std::ofstream(bundle) << "{}\n"; }
    REQUIRE(std::filesystem::is_regular_file(bundle));
    auto cfg = dv2_config("dv2_split_move", bundle.string());
    cfg["dv2_py_module"] = "wct_doctest_no_such_module";
    const auto msg = configure_error(cfg);
    std::filesystem::remove(bundle);
    CHECK(msg.find("cannot load DV2-MAX") != std::string::npos);
    CHECK(msg.find("wct_doctest_no_such_module") != std::string::npos);
    CHECK(msg.find("ModuleNotFoundError") != std::string::npos);
}
#endif
