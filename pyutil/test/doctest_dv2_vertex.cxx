// WCPPyUtil::DV2_Load / DV2_Vertex: the embedded-Python wrapper of the DV2-MAX
// vertexer (vertex_method dv2_*), against a fake module written to a temp
// directory -- no model, numpy or torch needed.  Covers the argument
// marshalling, "no prediction" (empty bytes), and that every Python error
// reaches the caller as std::runtime_error with the Python type and message
// (which TaggerCheckNeutrino turns into a configure or job failure).

#include <Python.h>

#include "WCPPyUtil/DV2_Vertex.h"
#include "WireCellUtil/Persist.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <filesystem>
#include <limits>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

using WCPPyUtil::FLOAT;

namespace {

    const char* kModule = "wct_doctest_dv2_fake";

    const char* kFake = R"PY(
EXPECT = (2, 0.5, -7.0, "float32")


def DV2_Load(weights):
    if weights == "raise":
        raise RuntimeError("injected load failure")
    return ("loaded " + weights).encode()


def DV2_Vertex(weights, x, y, z, q, n_vertex_rows, trad, q_scale, q_offset, dtype="float32"):
    if weights == "empty":
        return b""
    if weights == "raise":
        raise RuntimeError("injected inference failure")
    if weights == "ragged":
        return b"abcde"
    if weights == "notbytes":
        return 3
    got = (n_vertex_rows, q_scale, q_offset, dtype)
    if got != EXPECT:
        raise ValueError("arguments %r != %r" % (got, EXPECT))
    if not (len(x) == len(y) == len(z) == len(q) == 12 and len(trad) == 12):
        raise ValueError("sizes %d %d %d %d %d" % (len(x), len(y), len(z), len(q), len(trad)))
    # echo the second cloud point (x, y, z) and trad x as the score, as float32 bytes
    return x[4:8] + y[4:8] + z[4:8] + trad[0:4]
)PY";

    // A writable scratch directory: $TMPDIR may name a directory that does not
    // exist (e.g. a batch job's TMPDIR inside a container), so fall back to
    // /tmp, then the working directory.
    std::filesystem::path scratch_dir()
    {
        std::error_code ec;
        auto p = std::filesystem::temp_directory_path(ec);
        if (!ec && std::filesystem::is_directory(p, ec)) return p;
        if (std::filesystem::is_directory("/tmp", ec)) return "/tmp";
        return std::filesystem::current_path();
    }

    // Write the fake module once and put its directory first on sys.path.
    struct FakeModule {
        std::filesystem::path dir;
        FakeModule()
        {
            dir = scratch_dir() / ("wct_doctest_dv2_" + std::to_string(::getpid()));
            std::filesystem::create_directories(dir);
            std::ofstream(dir / (std::string(kModule) + ".py")) << kFake;
            REQUIRE(std::filesystem::is_regular_file(dir / (std::string(kModule) + ".py")));
            Py_Initialize();
            const std::string cmd = "import sys; sys.path.insert(0, r'" + dir.string() + "')";
            REQUIRE(PyRun_SimpleString(cmd.c_str()) == 0);
        }
        ~FakeModule()
        {
            std::error_code ec;
            std::filesystem::remove_all(dir, ec);
        }
    };

    void fake_module()
    {
        static FakeModule fm;
    }

    const std::vector<std::vector<FLOAT>> kCloud{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}};
    const std::vector<FLOAT> kTrad{0.25f, 0, 0};

    std::vector<FLOAT> call(const std::string& weights, const std::string& module = kModule,
                            const std::string& function = "DV2_Vertex")
    {
        fake_module();
        return WCPPyUtil::DV2_Vertex(module, function, weights, kCloud, 2, kTrad, 0.5, -7.0, "float32");
    }

    // The message of the runtime_error f() throws ("" if it does not throw).
    template <typename F>
    std::string error_of(F f)
    {
        try {
            f();
        }
        catch (const std::runtime_error& e) {
            return e.what();
        }
        return "";
    }

    bool has(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }
}

TEST_CASE("DV2_Vertex: arguments reach Python and the returned floats come back")
{
    const auto v = call("ok");
    REQUIRE(v.size() == 4);
    CHECK(v[0] == 2.0f);
    CHECK(v[1] == 5.0f);
    CHECK(v[2] == 8.0f);
    CHECK(v[3] == 0.25f);
}

TEST_CASE("DV2_Vertex: empty bytes are no prediction (an empty vector)")
{
    CHECK(call("empty").empty());
}

TEST_CASE("DV2_Vertex: a Python exception throws with its type and message")
{
    const auto msg = error_of([] { call("raise"); });
    CHECK(has(msg, "Python function call failed"));
    CHECK(has(msg, "RuntimeError: injected inference failure"));
}

TEST_CASE("DV2_Vertex: a malformed return throws")
{
    CHECK(has(error_of([] { call("ragged"); }), "not a whole number of floats"));
    CHECK(has(error_of([] { call("notbytes"); }), "did not return bytes"));
}

TEST_CASE("DV2_Vertex: a missing module or function throws")
{
    const auto msg = error_of([] { call("ok", "wct_doctest_no_such_module"); });
    CHECK(has(msg, "import failed for module"));
    CHECK(has(msg, "ModuleNotFoundError"));
    CHECK(has(error_of([] { call("ok", kModule, "no_such_function"); }), "function not found or not callable"));
}

TEST_CASE("DV2_Vertex: bad inputs throw before Python is called")
{
    fake_module();
    const std::vector<std::vector<FLOAT>> three{{1}, {2}, {3}};
    const std::vector<std::vector<FLOAT>> ragged{{1, 2}, {3}, {4}, {5}};
    CHECK(has(error_of([&] { WCPPyUtil::DV2_Vertex(kModule, "DV2_Vertex", "ok", three, 0, kTrad, 1, 0); }),
              "input.size() != 4"));
    CHECK(has(error_of([&] { WCPPyUtil::DV2_Vertex(kModule, "DV2_Vertex", "ok", ragged, 0, kTrad, 1, 0); }),
              "input size mismatch"));
    CHECK(has(error_of([&] { WCPPyUtil::DV2_Vertex(kModule, "DV2_Vertex", "ok", kCloud, 0, {1, 2}, 1, 0); }),
              "trad must have 3 values"));
    CHECK(has(error_of([&] { WCPPyUtil::DV2_Vertex(kModule, "DV2_Vertex", "ok", kCloud, 4, kTrad, 1, 0); }),
              "n_vertex_rows out of range"));
    CHECK(has(error_of([&] { WCPPyUtil::DV2_Vertex(kModule, "DV2_Vertex", "ok", kCloud, -1, kTrad, 1, 0); }),
              "n_vertex_rows out of range"));
}

TEST_CASE("DV2_Vertex: a failed call leaves no Python error behind")
{
    fake_module();
    REQUIRE(Py_IsInitialized());
    error_of([] { call("raise"); });
    CHECK_FALSE(PyErr_Occurred());
    CHECK(call("ok").size() == 4);
}

TEST_CASE("DV2_Load: returns the module's load report")
{
    fake_module();
    CHECK(WCPPyUtil::DV2_Load(kModule, "DV2_Load", "bundle.json") == "loaded bundle.json");
}

TEST_CASE("DV2_Load: every failure throws with the reason")
{
    fake_module();
    REQUIRE(Py_IsInitialized());
    {
        const auto msg = error_of([] { WCPPyUtil::DV2_Load("wct_doctest_no_such_module", "DV2_Load", "b"); });
        CHECK(has(msg, "import failed for module wct_doctest_no_such_module"));
        CHECK(has(msg, "ModuleNotFoundError"));
    }
    CHECK(has(error_of([] { WCPPyUtil::DV2_Load(kModule, "no_such_function", "b"); }),
              "function not found or not callable"));
    {
        const auto msg = error_of([] { WCPPyUtil::DV2_Load(kModule, "DV2_Load", "raise"); });
        CHECK(has(msg, "loading raise failed"));
        CHECK(has(msg, "RuntimeError: injected load failure"));
    }
    CHECK_FALSE(PyErr_Occurred());
}

// ---- the real DV2-MAX model through this wrapper ----------------------------
// Replays the first held-out calls of data/dv2max_golden.json (written by
// make_dv2max_golden.py; test_dv2max.py replays all of them in Python) with
// the real DV2_Vertex module and bundle.  Skipped, with the reason, when the
// fixture, the bundle (wire-cell-data dv2max/sbnd/r7_cloud/ on WIRECELL_PATH)
// or the model's Python packages are not available.  Run from the build
// directory, like the other fixture tests (../pyutil/...).

TEST_CASE("DV2_Vertex: the real DV2-MAX model reproduces the golden held-out calls")
{
    const std::string fixture = "../pyutil/test/data/dv2max_golden.json";
    if (!std::filesystem::is_regular_file(fixture)) {
        MESSAGE("Skipping: no " << fixture << " (run from the build directory)");
        return;
    }
    const std::string bundle = WireCell::Persist::resolve("dv2max/sbnd/r7_cloud/bundle.json");
    if (bundle.empty()) {
        MESSAGE("Skipping: dv2max/sbnd/r7_cloud/bundle.json not on WIRECELL_PATH");
        return;
    }
    fake_module();   // initialises the interpreter
    if (std::filesystem::is_directory("../pyutil/python")) {
        const std::string dir = std::filesystem::canonical("../pyutil/python").string();
        REQUIRE(PyRun_SimpleString(("import sys; sys.path.insert(0, r'" + dir + "')").c_str()) == 0);
    }
    if (PyRun_SimpleString("import numpy, scipy, torch, lightgbm, DV2_Vertex") != 0) {
        MESSAGE("Skipping: the DV2-MAX Python packages (numpy, scipy, torch, lightgbm, DV2_Vertex) do not import");
        return;
    }
    const auto g = WireCell::Persist::load(fixture);
    const double tol_cm = g["meta"]["tolerance"]["position_cm"].asDouble();
    const double tol_score = g["meta"]["tolerance"]["score"].asDouble();
    auto floats = [](const Json::Value& a) {
        std::vector<FLOAT> v;
        for (const auto& x : a) v.push_back(x.isNull() ? std::numeric_limits<FLOAT>::quiet_NaN() : x.asFloat());
        return v;
    };
    int n_real = 0, n_synth = 0;
    for (const auto& c : g["calls"]) {
        const bool synth = c["key"].asString().rfind("synthetic", 0) == 0;
        if ((synth && n_synth >= 1) || (!synth && n_real >= 4)) continue;
        (synth ? n_synth : n_real)++;
        const std::string key = c["key"].asString();
        CAPTURE(key);
        const std::vector<std::vector<FLOAT>> in{floats(c["x"]), floats(c["y"]), floats(c["z"]), floats(c["q"])};
        const auto out = WCPPyUtil::DV2_Vertex("DV2_Vertex", "DV2_Vertex", bundle, in, c["n_vertex_rows"].asInt(),
                                               floats(c["trad"]), c["q_scale"].asDouble(), c["q_offset"].asDouble());
        const auto expect = floats(c["expect"]);
        REQUIRE(out.size() == expect.size());
        if (expect.empty()) continue;
        const double d = std::sqrt(std::pow(out[0] - expect[0], 2) + std::pow(out[1] - expect[1], 2) +
                                   std::pow(out[2] - expect[2], 2));
        CAPTURE(d);
        CHECK(d < tol_cm);
        CHECK(std::abs(out[3] - expect[3]) < tol_score);
    }
    CHECK(n_real == 4);
    CHECK(n_synth == 1);
}
