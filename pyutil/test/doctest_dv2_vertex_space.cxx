// WCPPyUtil::DV2_Vertex with space points (the overload for a bundle of format
// dv2max-bundle-v2): the five extra arguments reach Python as bytes, in order,
// with the cluster ids as int32; size checks throw before Python is called; the
// cloud-only overload still passes none.  Against a fake module (no model,
// numpy or torch needed).

#include <Python.h>

#include "WCPPyUtil/DV2_Vertex.h"
#include "WCPPyUtil/Interpreter.h"
#include "WireCellUtil/doctest.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

using WCPPyUtil::FLOAT;

namespace {
    struct SpaceGIL {
        PyGILState_STATE s;
        SpaceGIL() : s(PyGILState_Ensure()) {}
        ~SpaceGIL() { PyGILState_Release(s); }
    };

    const char* kSpaceModule = "wct_doctest_dv2_fake_space";

    const char* kSpaceFake = R"PY(
import struct


def DV2_Vertex(weights, x, y, z, q, n_vertex_rows, trad, q_scale, q_offset, dtype="float32", *space):
    if weights == "cloud":
        if space:
            raise ValueError("the cloud-only call passed %d space arguments" % len(space))
        return struct.pack("<4f", 0.0, 0.0, 0.0, -1.0)
    if len(space) != 5:
        raise ValueError("%d space arguments" % len(space))
    sx, sy, sz, sq, sc = space
    n = len(sq) // 4
    if not (len(sx) == len(sy) == len(sz) == len(sq) == len(sc) == 4 * n):
        raise ValueError("space sizes %d %d %d %d %d" % (len(sx), len(sy), len(sz), len(sq), len(sc)))
    if n == 0:
        return b""
    ids = struct.unpack("<%di" % n, sc)
    # echo: last space x, first space q, sum of the int32 ids, number of points
    return struct.pack("<4f", struct.unpack("<f", sx[-4:])[0], struct.unpack("<f", sq[:4])[0], float(sum(ids)), float(n))
)PY";

    std::filesystem::path space_scratch_dir()
    {
        std::error_code ec;
        auto p = std::filesystem::temp_directory_path(ec);
        if (!ec && std::filesystem::is_directory(p, ec)) return p;
        if (std::filesystem::is_directory("/tmp", ec)) return "/tmp";
        return std::filesystem::current_path();
    }

    struct SpaceFakeModule {
        std::filesystem::path dir;
        SpaceFakeModule()
        {
            dir = space_scratch_dir() / ("wct_doctest_dv2_space_" + std::to_string(::getpid()));
            std::filesystem::create_directories(dir);
            std::ofstream(dir / (std::string(kSpaceModule) + ".py")) << kSpaceFake;
            REQUIRE(std::filesystem::is_regular_file(dir / (std::string(kSpaceModule) + ".py")));
            WCPPyUtil::ensure_interpreter();
            const std::string cmd = "import sys; sys.path.insert(0, r'" + dir.string() + "')";
            SpaceGIL g;
            REQUIRE(PyRun_SimpleString(cmd.c_str()) == 0);
        }
        ~SpaceFakeModule()
        {
            std::error_code ec;
            std::filesystem::remove_all(dir, ec);
        }
    };

    void space_fake_module() { static SpaceFakeModule fm; }

    const std::vector<std::vector<FLOAT>> kSCloud{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}};
    const std::vector<FLOAT> kSTrad{0.25f, 0, 0};

    std::vector<FLOAT> scall(const std::vector<std::vector<FLOAT>>& space, const std::vector<int>& ids,
                             const std::string& weights = "space")
    {
        space_fake_module();
        return WCPPyUtil::DV2_Vertex(kSpaceModule, "DV2_Vertex", weights, kSCloud, 2, kSTrad, 0.5, -7.0, "float32",
                                     space, ids);
    }

    template <typename F>
    std::string serror_of(F f)
    {
        try {
            f();
        }
        catch (const std::runtime_error& e) {
            return e.what();
        }
        return "";
    }

    bool shas(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }
}

TEST_CASE("DV2_Vertex space: the five space arrays reach Python in order, ids as int32")
{
    const std::vector<std::vector<FLOAT>> sp{{1.5f, 2.5f, 3.5f}, {0, 0, 0}, {0, 0, 0}, {10.f, 20.f, 30.f}};
    const auto v = scall(sp, {7, 7, 42});
    REQUIRE(v.size() == 4);
    CHECK(v[0] == 3.5f);      // last space x
    CHECK(v[1] == 10.f);      // first space q
    CHECK(v[2] == 56.f);      // 7 + 7 + 42
    CHECK(v[3] == 3.f);       // points
    const std::vector<int> neg{-3, 100000, 5};
    CHECK(scall(sp, neg)[2] == 100002.f);
}

TEST_CASE("DV2_Vertex space: an empty space set is passed as five empty arrays")
{
    const std::vector<std::vector<FLOAT>> none(4);
    CHECK(scall(none, {}).empty());      // the fake returns no prediction for zero points
}

TEST_CASE("DV2_Vertex space: the cloud-only overload passes no space arguments")
{
    space_fake_module();
    const auto v = WCPPyUtil::DV2_Vertex(kSpaceModule, "DV2_Vertex", "cloud", kSCloud, 2, kSTrad, 0.5, -7.0, "float32");
    REQUIRE(v.size() == 4);
    CHECK(v[3] == -1.f);
}

TEST_CASE("DV2_Vertex space: bad space sizes throw before Python is called")
{
    const std::vector<std::vector<FLOAT>> three{{1}, {2}, {3}};
    const std::vector<std::vector<FLOAT>> ragged{{1, 2}, {3}, {4}, {5}};
    const std::vector<std::vector<FLOAT>> ok{{1}, {2}, {3}, {4}};
    CHECK(shas(serror_of([&] { scall(three, {1}); }), "space.size() != 4"));
    CHECK(shas(serror_of([&] { scall(ragged, {1}); }), "space size mismatch"));
    CHECK(shas(serror_of([&] { scall(ok, {}); }), "space_cluster_id size mismatch"));
    CHECK(shas(serror_of([&] { scall(ok, {1, 2}); }), "space_cluster_id size mismatch"));
    CHECK(scall(ok, {9}).size() == 4);
}

TEST_CASE("DV2_Vertex space: a Python error in the space call throws with its message")
{
    const std::vector<std::vector<FLOAT>> ok{{1}, {2}, {3}, {4}};
    const auto msg = serror_of([&] { scall(ok, {9}, "cloud"); });     // the fake's cloud path refuses space args
    CHECK(shas(msg, "ValueError: the cloud-only call passed 5 space arguments"));
}
