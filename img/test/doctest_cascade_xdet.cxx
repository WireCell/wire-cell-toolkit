/** CascadeDeghosting's xdet_columns_file test path (wcfm doc 53 W7).

    - knob off: the forward gets the 9 tensors of the shipped models and the level dumps carry no xdet arrays;
    - knob on: the forward gets the 12 positional inputs of xdet/lean_x.py CascadeGNNLeanX (xb, xb_cols[N,8],
      wq, wplane, xw_cols[M,4], live[M], the six edge arrays) with the file's values in the level's node order, and
      the dumps carry the injected arrays; the output is the knob-off output (the fake forward ignores the columns);
    - a node missing from the file and a key given twice with different values are refused.

    The column file is written here from the knob-off run's own dumps (keys: face, slice ident, wire bounds; plane,
    first channel), with synthetic values that are a function of the key.
*/
#include "WireCellImg/CascadeDeghosting.h"
#include "WireCellImg/GridTiling.h"
#include "WireCellImg/BlobClustering.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/Testing.h"
#include "WireCellIface/IAnodeFace.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/cnpy.h"
#include "WireCellUtil/doctest.h"

#include <algorithm>
#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

using namespace WireCell;
using namespace WireCell::Img;
namespace fs = std::filesystem;

namespace {
    // per forward call of the knob-on run: the three injected tensors as received
    struct Seen {
        size_t N{0}, M{0};
        std::vector<float> xb_cols, xw_cols, live;
    };
    std::vector<Seen> g_seen;

    ITensorSet::pointer keep_all(const ITensorSet::pointer& in, size_t N)
    {
        std::vector<float> logit(N, 10.0f), qhat(N, 1.0f);
        ITensor::vector out{std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{N}, logit.data()),
                            std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{N}, qhat.data())};
        return std::make_shared<Aux::SimpleTensorSet>(in->ident(), Configuration(), std::make_shared<ITensor::vector>(out));
    }

    template <typename T>
    std::vector<T> copy_of(const ITensor::pointer& t)
    {
        const T* p = reinterpret_cast<const T*>(t->data());
        return std::vector<T>(p, p + t->size() / sizeof(T));
    }

    struct XdetOffForward : public ITensorForward {
        virtual ~XdetOffForward() {}
        virtual ITensorSet::pointer forward(const ITensorSet::pointer& in) const
        {
            auto tv = in->tensors();
            REQUIRE(tv->size() == 9);
            REQUIRE(tv->at(0)->shape()[1] == 15);
            return keep_all(in, tv->at(0)->shape()[0]);
        }
    };

    struct XdetOnForward : public ITensorForward {
        virtual ~XdetOnForward() {}
        virtual ITensorSet::pointer forward(const ITensorSet::pointer& in) const
        {
            auto tv = in->tensors();
            REQUIRE(tv->size() == 12);
            const auto& xb = tv->at(0);
            REQUIRE(xb->dtype() == "f4");
            REQUIRE(xb->shape().size() == 2);
            REQUIRE(xb->shape()[1] == 15);
            const size_t N = xb->shape()[0];
            const size_t M = tv->at(2)->shape()[0];
            // lean_x order: xb, xb_cols, wq, wplane, xw_cols, live, bw_src, bw_dst, bw_w, bb, bb_in, ww
            REQUIRE(tv->at(1)->dtype() == "f4");
            REQUIRE(tv->at(1)->shape() == ITensor::shape_t{N, 8});
            REQUIRE(tv->at(2)->dtype() == "f4");
            REQUIRE(tv->at(3)->dtype() == "i8");
            REQUIRE(tv->at(3)->shape() == ITensor::shape_t{M});
            REQUIRE(tv->at(4)->dtype() == "f4");
            REQUIRE(tv->at(4)->shape() == ITensor::shape_t{M, 4});
            REQUIRE(tv->at(5)->dtype() == "f4");
            REQUIRE(tv->at(5)->shape() == ITensor::shape_t{M});
            REQUIRE(tv->at(6)->dtype() == "i8");
            REQUIRE(tv->at(7)->dtype() == "i8");
            REQUIRE(tv->at(8)->dtype() == "f4");
            for (int k = 9; k < 12; ++k) REQUIRE(tv->at(k)->shape().size() == 2);
            Seen s;
            s.N = N;
            s.M = M;
            s.xb_cols = copy_of<float>(tv->at(1));
            s.xw_cols = copy_of<float>(tv->at(4));
            s.live = copy_of<float>(tv->at(5));
            g_seen.push_back(std::move(s));
            return keep_all(in, N);
        }
    };
}  // namespace

WIRECELL_FACTORY(XdetOffForward, XdetOffForward, WireCell::ITensorForward)
WIRECELL_FACTORY(XdetOnForward, XdetOnForward, WireCell::ITensorForward)

namespace {
    ISlice::pointer make_slice(IAnodeFace::pointer face, IFrame::pointer frame, int slice_ident, double start,
                               size_t nchan, double first = 0.5)
    {
        ISlice::map_t activity;
        for (const auto& plane : face->planes()) {
            const auto& chans = plane->channels();
            const size_t lo = (size_t) (chans.size() * first);
            for (size_t i = lo; i < std::min(lo + nchan, chans.size()); ++i) {
                activity[chans[i]] = ISlice::value_t(1000.0, 10.0);
            }
        }
        return std::make_shared<Aux::SimpleSlice>(frame, slice_ident, start, 2.0 * units::microsecond, activity);
    }

    ICluster::pointer make_cluster(size_t nchan)
    {
        auto anodes = Testing::anodes("uboone");
        REQUIRE(anodes.size() > 0);
        auto face = anodes[0]->face(0);
        Img::GridTiling gt;
        auto gcfg = gt.default_configuration();
        gcfg["anode"] = "AnodePlane:0";
        gcfg["face"] = 0;
        gt.configure(gcfg);
        Img::BlobClustering bcl;
        bcl.configure(bcl.default_configuration());
        auto frame = std::make_shared<Aux::SimpleFrame>(7, 0.0);
        Img::BlobClustering::output_queue q;
        for (int s = 0; s < 2; ++s) {
            IBlobSet::pointer bs;
            REQUIRE(gt(make_slice(face, frame, s, s * 2.0 * units::microsecond, nchan), bs));
            REQUIRE(bcl(bs, q));
        }
        REQUIRE(bcl(nullptr, q));
        REQUIRE(q.size() == 2);
        REQUIRE(q.front() != nullptr);
        return q.front();
    }

    std::vector<IBlob::pointer> blobs_of(const cluster_graph_t& gr)
    {
        std::vector<IBlob::pointer> v;
        for (auto vtx : boost::make_iterator_range(boost::vertices(gr))) {
            if (gr[vtx].code() == 'b') v.push_back(std::get<IBlob::pointer>(gr[vtx].ptr));
        }
        return v;
    }

    Configuration cascade_config(const std::string& fwd, const std::string& dump_dir)
    {
        Configuration cfg;
        const int ks[5] = {16, 8, 4, 2, 1};
        const int ws[5] = {0, 32, 16, 8, 4};
        for (int i = 0; i < 5; ++i) {
            Configuration l;
            l["width"] = ws[i];
            l["superwire"] = ks[i];
            l["forward"] = fwd;
            l["threshold"] = i < 4 ? -3.0 : -0.2;
            cfg["levels"].append(l);
        }
        cfg["charge_scale"] = 0.25;
        cfg["repair"] = false;
        cfg["dump_dir"] = dump_dir;
        return cfg;
    }

    void setup_components()
    {
        Testing::load_plugins({"WireCellImg"});
        make_XdetOffForward_factory();
        make_XdetOnForward_factory();
        Factory::lookup_tn<ITensorForward>("XdetOffForward");
        Factory::lookup_tn<ITensorForward>("XdetOnForward");
    }

    std::string fresh_dir(const std::string& name)
    {
        fs::remove_all(name);
        fs::create_directories(name);
        return name;
    }

    // dumps cascade-<ident>-L<i>.npz of one cluster, by level
    std::vector<cnpy::npz_t> load_dumps(const std::string& dir, int nlev)
    {
        std::vector<cnpy::npz_t> out;
        for (int li = 0; li < nlev; ++li) {
            std::string fn;
            for (const auto& e : fs::directory_iterator(dir)) {
                const auto n = e.path().filename().string();
                const std::string suf = "-L" + std::to_string(li) + ".npz";
                if (n.size() > suf.size() && n.compare(n.size() - suf.size(), suf.size(), suf) == 0) fn = e.path().string();
            }
            REQUIRE(!fn.empty());
            out.push_back(cnpy::npz_load(fn));
        }
        return out;
    }

    // synthetic columns, a function of the key
    float bcol(const int32_t* k, int j) { return 0.5f * j + k[2] + 0.01f * k[4] + 1e-4f * k[6] + 7.0f * k[0] + 0.125f * k[1]; }
    float wcol(int plane, int ch, int j) { return plane + 0.25f * j + 1e-3f * ch; }
    float wlive(int ch) { return (ch % 3 == 0) ? 0.0f : 1.0f; }

    struct LevelKeys {
        std::vector<int32_t> bkey;   // N x 8
        std::vector<int32_t> wkey;   // M x 2
    };

    LevelKeys keys_of(cnpy::npz_t& z)
    {
        LevelKeys lk;
        const size_t N = z["face"].shape[0];
        const int32_t* wip = z["wip"].data<int32_t>();
        const int32_t* face = z["face"].data<int32_t>();
        const int32_t* sl = z["slice_ident"].data<int32_t>();
        for (size_t i = 0; i < N; ++i) {
            lk.bkey.push_back(face[i]);
            lk.bkey.push_back(sl[i]);
            for (int j = 0; j < 6; ++j) lk.bkey.push_back(wip[6 * i + j]);
        }
        const size_t M = z["wplane"].shape[0];
        const int64_t* wp = z["wplane"].data<int64_t>();
        const int32_t* wc = z["wchan"].data<int32_t>();
        for (size_t m = 0; m < M; ++m) {
            lk.wkey.push_back((int32_t) wp[m]);
            lk.wkey.push_back(wc[m]);
        }
        return lk;
    }

    // write the column file from the dumps' keys; drop_level >= 0 leaves out that level's first blob row,
    // dup_level >= 0 repeats that level's first blob row with a different value
    void write_columns(const std::string& fn, std::vector<cnpy::npz_t>& dumps, int drop_level = -1, int dup_level = -1)
    {
        std::string mode = "w";
        for (size_t li = 0; li < dumps.size(); ++li) {
            auto lk = keys_of(dumps[li]);
            size_t N = lk.bkey.size() / 8, M = lk.wkey.size() / 2;
            std::vector<float> bc(N * 8), wc(M * 4), lv(M);
            for (size_t i = 0; i < N; ++i)
                for (int j = 0; j < 8; ++j) bc[8 * i + j] = bcol(&lk.bkey[8 * i], j);
            for (size_t m = 0; m < M; ++m) {
                for (int j = 0; j < 4; ++j) wc[4 * m + j] = wcol(lk.wkey[2 * m], lk.wkey[2 * m + 1], j);
                lv[m] = wlive(lk.wkey[2 * m + 1]);
            }
            if ((int) li == drop_level && N > 0) {
                lk.bkey.erase(lk.bkey.begin(), lk.bkey.begin() + 8);
                bc.erase(bc.begin(), bc.begin() + 8);
                --N;
            }
            if ((int) li == dup_level && N > 0) {
                std::vector<int32_t> k0(lk.bkey.begin(), lk.bkey.begin() + 8);
                std::vector<float> c0(bc.begin(), bc.begin() + 8);
                c0[0] += 1.0f;
                lk.bkey.insert(lk.bkey.end(), k0.begin(), k0.end());
                bc.insert(bc.end(), c0.begin(), c0.end());
                ++N;
            }
            const std::string pre = "L" + std::to_string(li) + "_";
            cnpy::npz_save(fn, pre + "bkey", lk.bkey.data(), {N, 8}, mode);
            mode = "a";
            cnpy::npz_save(fn, pre + "bcols", bc.data(), {N, 8}, mode);
            cnpy::npz_save(fn, pre + "wkey", lk.wkey.data(), {M, 2}, mode);
            cnpy::npz_save(fn, pre + "wcols", wc.data(), {M, 4}, mode);
            cnpy::npz_save(fn, pre + "live", lv.data(), {M}, mode);
        }
    }

    // face, slice ident, U / V / W strip bounds (the file's blob key)
    std::vector<int> key_of(const IBlob::pointer& b)
    {
        std::vector<int> k{b->face()->which(), b->slice()->ident(), 0, 0, 0, 0, 0, 0};
        for (const auto& s : b->shape().strips()) {
            const int ip = s.layer - 2;
            if (ip < 0 || ip > 2) continue;
            k[2 + 2 * ip] = s.bounds.first;
            k[3 + 2 * ip] = s.bounds.second;
        }
        return k;
    }

    ICluster::pointer run(const Configuration& cfg, const ICluster::pointer& cl)
    {
        Img::CascadeDeghosting cd;
        cd.configure(cfg);
        ICluster::pointer out;
        REQUIRE(cd(CascadeDeghosting::input_tuple_type{cl, nullptr}, out));
        REQUIRE(out);
        return out;
    }
}  // namespace

TEST_CASE("xdet_columns_file: off = 9 inputs and no xdet arrays; on = lean_x's 12 inputs with the file's values")
{
    setup_components();
    auto cl = make_cluster(60);

    // knob off
    const std::string doff = fresh_dir("doctest-cascade-xdet-off");
    auto o_off = run(cascade_config("XdetOffForward", doff), cl);
    auto dumps = load_dumps(doff, 5);
    for (auto& z : dumps) {
        CHECK(z.count("xb_cols") == 0);
        CHECK(z.count("xw_cols") == 0);
        CHECK(z.count("live") == 0);
    }

    // the file from the knob-off dumps (keep-everything forwards: the knob-on run builds the same level graphs)
    const std::string fn = "doctest-cascade-xdet-columns.npz";
    write_columns(fn, dumps);

    const std::string don = fresh_dir("doctest-cascade-xdet-on");
    auto cfg = cascade_config("XdetOnForward", don);
    cfg["xdet_columns_file"] = fn;
    g_seen.clear();
    auto o_on = run(cfg, cl);
    REQUIRE(g_seen.size() == 5);
    auto dumps_on = load_dumps(don, 5);
    for (size_t li = 0; li < 5; ++li) {
        auto& z = dumps_on[li];
        auto lk = keys_of(z);
        const size_t N = lk.bkey.size() / 8, M = lk.wkey.size() / 2;
        CHECK(N > 0);
        REQUIRE(g_seen[li].N == N);
        REQUIRE(g_seen[li].M == M);
        // what the forward received is the file's value of each node's key, in node order
        for (size_t i = 0; i < N; ++i)
            for (int j = 0; j < 8; ++j) CHECK(g_seen[li].xb_cols[8 * i + j] == bcol(&lk.bkey[8 * i], j));
        for (size_t m = 0; m < M; ++m) {
            for (int j = 0; j < 4; ++j) CHECK(g_seen[li].xw_cols[4 * m + j] == wcol(lk.wkey[2 * m], lk.wkey[2 * m + 1], j));
            CHECK(g_seen[li].live[m] == wlive(lk.wkey[2 * m + 1]));
        }
        // and the dump carries the same arrays
        REQUIRE(z.count("xb_cols") == 1);
        REQUIRE(z["xb_cols"].shape == std::vector<size_t>{N, 8});
        REQUIRE(z["xw_cols"].shape == std::vector<size_t>{M, 4});
        REQUIRE(z["live"].shape == std::vector<size_t>{M});
        CHECK(std::equal(g_seen[li].xb_cols.begin(), g_seen[li].xb_cols.end(), z["xb_cols"].data<float>()));
        CHECK(std::equal(g_seen[li].xw_cols.begin(), g_seen[li].xw_cols.end(), z["xw_cols"].data<float>()));
        CHECK(std::equal(g_seen[li].live.begin(), g_seen[li].live.end(), z["live"].data<float>()));
        // the knob changes nothing else in the level graph
        for (const char* name : {"xb", "wq", "bw_src", "bw_dst", "bb", "bb_in", "ww", "logit", "decision"}) {
            auto& a = dumps[li][name];
            auto& b = z[name];
            REQUIRE(a.shape == b.shape);
            CHECK(std::equal(a.data<char>(), a.data<char>() + a.num_bytes(), b.data<char>()));
        }
    }
    // same output blobs
    auto b_off = blobs_of(o_off->graph());
    auto b_on = blobs_of(o_on->graph());
    REQUIRE(b_off.size() == b_on.size());
    for (size_t i = 0; i < b_off.size(); ++i) CHECK(key_of(b_off[i]) == key_of(b_on[i]));

    // a node missing from the file is refused
    {
        const std::string fmiss = "doctest-cascade-xdet-missing.npz";
        write_columns(fmiss, dumps, 2, -1);
        auto c = cascade_config("XdetOnForward", "");
        c["xdet_columns_file"] = fmiss;
        Img::CascadeDeghosting cd;
        cd.configure(c);
        ICluster::pointer out;
        CHECK_THROWS(cd(CascadeDeghosting::input_tuple_type{cl, nullptr}, out));
    }
    // a key given twice with different values is refused when the file is read
    {
        const std::string fdup = "doctest-cascade-xdet-dup.npz";
        write_columns(fdup, dumps, -1, 3);
        auto c = cascade_config("XdetOnForward", "");
        c["xdet_columns_file"] = fdup;
        Img::CascadeDeghosting cd;
        CHECK_THROWS(cd.configure(c));
    }
    // a file with fewer levels than the cascade is refused at configure
    {
        const std::string fshort = "doctest-cascade-xdet-short.npz";
        std::vector<cnpy::npz_t> two(dumps.begin(), dumps.begin() + 2);
        write_columns(fshort, two);
        auto c = cascade_config("XdetOnForward", "");
        c["xdet_columns_file"] = fshort;
        Img::CascadeDeghosting cd;
        CHECK_THROWS(cd.configure(c));
    }
}
