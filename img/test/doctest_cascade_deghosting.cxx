/** CascadeDeghosting and its pieces (wcfm doc 14).

    - build_level: the node features, bw / bw_w edges and bins follow gnn_dataset / contract_level;
    - inslice_pairs: overlap or abut on all three planes, within a (face, slice) group only;
    - guard_prune: a candidate that is a wire's last explanation is kept;
    - steiner_repair: a gap of low-P cells between two confident fragments is bridged within the budget, a
      longer gap is not, a weak island is dropped;
    - the whole cascade with a keep-everything forward reproduces direct BlobCutting to 4 wires (same set of
      (slice, wire bounds)), and two runs give identical output graphs.

    The forward is an in-test ITensorForward (no model file, no torch): it checks the input schema and returns
    logits from a rule.
*/
#include "WireCellImg/CascadeDeghosting.h"
#include "WireCellImg/CascadeGraph.h"
#include "WireCellImg/CellSteiner.h"
#include "WireCellImg/GridTiling.h"
#include "WireCellImg/BlobClustering.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/Testing.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IWirePlane.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Testing.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <tuple>

using namespace WireCell;
using namespace WireCell::Img;

namespace {
    int g_mode = 0;   // 0: keep everything; 1: logit = +10 if the node has >= 3 U channels, else -10

    struct CascadeFakeForward : public ITensorForward {
        virtual ~CascadeFakeForward() {}
        virtual ITensorSet::pointer forward(const ITensorSet::pointer& in) const
        {
            auto tv = in->tensors();
            REQUIRE(tv->size() == 9);
            const auto& xb = tv->at(0);
            REQUIRE(xb->dtype() == "f4");
            REQUIRE(xb->shape().size() == 2);
            REQUIRE(xb->shape()[1] == 15);
            REQUIRE(tv->at(2)->dtype() == "i8");
            REQUIRE(tv->at(3)->dtype() == "i8");
            REQUIRE(tv->at(5)->dtype() == "f4");
            REQUIRE(tv->at(6)->shape().size() == 2);
            const size_t N = xb->shape()[0];
            const float* x = reinterpret_cast<const float*>(xb->data());
            std::vector<float> logit(N), qhat(N, 1.0f);
            for (size_t i = 0; i < N; ++i) {
                logit[i] = (g_mode == 0 || x[i * 15 + 0] >= 3) ? 10.0f : -10.0f;
            }
            ITensor::vector out{std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{N}, logit.data()),
                                std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{N}, qhat.data())};
            return std::make_shared<Aux::SimpleTensorSet>(in->ident(), Configuration(),
                                                          std::make_shared<ITensor::vector>(out));
        }
    };
}  // namespace

WIRECELL_FACTORY(CascadeFakeForward, CascadeFakeForward, WireCell::ITensorForward)

namespace {
    // nchan contiguous channels from `first` (fraction of the plane) in each plane
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

    // a cluster of two consecutive slices, tiled and clustered as in production
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

    using bkey = std::tuple<double, int, int, int, int, int, int>;
    bkey key_of(const IBlob::pointer& b)
    {
        int r[6] = {0, 0, 0, 0, 0, 0};
        for (const auto& s : b->shape().strips()) {
            const int ip = s.layer - 2;
            if (ip < 0 || ip > 2) continue;
            r[2 * ip] = s.bounds.first;
            r[2 * ip + 1] = s.bounds.second;
        }
        return {b->slice()->start(), r[0], r[1], r[2], r[3], r[4], r[5]};
    }

    void setup_components()
    {
        Testing::load_plugins({"WireCellImg"});
        make_CascadeFakeForward_factory();
        Factory::lookup_tn<ITensorForward>("CascadeFakeForward");   // create the instance configure() finds
    }

    Configuration cascade_config()
    {
        Configuration cfg;
        const int ks[5] = {16, 8, 4, 2, 1};
        const int ws[5] = {0, 32, 16, 8, 4};
        for (int i = 0; i < 5; ++i) {
            Configuration l;
            l["width"] = ws[i];
            l["superwire"] = ks[i];
            l["forward"] = "CascadeFakeForward";
            l["threshold"] = i < 4 ? -3.0 : -0.2;
            cfg["levels"].append(l);
        }
        cfg["charge_scale"] = 0.25;
        return cfg;
    }
}  // namespace

TEST_CASE("cascade level graph: features, bw weights, bins")
{
    auto cl = make_cluster(12);
    const auto& gr = cl->graph();
    auto sc = Cascade::make_slice_charge(gr, 0.25, 1e11);
    CHECK(sc.slice_of.size() == 2);
    auto blobs = blobs_of(gr);
    REQUIRE(blobs.size() > 0);
    auto uni = Cascade::make_universe(blobs, sc);
    for (int k : {1, 4}) {
        auto lev = Cascade::build_level(blobs, sc, uni, k, "uboone");
        REQUIRE(lev.nnodes() == blobs.size());
        std::vector<double> wsum(lev.nnodes() * 3, 0.0);
        for (size_t e = 0; e < lev.bw_src.size(); ++e) {
            wsum[lev.bw_src[e] * 3 + lev.wplane[lev.bw_dst[e]]] += lev.bw_w[e];
            if (k == 1) CHECK(lev.bw_w[e] == 1.0f);
        }
        for (size_t i = 0; i < lev.nnodes(); ++i) {
            const auto ch = Cascade::blob_channels(blobs[i]);
            double sums[3];
            for (int p = 0; p < 3; ++p) {
                const size_t nch = ch[p].size();
                // every channel of the test activity carries 1000 x 0.25
                const double sq = 250.0 * nch;
                CHECK(lev.xb[i * 15 + 4 * p + 0] == doctest::Approx(nch));
                CHECK(lev.xb[i * 15 + 4 * p + 1] == doctest::Approx(nch));
                CHECK(lev.xb[i * 15 + 4 * p + 2] == doctest::Approx(std::log1p(sq)));
                CHECK(lev.xb[i * 15 + 4 * p + 3] == doctest::Approx(nch ? std::log1p(sq / nch) : 0.0));
                CHECK(wsum[i * 3 + p] == doctest::Approx(nch));   // a node's channels, over its bins
                sums[p] = std::log1p(sq);
            }
            CHECK(lev.xb[i * 15 + 12] == doctest::Approx(sums[0] - sums[1]));
            CHECK(lev.xb[i * 15 + 14] == doctest::Approx(sums[1] - sums[2]));
        }
        // bins: wire charge sums to the universe total
        double tq = 0, tu = 0;
        for (float q : lev.wq) tq += q;
        for (float q : uni.q) tu += q;
        CHECK(tq == doctest::Approx(tu));
        if (k == 4) CHECK(lev.wq.size() < uni.wires.size());
        // the two slices overlap: some cross-slice edges, all between different slices
        REQUIRE(lev.bb.size() > 0);
        for (size_t e = 0; e < lev.bb.size(); e += 2) CHECK(lev.sidx[lev.bb[e]] != lev.sidx[lev.bb[e + 1]]);
        for (size_t e = 0; e < lev.bb_in.size(); e += 2) CHECK(lev.sidx[lev.bb_in[e]] == lev.sidx[lev.bb_in[e + 1]]);
    }
}

TEST_CASE("frame charge: scale x the sum over the slice's ticks, all ticks, 0 without a trace")
{
    auto cl = make_cluster(12);
    const auto& gr = cl->graph();
    // two traces: channel 100 from tick 0 with q[t] = t (negative ticks included), channel 101 from tick 6
    std::vector<float> q100(20), q101(4, 2.0f);
    for (int t = 0; t < 20; ++t) q100[t] = (t == 5) ? -3.0f : (float) t;
    ITrace::vector trs{std::make_shared<Aux::SimpleTrace>(100, 0, q100), std::make_shared<Aux::SimpleTrace>(101, 6, q101)};
    auto frame = std::make_shared<Aux::SimpleFrame>(7, 0.0, trs, 0.5 * units::microsecond);
    auto sc = Cascade::make_slice_charge_frame(gr, frame, "", 0.25);
    REQUIRE(sc.slice_of.size() == 2);            // slices at 0 and 2 us, 4 ticks each
    CHECK(sc.t0[0] == 0);
    CHECK(sc.t0[1] == 4);
    CHECK(sc.nt[0] == 4);
    CHECK(sc.charge(0, 100) == doctest::Approx(0.25 * (0 + 1 + 2 + 3)));
    CHECK(sc.charge(1, 100) == doctest::Approx(0.25 * (4 - 3 + 6 + 7)));
    CHECK(sc.charge(0, 101) == 0.0f);            // trace starts at tick 6
    CHECK(sc.charge(1, 101) == doctest::Approx(0.25 * 4.0));
    CHECK(sc.charge(1, 999) == 0.0f);
}

TEST_CASE("in-slice pairs: overlap or abut on all planes, within a group")
{
    std::vector<int64_t> g{0, 0, 0, 1};
    std::vector<std::array<int, 6>> r{{0, 4, 0, 4, 0, 4},     // a
                                      {4, 8, 2, 6, 3, 9},     // abuts a on U, overlaps V and W
                                      {9, 12, 0, 4, 0, 4},    // U gap to a (and to b: 8 < 9)
                                      {0, 4, 0, 4, 0, 4}};    // same ranges as a, other group
    auto p = Cascade::inslice_pairs(g, r);
    REQUIRE(p.size() == 1);
    CHECK(p[0][0] == 0);
    CHECK(p[0][1] == 1);
}

TEST_CASE("coverage guard keeps a wire's last explanation")
{
    Cascade::Level lev;
    lev.blobs.resize(3);
    lev.wq = {1, 1};
    lev.wplane = {0, 0};
    // node 0 -> wire 0; node 1 -> wires 0, 1; node 2 -> wire 1
    lev.bw_src = {0, 1, 1, 2};
    lev.bw_dst = {0, 0, 1, 1};
    size_t ng = 0;
    auto pr = Cascade::guard_prune(lev, {1, 0, 2}, ng);
    CHECK(pr[1]);          // first: both wires still have another node
    CHECK(!pr[0]);         // wire 0 would lose its last node
    CHECK(!pr[2]);
    CHECK(ng == 2);
}

TEST_CASE("steiner repair: bridge within budget, not beyond, weak island dropped")
{
    auto chain = [](size_t n) {
        std::vector<std::array<int64_t, 2>> e;
        for (size_t i = 0; i + 1 < n; ++i) e.push_back({(int64_t) i, (int64_t) i + 1});
        return e;
    };
    const float hi = 4.6f, lo = (float) std::log(0.2 / 0.8);   // P = 0.99, P = 0.2
    Cascade::SteinerParams par;   // p_term 0.8, q_floor 1e4, budget 5.5
    {
        // 0 1 [2 3 4] 5 6 : three P = 0.2 cells between two confident fragments; node 7 an isolated weak island
        std::vector<float> logit{hi, hi, lo, lo, lo, hi, hi, 0.4f};
        std::vector<float> qhat(8, 0.1f);
        std::vector<bool> keep{true, true, false, false, false, true, true, true};
        auto e = chain(7);   // node 7 unconnected
        auto res = Cascade::steiner_repair(8, e, logit, qhat, keep, par);
        CHECK(res.nbridges == 1);
        CHECK(res.nadded == 3);
        for (int i = 0; i < 7; ++i) CHECK(res.keep[i]);
        CHECK(!res.keep[7]);
        CHECK(res.nweak_cells == 1);
    }
    {
        // a gap of five P = 0.2 cells costs more than the budget
        std::vector<float> logit{hi, hi, lo, lo, lo, lo, lo, hi, hi};
        std::vector<float> qhat(9, 1.0f);
        std::vector<bool> keep{true, true, false, false, false, false, false, true, true};
        auto res = Cascade::steiner_repair(9, chain(9), logit, qhat, keep, par);
        CHECK(res.nbridges == 0);
        CHECK(res.nadded == 0);
    }
}

TEST_CASE("cascade with a keep-everything forward equals BlobCutting to 4 wires; runs are identical")
{
    setup_components();
    auto cl = make_cluster(60);
    auto blobs = blobs_of(cl->graph());
    REQUIRE(blobs.size() > 0);

    // reference: every uncut blob cut directly to 4 wires
    std::multiset<bkey> ref;
    {
        auto icfg = Factory::lookup_tn<IConfigurable>("BlobCutting:cdgref");
        auto cfg = icfg->default_configuration();
        cfg["length_threshold"] = 4;     // max_depth stays at BlobCutting's default 10: the doc 09 tier
        icfg->configure(cfg);
        auto cutter = Factory::find_tn<IFunctionNode<IBlobSet, IBlobSet>>("BlobCutting:cdgref");
        for (const auto& b : blobs) {
            IBlobSet::pointer bs = std::make_shared<Aux::SimpleBlobSet>(0, b->slice(), IBlob::vector{b});
            IBlobSet::pointer out;
            (*cutter)(bs, out);
            for (const auto& c : out->blobs()) ref.insert(key_of(c));
        }
    }
    REQUIRE(ref.size() > blobs.size());
    // the depth cap binds for some cells (still wider than 4 wires): the case a per-level depth restart breaks
    size_t nwide = 0;
    for (const auto& k : ref) {
        if (std::get<2>(k) - std::get<1>(k) > 4 || std::get<4>(k) - std::get<3>(k) > 4 || std::get<6>(k) - std::get<5>(k) > 4) ++nwide;
    }
    CHECK(nwide > 0);

    g_mode = 0;
    auto run = [&](Configuration cfg) {
        Img::CascadeDeghosting cd;
        cd.configure(cfg);
        ICluster::pointer out;
        REQUIRE(cd(CascadeDeghosting::input_tuple_type{cl, nullptr}, out));   // no frame: activity charge
        REQUIRE(out);
        return out;
    };
    auto cfg = cascade_config();
    cfg["repair"] = false;
    auto o1 = run(cfg);
    auto o2 = run(cfg);
    auto b1 = blobs_of(o1->graph());
    auto b2 = blobs_of(o2->graph());
    std::multiset<bkey> got;
    for (const auto& b : b1) got.insert(key_of(b));
    CHECK(got == ref);
    REQUIRE(b1.size() == b2.size());
    for (size_t i = 0; i < b1.size(); ++i) {
        CHECK(b1[i]->ident() == b2[i]->ident());
        CHECK(key_of(b1[i]) == key_of(b2[i]));
    }
    CHECK(boost::num_edges(o1->graph()) == boost::num_edges(o2->graph()));
    CHECK(boost::num_vertices(o1->graph()) == boost::num_vertices(o2->graph()));
    CHECK(b1.front()->ident() == (1 << 20));

    // a pruning forward keeps fewer blobs, and never more than the cut-everything reference
    g_mode = 1;
    auto cfg2 = cascade_config();
    auto o3 = run(cfg2);
    const size_t n3 = blobs_of(o3->graph()).size();
    CHECK(n3 <= ref.size());
    g_mode = 0;
}
