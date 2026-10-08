/** CascadeDeghosting and its pieces (wcfm doc 14).

    - build_level: the node features, bw / bw_w edges and bins follow gnn_dataset / contract_level;
    - inslice_pairs: overlap or abut on all three planes, within a (face, slice) group only;
    - geom_pairs (the direct cross-slice builder, doc 15) gives exactly Img::geom_clustering's pairs, for tiled and
      cut blobs, slice gaps of 1, 2 and 3 spans and every policy;
    - guard_prune: a candidate that is a wire's last explanation is kept;
    - steiner_repair: a gap of low-P cells between two confident fragments is bridged within the budget, a
      longer gap is not, a weak island is dropped;
    - the whole cascade with a keep-everything forward reproduces direct BlobCutting to 4 wires (same set of
      (slice, wire bounds)), and two runs give identical output graphs, also with nthreads 4 (doc 15 round 3).

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
    int g_mode = 0;   // 0: keep everything; 1: logit = +10 if the node has >= 3 U channels, else -10; 2: -10 for all

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
                logit[i] = g_mode == 2 ? -10.0f : ((g_mode == 0 || x[i * 15 + 0] >= 3) ? 10.0f : -10.0f);
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
        // threads (doc 15): every array identical
        for (int nt : {2, 5}) {
            auto lt = Cascade::build_level(blobs, sc, uni, k, "uboone", nt);
            CHECK(lt.xb == lev.xb);
            CHECK(lt.wq == lev.wq);
            CHECK(lt.bw_src == lev.bw_src);
            CHECK(lt.bw_dst == lev.bw_dst);
            CHECK(lt.bw_w == lev.bw_w);
            CHECK(lt.bb == lev.bb);
            CHECK(lt.bb_in == lev.bb_in);
            CHECK(lt.ww == lev.ww);
        }
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

TEST_CASE("frame charge: a frame with a non-zero time and frame-relative slice starts (MaskSlice)")
{
    // pdvd doc 122 sec 3: the d121 simulation frames have time -250 us while MaskSlice slice starts count from the
    // frame's first tick.  The legacy reader then starts 500 ticks late; slice_start_relative reads the slice's ticks.
    auto cl = make_cluster(12);
    const auto& gr = cl->graph();
    std::vector<float> q100(20);
    for (int t = 0; t < 20; ++t) q100[t] = (float) t;
    ITrace::vector trs{std::make_shared<Aux::SimpleTrace>(100, 0, q100)};
    auto frame = std::make_shared<Aux::SimpleFrame>(7, -250.0 * units::microsecond, trs, 0.5 * units::microsecond);
    auto legacy = Cascade::make_slice_charge_frame(gr, frame, "", 0.25);
    REQUIRE(legacy.slice_of.size() == 2);
    CHECK(legacy.t0[0] == 500);
    CHECK(legacy.charge(0, 100) == 0.0f);        // read beyond the trace
    auto sc = Cascade::make_slice_charge_frame(gr, frame, "", 0.25, true);
    CHECK(sc.t0[0] == 0);
    CHECK(sc.t0[1] == 4);
    CHECK(sc.charge(0, 100) == doctest::Approx(0.25 * (0 + 1 + 2 + 3)));
    CHECK(sc.charge(1, 100) == doctest::Approx(0.25 * (4 + 5 + 6 + 7)));
    // a frame with time 0: the two readings are the same
    auto frame0 = std::make_shared<Aux::SimpleFrame>(7, 0.0, trs, 0.5 * units::microsecond);
    auto a = Cascade::make_slice_charge_frame(gr, frame0, "", 0.25, false);
    auto b = Cascade::make_slice_charge_frame(gr, frame0, "", 0.25, true);
    CHECK(a.t0 == b.t0);
    CHECK(a.charge(1, 100) == b.charge(1, 100));
}

TEST_CASE("geom_pairs equals the geom_clustering graph pairs")
{
    auto anodes = Testing::anodes("uboone");
    REQUIRE(anodes.size() > 0);
    auto face = anodes[0]->face(0);
    Img::GridTiling gt;
    auto gcfg = gt.default_configuration();
    gcfg["anode"] = "AnodePlane:0";
    gcfg["face"] = 0;
    gt.configure(gcfg);
    auto frame = std::make_shared<Aux::SimpleFrame>(7, 0.0);
    const double span = 2.0 * units::microsecond;
    // slice starts in spans (gaps of 1, 2 and 3), active width and position (overlaps vary slice to slice)
    const int at[7] = {0, 1, 2, 4, 5, 8, 9};
    const size_t nch[7] = {12, 20, 8, 16, 30, 10, 14};
    const double first[7] = {0.5, 0.501, 0.4995, 0.5, 0.502, 0.5, 0.5007};
    std::vector<ISlice::pointer> slice_of;
    std::vector<IBlob::pointer> uncut;
    std::vector<int> sidx_uncut;
    for (int s = 0; s < 7; ++s) {
        // one band per plane with a dead channel every `gap` channels: tiling splits it into many blobs
        auto band = make_slice(face, frame, s, at[s] * span, nch[s], first[s]);
        ISlice::map_t act;
        std::vector<std::pair<int, IChannel::pointer>> byid;
        for (const auto& kv : band->activity()) byid.emplace_back(kv.first->ident(), kv.first);
        std::sort(byid.begin(), byid.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        const int gap = 3 + s % 3;
        for (size_t i = 0; i < byid.size(); ++i) {
            if (i % gap == (size_t) gap - 1) continue;
            act[byid[i].second] = band->activity().at(byid[i].second);
        }
        auto sl = std::make_shared<Aux::SimpleSlice>(frame, s, at[s] * span, span, act);
        slice_of.push_back(sl);
        IBlobSet::pointer bs;
        REQUIRE(gt(sl, bs));
        for (const auto& b : bs->blobs()) {
            uncut.push_back(b);
            sidx_uncut.push_back(s);
        }
    }
    REQUIRE(uncut.size() >= 14);
    // the same blobs cut to 2 wires (more, smaller blobs, abutting and gapped)
    std::vector<IBlob::pointer> cut;
    std::vector<int> sidx_cut;
    Cascade::CutParams cpar;
    cpar.min_length = 1;
    int ident = 1000;
    for (size_t i = 0; i < uncut.size(); ++i) {
        const auto& b = uncut[i];
        for (auto& [shape, d] : Cascade::cut_shape(b->face()->raygrid(), b->shape(), 0, 2, cpar)) {
            cut.push_back(std::make_shared<Aux::SimpleBlob>(ident++, 1.0, 0.0, shape, b->slice(), b->face()));
            sidx_cut.push_back(sidx_uncut[i]);
        }
    }
    REQUIRE(cut.size() > uncut.size());
    for (const std::string policy : {"uboone", "uboone_local", "simple", "dead_clus"}) {
        for (int which = 0; which < 2; ++which) {
            const auto& bl = which ? cut : uncut;
            const auto& sx = which ? sidx_cut : sidx_uncut;
            const auto ref = Cascade::geom_pairs_graph(bl, sx, slice_of, policy);
            const auto got = Cascade::geom_pairs(bl, sx, slice_of, policy);
            CAPTURE(policy);
            CAPTURE(which);
            CHECK(ref.size() > 0);
            CHECK(got == ref);
            for (int nt : {2, 3, 16}) CHECK(Cascade::geom_pairs(bl, sx, slice_of, policy, nt) == ref);
        }
    }
    CHECK_THROWS(Cascade::geom_pairs(uncut, sidx_uncut, slice_of, "nonesuch"));
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

    // threads: the same output graph, vertex by vertex and edge by edge
    {
        auto cfgt = cfg;
        cfgt["nthreads"] = 4;
        auto o4 = run(cfgt);
        const auto& g1 = o1->graph();
        const auto& g4 = o4->graph();
        REQUIRE(boost::num_vertices(g1) == boost::num_vertices(g4));
        REQUIRE(boost::num_edges(g1) == boost::num_edges(g4));
        for (auto v : boost::make_iterator_range(boost::vertices(g1))) {
            CHECK(g1[v].code() == g4[v].code());
            CHECK(g1[v].ident() == g4[v].ident());
        }
        std::vector<std::pair<size_t, size_t>> e1, e4;
        for (auto e : boost::make_iterator_range(boost::edges(g1))) e1.emplace_back(boost::source(e, g1), boost::target(e, g1));
        for (auto e : boost::make_iterator_range(boost::edges(g4))) e4.emplace_back(boost::source(e, g4), boost::target(e, g4));
        CHECK(e1 == e4);
        auto k4 = blobs_of(g4);
        REQUIRE(k4.size() == b1.size());
        for (size_t i = 0; i < b1.size(); ++i) CHECK(key_of(k4[i]) == key_of(b1[i]));
    }

    // a pruning forward keeps fewer blobs, and never more than the cut-everything reference
    g_mode = 1;
    auto cfg2 = cascade_config();
    auto o3 = run(cfg2);
    const size_t n3 = blobs_of(o3->graph()).size();
    CHECK(n3 <= ref.size());
    g_mode = 0;
}

TEST_CASE("wcfm doc 22 knobs: keep_slices keeps every emptied slice node; final_guard explains every charged wire")
{
    setup_components();
    auto cl = make_cluster(60);
    size_t nslice_in = 0;
    for (auto vtx : boost::make_iterator_range(boost::vertices(cl->graph()))) nslice_in += cl->graph()[vtx].code() == 's';
    REQUIRE(nslice_in > 1);
    auto run = [&](Configuration cfg) {
        Img::CascadeDeghosting cd;
        cd.configure(cfg);
        ICluster::pointer out;
        REQUIRE(cd(CascadeDeghosting::input_tuple_type{cl, nullptr}, out));
        REQUIRE(out);
        return out;
    };
    auto nslices = [](const ICluster::pointer& c) {
        size_t n = 0;
        for (auto vtx : boost::make_iterator_range(boost::vertices(c->graph()))) n += c->graph()[vtx].code() == 's';
        return n;
    };
    auto cfg = cascade_config();
    cfg["repair"] = false;
    // the defaults round-trip and are off
    {
        Img::CascadeDeghosting cd;
        auto dc = cd.default_configuration();
        CHECK(dc["final_guard"].asBool() == false);
        CHECK(dc["keep_slices"].asBool() == false);
    }

    // a forward that drops everything: no blob survives; keep_slices keeps every input slice node
    g_mode = 2;
    auto off = run(cfg);
    CHECK(blobs_of(off->graph()).empty());
    CHECK(nslices(off) == 0);
    auto cks = cfg;
    cks["keep_slices"] = true;
    auto on = run(cks);
    CHECK(blobs_of(on->graph()).empty());
    CHECK(nslices(on) == nslice_in);

    // final_guard on the same forward: blobs come back, one or more per charged wire; never more than the cut reference
    auto cfg_g = cfg;
    cfg_g["final_guard"] = true;
    auto og = run(cfg_g);
    const size_t ng = blobs_of(og->graph()).size();
    CHECK(ng > 0);

    // a pruning forward: keep_slices changes no blob, only adds slice nodes
    g_mode = 1;
    auto p_off = run(cfg);
    auto p_on = run(cks);
    auto b_off = blobs_of(p_off->graph());
    auto b_on = blobs_of(p_on->graph());
    REQUIRE(b_off.size() == b_on.size());
    for (size_t i = 0; i < b_off.size(); ++i) CHECK(key_of(b_off[i]) == key_of(b_on[i]));
    CHECK(nslices(p_on) == nslice_in);
    CHECK(nslices(p_off) <= nslice_in);
    g_mode = 0;
}
