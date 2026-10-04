/** CascadeDeghostingFM (wcfm doc 30): CascadeDeghosting + the head-score columns at the final level.

    - head off: the node's output equals CascadeDeghosting's on the same cluster (same kept cells);
    - head on: the head forward receives, per final-level node, has_p = 1 iff the cell has a channel with an FM pixel
      in plane p and e_p = the mean of those pixels' features (half-rounded); the final level's forward receives 18
      columns = the 15 charge columns + [score, ls_mean, ls_min] with score = the head's output, ls_mean >= ls_min,
      both <= 0 (log-softmax), and every active wire's candidates' exp(ls) sum to 1;
    - head on without an FM set on port 2: no active pixel, every has = 0, the ls columns take the fill value.
    The forwards are in-test ITensorForwards (no model file, no torch); the fixtures duplicate doctest_cascade_deghosting's.
*/
#include "WireCellImg/CascadeDeghosting.h"
#include "WireCellImg/CascadeDeghostingFM.h"
#include "WireCellImg/CascadeGraph.h"
#include "WireCellImg/GridTiling.h"
#include "WireCellImg/BlobClustering.h"

#include "WireCellAux/SimpleBlob.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleSlice.h"
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
#include <map>
#include <set>
#include <tuple>

using namespace WireCell;
using namespace WireCell::Img;

namespace {
    const int FMDIM = 4;
    // what the fakes saw in their last call
    std::vector<float> g_xb;          // final-level xb (N x ncol)
    size_t g_ncol = 0;
    std::vector<float> g_eu, g_ev, g_ew, g_has;   // the head's inputs
    std::vector<float> g_score;

    float fake_score(const float* has) { return 10.0f * has[0] + 5.0f * has[1] + has[2] - 3.0f; }

    struct FMFakeLevelForward : public ITensorForward {   // keep everything; records xb
        virtual ~FMFakeLevelForward() {}
        virtual ITensorSet::pointer forward(const ITensorSet::pointer& in) const
        {
            auto tv = in->tensors();
            REQUIRE(tv->size() == 9);
            const auto& xb = tv->at(0);
            REQUIRE(xb->dtype() == "f4");
            REQUIRE(xb->shape().size() == 2);
            const size_t N = xb->shape()[0];
            g_ncol = xb->shape()[1];
            const float* x = reinterpret_cast<const float*>(xb->data());
            g_xb.assign(x, x + N * g_ncol);
            std::vector<float> logit(N, 10.0f), qhat(N, 1.0f);
            ITensor::vector out{std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{N}, logit.data()),
                                std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{N}, qhat.data())};
            return std::make_shared<Aux::SimpleTensorSet>(in->ident(), Configuration(), std::make_shared<ITensor::vector>(out));
        }
    };
    struct FMFakeHeadForward : public ITensorForward {
        virtual ~FMFakeHeadForward() {}
        virtual ITensorSet::pointer forward(const ITensorSet::pointer& in) const
        {
            auto tv = in->tensors();
            REQUIRE(tv->size() == 4);
            for (int k = 0; k < 4; ++k) REQUIRE(tv->at(k)->dtype() == "f4");
            const size_t N = tv->at(3)->shape()[0];
            REQUIRE(tv->at(3)->shape()[1] == 3);
            REQUIRE(tv->at(0)->shape()[1] == (size_t) FMDIM);
            auto grab = [&](int k, std::vector<float>& v) {
                const float* p = reinterpret_cast<const float*>(tv->at(k)->data());
                v.assign(p, p + tv->at(k)->size() / sizeof(float));
            };
            grab(0, g_eu); grab(1, g_ev); grab(2, g_ew); grab(3, g_has);
            g_score.resize(N);
            for (size_t i = 0; i < N; ++i) g_score[i] = fake_score(&g_has[i * 3]);
            ITensor::vector out{std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{N}, g_score.data())};
            return std::make_shared<Aux::SimpleTensorSet>(in->ident(), Configuration(), std::make_shared<ITensor::vector>(out));
        }
    };
}  // namespace

WIRECELL_FACTORY(FMFakeLevelForward, FMFakeLevelForward, WireCell::ITensorForward)
WIRECELL_FACTORY(FMFakeHeadForward, FMFakeHeadForward, WireCell::ITensorForward)

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
    std::set<bkey> keys_of(const ICluster::pointer& cl)
    {
        std::set<bkey> out;
        for (const auto& b : blobs_of(cl->graph())) out.insert(key_of(b));
        return out;
    }

    void setup_components()
    {
        Testing::load_plugins({"WireCellImg"});
        make_FMFakeLevelForward_factory();
        make_FMFakeHeadForward_factory();
        Factory::lookup_tn<ITensorForward>("FMFakeLevelForward");
        Factory::lookup_tn<ITensorForward>("FMFakeHeadForward");
    }

    Configuration cascade_config(bool head)
    {
        Configuration cfg;
        const int ks[5] = {16, 8, 4, 2, 1};
        const int ws[5] = {0, 32, 16, 8, 4};
        for (int i = 0; i < 5; ++i) {
            Configuration l;
            l["width"] = ws[i];
            l["superwire"] = ks[i];
            l["forward"] = "FMFakeLevelForward";
            l["threshold"] = i < 4 ? -3.0 : -0.2;
            cfg["levels"].append(l);
        }
        cfg["charge_scale"] = 0.25;
        if (head) {
            cfg["head"]["forward"] = "FMFakeHeadForward";
            cfg["head"]["fm_dim"] = FMDIM;
        }
        return cfg;
    }

    // feature value of a pixel: channel ident / 1000, the same in every dimension
    float pix_value(int ch) { return ch / 1000.0f; }

    // an FM tensor set: one pixel per (plane, channel, slice ident) of the cluster's blobs, planes U and V only
    // when `uv_only` (so plane W has no pixel anywhere -> has_W = 0 for every cell)
    ITensorSet::pointer make_fm(const ICluster::pointer& cl, bool uv_only)
    {
        std::array<std::set<std::pair<int, int>>, 3> pix;
        for (const auto& b : blobs_of(cl->graph())) {
            const auto ch = Cascade::blob_channels(b);
            for (int p = 0; p < (uv_only ? 2 : 3); ++p) {
                for (int c : ch[p]) pix[p].insert({c, b->slice()->ident()});
            }
        }
        auto itv = std::make_shared<ITensor::vector>();
        for (int p = 0; p < 3; ++p) {
            std::vector<int32_t> coords;
            std::vector<float> feat;
            for (const auto& [c, k] : pix[p]) {
                coords.push_back(c);
                coords.push_back(k);
                for (int d = 0; d < FMDIM; ++d) feat.push_back(pix_value(c) * (d == 0 ? 1.0f : 0.5f));
            }
            const size_t N = pix[p].size();
            Configuration cmd;
            cmd["name"] = "coords";
            cmd["plane"] = p;
            itv->push_back(std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{N, 2}, coords.data(), cmd));
            Configuration fmd;
            fmd["name"] = "feat";
            fmd["plane"] = p;
            itv->push_back(std::make_shared<Aux::SimpleTensor>(ITensor::shape_t{N, (size_t) FMDIM}, feat.data(), fmd));
        }
        return std::make_shared<Aux::SimpleTensorSet>(7, Configuration(), itv);
    }
}  // namespace

TEST_CASE("CascadeDeghostingFM with the head off equals CascadeDeghosting")
{
    setup_components();
    auto cl = make_cluster(24);
    Img::CascadeDeghosting ref;
    ref.configure(cascade_config(false));
    Img::CascadeDeghostingFM fm;
    fm.configure(cascade_config(false));
    ICluster::pointer o1, o2;
    REQUIRE(ref(std::make_tuple(cl, IFrame::pointer(nullptr)), o1));
    REQUIRE(fm(std::make_tuple(cl, IFrame::pointer(nullptr), ITensorSet::pointer(nullptr)), o2));
    REQUIRE(o1);
    REQUIRE(o2);
    CHECK(g_ncol == 15);
    CHECK(keys_of(o1) == keys_of(o2));
    CHECK(boost::num_edges(o1->graph()) == boost::num_edges(o2->graph()));
}

TEST_CASE("CascadeDeghostingFM head columns: has, per-view means, score, per-wire log-softmax")
{
    setup_components();
    auto cl = make_cluster(24);
    auto fmset = make_fm(cl, true);
    Img::CascadeDeghostingFM fm;
    fm.configure(cascade_config(true));
    ICluster::pointer out;
    REQUIRE(fm(std::make_tuple(cl, IFrame::pointer(nullptr), fmset), out));
    REQUIRE(out);
    REQUIRE(g_ncol == 18);
    const size_t N = g_xb.size() / 18;
    REQUIRE(N > 0);
    REQUIRE(g_has.size() == N * 3);
    // the final-level cells are the 4-wire cut of the input; every cell has channels in all three planes,
    // the FM set has pixels in U and V only
    size_t nu = 0, nv = 0, nw = 0;
    for (size_t i = 0; i < N; ++i) {
        nu += g_has[i * 3] > 0;
        nv += g_has[i * 3 + 1] > 0;
        nw += g_has[i * 3 + 2] > 0;
        // has_p == (the cell has >= 1 channel with charge in plane p) for U, V; W has no pixel
        CHECK(g_has[i * 3] == ((g_xb[i * 18 + 1] > 0) ? 1.0f : 0.0f));
        CHECK(g_has[i * 3 + 1] == ((g_xb[i * 18 + 5] > 0) ? 1.0f : 0.0f));
        CHECK(g_has[i * 3 + 2] == 0.0f);
        for (int d = 0; d < FMDIM; ++d) CHECK(g_ew[i * FMDIM + d] == 0.0f);
        // e_U dimension 1 = half of dimension 0 (the pattern), e_U in the channel range of the plane / 1000
        if (g_has[i * 3] > 0) {
            CHECK(g_eu[i * FMDIM + 1] == doctest::Approx(0.5 * g_eu[i * FMDIM]).epsilon(2e-3));
            CHECK(g_eu[i * FMDIM] > 0.0f);
        }
        // the three columns
        CHECK(g_xb[i * 18 + 15] == doctest::Approx(fake_score(&g_has[i * 3])));
        CHECK(g_xb[i * 18 + 16] <= 1e-6f);
        CHECK(g_xb[i * 18 + 17] <= 1e-6f);
        CHECK(g_xb[i * 18 + 16] >= g_xb[i * 18 + 17] - 1e-6f);
    }
    CHECK(nu == N);
    CHECK(nv == N);
    CHECK(nw == 0);
    // every cell has every pixel: the score is the same for all, so ls = -log(number of cells on the wire):
    // the strictest cell (ls_min) sits on the busiest of its wires, and exp(ls_min) * n_cells_on_that_wire = 1.
    // Rebuild the final level to count cells per wire.
    const auto& gr = cl->graph();
    auto sc = Cascade::make_slice_charge(gr, 0.25, 1e11);
    // the output cells are the final-level cells (keep-everything forward); build the k = 1 level on them
    auto cells = blobs_of(out->graph());
    REQUIRE(cells.size() == N);
    auto uni = Cascade::make_universe(blobs_of(gr), sc);
    auto lev = Cascade::build_level(cells, sc, uni, 1, "uboone");
    std::map<int64_t, int> ncell_on_wire;
    for (size_t e = 0; e < lev.bw_src.size(); ++e) ncell_on_wire[lev.bw_dst[e]] += 1;
    // the output cells are renumbered in the order of the final level's kept list, which is node order
    for (size_t i = 0; i < N; ++i) {
        int busiest = 0;
        double tot = 0;
        int nact = 0;
        for (size_t e = 0; e < lev.bw_src.size(); ++e) {
            if ((size_t) lev.bw_src[e] != i) continue;
            if (lev.wplane[lev.bw_dst[e]] == 2) continue;   // W wires are not active (no pixel)
            const int n = ncell_on_wire[lev.bw_dst[e]];
            busiest = std::max(busiest, n);
            tot += -std::log((double) n);
            ++nact;
        }
        REQUIRE(nact > 0);
        CHECK(g_xb[i * 18 + 17] == doctest::Approx(-std::log((double) busiest)).epsilon(1e-4));
        CHECK(g_xb[i * 18 + 16] == doctest::Approx(tot / nact).epsilon(1e-4));
    }
}

TEST_CASE("CascadeDeghostingFM head on without an FM set: no active pixel, filled columns")
{
    setup_components();
    auto cl = make_cluster(12);
    Img::CascadeDeghostingFM fm;
    fm.configure(cascade_config(true));
    ICluster::pointer out;
    REQUIRE(fm(std::make_tuple(cl, IFrame::pointer(nullptr), ITensorSet::pointer(nullptr)), out));
    REQUIRE(out);
    REQUIRE(g_ncol == 18);
    const size_t N = g_xb.size() / 18;
    for (size_t i = 0; i < N; ++i) {
        CHECK(g_has[i * 3] == 0.0f);
        CHECK(g_has[i * 3 + 1] == 0.0f);
        CHECK(g_has[i * 3 + 2] == 0.0f);
        CHECK(g_xb[i * 18 + 15] == doctest::Approx(-3.0f));
        CHECK(g_xb[i * 18 + 16] == 0.0f);
        CHECK(g_xb[i * 18 + 17] == 0.0f);
    }
}
