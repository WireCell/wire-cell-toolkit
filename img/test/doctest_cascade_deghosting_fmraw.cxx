/** CascadeDeghostingFMRaw (wcfm doc 53 part 43, pr_ml doc 30 sec 9): CascadeDeghostingFM with the head call replaced by
    a PCA projection of the per-view FM means at the final level.  Fork by duplication of doctest_cascade_deghosting_fm.cxx.
    - fmraw off: the node's output equals CascadeDeghosting's on the same cluster (same kept cells);
    - fmraw on, identity basis: the final-level xb = 15 charge columns + [x_U (k), x_V (k), x_W (k), has_U, has_V,
      has_W] with x_p = e_p (the head stage's per-view mean: dimension 1 = half of dimension 0), x_W = 0 where W has no
      pixel, has_p as CascadeDeghostingFM;
    - fmraw on, a general mean and basis: x_p = ((e_p - mean) @ basis) * has_p against e_p from the identity run;
    - fmraw on without an FM set on port 2: every has and every x = 0.
 */
#include "WireCellImg/CascadeDeghosting.h"
#include "WireCellImg/CascadeDeghostingFMRaw.h"
#include "WireCellUtil/Persist.h"
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
#include <fstream>
#include <tuple>

using namespace WireCell;
using namespace WireCell::Img;

namespace {
    const int FMDIM = 4;
    // what the fakes saw in their last call
    std::vector<float> g_xb;          // final-level xb (N x ncol)
    size_t g_ncol = 0;

    struct RawFakeLevelForward : public ITensorForward {   // keep everything; records xb
        virtual ~RawFakeLevelForward() {}
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
}  // namespace

WIRECELL_FACTORY(RawFakeLevelForward, RawFakeLevelForward, WireCell::ITensorForward)

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
        make_RawFakeLevelForward_factory();
        Factory::lookup_tn<ITensorForward>("RawFakeLevelForward");
    }

    Configuration cascade_config(const std::string& pca)
    {
        Configuration cfg;
        const int ks[5] = {16, 8, 4, 2, 1};
        const int ws[5] = {0, 32, 16, 8, 4};
        for (int i = 0; i < 5; ++i) {
            Configuration l;
            l["width"] = ws[i];
            l["superwire"] = ks[i];
            l["forward"] = "RawFakeLevelForward";
            l["threshold"] = i < 4 ? -3.0 : -0.2;
            cfg["levels"].append(l);
        }
        cfg["charge_scale"] = 0.25;
        if (!pca.empty()) {
            cfg["fmraw"]["pca"] = pca;
            cfg["fmraw"]["fm_dim"] = FMDIM;
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


namespace {
    // a pca json {mean[dim], basis[dim][k], k, dim} in the working directory
    std::string write_pca(const std::string& name, const std::vector<float>& mean, const std::vector<std::vector<float>>& basis)
    {
        Configuration j;
        j["dim"] = FMDIM;
        j["k"] = (int) basis[0].size();
        for (float m : mean) j["mean"].append(m);
        for (const auto& row : basis) {
            Configuration r;
            for (float v : row) r.append(v);
            j["basis"].append(r);
        }
        Persist::dump(name, j);
        return name;
    }
    std::string identity_pca()
    {
        std::vector<std::vector<float>> b(FMDIM, std::vector<float>(FMDIM, 0.0f));
        for (int d = 0; d < FMDIM; ++d) b[d][d] = 1.0f;
        return write_pca("doctest_fmraw_pca_identity.json", std::vector<float>(FMDIM, 0.0f), b);
    }
}  // namespace

TEST_CASE("CascadeDeghostingFMRaw with fmraw off equals CascadeDeghosting")
{
    setup_components();
    auto cl = make_cluster(24);
    Img::CascadeDeghosting ref;
    ref.configure(cascade_config(""));
    Img::CascadeDeghostingFMRaw raw;
    raw.configure(cascade_config(""));
    ICluster::pointer o1, o2;
    REQUIRE(ref(std::make_tuple(cl, IFrame::pointer(nullptr)), o1));
    REQUIRE(raw(std::make_tuple(cl, IFrame::pointer(nullptr), ITensorSet::pointer(nullptr)), o2));
    REQUIRE(o1);
    REQUIRE(o2);
    CHECK(g_ncol == 15);
    CHECK(keys_of(o1) == keys_of(o2));
    CHECK(boost::num_edges(o1->graph()) == boost::num_edges(o2->graph()));
}

TEST_CASE("CascadeDeghostingFMRaw columns: identity basis gives the per-view means, has flags, zeros without pixels")
{
    setup_components();
    auto cl = make_cluster(24);
    auto fmset = make_fm(cl, true);
    Img::CascadeDeghostingFMRaw raw;
    raw.configure(cascade_config(identity_pca()));
    ICluster::pointer out;
    REQUIRE(raw(std::make_tuple(cl, IFrame::pointer(nullptr), fmset), out));
    REQUIRE(out);
    const size_t K = FMDIM, NC = 15 + 3 * K + 3;
    REQUIRE(g_ncol == NC);
    const size_t N = g_xb.size() / NC;
    REQUIRE(N > 0);
    for (size_t i = 0; i < N; ++i) {
        const float* x = &g_xb[i * NC];
        const float hu = x[15 + 3 * K], hv = x[15 + 3 * K + 1], hw = x[15 + 3 * K + 2];
        // has_p as CascadeDeghostingFM: U and V pixels for every charged channel, no W pixel
        CHECK(hu == ((x[1] > 0) ? 1.0f : 0.0f));
        CHECK(hv == ((x[5] > 0) ? 1.0f : 0.0f));
        CHECK(hw == 0.0f);
        CHECK(hu == 1.0f);
        for (size_t d = 0; d < K; ++d) CHECK(x[15 + 2 * K + d] == 0.0f);   // x_W
        // identity basis, zero mean: x_U = e_U, whose dimension 1 is half of dimension 0 (the make_fm pattern)
        CHECK(x[15] > 0.0f);
        CHECK(x[15 + 1] == doctest::Approx(0.5 * x[15]).epsilon(2e-3));
        CHECK(x[15 + K + 1] == doctest::Approx(0.5 * x[15 + K]).epsilon(2e-3));
    }
}

TEST_CASE("CascadeDeghostingFMRaw columns: x_p = ((e_p - mean) @ basis) * has_p")
{
    setup_components();
    auto cl = make_cluster(24);
    auto fmset = make_fm(cl, true);
    // e_p from the identity run
    Img::CascadeDeghostingFMRaw id;
    id.configure(cascade_config(identity_pca()));
    ICluster::pointer o1;
    REQUIRE(id(std::make_tuple(cl, IFrame::pointer(nullptr), fmset), o1));
    const size_t D = FMDIM, NC1 = 15 + 3 * D + 3;
    const std::vector<float> e = g_xb;
    const size_t N = e.size() / NC1;
    // a 2-component basis and a non-zero mean
    const std::vector<float> mean{0.1f, -0.2f, 0.05f, 0.3f};
    const std::vector<std::vector<float>> basis{{1.0f, 0.5f}, {-2.0f, 0.0f}, {0.25f, 1.0f}, {0.0f, -1.5f}};
    Img::CascadeDeghostingFMRaw raw;
    raw.configure(cascade_config(write_pca("doctest_fmraw_pca_k2.json", mean, basis)));
    ICluster::pointer o2;
    REQUIRE(raw(std::make_tuple(cl, IFrame::pointer(nullptr), fmset), o2));
    const size_t K = 2, NC = 15 + 3 * K + 3;
    REQUIRE(g_ncol == NC);
    REQUIRE(g_xb.size() == N * NC);
    for (size_t i = 0; i < N; ++i) {
        for (size_t c = 0; c < 15; ++c) CHECK(g_xb[i * NC + c] == e[i * NC1 + c]);
        for (size_t p = 0; p < 3; ++p) {
            const float has = e[i * NC1 + 15 + 3 * D + p];
            CHECK(g_xb[i * NC + 15 + 3 * K + p] == has);
            for (size_t j = 0; j < K; ++j) {
                double want = 0;
                for (size_t d = 0; d < D; ++d) want += (e[i * NC1 + 15 + p * D + d] - mean[d]) * basis[d][j];
                CHECK(g_xb[i * NC + 15 + p * K + j] == doctest::Approx(has * want).epsilon(1e-5));
            }
        }
    }
}

TEST_CASE("CascadeDeghostingFMRaw on without an FM set: no active pixel, zero columns")
{
    setup_components();
    auto cl = make_cluster(12);
    Img::CascadeDeghostingFMRaw raw;
    raw.configure(cascade_config(identity_pca()));
    ICluster::pointer out;
    REQUIRE(raw(std::make_tuple(cl, IFrame::pointer(nullptr), ITensorSet::pointer(nullptr)), out));
    REQUIRE(out);
    const size_t K = FMDIM, NC = 15 + 3 * K + 3;
    REQUIRE(g_ncol == NC);
    const size_t N = g_xb.size() / NC;
    for (size_t i = 0; i < N; ++i) {
        for (size_t c = 15; c < NC; ++c) CHECK(g_xb[i * NC + c] == 0.0f);
    }
}
