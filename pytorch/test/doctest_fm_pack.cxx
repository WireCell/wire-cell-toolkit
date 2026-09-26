/** FMFeatureExtract packing test (wcfm doc 04, doc 01 sec 7 F3): rows, tick sum x scale,
    log map, bbox floor, gather alignment, halo tiling seamlessness, half storage, EOS.

    The forward is a POINTWISE fake (feature c of a pixel is a function of that pixel's two
    input values only), so tiling must reproduce the untiled features exactly and the gather
    alignment is checked through the feature values themselves.  No model file, no torch call.

    Revert-proven: with the tick sum replaced by the first tick, the q check fails; with the
    gather row/col swapped, the feature-0 check fails; with the halo columns gathered, the
    tiled-vs-untiled equality fails.
*/
#include "WireCellPytorch/FMFeatureExtract.h"

#include "WireCellIface/ITensorForward.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/Testing.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/doctest.h"

#include <cmath>
#include <cstdint>
#include <map>

using namespace WireCell;

namespace {
    const int FAKE_C = 6;

    // out[c] = x0 (c=0), x1 (c=1), (c+1)*x0 + 0.5*x1 (c>=2): pointwise in the pixel.
    class FMFakeForward : public ITensorForward {
      public:
        virtual ~FMFakeForward() {}
        virtual ITensorSet::pointer forward(const ITensorSet::pointer& in) const
        {
            auto iten = in->tensors()->front();
            const auto shape = iten->shape();
            REQUIRE(shape.size() == 4);
            REQUIRE(shape[0] == 1);
            REQUIRE(shape[1] == 2);
            REQUIRE(shape[2] >= 64);
            REQUIRE(shape[3] >= 64);
            const size_t h = shape[2], w = shape[3];
            const float* x = reinterpret_cast<const float*>(iten->data());
            std::vector<float> y((size_t) FAKE_C * h * w, 0.0f);
            for (size_t p = 0; p < h * w; ++p) {
                const float x0 = x[p], x1 = x[h * w + p];
                y[p] = x0;
                y[h * w + p] = x1;
                for (int c = 2; c < FAKE_C; ++c) {
                    y[(size_t) c * h * w + p] = (c + 1) * x0 + 0.5f * x1;
                }
            }
            ITensor::shape_t oshape = {1, (size_t) FAKE_C, h, w};
            auto oten = std::make_shared<Aux::SimpleTensor>(oshape, y.data());
            return std::make_shared<Aux::SimpleTensorSet>(in->ident(), Configuration(),
                                                          std::make_shared<ITensor::vector>(ITensor::vector{oten}));
        }
    };
}
WIRECELL_FACTORY(FMFakeForward, FMFakeForward, WireCell::ITensorForward)

namespace {
    float logmap(float x, double m, double M)
    {
        const float y0 = (float) std::log10(m);
        const float sc = (float) (2.0 / (std::log10(M + m) - std::log10(m)));
        return (std::log10(x + (float) m) - y0) * sc - 1.0f;
    }

    // frame: channel ch gets charges q[0..n) starting at tbin.
    IFrame::pointer make_frame(int ident, const std::vector<std::tuple<int, int, std::vector<float>>>& specs)
    {
        ITrace::vector traces;
        IFrame::trace_list_t idx;
        for (const auto& [ch, tbin, q] : specs) {
            traces.push_back(std::make_shared<Aux::SimpleTrace>(ch, tbin, q));
            idx.push_back(idx.size());
        }
        auto sf = new Aux::SimpleFrame(ident, 0.0, traces, 0.5 * units::microsecond);
        sf->tag_traces("gauss0", idx);
        return IFrame::pointer(sf);
    }

    Configuration base_cfg(Pytorch::FMFeatureExtract& fm)
    {
        auto cfg = fm.default_configuration();
        cfg["anode"] = "AnodePlane:0";
        cfg["forward"] = "FMFakeForward";
        cfg["input_tag"] = "gauss0";
        cfg["planes"] = Json::arrayValue;
        cfg["planes"].append(2);   // uboone Y: channels 4800..8255
        cfg["feature_dim"] = FAKE_C;
        cfg["nticks"] = 2000;
        cfg["store_half"] = false;
        return cfg;
    }

    struct Out {
        std::vector<std::array<int, 2>> coords;
        std::vector<float> feat;
        Configuration md;
        std::string feat_name;
        std::string feat_dtype;
    };

    Out unpack(ITensorSet::pointer set)
    {
        Out o;
        o.md = set->metadata();
        REQUIRE(set->tensors()->size() == 2);
        auto ct = set->tensors()->at(0);
        auto ft = set->tensors()->at(1);
        REQUIRE(ct->metadata()["name"].asString() == "coords");
        REQUIRE(ct->dtype() == "i4");
        const size_t N = ct->shape()[0];
        REQUIRE(ct->shape()[1] == 2);
        const int32_t* cd = reinterpret_cast<const int32_t*>(ct->data());
        for (size_t i = 0; i < N; ++i) o.coords.push_back({cd[2 * i], cd[2 * i + 1]});
        o.feat_name = ft->metadata()["name"].asString();
        o.feat_dtype = ft->dtype();
        REQUIRE(ft->shape()[0] == N);
        REQUIRE(ft->shape()[1] == (size_t) FAKE_C);
        if (o.feat_dtype == "f4") {
            const float* fd = reinterpret_cast<const float*>(ft->data());
            o.feat.assign(fd, fd + N * FAKE_C);
        }
        else {
            REQUIRE(o.feat_dtype == "u2");
            const uint16_t* hd = reinterpret_cast<const uint16_t*>(ft->data());
            for (size_t i = 0; i < N * FAKE_C; ++i) {
                // IEEE half -> float
                const uint16_t b = hd[i];
                const int sign = (b >> 15) & 1, exp = (b >> 10) & 0x1f, man = b & 0x3ff;
                float v;
                if (exp == 0) v = std::ldexp((float) man, -24);
                else v = std::ldexp((float) (man | 0x400), exp - 25);
                o.feat.push_back(sign ? -v : v);
            }
        }
        return o;
    }
}

TEST_CASE("fm pack: rows, tick sum, log map, bbox floor, gather")
{
    auto anodes = Testing::anodes("uboone");
    REQUIRE(anodes.size() > 0);
    make_FMFakeForward_factory();  // the factory lives in this binary, not in a plugin (util doctest precedent)
    Factory::lookup_tn<ITensorForward>("FMFakeForward");

    Pytorch::FMFeatureExtract fm;
    auto cfg = base_cfg(fm);
    CHECK_EQ(cfg["tick_span"].asInt(), 4);
    CHECK_EQ(cfg["input_scale"].asDouble(), 0.25);
    CHECK_EQ(cfg["min_canvas"].asInt(), 64);
    CHECK_EQ(cfg["feature_dim"].asInt(), FAKE_C);
    CHECK_EQ(cfg["view_norm"][2][0].asDouble(), 3.75);
    fm.configure(cfg);

    // channel 5000: ticks 398..403 = [1,2,3,4,5,6] -> slice 99 = 0.25*(1+2) = 0.75, slice 100 = 0.25*(3+4+5+6) = 4.5
    // channel 5001: ticks 400..403 = [4,8,12,16] -> slice 100 = 10
    // channel 5010: one sample at tick 1200 = 100 -> slice 300 = 25
    // channel 4700 (plane V, not packed): ignored
    auto frame = make_frame(7, {{5000, 398, {1, 2, 3, 4, 5, 6}}, {5001, 400, {4, 8, 12, 16}},
                                {5010, 1200, {100}}, {4700, 400, {50, 50, 50, 50}}});
    ITensorSet::pointer set;
    REQUIRE(fm(frame, set));
    REQUIRE(set != nullptr);
    CHECK_EQ(set->ident(), 7);
    auto o = unpack(set);

    // sorted by (channel, slice)
    std::vector<std::array<int, 2>> want = {{5000, 99}, {5000, 100}, {5001, 100}, {5010, 300}};
    CHECK(o.coords == want);
    std::vector<float> qwant = {0.75f, 4.5f, 10.0f, 25.0f};
    const double m = 3.75, M = 83861.2;
    for (size_t i = 0; i < want.size(); ++i) {
        CHECK(o.feat[i * FAKE_C + 0] == doctest::Approx(logmap(qwant[i], m, M)).epsilon(1e-6));
        CHECK(o.feat[i * FAKE_C + 1] == 1.0f);
        CHECK(o.feat[i * FAKE_C + 3] == doctest::Approx(4 * logmap(qwant[i], m, M) + 0.5f).epsilon(1e-6));
    }
    // metadata: bbox rows 200..210 (5000-4800), slices 99..300; canvas floored in rows, tight in slices
    const auto& p = o.md["planes"][0];
    CHECK_EQ(p["plane"].asInt(), 2);
    CHECK_EQ(p["base_channel"].asInt(), 4800);
    CHECK_EQ(p["n_active"].asInt(), 4);
    CHECK_EQ(p["bbox"][0].asInt(), 200);
    CHECK_EQ(p["bbox"][1].asInt(), 210);
    CHECK_EQ(p["bbox"][2].asInt(), 99);
    CHECK_EQ(p["bbox"][3].asInt(), 300);
    CHECK_EQ(p["canvas"][0].asInt(), 64);
    CHECK_EQ(p["canvas"][1].asInt(), 202);
    CHECK_EQ(p["tiles"].asInt(), 1);
    CHECK_EQ(o.md["tick_span"].asInt(), 4);
    CHECK_EQ(o.md["input_tag"].asString(), "gauss0");
    CHECK_EQ(o.md["frame_ident"].asInt(), 7);

    // EOS
    ITensorSet::pointer eos;
    REQUIRE(fm(nullptr, eos));
    CHECK(eos == nullptr);
}

TEST_CASE("fm pack: halo tiling is seamless and half storage round-trips")
{
    auto anodes = Testing::anodes("uboone");
    REQUIRE(anodes.size() > 0);
    make_FMFakeForward_factory();  // the factory lives in this binary, not in a plugin (util doctest precedent)
    Factory::lookup_tn<ITensorForward>("FMFakeForward");

    // a track-like diagonal: 120 channels x 300 slices
    std::vector<std::tuple<int, int, std::vector<float>>> specs;
    for (int i = 0; i < 120; ++i) {
        specs.emplace_back(5000 + i, 40 + 10 * i, std::vector<float>{10.0f + i, 20.0f, 5.0f, 1.0f, 3.0f});
    }
    auto frame = make_frame(11, specs);

    Pytorch::FMFeatureExtract plain;
    auto cfg = base_cfg(plain);
    plain.configure(cfg);
    ITensorSet::pointer sp;
    REQUIRE(plain(frame, sp));
    auto op = unpack(sp);
    CHECK_EQ(op.md["planes"][0]["tiles"].asInt(), 1);
    CHECK(op.coords.size() >= 120);

    Pytorch::FMFeatureExtract tiled;
    auto tcfg = base_cfg(tiled);
    tcfg["max_dense_pixels"] = 120 * 90;  // canvas 120 x ~303 -> several tiles of core 90-2*16
    tcfg["halo"] = 16;
    tiled.configure(tcfg);
    ITensorSet::pointer st;
    REQUIRE(tiled(frame, st));
    auto ot = unpack(st);
    CHECK(ot.md["planes"][0]["tiles"].asInt() > 2);
    CHECK(ot.coords == op.coords);
    REQUIRE(ot.feat.size() == op.feat.size());
    for (size_t i = 0; i < op.feat.size(); ++i) {
        CHECK(ot.feat[i] == op.feat[i]);
    }

    Pytorch::FMFeatureExtract half;
    auto hcfg = base_cfg(half);
    hcfg["store_half"] = true;
    half.configure(hcfg);
    ITensorSet::pointer sh;
    REQUIRE(half(frame, sh));
    auto oh = unpack(sh);
    CHECK_EQ(oh.feat_name, "feat_half");
    CHECK_EQ(oh.feat_dtype, "u2");
    CHECK(oh.coords == op.coords);
    REQUIRE(oh.feat.size() == op.feat.size());
    for (size_t i = 0; i < op.feat.size(); ++i) {
        CHECK(std::abs(oh.feat[i] - op.feat[i]) <= 1e-3 * std::max(1.0f, std::abs(op.feat[i])));
    }
}
