/** FMFeatureExtract direct-model test (wcfm doc 34): with `model` set the node runs the TorchScript file itself
    and reads the output tensor in place.  The features must equal, bit for bit, those of the same file run
    through a TorchService (the ITensorForward path of doctest_fm_pack.cxx), untiled and tiled.

    The model is a pointwise script written here: out[c] = x0, x1, (c+1)*x0 + 0.5*x1.
*/
#include "WireCellPytorch/FMFeatureExtract.h"
#include "WireCellPytorch/Torch.h"

#include "WireCellIface/IConfigurable.h"
#include "WireCellAux/SimpleFrame.h"
#include "WireCellAux/SimpleTrace.h"
#include "WireCellAux/Testing.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/PluginManager.h"
#include "WireCellUtil/Units.h"
#include "WireCellUtil/doctest.h"

#include <cstdio>
#include <cstring>

using namespace WireCell;

namespace {
    const int FAKE_C = 6;

    IFrame::pointer make_frame(int ident)
    {
        // a track-like diagonal: 120 channels x 300 slices of the uboone Y plane (channels 4800..8255)
        ITrace::vector traces;
        IFrame::trace_list_t idx;
        for (int i = 0; i < 120; ++i) {
            traces.push_back(std::make_shared<Aux::SimpleTrace>(5000 + i, 40 + 10 * i,
                                                                std::vector<float>{10.0f + i, 20.0f, 5.0f, 1.0f, 3.0f}));
            idx.push_back(idx.size());
        }
        auto sf = new Aux::SimpleFrame(ident, 0.0, traces, 0.5 * units::microsecond);
        sf->tag_traces("gauss0", idx);
        return IFrame::pointer(sf);
    }

    // the raw bytes of every tensor of the set, in order
    std::vector<std::vector<char>> bytes_of(const ITensorSet::pointer& set)
    {
        std::vector<std::vector<char>> out;
        for (const auto& t : *set->tensors()) {
            const char* p = reinterpret_cast<const char*>(t->data());
            out.emplace_back(p, p + t->size());
        }
        return out;
    }
}

TEST_CASE("fm direct: the model run by the node equals the same model run by a TorchService")
{
    auto anodes = Testing::anodes("uboone");
    REQUIRE(anodes.size() > 0);
    const std::string path = "doctest_fm_direct.ts";
    {
        torch::jit::Module m("fmfake");
        m.define(R"(
def forward(self, x):
    v = x[:, 0:1]
    k = x[:, 1:2]
    return torch.cat([v, k, 3.0 * v + 0.5 * k, 4.0 * v + 0.5 * k, 5.0 * v + 0.5 * k, 6.0 * v + 0.5 * k], 1)
)");
        m.save(path);
    }
    PluginManager::instance().add("WireCellPytorch");
    {
        auto icfg = Factory::lookup_tn<IConfigurable>("TorchService:fmdirect");
        auto cfg = icfg->default_configuration();
        cfg["model"] = path;
        cfg["device"] = "cpu";
        icfg->configure(cfg);
    }
    auto frame = make_frame(12);
    for (int tiled = 0; tiled < 2; ++tiled) {
        Pytorch::FMFeatureExtract svc;
        auto scfg = svc.default_configuration();
        scfg["anode"] = "AnodePlane:0";
        scfg["forward"] = "TorchService:fmdirect";
        scfg["input_tag"] = "gauss0";
        scfg["planes"] = Json::arrayValue;
        scfg["planes"].append(2);
        scfg["feature_dim"] = FAKE_C;
        scfg["nticks"] = 2000;
        scfg["store_half"] = false;
        if (tiled) {
            scfg["max_dense_pixels"] = 120 * 90;
            scfg["halo"] = 16;
        }
        svc.configure(scfg);
        ITensorSet::pointer ss;
        REQUIRE(svc(frame, ss));
        REQUIRE(ss);

        Pytorch::FMFeatureExtract direct;
        auto dcfg = scfg;
        dcfg["model"] = path;
        dcfg["forward"] = "NoSuchForward";   // not looked up when model is set
        direct.configure(dcfg);
        ITensorSet::pointer sd;
        REQUIRE(direct(frame, sd));
        REQUIRE(sd);

        const auto bs = bytes_of(ss), bd = bytes_of(sd);
        REQUIRE(bs.size() == 2);
        REQUIRE(bd.size() == 2);
        const size_t npix = ss->tensors()->at(0)->shape()[0];
        CHECK(npix >= 120);
        const bool same_coords = bs[0] == bd[0], same_feat = bs[1] == bd[1];
        CHECK(same_coords);
        CHECK(same_feat);
        const int ntiles = sd->metadata()["planes"][0]["tiles"].asInt();
        CHECK((tiled ? ntiles > 2 : ntiles == 1));
        // the features are the model's: channel 0 of the first pixel is its log-mapped value, channel 1 its mask
        const float* f = reinterpret_cast<const float*>(sd->tensors()->at(1)->data());
        CHECK(f[1] == 1.0f);
        CHECK(f[2] == doctest::Approx(3.0f * f[0] + 0.5f));
    }
    std::remove(path.c_str());
}
