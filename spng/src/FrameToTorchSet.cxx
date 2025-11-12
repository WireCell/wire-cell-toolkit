#include "WireCellSpng/FrameToTorchSet.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/SimpleFrame.h"

WIRECELL_FACTORY(FrameToTorchSet, WireCell::SPNG::FrameToTorchSet,
                 WireCell::INamed,
                 WireCell::SPNG::IFrameToTorchSet)

using namespace WireCell;
using namespace WireCell::SPNG;

FrameToTorchSet::FrameToTorchSet()
    : Aux::Logger("FrameToTorchSet", "spng")
{
}


WireCell::Configuration FrameToTorchSet::default_configuration() const
{
    Configuration cfg;
    cfg["intags"] = Json::arrayValue;
    return cfg;
}

void FrameToTorchSet::configure(const WireCell::Configuration& cfg)
{
    m_intags.clear();
    if (cfg.isMember("intags")) {
        for (const auto& one : cfg["intags"]) {
            m_intags.push_back(one.asString());
        }
    }
}

bool FrameToTorchSet::operator()(const input_pointer& in, output_pointer& out) {

    out = nullptr;
    if (!in) return true;
    auto ntraces = in->traces()->size();
    
    log->debug("Ntraces: {}", ntraces);

    if (ntraces == 0) {
        log->debug("No traces, exiting");
        return true;
    }

    auto shared_vec = std::make_shared<ITorchTensor::vector>();
    for (const auto& tag : m_intags) {
        auto traces = Aux::tagged_traces(in, tag);
        auto nticks = traces.empty() ? 0 : traces[0]->charge().size();
        log->debug("Making tensor for tag '{}' of shape: {} {}", tag, traces.size(), nticks);
        //make tensor 3D with batch=1
        torch::Tensor tensor = torch::zeros({1, (int64_t)traces.size(), (int64_t)nticks});
        //torch::Tensor tensor = torch::zeros({(int64_t)traces.size(), (int64_t)nticks});
        for (size_t i = 0; i < traces.size(); ++i) {
            for (size_t j = 0; j < nticks; ++j) {
                //tensor.index_put_({static_cast<int64_t>(i), static_cast<int64_t>(j)}, traces[i]->charge()[j]);
                tensor.index_put_({0, static_cast<int64_t>(i), static_cast<int64_t>(j)}, traces[i]->charge()[j]);
            }    
        }
        Configuration md;
        md["tag"] = tag;
        auto out_ptr = std::make_shared<SimpleTorchTensor>(tensor, md);
        shared_vec->push_back(out_ptr);
        log->debug("Put for tag '{}'", tag);
    }
    Configuration md;
    out = std::make_shared<SimpleTorchTensorSet>(
        in->ident(), md,
        shared_vec
    );
    return true;
}
