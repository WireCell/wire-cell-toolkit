#include "WireCellSpng/FrameToTorchSet.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/Util.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/ITrace.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellAux/FrameTools.h"
#include "WireCellAux/PlaneTools.h"
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
    for(const auto& one : m_cfg.intags){
        cfg["intags"].append(one);
    }
    cfg["nticks"] = m_cfg.nticks;
    cfg["anode"] = m_cfg.anode;
    cfg["plane"] = m_cfg.plane;
    cfg["tick0"] = m_cfg.tick0;
    cfg["nticks"] = m_cfg.nticks;
    cfg["sort_chanids"] = m_cfg.sort_chanids;
    return cfg;
}

void FrameToTorchSet::configure(const WireCell::Configuration& cfg)
{
    m_cfg.intags.clear();
    if (cfg.isMember("intags")) {
        for (const auto& one : cfg["intags"]) {
            m_cfg.intags.push_back(one.asString());
        }
    }

    m_cfg.nticks = get(cfg, "nticks", m_cfg.nticks);
    m_cfg.tick0 = get(cfg, "tick0", m_cfg.tick0);

    m_cfg.anode = get(cfg, "anode", m_cfg.anode);
    m_cfg.plane = get(cfg, "plane", m_cfg.plane);
    auto anode = Factory::find_tn<IAnodePlane>(m_cfg.anode);
    auto ichans = Aux::plane_channels(anode, m_cfg.plane);
    //size of selected channels
    log->debug("FrameToTorchSet: ichans size: {}", ichans.size());
    m_chset.clear();
    m_chlist.clear();
    for (const auto& ichan : ichans) {
        auto chid = ichan->ident();
        m_chset.insert(chid);
        m_chlist.push_back(chid);
    }
    m_cfg.sort_chanids = get(cfg, "sort_chanids", m_cfg.sort_chanids);
    if (m_cfg.sort_chanids) {
        std::sort(m_chlist.begin(), m_chlist.end());
    }
    m_nrows = m_chlist.size();
    m_ncols = m_cfg.nticks;

}

ITrace::vector FrameToTorchSet::select(ITrace::vector traces){
    auto end = std::remove_if(traces.begin(), traces.end(),
                              [&](const ITrace::pointer& t) {
                                  return m_chset.find(t->channel()) == m_chset.end();
                              });
   traces.resize(end - traces.begin());
   log->debug("Selected {} traces after channel filtering", traces.size());
   return traces;
}

torch::Tensor FrameToTorchSet::traces_to_tensor(ITrace::vector traces)
{
    Array::array_xxf arr = Array::array_xxf::Zero(m_nrows, m_ncols);
    torch::Tensor tensor = torch::zeros({static_cast<int64_t>(m_ncols), static_cast<int64_t>(m_nrows)});
    traces = select(traces);
    if (traces.empty()) {
        log->debug("No traces selected, returning zero tensor");
        return tensor;
    }
    Aux::fill(arr, traces, m_chlist.begin(), m_chlist.end(), m_cfg.tick0);
    if(arr.sum() == 0.0){
        log->debug("All-zero array after filling traces, returning zero tensor");
        return tensor;
    }
    //print the sum of tensor
    log->debug("Array sum after filling traces: {}", arr.sum());
    //now convert to tensor
    tensor = torch::from_blob(arr.data(), {static_cast<int64_t>(arr.cols()), static_cast<int64_t>(arr.rows())}).clone();
    //print the sum of tensor
    log->debug("Tensor sum after conversion: {}", tensor.sum().item<float>());
    return tensor;
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
    shared_vec->reserve(m_cfg.intags.size());
    for (const auto& tag : m_cfg.intags) {
        auto traces = Aux::tagged_traces(in, tag);
        //auto nticks = traces.empty() ? 0 : traces[0]->charge().size();
        //log->debug("Making tensor for tag '{}' of shape: {} {}", tag, traces.size(), nticks);
        //make tensor 3D with batch=1
        //torch::Tensor tensor = torch::zeros({1, (int64_t)traces.size(), (int64_t)nticks});
        //torch::Tensor tensor = torch::zeros({(int64_t)traces.size(), (int64_t)nticks});
        torch::Tensor tensor = traces_to_tensor(traces);
        SPNG::write_torch_to_npy(tensor, fmt::format("FrameToTorchSet_tensor_tag_{}.pt", tag));
        log->debug("Tensor shape and sum for tag '{}': {} , {}", tag, tensor_shape_string(tensor), tensor.sum().item<float>());
        //These tensors have shape (nticks, nchannels), need to make it (nchannels, nticks)
        tensor = torch::transpose(tensor, 0, 1);
        //add batch dimension
        tensor = tensor.unsqueeze(0); //add batch dimension
        //write the tensor to file for debugging
        log->debug("Made tensor for tag '{}' of shape: {}", tag, tensor_shape_string(tensor));
        Configuration md;
        md["tag"] = /*if mp3 in tag, name it mp3*/ (tag.find("mp3") != std::string::npos) ? "mp3" :
                    (tag.find("mp2") != std::string::npos) ? "mp2" : "target";
        auto out_ptr = std::make_shared<SimpleTorchTensor>(tensor, md);
        shared_vec->push_back(out_ptr);
    }    
    Configuration md;
    out = std::make_shared<SimpleTorchTensorSet>(
        in->ident(), md,
        shared_vec
    );
    return true;
}
