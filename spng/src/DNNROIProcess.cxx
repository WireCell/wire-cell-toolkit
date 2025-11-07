#include "WireCellSpng/DNNROIProcess.h"
#include "WireCellSpng/DNNROIPreProcess.h"
#include "WireCellSpng/DNNROIPostProcess.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/Util.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellSpng/TorchnvTools.h"

WIRECELL_FACTORY(SPNGDNNROIProcess,
     WireCell::SPNG::DNNROIProcess, 
     WireCell::INamed,
     WireCell::SPNG::ITorchTensorSetFilter,
     WireCell::IConfigurable
    )

using namespace WireCell;
using namespace WireCell::SPNG;

DNNROIProcess::DNNROIProcess()
: Aux::Logger("SPNGDNNROIProcess", "spng")
{
}

DNNROIProcess::~DNNROIProcess()
{
}

//Configuration
void DNNROIProcess::configure(const WireCell::Configuration& cfg)
{
    m_cfg.plane = get(cfg, "plane", m_cfg.plane);
    log->debug("DNNROIProcess: Configured for plane {}", m_cfg.plane);

    m_is_gpu = get(cfg, "is_gpu", m_is_gpu);
    m_cfg.forward = get(cfg, "forward", m_cfg.forward);
    m_cfg.preprocess = get(cfg, "preprocess", m_cfg.preprocess);
    m_cfg.postprocess = get(cfg, "postprocess", m_cfg.postprocess);
    log->debug("DNNROIProcess: Using forwarder {}", m_cfg.forward);

    //now the pre-processor initialization
    //m_preprocess = std::make_shared<DNNROIPreProcess>();
    m_preprocess = Factory::find_tn<IDNNROIPreProcess>(m_cfg.preprocess);
    Configuration pre_cfg;
    pre_cfg["input_scale"] = get(cfg, "input_scale", 1.0);
    pre_cfg["input_offset"] = get(cfg, "input_offset", 0.0);
    pre_cfg["nchunks"] = get(cfg, "nchunks", 1);
    //m_preprocess->configure(pre_cfg);
    log->debug("DNNROIProcess: Preprocess configured");
    
    //now the post-processor initialization
    //m_postprocess = std::make_shared<DNNROIPostProcess>();
    m_postprocess = Factory::find_tn<IDNNROIPostProcess>(m_cfg.postprocess);
    Configuration post_cfg;
    post_cfg["output_scale"] = get(cfg, "output_scale", 1.0);
    post_cfg["output_offset"] = get(cfg, "output_offset", 0.0);
    //m_postprocess->configure(post_cfg);
    log->debug("DNNROIProcess: Postprocess configured");

    //find the TorchForward instance
    try{
        m_forward = Factory::find_tn<ITorchForward>(m_cfg.forward);
    }
    catch(const std::exception& e){
        log->error("DNNROIProcess: Failed to find ITorchForward instance named {}: {}", m_cfg.forward, e.what());
        throw std::runtime_error("DNNROIProcess: Failed to find ITorchForward instance");
    }
}

bool DNNROIProcess::operator()(const input_pointer& in, output_pointer& out)
{
    NVTX_SCOPED_RANGE("DNNROIProcess::operator()");
    out = nullptr;
    log->debug("DNNROIProcess: Calling operator()");
    if (!in) {
        log->debug("DNNROIProcess: EOS ");
        return true;
    }

    auto tensors = in->tensors();
    if (tensors->size()!=3){
        log->error("DNNROIProcess: Expected 3 tensors in input ITorchTensorSet, got {}", tensors->size());
        return false;
    }
    
    //now the pre-processing
    std::vector<torch::Tensor> chunks;
    Configuration preprocess_metadata;

    chunks = m_preprocess->preprocess(in);
    log->debug("DNNROIProcess: Preprocess done, got {} chunks", chunks.size());
    //now the inference
    std::vector<torch::Tensor> outputs;
    for(auto& chunk : chunks){
        log->debug("DNNROIProcess: Sending chunk of size {} to forwarder {}", tensor_shape_string(chunk), m_cfg.forward);
        auto out_tensor = m_forward->forward(chunk);
        outputs.push_back(out_tensor);
    }
    log->debug("DNNROIProcess: Inference done, got {} output chunks", outputs.size());
    SPNG::write_torch_to_npy(outputs[0], "DNNROIProcess_inference_output_chunk0.pt");
    //now the post-processing
    std::vector<torch::Tensor> postprocessed = m_postprocess->postprocess(outputs, preprocess_metadata);
    log->debug("DNNROIProcess: Postprocess done, output ITorchTensorSet has {} tensors", postprocessed.size());
    SPNG::write_torch_to_npy(postprocessed[0], "DNNROIProcess_postprocess_output_tensor0.pt");
    //now convert to ITorchTensorSet
    //loop over postprocessed and make a SimpleTorchTensor for each
    //then make a SimpleTorchTensorSet from the vector of SimpleTorchTensor
    auto shared_vec = std::make_shared<ITorchTensor::vector>();
    for(auto& t : postprocessed){
        auto out_ptr = std::make_shared<SimpleTorchTensor>(t, tensors->at(0)->metadata());
        shared_vec->push_back(out_ptr);
    }
    out = std::make_shared<SimpleTorchTensorSet>(in->ident(), in->metadata(), shared_vec);
    log->debug("DNNROIProcess: Finished operator()");
    return true;
}


void DNNROIProcess::finalize()
{
}
