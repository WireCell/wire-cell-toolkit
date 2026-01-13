#include "WireCellSpng/DNNROIPostProcess.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/Util.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellSpng/TorchnvTools.h"



WIRECELL_FACTORY(SPNGDNNROIPostProcess,
    WireCell::SPNG::DNNROIPostProcess,
    WireCell::INamed,
    WireCell::SPNG::ITorchTensorSetFilter,
    WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::SPNG;

DNNROIPostProcess::DNNROIPostProcess()
    : Aux::Logger("SPNGDNNROIPostProcess", "spng")
{
}

DNNROIPostProcess::~DNNROIPostProcess()
{
}

Configuration DNNROIPostProcess::default_configuration() const
{
    Configuration cfg;
    cfg["output_scale"] = 1.0;
    cfg["output_offset"] = 0.0;
    cfg["ntick"] = 6000;
    cfg["nchunks"] = 1;
    cfg["tick_per_slice"] = 4;
    return cfg;
}

void DNNROIPostProcess::configure(const Configuration& cfg)
{

    #ifdef HAVE_NVTX
    log->debug("NVTX support enabled");
    #else
    log->debug("NVTX support not enabled");
    #endif
    
    m_cfg.output_scale = get(cfg, "output_scale", m_cfg.output_scale);
    m_cfg.output_offset = get(cfg, "output_offset", m_cfg.output_offset);
    m_cfg.ntick = get(cfg, "ntick", m_cfg.ntick);
    m_cfg.tick_per_slice = get(cfg, "tick_per_slice", m_cfg.tick_per_slice);
    m_cfg.nchunk = get(cfg, "nchunks", m_cfg.nchunk);
}
//TODO: Should I use std::vector<torch::Tensor>& as input instead?
bool DNNROIPostProcess::operator()(const input_pointer &in, output_pointer& out)
{
    NVTX_SCOPED_RANGE("DNNROIPostProcess::postprocess");
    out = nullptr;
    if(!in){
        log->debug("DNNROIPostProcess: EOS ");
        return true;
    }
    //int nchunks = m_cfg.nchunk;
    //Needs proper implementation to account for the chunking....
    auto tensors = in->tensors();
    std::vector<torch::Tensor> merged_dnn_output;
    std::vector<torch::Tensor> merged_chunks;
    int nchunks = m_cfg.nchunk;
    merged_chunks.reserve(nchunks);
    for (int i = 0; i < nchunks; ++i) {
        merged_chunks.push_back(tensors->at(i)->tensor());
    }
    // concatenate along tick dimension (dim=3) as in preprocessing
    torch::Tensor merged = torch::cat(merged_chunks, /*dim=*/3);
    merged_dnn_output = {merged};
     

    // Apply output scaling and upsample along the tick dimension
    using torch::indexing::Slice;
    //TODO: Do we really need the for loop here?
    auto shared_vec = std::make_shared<ITorchTensor::vector>();
    for (const auto& tensor : merged_dnn_output) {
        if (!tensor.defined()) {
            log->warn("DNNROIPostProcess: encountered undefined tensor; skipping");
            continue;
        }
        torch::Tensor output = tensor.to(torch::kFloat32) * m_cfg.output_scale + m_cfg.output_offset;
        //shape of output
        log->debug("Postprocessed output before upsampling shape: {}", tensor_shape_string(output));

        //output = torch::transpose(output,2,3); //transpose back to [B,C,T,W]
        log->debug("Postprocessed output after transpose shape: {}", tensor_shape_string(output));
        //remove the batch dimension if exists
        if (output.dim() == 4 && output.size(0) == 1) {
            output = output.squeeze(0); // [C,T,W]
        }
        else if (output.dim() == 3) {
            // do nothing
        }
        else {
            log->error("DNNROIPostProcess: unexpected tensor shape after transpose: {}", tensor_shape_string(output));
            throw std::runtime_error("DNNROIPostProcess: unexpected tensor shape");
        }
        //now the upsampling
        //*******NOTE: How you upsample depends upon how you downsample...be consistent*****
        //This assumes you used simple mean downsampling (x.mean(-1))
        output = output.repeat_interleave(m_cfg.tick_per_slice, /*dim=*/2); //upsample along tick dimension
        //write the output after upsample
        SPNG::write_torch_to_npy(output, "DNNROIPostProcess_output_after_upsample.pt");
        shared_vec->push_back(std::make_shared<SimpleTorchTensor>(output, in->tensors()->at(0)->metadata()));
        log->debug("Postprocessed output shape: {}", tensor_shape_string(output));
    }
    out = std::make_shared<SimpleTorchTensorSet>(in->ident(), in->metadata(), shared_vec);
    return true;
}