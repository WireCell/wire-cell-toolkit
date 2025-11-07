#include "WireCellSpng/DNNROIPostProcess.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/Util.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellSpng/TorchnvTools.h"



WIRECELL_FACTORY(SPNGDNNROIPostProcess,
    WireCell::SPNG::DNNROIPostProcess,
    WireCell::SPNG::IDNNROIPostProcess,
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
std::vector<torch::Tensor> DNNROIPostProcess::postprocess(
    const std::vector<torch::Tensor>& dnn_output,
    const Configuration& preprocess_metadata)
{
    NVTX_SCOPED_RANGE("DNNROIPostProcess::postprocess");

    if (dnn_output.empty()) {
        log->warn("DNNROIPostProcess: empty DNN output");
        return {};
    }
    int nchunks = m_cfg.nchunk;
    if (nchunks > static_cast<int>(dnn_output.size())) {
        log->warn("DNNROIPostProcess: metadata nchunks ({}) exceeds provided outputs ({}); clamping", nchunks, dnn_output.size());
        nchunks = static_cast<int>(dnn_output.size());
    }

    // chunks should be merged before other operations
    std::vector<torch::Tensor> merged_dnn_output;
    if (nchunks > 1) {
        std::vector<torch::Tensor> merged_chunks;
        merged_chunks.reserve(nchunks);
        for (int i = 0; i < nchunks; ++i) {
            merged_chunks.push_back(dnn_output[i]);
        }
        // concatenate along tick dimension (dim=3) as in preprocessing
        torch::Tensor merged = torch::cat(merged_chunks, /*dim=*/3);
        merged_dnn_output = {merged};
    } else {
        merged_dnn_output = dnn_output; // pass-through
    }

    // Apply output scaling and upsample along the tick dimension
    using torch::indexing::Slice;
    std::vector<torch::Tensor> processed_outputs;
    processed_outputs.reserve(merged_dnn_output.size());

    for (const auto& tensor : merged_dnn_output) {
        if (!tensor.defined()) {
            log->warn("DNNROIPostProcess: encountered undefined tensor; skipping");
            continue;
        }
        torch::Tensor output = tensor.to(torch::kFloat32) * m_cfg.output_scale + m_cfg.output_offset;
        //shape of output
        log->debug("Postprocessed output before upsampling shape: {}", tensor_shape_string(output));

        output = torch::transpose(output,2,3); //transpose back to [B,C,T,W]
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
        //before upsampling shape
        SPNG::write_torch_to_npy(output, "DNNROIPostProcess_output_before_upsample.pt");
        output = output.repeat_interleave(m_cfg.tick_per_slice, /*dim=*/2); //upsample along tick dimension

        log->debug("Postprocessed output shape: {}", tensor_shape_string(output));
        processed_outputs.push_back(output);
    }
    return processed_outputs;
}