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
    m_cfg.save_debug = get(cfg, "save_debug", m_cfg.save_debug);

    log->debug("DNNROIPostProcess configured with scale: {}, offset: {}, save_debug: {}",
               m_cfg.output_scale, m_cfg.output_offset, m_cfg.save_debug);
}

std::vector<torch::Tensor> DNNROIPostProcess::postprocess(
    const std::vector<torch::Tensor>& dnn_output,
    const Configuration& preprocess_metadata)
{
    // Extract metadata
    NVTX_SCOPED_RANGE("DNNROIPostProcess::postprocess");
    int tick_per_slice = get(preprocess_metadata, "tick_per_slice", 1);
    int nticks_a = get(preprocess_metadata, "nticks_a", 6000);
    int nchunks = get(preprocess_metadata, "nchunks", 1);
    
    // chunks should be merged before other operations
    std::vector<torch::Tensor> merged_dnn_output;
    if (nchunks > 1) {
        std::vector<torch::Tensor> merged_chunks;
        for (int i = 0; i < nchunks; ++i) {
            merged_chunks.push_back(dnn_output[i]);
        }
        torch::Tensor merged = torch::cat(merged_chunks, 2);
        merged_dnn_output = {merged};
    }

    // Apply output scaling
    std::vector<torch::Tensor> processed_outputs;
    for (const auto& tensor : merged_dnn_output) {
        torch::Tensor output = tensor * m_cfg.output_scale + m_cfg.output_offset;
        // Reshape and transpose

        output = output.squeeze(); // Remove batch and channel dimensions
        output = torch::transpose(output, 0, 1);
        output = output.unsqueeze(0); // Add batch dimension back
        
        // Upsample to original time resolution
        output = output.repeat({1, 1, tick_per_slice});
        output = output.index({
            torch::indexing::Slice(), 
            torch::indexing::Slice(), 
            torch::indexing::Slice(0, nticks_a)
        });
        
        log->debug("Postprocessed output shape: {}", tensor_shape_string(output));
        processed_outputs.push_back(output);
    }
    return processed_outputs;
}