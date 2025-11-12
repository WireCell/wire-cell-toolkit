#include "WireCellSpng/DNNROIPreProcess.h"
#include "WireCellSpng/Util.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/TorchnvTools.h"

WIRECELL_FACTORY(SPNGDNNROIPreProcess,
     WireCell::SPNG::DNNROIPreProcess, 
    WireCell::INamed,
    WireCell::SPNG::ITorchTensorSetFilter,
     WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::SPNG;

//Implementation of inc/WireCellSpng/DNNROIPreProcess.h
DNNROIPreProcess::DNNROIPreProcess()
: Aux::Logger("SPNGDNNROIPreProcess", "spng")
{
}

DNNROIPreProcess::~DNNROIPreProcess()
{
}

Configuration DNNROIPreProcess::default_configuration() const
{
    Configuration cfg;
    cfg["input_scale"] = 1.0;
    cfg["input_offset"] = 0.0;
    cfg["nchunks"] = 1;
    return cfg;
}

void DNNROIPreProcess::configure(const WireCell::Configuration& cfg)
{
    m_cfg.input_scale = get(cfg, "input_scale", m_cfg.input_scale);
    m_cfg.input_offset = get(cfg, "input_offset", m_cfg.input_offset);
    m_cfg.nchunks = get(cfg, "nchunks", m_cfg.nchunks);
    m_cfg.nticks = get(cfg, "nticks", m_cfg.nticks);
    m_cfg.tick_per_slice = get(cfg, "tick_per_slice", m_cfg.tick_per_slice);
    m_cfg.preprocess_all = get(cfg, "preprocess_all", m_cfg.preprocess_all);
}

//preprocessing
bool DNNROIPreProcess::operator()(const input_pointer& in, output_pointer& out)
{
    NVTX_SCOPED_RANGE("DNNROIPreProcess::preprocess");
    out = nullptr;
    if(!in){
        log->debug("DNNROIPreProcess: EOS ");
        return true;
    }

    auto tensors = in->tensors();
    if (!tensors || tensors->size() != 3) {
        log->error("DNNROIPreProcess: Expecting 3 input tensors, got {}", tensors ? tensors->size() : 0);
        throw std::runtime_error("DNNROIPreProcess: Invalid number of input tensors");
    }

    torch::Tensor ten_target = tensors->at(2)->tensor().clone(); // target plane
    torch::Tensor ten_mp2    = tensors->at(1)->tensor().clone(); // mp2 plane
    torch::Tensor ten_mp3    = tensors->at(0)->tensor().clone(); // mp3 plane
    //print shape of each tensor
    
    log->debug("DNNROIPreProcess: Input tensor shapes - target: {}, mp2: {}, mp3: {}",
               tensor_shape_string(ten_target),
               tensor_shape_string(ten_mp2),
               tensor_shape_string(ten_mp3));
    
    if (ten_target.dim() < 3 || ten_mp2.dim() < 3 || ten_mp3.dim() < 3) {
        log->error("DNNROIPreProcess: Each input tensor must have at least 3 dims, got target={}, mp2={}, mp3={}",
                   ten_target.dim(), ten_mp2.dim(), ten_mp3.dim());
        throw std::runtime_error("DNNROIPreProcess: Invalid tensor rank");
    }

    if (ten_mp2.size(2) != ten_mp3.size(2)) {
        log->error("DNNROIPreProcess: mp2 and mp3 must have same tick dimension: mp2={}, mp3={}",
                   ten_mp2.size(2), ten_mp3.size(2));
        throw std::runtime_error("DNNROIPreProcess: Inconsistent tick dimension between mp2 and mp3");
    }
    SPNG::write_torch_to_npy(ten_target, "DNNROIPreProcess_target_before.pt");
    SPNG::write_torch_to_npy(ten_mp2, "DNNROIPreProcess_mp2_before.pt");
    SPNG::write_torch_to_npy(ten_mp3, "DNNROIPreProcess_mp3_before.pt");
    // Downsample the target tensor along tick dimension to match mp2/mp3
    //Save before downsampling
    SPNG::write_torch_to_npy(ten_target, "DNNROIPreProcess_target_before_downsample.pt");
    ten_target = ten_target.view({ten_target.size(0), ten_target.size(1), m_cfg.nticks/m_cfg.tick_per_slice, m_cfg.tick_per_slice});
    //average it...
    ten_target = ten_target.mean(-1);
    log->debug("DNNROIPreProcess: Target tensor shape after downsampling: {}", tensor_shape_string(ten_target));
    SPNG::write_torch_to_npy(ten_target, "DNNROIPreProcess_target_after_downsample.pt");
    //if all preprocess is true, then preprocess all tensors
    if(m_cfg.preprocess_all){
        //downsample mp2
        ten_mp2 = ten_mp2.view({ten_mp2.size(0), ten_mp2.size(1), m_cfg.nticks/m_cfg.tick_per_slice, m_cfg.tick_per_slice});
        ten_mp2 = ten_mp2.mean(-1);
        log->debug("DNNROIPreProcess: mp2 tensor shape after downsampling: {}", tensor_shape_string(ten_mp2));
        SPNG::write_torch_to_npy(ten_mp2, "DNNROIPreProcess_mp2_after_downsample.pt");
        //downsample mp3
        ten_mp3 = ten_mp3.view({ten_mp3.size(0), ten_mp3.size(1), m_cfg.nticks/m_cfg.tick_per_slice, m_cfg.tick_per_slice});
        ten_mp3 = ten_mp3.mean(-1);
        log->debug("DNNROIPreProcess: mp3 tensor shape after downsampling: {}", tensor_shape_string(ten_mp3));
        SPNG::write_torch_to_npy(ten_mp3, "DNNROIPreProcess_mp3_after_downsample.pt");
    }
    // Apply scaling and offset
    
    ten_target = ten_target * m_cfg.input_scale + m_cfg.input_offset;
    ten_mp2    = ten_mp2    * m_cfg.input_scale + m_cfg.input_offset;
    ten_mp3    = ten_mp3    * m_cfg.input_scale + m_cfg.input_offset;

    // Stack channels: [B, 3, C, T]
    torch::Tensor stacked = torch::stack({ten_target, ten_mp2, ten_mp3}, /*dim=*/1);
    log->debug("DNNROIPreProcess: Stacked tensor shape: {}", tensor_shape_string(stacked));

    // Remove the batch dimension if it's singleton -> [3, C, T]
    stacked = stacked.squeeze(0);

    // Transpose to [3, T, C] and re-introduce batch dim -> [1, 3, T, C]
    auto transposed = torch::transpose(stacked, /*dim0=*/1, /*dim1=*/2).unsqueeze(0);
    log->debug("DNNROIPreProcess: Transposed tensor shape: {}", tensor_shape_string(transposed));

    if (transposed.dtype() != torch::kFloat32) {
        transposed = transposed.to(torch::kFloat32);
        log->debug("DNNROIPreProcess: Converted tensor to float32");
    }

    auto shared_vec = std::make_shared<ITorchTensor::vector>();
    shared_vec->push_back(std::make_shared<SimpleTorchTensor>(transposed, tensors->at(0)->metadata()));
    shared_vec->push_back(std::make_shared<SimpleTorchTensor>(ten_mp2, tensors->at(1)->metadata()));
    shared_vec->push_back(std::make_shared<SimpleTorchTensor>(ten_mp3, tensors->at(2)->metadata()));
    out = std::make_shared<SimpleTorchTensorSet>(in->ident(), in->metadata(), shared_vec);
    log->debug("DNNROIPreProcess: Output ITorchTensorSet created with {} tensors", shared_vec->size());
    return true;
}

