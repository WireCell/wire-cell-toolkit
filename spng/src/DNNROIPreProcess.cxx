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

void DNNROIPreProcess::configure(const WireCell::Configuration& cfg)
{
    #ifdef HAVE_NVTX
    log->debug("NVTX support enabled");
    #else
    log->debug("NVTX support not enabled");
    #endif
    m_cfg.input_scale = get(cfg, "input_scale", m_cfg.input_scale);
    m_cfg.input_offset = get(cfg, "input_offset", m_cfg.input_offset);

    m_cfg.nchunks = get(cfg, "nchunks", m_cfg.nchunks);
}

//preprocessing
std::vector<torch::Tensor> DNNROIPreProcess::preprocess(const ITorchTensorSet::pointer& input)
{
    NVTX_SCOPED_RANGE("DNNROIPreProcess::preprocess");
    std::vector<torch::Tensor> to_save = {};
    auto tensors = input->tensors();
    //make sure there is target and mp2/mp3 information
    if(tensors->size()!=3) {
        log->error("DNNROIPreProcess: Expecting 3 input tensors, got {}", tensors->size());
        throw std::runtime_error("DNNROIPreProcess: Invalid number of input tensors");
    }
    torch::Tensor ten_target = tensors->at(2)->tensor().clone(); //target plane
    torch::Tensor ten_mp2 = tensors->at(1)->tensor().clone(); //mp2 plane
    torch::Tensor ten_mp3 = tensors->at(0)->tensor().clone(); //mp3 plane

    //Apply scaling and offset
    ten_target = ten_target*m_cfg.input_scale + m_cfg.input_offset;
    ten_mp2 = ten_mp2*m_cfg.input_scale + m_cfg.input_offset;
    ten_mp3 = ten_mp3*m_cfg.input_scale + m_cfg.input_offset;

    //downsample target tensor

    int nticks_a = ten_target.size(2); //number of ticks in target plane
    int nticks_mp2 = ten_mp2.size(2); //number of ticks in mp2 plane
    int tick_per_slice = nticks_a/nticks_mp2; //downsampling factor
    int nticks_ds_a = nticks_a / tick_per_slice; //number of ticks after downsampling
    //make sure that the target plane is downsampled to the mp2/mp3 plane size
    torch::Tensor a_trimmed = ten_target.index({torch::indexing::Slice(), torch::indexing::Slice(), torch::indexing::Slice(0, nticks_ds_a*tick_per_slice)});
    torch::Tensor a_reshaped = a_trimmed.view({a_trimmed.size(0), a_trimmed.size(1), nticks_ds_a, tick_per_slice});
    torch::Tensor a_downsampled = torch::mean(a_reshaped, /*dim=*/3);
    to_save.push_back(a_downsampled);
    to_save.push_back(ten_target);
    to_save.push_back(ten_mp2);
    to_save.push_back(ten_mp3);

    //now stack tensors along the new dimension [3, 800, 1500]
    torch::Tensor stacked = torch::stack({a_downsampled, ten_mp2, ten_mp3}, /*dim=*/1);
    log->debug("DNNROIPreProcess: Stacked tensor shape: {}", tensor_shape_string(stacked));
    //remove the batch dimension
    stacked = stacked.squeeze(0);
    //now transpose the tensor
    auto transposed = torch::stack({torch::transpose(stacked,1,2)},0);
    log->debug("DNNROIPreProcess: Transposed tensor shape: {}", tensor_shape_string(transposed));
    if(transposed.dtype() != torch::kFloat32){
        transposed = transposed.to(torch::kFloat32);
        log->debug("DNNROIPreProcess: Converted tensor to float32");
    }

    //now chunk along the channel dimension (2)
    int nchunks = m_cfg.nchunks;
    std::vector<torch::Tensor> chunk = transposed.chunk(nchunks, /*dim=*/2);

    //store the meta-data for post-processing if needed
    m_metadata = Configuration{};
    m_metadata["tick_per_slice"] = tick_per_slice;
    m_metadata["nticks_target"] = nticks_a;
    m_metadata["nticks_mp2"] = nticks_mp2;
    m_metadata["nticks_ds_a"] = nticks_ds_a;
    m_metadata["nchunks"] = nchunks;
    return chunk;
}

Configuration DNNROIPreProcess::get_metadata() const
{
    return m_metadata;
}