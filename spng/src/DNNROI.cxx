#include "WireCellSpng/DNNROI.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/Util.h"
#include "WireCellUtil/NamedFactory.h"

#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/ITrace.h"
#include "WireCellAux/PlaneTools.h"

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>


//read the tag to get the apa, filter-name and filter-type information
std::tuple<int, std::string, std::string> parseTag(const std::string& tag) {
    std::istringstream ss(tag);
    std::string segment;
    std::string segments[3];
    int i = 0;

    while (std::getline(ss, segment, ':') && i < 3) {
        segments[i++] = segment;
    }

    int plane = -1;
    if (segments[0].find("collated_") == 0) {
        plane = std::stoi(segments[0].substr(9));  // "collated_" is 9 chars
    }

    return {plane, segments[1], segments[2]};
}

//register this (DNNROI) as a factory
WIRECELL_FACTORY(SPNGDNNROI,// name of the factory
    WireCell::SPNG::DNNROI, // name of the class
    WireCell::INamed, // name of the interface 1 (allows object to have unique name)
    WireCell::ITorchTensorSetFilter, // interface 2 (process ITorchTensorSet)
    WireCell::SPNG::ITorchForward, //interface 2 (process ITorchForward)
    WireCell::IConfigurable // interface 3 (allows configuration)
    )

using namespace WireCell;
using namespace WireCell::SPNG;




//First step. Get the output of the Collator here as input. 
DNNROI::DNNROI()
    : Aux::Logger("SPNGDNNROI", "spng")
{
}
DNNROI::~DNNROI()
{
}

void DNNROI::configure(const WireCell::Configuration& cfg)
{
   m_cfg.apa = get(cfg, "apa",m_cfg.apa);
   m_cfg.plane = get(cfg, "plane", m_cfg.plane); 
   log->debug("DNNROI: Configuring with apa: {}, plane: {}", m_cfg.apa, m_cfg.plane);
   // Is it implemented already>
   /*
   auto apa = Factory::find_tn<IAnodePlane>(m_cfg.apa);
   auto ichans = Aux::plane_channels(apa,m_cfg.plane); //aux/src/PlaneTools.cxx

   //channel information

   for(const auto & ichan : ichans) {
    auto chid = ichan->ident();
    m_chset.insert(chid);
    m_chlist.push_back(chid);
   }

   //sort the channels
   m_cfg.sort_chanids = get(cfg, "sort_chanids", m_cfg.sort_chanids);
   if(m_cfg.sort_chanids) {
       std::sort(m_chlist.begin(), m_chlist.end());
   }
    */
   m_cfg.input_scale = get(cfg, "input_scale", m_cfg.input_scale);
   m_cfg.input_offset = get(cfg, "input_offset", m_cfg.input_offset);
   m_cfg.output_scale = get(cfg, "output_scale", m_cfg.output_scale);
   m_cfg.output_offset = get(cfg, "output_offset", m_cfg.output_offset);
   if (m_cfg.output_scale != 1.0) {
       log->debug("using output scale: {}", m_cfg.output_scale);
   }
   m_is_gpu = get(cfg, "is_gpu", m_is_gpu);
   m_cfg.forward = get(cfg, "forward", m_cfg.forward);
   log->debug("DNNROI: Configured with input_scale: {}, input_offset: {}, output_scale: {}, output_offset: {}, is_gpu: {}, forward: {}",
                m_cfg.input_scale, m_cfg.input_offset, m_cfg.output_scale, m_cfg.output_offset, m_is_gpu, m_cfg.forward);
  
   try{
    m_forward = Factory::find_tn<ITorchForward>(get(cfg, "forward", m_cfg.forward));
   }
   catch (const std::exception& e) {
       log->debug("DNNROI: Failed to find TorchForward instance with name: {}", m_cfg.forward);
       THROW(ValueError() << errmsg{"Failed to find TorchForward instance"});
   }
}

void DNNROI::finalize()
{

}

std::string tensor_shape_string(const at::Tensor& t) {
    std::ostringstream oss;
    oss << "[";
    for (size_t i = 0; i < t.sizes().size(); ++i) {
        oss << t.sizes()[i];
        if (i != t.sizes().size() - 1)
            oss << ", ";
    }
    oss << "]";
    return oss.str();
}


bool DNNROI::operator()(const input_pointer& in, output_pointer& out)
{
    out = nullptr;
    log->debug("Calling DNNROI operator()");
    if (!in) {
        log->debug("DNNROI: EOS ");
        return true;
    }
    log->debug("Running DNNROI");
    
     //TODO -- Loop over input tensors
    auto tensors = in->tensors();
    if (tensors->empty()) {
        log->debug("DNNROI: No tensors in input set");
        return false;
    }
    //TimeKeeper tk(fmt::format("call={}",m_save_count));
    //Process Each TorchTensors in input set to torch::Tensor
    //std::vector<torch::Tensor> ch_tensors; //induction planes
    ROIData ch_tensors; //induction planes
    ROIData coll_tensors; //collection planes
    //save the metadata tags from the input tensors
    log->debug("DNNROI: Processing {} tensors", tensors->size());
    for (size_t i = 0; i < tensors->size(); ++i) {
        auto tensor = tensors->at(i)->tensor();
        // process tensor to torch_tensor to write to ch_tensors

        // Check the tensor tags (metadata)
        auto metadata = tensors->at(i)->metadata();
        if (metadata.isMember("tag") && !metadata["tag"].asString().empty()) {
            log->debug("DNNROI: Found tag in metadata: {}", metadata["tag"].asString());
        } else {
            log->warn("DNNROI: No tag found in metadata for tensor {}", i);
        }

        torch::Tensor scaled = tensor*m_cfg.input_scale + m_cfg.input_offset;

        if(i==0){
            //save scaled tensor
            save_torchtensor_data(scaled, "scaled_tensor_0.pt");
        }

        //  ch_eigen.push_back(Array::downsample(arr, m_cfg.tick_per_slice, 1));
        auto tick_per_slice = m_cfg.tick_per_slice;
        auto nticks = tensor.size(1); //0 is channel, 1 is time
        int nticks_ds = nticks / tick_per_slice;
        //keep all the dimensions unchanged, select range from 0, nticks_ds* tick_per_slice
        //This will downsample the tensor by tick_per_slice
        auto trimmed = tensor.index({"...",torch::indexing::Slice(0, nticks_ds * tick_per_slice)}); 
        //save the trimmed tensor
        if(i==0){
            save_torchtensor_data(trimmed, "trimmed_tensor_0.pt");
        }

        //now reshape the tensor
        auto reshaped = trimmed.view({tensor.size(0), nticks_ds, tick_per_slice});
        //save the reshaped tensor
        if(i==0){
            save_torchtensor_data(reshaped, "reshaped_tensor_0.pt");
        }
        //reshaped tensor has the dimensions [channels, downsampled_time, tick_per_slice]
        //now take the mean along the last dimension (tick_per_slice)
        auto downsampled = reshaped.mean(2);
        //save the downsampled tensor
        if(i==0){   
            save_torchtensor_data(downsampled, "downsampled_tensor_0.pt");
        }

        //now cscale and offset the downsampled tensor
        auto dscaled = downsampled * m_cfg.input_scale + m_cfg.input_offset;
        auto nchannels = dscaled.size(0); //number of channels 
        if(nchannels == 800){
            ch_tensors.tensors.push_back(dscaled); //induction plane
            ch_tensors.r_tags.push_back(metadata["tag"].asString()); //store the tag
        }
        else {
            coll_tensors.tensors.push_back(dscaled); //collection plane
            coll_tensors.r_tags.push_back(metadata["tag"].asString()); //store the tag
        }
    }
    for (size_t i = 0; i < ch_tensors.tensors.size(); ++i) {
        log->debug("DNNROI: Tensor {} shape: {} tags: {} ", i, tensor_shape_string(ch_tensors.tensors[i]),ch_tensors.r_tags[i]);
    }
    
    //stack vector of tensors along a new dimension (0)
    auto img = torch::stack(ch_tensors.tensors, 0); //stack all the tensors along a new dimension (0)
    //img is now a 3D stacked tensors 
    auto transposed = img.transpose(1,2); //transposition into ntags, nticks_ds, nchannels
    auto batch = torch::unsqueeze(transposed, 0); //add a batch dimension at the start
    //chunking -->Check the edging effects.
    auto chunks = batch.chunk(m_cfg.nchunks, 2); // chunk the batch into m_cfg.nchunks along the time dimension (2)
    std::vector<torch::Tensor> outputs;
    for (auto& chunk : chunks) {
        log->debug("DNNROI: Chunk shape: {}, nchunks {}", chunk.sizes()[1], m_cfg.nchunks);
        std::vector<torch::IValue> inputs = {chunk};
        log->debug("DNNROI: Chunk shape: {}", tensor_shape_string(chunk));
        auto iitens = to_itensor(inputs); // convert inputs to ITorchTensorSet        
        log->debug("DNNROI: ITorchTensorSet shape: {}", tensor_shape_string(iitens->tensors()->at(0)->tensor()));
        log->debug("DNNROI: Forwarding chunk with shape: {} from inputs of shape {}", tensor_shape_string(iitens->tensors()->at(0)->tensor()),tensor_shape_string(chunk));
        auto oitens = m_forward->forward(iitens);
        //log->debug("DNNROI: Output chunk shape: {}", tensor_shape_string(oitens->tensors()->at(0)->tensor()));
        //torch::Tensor out_chunk = oitens.toTensor().to(torch::kCUDA); // keep the data in gpu if needed for other stuff.
        torch::Tensor out_chunk = from_itensor(oitens, m_is_gpu)[0].toTensor().cpu(); // convert ITorchTensorSet to torch::Tensor
        //AB: Mostly for the debugging purposes, we can remove it later.
       // torch::Tensor out_chunk = from_itensor(iitens, m_is_gpu)[0].toTensor().cpu(); // convert ITorchTensorSet to torch::Tensor
        outputs.push_back(out_chunk);
    }
    
    torch::Tensor out_tensor = torch::cat(outputs, 2); // concatenate the output chunks along the time dimension (2)
    //save the out_tensor
    save_torchtensor_data(out_tensor, "out_tensor.pt");
    auto mask = out_tensor.gt(m_cfg.mask_thresh);
    auto finalized_tensor = out_tensor*mask.to(out_tensor.dtype()); // apply the mask to the output tensor  
    finalized_tensor = finalized_tensor * m_cfg.output_scale + m_cfg.output_offset; // scale and offset the output tensor
    //save the finalized tensor
    save_torchtensor_data(finalized_tensor, "finalized_tensor.pt");
    //std::vector<ITorchTensor::pointer>processed_tensors;
    auto shared_vec = std::make_shared<ITorchTensor::vector>();
    for(int i = 0; i < finalized_tensor.size(0); ++i) {
        auto single_tensor = finalized_tensor[i].detach().clone(); // detach and clone the tensor to avoid modifying the original tensor
        auto meta = tensors->at(i)->metadata(); // get the metadata from the original tensors
        shared_vec->push_back(
            std::make_shared<SimpleTorchTensor>(single_tensor,meta));
    }
      
    out = std::make_shared<SimpleTorchTensorSet>(
        in->ident(), in->metadata(), shared_vec
    );
    //save out tensor set
    save_simpletensor_data(out, "out_tensor_set.pt");
    
    return true;
}
