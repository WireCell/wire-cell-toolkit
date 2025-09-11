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
   //m_cfg.apa = get(cfg, "apa",m_cfg.apa);
   m_cfg.plane = get(cfg, "plane", m_cfg.plane); 
   log->debug("DNNROI: Configuring with plane: {}", m_cfg.plane);
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
    if(tensors->size()!=3) {
        log->error("DNNROI: Expecting 3 input tensors, got {}", tensors->size());
        return false;   
    } 
    torch::Tensor a = tensors->at(2)->tensor().clone(); //target plane
    torch::Tensor b = tensors->at(1)->tensor().clone(); //MP2
    torch::Tensor c = tensors->at(0)->tensor().clone(); //MP3
    //print the shapes of the tensors
    log->debug("DNNROI: Input tensor a shape: {}", tensor_shape_string(a));
    log->debug("DNNROI: Input tensor b shape: {}", tensor_shape_string(b));
    log->debug("DNNROI: Input tensor c shape: {}", tensor_shape_string(c));

    a = a*m_cfg.input_scale + m_cfg.input_offset;
    b = b*m_cfg.input_scale + m_cfg.input_offset;
    c = c*m_cfg.input_scale + m_cfg.input_offset;

    // b and c have 1500 ticks and a has 6000 ticks
    //assert that both b and c have the same number of ticks
    if(b.size(2) != c.size(2)) {
        log->error("DNNROI: Expecting aux_tensor_l and aux_tensor_m to have the same number of ticks, got {} and {}", b.size(2), c.size(2));
        return false;   
    }
    int tick_per_slice = a.size(2) / b.size(2);
    log->debug("DNNROI: tick_per_slice: {}", tick_per_slice);
    int nticks_a = a.size(2); 
    log->debug("DNNROI: nticks_a: {}", nticks_a);
    int nticks_ds_a = nticks_a / tick_per_slice;
    log->debug("DNNROI: nticks_ds_a: {}", nticks_ds_a);
    //preprocessing of the target tensor
    //shape of the each tensor is [1, 800, 6000] for a and [1, 800, 1500] for b and c
    //downsample a to have the same number of ticks as b and c
    torch::Tensor a_trimmed = a.index({torch::indexing::Slice(), torch::indexing::Slice(), torch::indexing::Slice(0, nticks_ds_a * tick_per_slice)});
    log->debug("DNNROI: a_trimmed shape: {}", tensor_shape_string(a_trimmed));
    torch::Tensor a_reshaped = a_trimmed.view({a_trimmed.size(0), a_trimmed.size(1), nticks_ds_a, tick_per_slice});
    log->debug("DNNROI: a_reshaped shape: {}", tensor_shape_string(a_reshaped));
    torch::Tensor a_ds = torch::mean(a_reshaped, 3);
    log->debug("DNNROI: a_ds shape: {}", tensor_shape_string(a_ds));
    std::vector<torch::Tensor> to_save = {a_ds};
    //save the tensors to a file for debugging
   // torch::save(to_save, "DNNROI_debug.pt");
    
    
    //Do the stacking along a new dimension [3, 800, 1500]
    //now a_ds, b and c have the same shape [1, 800, 1500]
    torch::Tensor stacked = torch::stack({a_ds, b, c}, 1);
    //torch::Tensor stacked = torch::stack({a_ds, b, c}, 0);
    log->debug("DNNROI: Stacked shape: {}", tensor_shape_string(stacked));
    stacked = stacked.squeeze(0); //remove the batch dimension if needed [3, 800, 1500]
    log->debug("DNNROI: Stacked shape after squeeze: {}", tensor_shape_string(stacked));
    //transpose the stacked tensor
    auto batch = torch::stack({torch::transpose(stacked, 1, 2)},0);

    log->debug("DNNROI: Batch shape: {}", tensor_shape_string(batch));
    
    //now the chunking
    //now chunk along the channel dimension (2)
    int nchunks = m_cfg.nchunks;
    auto chunks = batch.chunk(nchunks, 2); // chunk the batch into

    //TO DO: Later on we want to forward the chunks in parallel without data having to leave the GPU
    std::vector<torch::Tensor> outputs;
    for (auto& chunk : chunks) {
        log->debug("DNNROI: Chunk shape: {}, nchunks {}", chunk.sizes()[3], m_cfg.nchunks);
        std::vector<torch::IValue> inputs = {chunk};
        log->debug("DNNROI: Chunk shape: {}", tensor_shape_string(chunk));
        auto iitens = to_itensor(inputs); // convert inputs to ITorchTensorSet        
        log->debug("DNNROI: ITorchTensorSet shape: {}", tensor_shape_string(iitens->tensors()->at(0)->tensor()));
        log->debug("DNNROI: Forwarding chunk with shape: {} from inputs of shape {}", tensor_shape_string(iitens->tensors()->at(0)->tensor()),tensor_shape_string(chunk));
        auto oitens = m_forward->forward(iitens);
        //log->debug("DNNROI: Output chunk shape: {}", tensor_shape_string(oitens->tensors()->at(0)->tensor()));
        //torch::Tensor out_chunk = oitens.toTensor().to(torch::kCUDA); // keep the data in gpu if needed for other stuff.
        torch::Tensor out_chunk = from_itensor(oitens, m_is_gpu)[0].toTensor(); // convert ITorchTensorSet to torch::Tensor
        //log->debug("DNNROI: Output chunk shape: {}", tensor_shape_string(out_chunk));
        outputs.push_back(out_chunk.clone());
    }

    //now concatenate along the time dimensions (3)
    torch::Tensor output = torch::cat(outputs, 3);
    log->debug("DNNROI: Output shape: {}", tensor_shape_string(output));
    //shape of the output is [1, 1, 1500, 800] 
    //we want to reshape it to [1, 3, 800, 1500] to match the input shape
    //so that we can replace the target plane with the output tensor
    //output = output.view({1, 3, 800, nticks_ds_a * tick_per_slice});
    //now create the output tensor with the same information as input tensor but replacing target plane with output tensor
    auto shared_vec = std::make_shared<ITorchTensor::vector>();
    //now post processing
    output = output * m_cfg.output_scale + m_cfg.output_offset;
    output = output.squeeze(); //remove the batch and channel dimensions [1,1,1500,600] --> [1500,600]
    log->debug("DNNROI: Output shape after squeeze: {}", tensor_shape_string(output));
    //transpose the output to have the same shape as input [1500,600] --> [600, 1500]
    output = torch::transpose(output, 0, 1);
    log->debug("DNNROI: Output shape after transpose: {}", tensor_shape_string(output));
    output = output.unsqueeze(0); //add the batch dimension back [1, 800
    //now upsampe the output to have the same number of ticks as input
    output = output.repeat({1, 1, tick_per_slice}); //repeat along the time dimension
    log->debug("DNNROI: Output shape after upsampling: {}", tensor_shape_string(output));
    //trim the output to have the same number of ticks as input
    output = output.index({torch::indexing::Slice(), torch::indexing::Slice(), torch::indexing::Slice(0, nticks_a)});
    log->debug("DNNROI: Output shape after trimming: {}", tensor_shape_string(output));
    //now the output tensor has the same shape as input tensor a [1, 800, 6000]
    //save the output tensor for debugging
    to_save.push_back(output);
    torch::save(to_save, "DNNROI_debug_output.pt");
    //split the output into 3 tensors
    log->debug("DNNROI: Final Output shape after unsqueeze: {}", tensor_shape_string(output));
    auto out_ptr = std::make_shared<SimpleTorchTensor>(output, tensors->at(0)->metadata());
    //first make Simple TorchTensor for each output tensor
    shared_vec->push_back(out_ptr);

    out = std::make_shared<SimpleTorchTensorSet>(
        in->ident(), in->metadata(),
        shared_vec
    );
    
    log->debug("DNNROI: Finished processing, returning output with {} tensors", out->tensors()->size());
    return true;
}
