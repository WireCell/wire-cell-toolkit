#include "WireCellIface/ITensorSet.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/SimpleTensor.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellPytorch/Torch.h" 

#include "grpc_client.h"           

#include <iostream>
#include <memory>
#include <vector>

using namespace WireCell;       
using namespace WireCell::Aux;  
using namespace triton::client; 


// -----------------------------------------------------------------------------
// Helper: convert a vector of Torch IValues into an ITensorSet
// -----------------------------------------------------------------------------
static ITensorSet::pointer to_itensor(const std::vector<torch::IValue>& inputs) {
  //expect exactly one 4‐D tensor per IValue
  if (inputs.size() != 1) {
    THROW(ValueError() << errmsg{"to_itensor: need exactly one input IValue"});
  }
  auto ten = inputs[0].toTensor().cpu();
  if (ten.dim() != 4) {
    THROW(ValueError() << errmsg{"to_itensor: tensor must be 4‐D"});
  }
  // shape = [N, C, H, W]
  std::vector<size_t> shape = {
    (size_t)ten.size(0), (size_t)ten.size(1),
    (size_t)ten.size(2), (size_t)ten.size(3)
  };
  // allocate a SimpleTensor of the right size
  auto st = std::make_shared<WireCell::Aux::SimpleTensor>(shape, (float*)nullptr);
  // copy raw data out of the torch tensor
  const auto* src = ten.data_ptr<float>();
  size_t nbyte = sizeof(float);
  for (auto d : shape) nbyte *= d;
  std::memcpy((void*)st->data(), src, nbyte);
  // wrap into a SimpleTensorSet
  auto seqno = 0;
  WireCell::Configuration md; 
  ITensor::shared_vector tv(new ITensor::vector{st});
  return std::make_shared<WireCell::Aux::SimpleTensorSet>(seqno, md, tv);
}

// -----------------------------------------------------------------------------
// Helper: convert an ITensorSet (with one 4‐D tensor) back into a torch::IValue
// -----------------------------------------------------------------------------
static torch::IValue from_itensor(const ITensorSet::pointer& its) {
  auto tens = its->tensors()->front();
  auto shape = tens->shape();
  if (shape.size() != 4) {
    THROW(ValueError() << errmsg{"from_itensor: tensor must be 4‐D"});
  }
  // build a torch::Tensor that wraps the same memory:
  auto iv = torch::from_blob((float*)tens->data(),
			     { (long)shape[0], (long)shape[1],
                                 (long)shape[2], (long)shape[3] },
			     torch::kFloat32);
  return iv;
}

// -----------------------------------------------------------------------------
// Main: build a dummy input, do Triton inference, print result shapes
// -----------------------------------------------------------------------------
int main() {
  using namespace triton::client;

  
  std::vector<torch::IValue> inputs;
  inputs.push_back(
		   torch::ones({1, 3, 800, 600},
			       torch::dtype(torch::kFloat32).device(torch::kCPU))
		   );

  
  auto iit = to_itensor(inputs);
  std::cout << "Converted to ITensorSet with shape: [";
    for (auto d : iit->tensors()->front()->shape()) std::cout << d << " ";
    std::cout << "]\n";

    
    auto tens = iit->tensors()->front();
    const auto* raw = (const uint8_t*)tens->data();
    size_t byte_size = tens->size();
    std::vector<int64_t> triton_shape = {
        (int64_t)tens->shape()[0],
        (int64_t)tens->shape()[1],
        (int64_t)tens->shape()[2],
        (int64_t)tens->shape()[3]
    };

    //create Triton gRPC client
    std::unique_ptr<InferenceServerGrpcClient> client;
    InferenceServerGrpcClient::Create(&client, "ailab01.fnal.gov:8001");
    
    // prepare InferInput
    InferInput* infer_input = nullptr;
    InferInput::Create(&infer_input,
                       "INPUT__0",      // must match triton deployed model's config.pbtxt
                       triton_shape,
                       "FP32");
    infer_input->AppendRaw(raw, byte_size);

    //prepare InferRequestedOutput
    InferRequestedOutput* infer_output = nullptr;
    InferRequestedOutput::Create(&infer_output, "OUTPUT__0");

    std::vector<InferInput*>  inputs_ptr  = { infer_input };
    std::vector<const InferRequestedOutput*> outputs = { infer_output };
    Headers headers;  

    //run inference
    InferResult* result = nullptr;
    InferOptions options("dnn");  // your model name
    client->Infer(&result, options, inputs_ptr, outputs, headers);

    
    const uint8_t* out_raw = nullptr;
    size_t out_byte_size = 0;
    result->RawData("OUTPUT__0", &out_raw, &out_byte_size);

    // Print mask output
    const float* mask_out = reinterpret_cast<const float*>(out_raw);
    size_t numel = out_byte_size / sizeof(float);

    std::cout << "Output mask size: " << numel << std::endl;
    std::cout << "First 10 mask values: ";
    for (size_t i = 0; i < std::min(numel, size_t(10)); ++i)
      std::cout << mask_out[i] << " ";
    std::cout << std::endl;
    
    
    std::vector<size_t> out_shape = {1, 1, 800, 600};
    auto out_tensor = std::make_shared<WireCell::Aux::SimpleTensor>(
                          out_shape,
                          reinterpret_cast<const float*>(out_raw)
                      );
    ITensor::shared_vector otv(new ITensor::vector{out_tensor});
    auto oit = std::make_shared<WireCell::Aux::SimpleTensorSet>(0, WireCell::Configuration(), otv);

    std::cout << "Received ITensorSet with shape: [";
    for (auto d : oit->tensors()->front()->shape()) std::cout << d << " ";
    std::cout << "]\n";

    
    auto out_ival = from_itensor(oit);
    std::cout << "Back to torch::IValue with dims: ";
    for (auto s : out_ival.toTensor().sizes()) std::cout << s << " ";
    std::cout << "\n";

    // clean up
    delete infer_input;
    delete infer_output;
    delete result;

    return 0;
}


