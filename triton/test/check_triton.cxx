#include <iostream>
#include <vector>
#include <thread>
#include <chrono>
#include "grpc_client.h"
#include "WireCellPytorch/Torch.h"

using namespace triton::client;

// Helper function to convert a torch::Tensor to std::vector<float>
std::vector<float> tensor_to_vector(const torch::Tensor& tensor) {
  torch::Tensor cpu_tensor = tensor.contiguous().to(torch::kCPU);
  return std::vector<float>(cpu_tensor.data_ptr<float>(), cpu_tensor.data_ptr<float>() + cpu_tensor.numel());
}

void triton_forward(int threadid, const std::string& url, const std::string& model_name) {
  std::vector<int64_t> input_shape = {1, 3, 800, 600};
  std::string input_name = "INPUT__0";
  std::string output_name = "OUTPUT__0";

  torch::Tensor tensor = torch::rand(input_shape, torch::dtype(torch::kFloat32).device(torch::kCPU));
  std::vector<float> input_data = tensor_to_vector(tensor);

  // 1. Prepare input object
  triton::client::InferInput* input = nullptr;
  triton::client::InferInput::Create(&input, input_name, input_shape, "FP32");
  input->AppendRaw(reinterpret_cast<uint8_t*>(input_data.data()), input_data.size() * sizeof(float));

  // 2. Prepare output request
  triton::client::InferRequestedOutput* output = nullptr;
  triton::client::InferRequestedOutput::Create(&output, output_name);
  std::vector<const triton::client::InferRequestedOutput*> outputs = {output};

  // 3. Prepare the client
  std::unique_ptr<triton::client::InferenceServerGrpcClient> client;
  triton::client::InferenceServerGrpcClient::Create(&client, url);

  // 4. Set inference options and inputs
  triton::client::InferOptions options(model_name);
  std::vector<triton::client::InferInput*> inputs = {input};

  // 5. Run inference
  triton::client::InferResult* result = nullptr;
  triton::client::Headers headers; 
  
  client->Infer(&result, options, inputs, outputs, headers, GRPC_COMPRESS_NONE);
  // 6. Read and print output
  const uint8_t* output_buf;
  size_t output_byte_size;
  result->RawData(output_name, &output_buf, &output_byte_size);
  std::cout << "Thread " << threadid << " got output of size: " << (output_byte_size / sizeof(float)) << std::endl;

  // 7. Clean up
  delete input;
  delete output;
  delete result;
}

int main() {
  const std::string url = "ailab01.fnal.gov:8001";
  const std::string model_name = "dnn";

  const int nthreads = 1;
  std::vector<std::thread> threads;
  for (int i = 0; i < nthreads; i++) {
    threads.emplace_back(triton_forward, i, url, model_name);
  }
  for (int i = 0; i < nthreads; i++) {
    threads[i].join();
  }
  return 0;
}
