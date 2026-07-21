// #include "WireCellSpng/TorchTensorHandle.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>


int main(int argc, char * argv[]) {

  torch::Device device = torch::kCPU;
  if (argc > 1) {
    device = torch::kCUDA;
  }

  std::vector<float> data = {
      0,  1,  2,  3,
      4,  5,  6,  7,
      8,  9, 10, 11,
      12, 13, 14, 15
  };
  float * device_ptr;
  if (device == torch::kCUDA) {
    cudaMalloc((void **)&device_ptr, data.size() * sizeof(float));
    cudaMemcpy(device_ptr, &data[0], data.size() * sizeof(float),
               cudaMemcpyHostToDevice);
  }

  auto options = torch::TensorOptions()
                 .dtype(at::kFloat)
                 .requires_grad(false)
                 .device(device);


  torch::Tensor the_tensor = torch::from_blob(
    (device == torch::kCPU ? &data[0] : device_ptr),
    {4, 4},
    options
  );

  std::cout << the_tensor << std::endl;

  WireCell::SimpleTorchTensor simple_tensor(the_tensor);
  auto cloned1 = simple_tensor.tensor();
  the_tensor[0][0] = -999.;
  auto cloned2 = simple_tensor.tensor();
  std::cout << "cloned1\n" << cloned1 << std::endl;
  std::cout << "cloned2\n" << cloned2 << std::endl;
  
  std::cout << "Device" << simple_tensor.device()<< std::endl;
  std::cout << "Dtype"  << simple_tensor.dtype()<< std::endl;
  std::cout << "Shape: ";
  for (const auto & i : simple_tensor.shape()) std::cout << i << " ";
  std::cout << std::endl;
  
  return 0;

}