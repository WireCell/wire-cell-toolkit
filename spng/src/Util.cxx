#include "WireCellSpng/Util.h"
#include <cmath>
#include <sstream>


namespace WireCell::Torch {

    std::string to_string(const torch::Device& device)
    {
        std::stringstream ss;
        ss << device;
        return ss.str();        
    }

}

using namespace WireCell;

torch::Tensor Torch::gaussian1d(double mean, double sigma,
                                int64_t npoints, double xmin, double xmax,
                                torch::TensorOptions options)
{
    auto x = torch::linspace(xmin, xmax, npoints, options);
    auto rel = (x - mean)/sigma;
    const double norm = sqrt(2*M_PI*sigma*sigma);
    return norm * torch::exp(-0.5*rel*rel);
}


std::vector<int64_t> Torch::linear_shape(const std::vector<torch::Tensor>& tens, 
                                         torch::IntArrayRef extra_shape)
{
    // Find shape that assures linear convolution.
    std::vector<int64_t> shape = {extra_shape[0], extra_shape[1]};
    for (const auto& ten : tens) {
        auto sizes = ten.sizes();
        shape[0] += sizes[0] - 1;
        shape[1] += sizes[1] - 1;
    }
    return shape;
}


torch::Tensor Torch::pad(torch::Tensor ten, double value, torch::IntArrayRef shape)
{
    using torch::indexing::Slice;

    torch::Tensor padded = torch::zeros(shape, ten.options()) + value;
    auto s = ten.sizes();
    padded.index_put_({
            Slice(0,std::min(s[0], shape[0])),
            Slice(0,std::min(s[1], shape[1]))
        }, ten);
    return padded;    
}

torch::Tensor Torch::convo_spec(const std::vector<torch::Tensor>& tens, 
                                torch::IntArrayRef shape)
{
    using torch::indexing::Slice;

    // Return value will be complex spectrum in Fourier domain.
    torch::Tensor fourier;
    // Allocate working array in interval domain.
    torch::Tensor interval = torch::zeros(shape, tens[0].options());

    const size_t ntens = tens.size();

    // First, accumulate denominator.
    for (size_t ind=0; ind<ntens; ++ind) {
        // Caveat: zero-padding is not always appropriate for every ten in tens.
        interval.zero_();

        const auto& ten = tens[ind];
        auto s = ten.sizes();
        interval.index_put_({
                Slice(0,std::min(s[0], shape[0])),
                Slice(0,std::min(s[1], shape[1]))
            }, ten);

        if (ind == 0) {
            fourier = torch::fft::fft2(interval);
        }
        else {
            fourier *= torch::fft::fft2(interval);
        }
    }

    return fourier;
}


torch::Tensor Torch::filtered_decon_2d(const std::vector<torch::Tensor>& numerator,
                                       const std::vector<torch::Tensor>& denominator,
                                       torch::IntArrayRef shape)
{
    if (denominator.empty()) {  // graceful degradation 
        return convo_spec(numerator, shape);
    }

    // Note: this suffers holding an extra array which could be avoided by
    // accumulating denominator convolutions, inverting that result and
    // continuing accumulating numerator convolutions.
    auto num = convo_spec(numerator, shape);
    auto den = convo_spec(denominator, shape);
    return torch::divide(num, den); // fixme, divide-by-zero?
}


torch::Tensor Torch::filtered_decon_2d_auto(const std::vector<torch::Tensor>& numerator,
                                          const std::vector<torch::Tensor>& denominator,
                                          torch::IntArrayRef extra_shape)
{
    std::vector<torch::Tensor> all_in(numerator.begin(), numerator.end());
    all_in.insert(all_in.end(), denominator.begin(), denominator.end());

    if (all_in.empty()) {
        return torch::Tensor();
    }

    auto shape = linear_shape(all_in, extra_shape);

    return filtered_decon_2d(numerator, denominator, shape);
}

void WireCell::SPNG::metadata_passthrough(
    const WireCell::Configuration & metadata_in,
    WireCell::Configuration & metadata_out,
    const Json::Value & passing_values) {
        // Throw
        // if (!passing_values.isArray()) {
        // }

        // std::cout << "Passing values:" << passing_values << std::endl;
        for (Json::ArrayIndex i = 0; i < passing_values.size(); ++i) {
            const auto & name = passing_values[i].asString();
            // std::cout << "Passing " << name << std::endl;
            // std::cout << "In input? " << metadata_in.isMember(name) << std::endl;
            metadata_out[name] = metadata_in[name];
        }
}

    void WireCell::SPNG::save_torchtensor_data(const  torch::Tensor& tensor, const std::string& filename) {
        //write the torch::Tensor in a file using torch::save
        try {
            torch::serialize::OutputArchive archive;
            archive.write("tensor", tensor.cpu());
            archive.save_to(filename);
        } catch (const c10::Error& e) {
            std::cerr << "Util::save_torchtensor_data  Error saving tensor: " << e.what() << std::endl;
        }
    }  
    
    void WireCell::SPNG::save_simpletensor_data(const ITorchTensorSet::pointer& in, const std::string& filename) {
        // Create a map<string, torch::Tensor> to hold all tensors with unique keys
        torch::serialize::OutputArchive archive;
    
        int idx = 0;
        for (const auto& tensor : *in->tensors()) {
            auto ten = tensor->tensor();
    

    
            // Move tensor to CPU if needed for portability
            torch::Tensor cpu_tensor = ten.device().is_cpu() ? ten : ten.to(torch::kCPU);
    
            // Store with a unique key
            std::string key = "tensor_" + std::to_string(idx++);
            archive.write(key, cpu_tensor);
        }
    
        try {
            // Save the entire map of tensors in one file
            archive.save_to(filename);
        } catch (const c10::Error& e) {
            std::cerr << "Util::save_simpletensor_data  Error saving tensors: " << e.what() << std::endl;
        }
    }
    

    
    

ITorchTensorSet::pointer WireCell::SPNG::to_itensor(const std::vector<torch::IValue>& inputs) {
    auto itv = std::make_shared<ITorchTensor::vector>();

    for (size_t i = 0; i < inputs.size(); ++i) {
        try {
            const auto& ivalue = inputs[i];
            if (!ivalue.isTensor()) {
                std::cerr << "Error: Expected torch::IValue at index " << i << " to be a Tensor\n";
                continue;
            }
            torch::Tensor ten = ivalue.toTensor();
            if (ten.dim() != 4) {
                std::cerr << "Error: Tensor at index " << i << " must be 4D, got " << ten.dim() << std::endl;
                continue;
            }

            if(ten.scalar_type() != torch::kFloat32) {
                ten = ten.to(torch::kFloat32);
                std::cout << "Converted tensor " << i << " to float32." << std::endl;
            }

            std::cout << "Tensor " << i << ": shape=" << ten.sizes() 
            << ", dtype=" << ten.dtype() 
            << ", device=" << ten.device() << std::endl;

            // No forced dtype or device conversion here to preserve input tensors as-is
            auto stp = std::make_shared<SimpleTorchTensor>(ten);
            itv->emplace_back(stp);

        } catch (const std::exception& e) {
            std::cerr << "Exception caught while processing tensor " << i << ": " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Unknown exception caught while processing tensor " << i << std::endl;
        }
    }
    return std::make_shared<SimpleTorchTensorSet>(0, Json::nullValue, itv);
}

//ITorchTensor --> torch::IValue
std::vector<torch::IValue> WireCell::SPNG::from_itensor(const ITorchTensorSet::pointer& in, bool is_gpu)
{
    // Create a new SimpleTorchTensorSet to hold the converted tensors
    std::vector<torch::IValue> ret;
    //Populate this function as needed...

    for(auto iten: *in->tensors()) {
        // Convert each tensor to IValue
        torch::Tensor ten = iten->tensor();
        //why casting to float?
        if(ten.scalar_type() != torch::kFloat32) {
            ten = ten.to(torch::kFloat32);
        }
        if (is_gpu) {
            ten = ten.to(torch::Device(torch::kCUDA, 0));
            assert(ten.device().type() == torch::kCUDA);
        } 
        ret.emplace_back(ten);  
    }
    return ret;
}

std::string WireCell::SPNG::tensor_shape_string(const torch::Tensor& t) {
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
