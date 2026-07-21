#include "WireCellSpng/TensorTools.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"



namespace WireCell::SPNG {

    void metadata_passthrough(
        const WireCell::Configuration & metadata_in,
        WireCell::Configuration & metadata_out,
        const Json::Value & passing_values) {
        // Throw
        // if (!passing_values.isArray()) {
        // }

        for (Json::ArrayIndex i = 0; i < passing_values.size(); ++i) {
            const auto & name = passing_values[i].asString();
            metadata_out[name] = metadata_in[name];
        }
    }

    

    ITorchTensorSet::pointer to_itensor(const std::vector<torch::IValue>& inputs) {
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
                }

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
    std::vector<torch::IValue> from_itensor(const ITorchTensorSet::pointer& in, bool is_gpu)
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

}
