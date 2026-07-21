#ifndef WIRECELL_SPNG_TORCHTOOLS
#define WIRECELL_SPNG_TORCHTOOLS

#include "WireCellSpng/ITorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"

#include "WireCellUtil/Exceptions.h"

namespace WireCell::SPNG {

    /// Assure tensors are on the device.
    template<typename Type>
    Type to_device(const Type& data, torch::Device device) {
        raise<TypeError>("to_device not implemented for this type");
        return data;            // quell compiler warning.
    }
    
    template<> inline
    torch::Tensor to_device<torch::Tensor>(const torch::Tensor& tensor, torch::Device device)
    {
        if (tensor.device() == device) {
            return tensor;
        }
        return tensor.to(device);
    }

    template<> inline
    std::vector<torch::Tensor> to_device<std::vector<torch::Tensor>>(const std::vector<torch::Tensor>& tenvec, torch::Device device)
    {
        std::vector<torch::Tensor> ret;
        for (const auto& ten : tenvec) {
            ret.push_back(to_device(ten, device));
        }
        return ret;
    }


    template<> inline
    ITorchTensor::pointer to_device<ITorchTensor::pointer>(const ITorchTensor::pointer& itensor, torch::Device device)
    {
        auto ten = itensor->tensor();
        if (ten.device() == device) {
            return itensor;
        }
        return std::make_shared<SimpleTorchTensor>(to_device(ten, device), itensor->metadata());
    }
    
    template<> inline
    ITorchTensor::vector to_device<ITorchTensor::vector>(const ITorchTensor::vector& itenvec, torch::Device device)
    {
        ITorchTensor::vector ret;
        for (const auto& iten : itenvec) {
            ret.push_back(to_device(iten, device));
        }
        return ret;
    }
    

    template<> inline
    ITorchTensorSet::pointer to_device<ITorchTensorSet::pointer>(const ITorchTensorSet::pointer& itenset, torch::Device device)
    {
        ITorchTensor::vector tenvec;
        for (const auto& iten : *itenset->tensors()) {
            tenvec.push_back(to_device(iten, device));
        }
        return std::make_shared<SimpleTorchTensorSet>(itenset->ident(), itenset->metadata(), tenvec);
    }
    
    template<> inline
    ITorchTensorSet::vector to_device<ITorchTensorSet::vector>(const ITorchTensorSet::vector& itensetvec, torch::Device device)
    {
        ITorchTensorSet::vector ret;
        for (const auto& itenset : itensetvec) {
            ret.push_back(to_device(itenset, device));
        }
        return ret;
    }
    

}


#include "WireCellUtil/Configuration.h"


namespace WireCell::SPNG {

    // FIXME: Remove.  This already exists Configuration.h as update().  It also
    // is not generic but defines specific data model policy.
  void metadata_passthrough(
    const WireCell::Configuration & metadata_in,
    WireCell::Configuration & metadata_out,
    const Json::Value & passing_values);

    // FIXME:  Remove, eventually.  Not generic, obsoleted by TDM
    std::vector<torch::IValue> from_itensor(const ITorchTensorSet::pointer& in, bool is_gpu = false);
    ITorchTensorSet::pointer to_itensor(const std::vector<torch::IValue>& inputs);
}

#endif


