#include "WireCellSpng/TorchFans.h"

#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"

#include "WireCellUtil/String.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGFaninTensorSets,
                 WireCell::SPNG::FaninTensorSets,
                 WireCell::SPNG::FaninTensorSets::fan_type,
                 WireCell::IConfigurable,
                 WireCell::INamed);
WIRECELL_FACTORY(SPNGFanoutTensorSets,
                 WireCell::SPNG::FanoutTensorSets,
                 WireCell::SPNG::FanoutTensorSets::fan_type,
                 WireCell::IConfigurable,
                 WireCell::INamed);

WIRECELL_FACTORY(SPNGFaninTensors,
                 WireCell::SPNG::FaninTensors,
                 WireCell::SPNG::FaninTensors::fan_type,
                 WireCell::IConfigurable,
                 WireCell::INamed);
WIRECELL_FACTORY(SPNGFanoutTensors,
                 WireCell::SPNG::FanoutTensors,
                 WireCell::SPNG::FanoutTensors::fan_type,
                 WireCell::IConfigurable,
                 WireCell::INamed);

namespace WireCell::SPNG {

    FaninTensorSets::FaninTensorSets(const std::string& type_name, const std::string& group_name)
        : FaninBase<ITorchTensorSet>(type_name, group_name)
    {
    }

    void FaninTensorSets::fanin_combine(const input_vector& inv, output_pointer& out)
    {
        const int ident = inv[0]->ident();
        Configuration md;
        auto tens = std::make_shared<ITorchTensor::vector>();
        for (const auto& its : inv) {
            auto new_md = its->metadata();
            update(md, new_md);
            for (const auto& ten : *its->tensors()) {
                tens->push_back(ten);
            }
        }
        out = std::make_shared<SimpleTorchTensorSet>(ident, md, tens);
        logit(out, "combined tensor set");
    }

    FanoutTensorSets::FanoutTensorSets(const std::string& type_name, const std::string& group_name)
        : FanoutBase<ITorchTensorSet>(type_name, group_name)
    {
    }

    void FanoutTensorSets::fanout_separate(const input_pointer& in, output_vector& outv)
    {
        logit(in, "separating tensor set");
        // Simply duplicate the pointer.
        for (size_t ind=0; ind<outv.size(); ++ind) {
            outv[ind] = in;
        }
    }


    // scalar


    FaninTensors::FaninTensors(const std::string& type_name, const std::string& group_name)
        : FaninBase<ITorchTensor>(type_name, group_name)
    {
    }

    void FaninTensors::fanin_combine(const input_vector& inv, output_pointer& out)
    {
        Configuration md;
        std::vector<torch::Tensor> tensors;
        for (const auto& in : inv) {
            update(md, in->metadata());
            tensors.push_back(in->tensor());
        }
        auto batched = torch::stack(tensors, 0);
        md["batches"] = (int) tensors.size();
        out = std::make_shared<SimpleTorchTensor>(batched, md);
        logit(out, "batched tensor");
    }

    FanoutTensors::FanoutTensors(const std::string& type_name, const std::string& group_name)
        : FanoutBase<ITorchTensor>(type_name, group_name)
    {
    }

    void FanoutTensors::fanout_separate(const input_pointer& in, output_vector& outv)
    {
        logit(in, "separating tensor");
        // Simply duplicate the pointer.
        for (size_t ind=0; ind<outv.size(); ++ind) {
            outv[ind] = in;
        }
    }
    
}
