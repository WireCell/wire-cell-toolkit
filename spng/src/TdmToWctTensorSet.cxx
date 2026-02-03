#include "WireCellSpng/TdmToWctTensorSet.h"
#include "WireCellSpng/TdmTools.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellAux/SimpleTensorSet.h"

namespace WireCell::SPNG {

    using HanaJsonCPP::to_json;
    using HanaJsonCPP::from_json;

    TdmToWctTensorSet::TdmToWctTensorSet()
        : Logger("TdmToWctTensorSet", "spng")
    {
    }

    TdmToWctTensorSet::~TdmToWctTensorSet()
    {
    }
    
    // IFunction
    bool TdmToWctTensorSet::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }

        logit(in, "input");
        
        ITorchTensor::vector included_tensors;
        if (m_config.include_rules.empty()) {
            included_tensors = *in->tensors();
        }
        else {
            std::unordered_set<ITorchTensor::pointer> seen;
            for (const auto& ir : m_config.include_rules) {
                auto some = TDM::select_tensors(*in->tensors(), ir);
                for (const auto& one : some) {
                    if (seen.find(one) == seen.end()) {
                        continue;
                    }
                    seen.insert(one);
                    included_tensors.push_back(one);
                }
            }
        }


        std::unordered_set<ITorchTensor::pointer> excluded_tensors;
        for (const auto& er : m_config.exclude_rules) {
            auto some = TDM::select_tensors(included_tensors, er);
            excluded_tensors.insert(some.begin(), some.end());
        }
        

        ITorchTensor::vector final_tensors;
        for (const auto& one : included_tensors) {
            if (excluded_tensors.find(one) == excluded_tensors.end()) {
                continue;
            }
            final_tensors.push_back(one);
        }

        // now convert....

        auto wct_tensors = std::make_shared<ITensor::vector>();
        for (const auto & one : final_tensors) {

            ITensor::pointer wct = TDM::tdm_to_wct(one);
            wct_tensors->push_back(wct);
        }

        out = std::make_shared<WireCell::Aux::SimpleTensorSet>(in->ident(), in->metadata(), wct_tensors);

        return true;
    }
    
    // IConfigurable
    void TdmToWctTensorSet::configure(const WireCell::Configuration& cfg)
    {
        this->ContextBase::configure(cfg);
        this->Logger::configure(cfg);
        from_json(m_config, cfg);
    }
    WireCell::Configuration TdmToWctTensorSet::default_configuration() const
    {
        auto cfg = this->ContextBase::default_configuration();
        auto cfg2 = this->Logger::default_configuration();
        update(cfg, cfg2);
        cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }

}
