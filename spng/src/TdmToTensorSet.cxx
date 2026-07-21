#include "WireCellSpng/TdmToTensorSet.h"
#include "WireCellSpng/TdmTools.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellAux/SimpleTensorSet.h"

#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGTdmToTensorSet,
                 WireCell::SPNG::TdmToTensorSet,
                 WireCell::SPNG::ITorchToTensorSet,
                 WireCell::IConfigurable,
                 WireCell::INamed);


namespace WireCell::SPNG {

    using HanaJsonCPP::to_json;
    using HanaJsonCPP::from_json;

    TdmToTensorSet::TdmToTensorSet()
        : Logger("TdmToTensorSet", "spng")
    {
    }

    TdmToTensorSet::~TdmToTensorSet()
    {
    }
    
    // IFunction
    bool TdmToTensorSet::operator()(const input_pointer& in, output_pointer& out)
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
            log->debug("including all {} tensors", included_tensors.size());
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
            log->debug("including select {} tensors", included_tensors.size());
        }


        std::unordered_set<ITorchTensor::pointer> excluded_tensors;
        for (const auto& er : m_config.exclude_rules) {
            auto some = TDM::select_tensors(included_tensors, er);
            const auto nexcluded = some.size();
            if (nexcluded) {
                excluded_tensors.insert(some.begin(), some.end());
                log->debug("excluding {} tensors with rule {}", nexcluded, er);
            }
        }
        log->debug("excluding {} tensors", excluded_tensors.size());
        

        ITorchTensor::vector final_tensors;
        for (const auto& one : included_tensors) {
            if (excluded_tensors.find(one) != excluded_tensors.end()) {
                logit(one, "excluding");
                continue;
            }
            logit(one, "including");
            final_tensors.push_back(one);
        }

        if (final_tensors.empty()) {
            log->warn("warning: no tensors survive");
        }

        // now convert....

        auto wct_tensors = std::make_shared<ITensor::vector>();
        for (const auto & one : final_tensors) {

            logit(one, "output");
            ITensor::pointer wct = TDM::tdm_to_wct(one);
            wct_tensors->push_back(wct);
        }

        out = std::make_shared<WireCell::Aux::SimpleTensorSet>(in->ident(), in->metadata(), wct_tensors);

        return true;
    }
    
    // IConfigurable
    void TdmToTensorSet::configure(const WireCell::Configuration& cfg)
    {
        this->ContextBase::configure(cfg);
        this->Logger::configure(cfg);
        from_json(m_config, cfg);
    }
    WireCell::Configuration TdmToTensorSet::default_configuration() const
    {
        auto cfg = this->ContextBase::default_configuration();
        auto cfg2 = this->Logger::default_configuration();
        update(cfg, cfg2);
        cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }

}
