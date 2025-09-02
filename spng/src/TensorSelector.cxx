#include "WireCellSpng/TensorSelector.h"
#include "WireCellSpng/SimpleTorchTensor.h"

namespace WireCell::SPNG {

    WireCell::Configuration TensorSelector::default_configuration() const
    {
        Configuration cfg;
        cfg["tensor_selection"] = Json::arrayValue;
        cfg["tensor_renaming"] = Json::arrayValue;
        return cfg;
    }

    void TensorSelector::configure(const WireCell::Configuration& cfg)
    {
        // pack into regex
        m_selection_rules.clear();
        for (auto jrule : cfg["tensor_selection"]) {
            SelectionRule rule;
            int count = 0;
            {
                auto jaccept = jrule["accept"];
                if (jaccept.isString()) {
                    rule.accept = jaccept.asString();
                    ++count;
                }
            }
            {
                auto jreject = jrule["reject"];
                if (jreject.isString()) {
                    rule.reject = jreject.asString();
                    ++count;
                }
            }
            if (! count) { continue; } // skip empty
            m_selection_rules.push_back(rule);
        }

        m_renaming_rules.clear();
        for (auto jrule : cfg["tensor_renaming"]) {
            RenamingRule rule;
            {
                auto jmatch = jrule["match"];
                if (! jmatch.isString()) { continue; }
                rule.match = jmatch.asString();
            }
            {
                auto jreplace = jrule["replace"];
                if (! jreplace.isString()) { continue; }
                rule.replace = jreplace.asString();
            }
            m_renaming_rules.push_back(rule);
        }
    }

    TensorSelector::SelectionResult TensorSelector::select_tensor(const ITorchTensor::pointer ten) const
    {
        const std::string datapath = ten->metadata()["datapath"].asString();

        for (const auto& rule : m_selection_rules) {
            if (std::regex_match(datapath, rule.accept)) {
                return SelectionResult::kAccept;
            }
            if (std::regex_match(datapath, rule.reject)) {
                return SelectionResult::kReject;
            }
        }
        return SelectionResult::kNoMatch;
    }

    ITorchTensor::pointer TensorSelector::rename_tensor(const ITorchTensor::pointer ten) const
    {
        std::string datapath = ten->metadata()["datapath"].asString();

        for (const auto& rule : m_renaming_rules) {

            // Fixme, this is a "try then do" pattern.  It is maybe better to
            // just "do" and then check for a change in the resulting string?

            if (! std::regex_match(datapath, rule.match)) {
                continue;
            }

            auto md = ten->metadata();
            md["datapath"] = std::regex_replace(datapath, rule.match, rule.replace);
            return std::make_shared<SimpleTorchTensor>(ten->tensor(), md);
        }
        return nullptr;
    }
}

