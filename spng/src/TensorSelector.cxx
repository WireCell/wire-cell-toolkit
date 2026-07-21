#include "WireCellSpng/TensorSelector.h"
#include "WireCellSpng/SimpleTorchTensor.h"

namespace WireCell::SPNG {

    WireCell::Configuration TensorSelector::default_configuration() const
    {
        Configuration cfg;
        cfg["tensor_selection"] = Json::arrayValue;
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
                    auto accept = jaccept.asString();
                    rule.accept = std::regex(accept);
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

    }

    TensorSelector::SelectionResult TensorSelector::select_tensor(ITorchTensor::pointer ten) const
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

    TensorIndex TensorSelector::apply(const TensorIndex& index,
                                      bool keep_unselected, bool select_parents) const
    {
        if (select_parents) return apply_parents(index, keep_unselected);
        return apply_all(index, keep_unselected);
    }

    bool TensorSelector::apply_one(ITorchTensor::pointer iten,
                                   bool keep_unselected, const std::string& kind) const
    {
        auto res = select_tensor(iten);

        if (res == SelectionResult::kReject) {
            return false;
        }

        if (res == SelectionResult::kAccept) {
            return true;
        }

        if (res == SelectionResult::kNoMatch && keep_unselected) {
            return true;
        }
        return false;
    }

    TensorIndex TensorSelector::apply_parents(const TensorIndex& index, bool keep_unselected) const
    {
        TensorIndex ti(index.ident(), index.metadata());
        for (auto iten : index.tree().child_values()) { // top level parents

            if (apply_one(iten, keep_unselected, "parent")) {
                const auto* node = index.tree_node(iten);
                ti.add(*node);
            }

        }        
        return ti;
    }

    TensorIndex TensorSelector::apply_all(const TensorIndex& index, bool keep_unselected) const
    {

        TensorIndex ti(index.ident(), index.metadata());
        for (const auto& node : index.tree().depth()) {
            auto iten = node.value;

            if (! iten) { // skip empty root node
                continue;
            }
            
            if (apply_one(iten, keep_unselected, "tensor")) {
                ti.add(iten);
            }

        }        
        return ti;
    }
}

