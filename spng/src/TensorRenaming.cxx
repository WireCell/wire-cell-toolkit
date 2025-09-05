#include "WireCellSpng/TensorRenaming.h"
#include "WireCellSpng/SimpleTorchTensor.h"

namespace WireCell::SPNG {

    WireCell::Configuration TensorRenaming::default_configuration() const
    {
        Configuration cfg;
        cfg["tensor_renaming"] = Json::arrayValue;
        return cfg;
    }

    void TensorRenaming::configure(const WireCell::Configuration& cfg)
    {
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

    std::string TensorRenaming::rename(const std::string& datapath) const
    {
        for (const auto& rule : m_renaming_rules) {

            if (! std::regex_match(datapath, rule.match)) {
                continue;
            }

            return std::regex_replace(datapath, rule.match, rule.replace);
        }
        return datapath;
    }


    ITorchTensor::pointer TensorRenaming::rename_tensor(const ITorchTensor::pointer ten) const
    {
        std::string old_datapath = ten->metadata()["datapath"].asString();
        std::string new_datapath = rename(old_datapath);

        if (new_datapath == old_datapath) {
            return ten;
        }

        auto md = ten->metadata();
        md["datapath"] = new_datapath;
        return std::make_shared<SimpleTorchTensor>(ten->tensor(), md);
    }

    TensorIndex TensorRenaming::apply(const TensorIndex& index) const
    {
        TensorIndex ti(index.ident(), index.metadata());

        // Keep track of datapath changes to adjust any "parent" attributes.
        std::unordered_map<std::string, std::string> old_new;

        for (const auto& node : index.tree().depth()) {
            if (! node.value) { // skip empty root node
                continue;
            }

            ITorchTensor::pointer iten = node.value;
            auto md = iten->metadata();

            const std::string old_datapath = get<std::string>(md, "datapath");
            const std::string new_datapath = rename(old_datapath);

            old_new[old_datapath] = new_datapath;

            const std::string old_parent = get<std::string>(md, "parent", "");
            std::string new_parent = old_parent.empty() ? "" : old_new[old_parent];

            // replace if anything differs
            if (old_datapath == new_datapath && old_parent == new_parent) {
                ti.add(iten);
                continue;
            }

            md["datapath"] = new_datapath;
            if (! new_parent.empty()) {
                md["parent"] = new_parent;
            }
            
            iten = std::make_shared<SimpleTorchTensor>(iten->tensor(), md);
            ti.add(iten);
        }        
        return ti;
    }

}

