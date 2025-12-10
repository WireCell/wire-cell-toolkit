#include "WireCellSpng/Transform.h"
#include "WireCellSpng/SimpleTorchTensor.h"

using namespace WireCell::HanaJsonCPP;                         

namespace WireCell::SPNG {


    Transform::Transform()
        : Logger("Transform", "spng")
    {
    }

    bool Transform::operator()(const input_pointer& in, output_pointer& out)
    {
        out = nullptr;
        if (!in) {
            logit("EOS");
            next_count();
            return true;
        }

        torch::Tensor tensor = in->tensor();
        for (auto op : m_ops) {
            tensor = op(tensor);
        }
        out = std::make_shared<SimpleTorchTensor>(tensor, in->metadata());
        next_count();
        return true;
    }

    void Transform::configure(const WireCell::Configuration& jconfig)
    {
        WireCell::configure_bases<Transform, ContextBase, Logger>(this, jconfig);
        from_json(m_config, jconfig);

        m_ops.clear();
        for (const auto& opcfg: m_config.operations) {

            std::vector<int64_t> dims(opcfg.dims.begin(), opcfg.dims.end());
            const float scalar = opcfg.scalar;
        
            if (opcfg.operation == "permute") {
                // dims size must match tensor shape size, leave it to torch to check
                if (dims.size() < 2) {
                    raise<ValueError>("nonsensical to permute less than two dimensions");
                }
                m_ops.push_back([dims](const torch::Tensor& tensor) -> torch::Tensor {
                    return tensor.permute(dims);
                });
                continue;
            }
            if (opcfg.operation == "transpose") {
                if (dims.size() != 2) {
                    raise<ValueError>("wrong size sims for transpose");
                }
                m_ops.push_back([dims](const torch::Tensor& tensor) -> torch::Tensor {
                    return torch::transpose(tensor, dims[0], dims[1]);
                });
                
                continue;
            }
            if (opcfg.operation == "scale") {
                m_ops.push_back([scalar](const torch::Tensor& tensor) -> torch::Tensor {
                    return torch::mul(tensor, scalar);
                });
                continue;
            }
            if (opcfg.operation == "offset") {
                m_ops.push_back([scalar](const torch::Tensor& tensor) -> torch::Tensor {
                    return tensor.add(scalar);
                });
                continue;
            }
            raise<ValueError>("unknown operation: %s", opcfg.operation);
        }

    }

    WireCell::Configuration Transform::default_configuration() const
    {
        auto cfg = WireCell::default_configuration_bases<Transform, ContextBase, Logger>(this);
        auto cfg2 = to_json(m_config);
        update(cfg, cfg2);
        return cfg;
    }

}

