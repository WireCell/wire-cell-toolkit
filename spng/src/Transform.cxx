#include "WireCellSpng/Transform.h"
#include "WireCellSpng/SimpleTorchTensor.h"


#include "WireCellUtil/NamedFactory.h"

WIRECELL_FACTORY(SPNGTransform,
                 WireCell::SPNG::Transform,
                 WireCell::ITorchTensorFilter,
                 WireCell::IConfigurable,
                 WireCell::INamed)

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

        logit(in, "input");

        torch::Tensor tensor = in->tensor();
        for (auto op : m_ops) {
            tensor = op(tensor);
        }
        // Fixme: TDM MD handling still needs thought
        auto md = in->metadata();
        md["datapath"] = md["datapath"].asString() + "/Transform/" + get_name();

        out = std::make_shared<SimpleTorchTensor>(tensor, md);

        logit(out, "output");

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
            if (opcfg.operation == "normalize") {
                m_ops.push_back([](const torch::Tensor& tensor) -> torch::Tensor {
                    torch::Tensor vmin = tensor.min();
                    torch::Tensor vmax = tensor.max();

                    // Special case that all values are same
                    if (vmax.item<float>() == vmin.item<float>()) {
                        return torch::full_like(tensor, 0.5); 
                    }

                    return (tensor - vmin) / (vmax - vmin);
                });
                continue;
            }
            if (opcfg.operation == "medsub") {
                m_ops.push_back([dims](const torch::Tensor& tensor) -> torch::Tensor {
                    torch::Tensor medians = std::get<0>(torch::median(tensor, dims[0], true));
                    return tensor.sub(medians);
                });
                continue;
            }
            if (opcfg.operation == "lowmedsub") {
                m_ops.push_back([dims, scalar](const torch::Tensor& tensor) -> torch::Tensor {
                    auto temp = tensor.clone();
                    auto mask = temp.gt(scalar);
                    temp.masked_fill_(mask, std::numeric_limits<float>::quiet_NaN());
                    temp = std::get<0>(torch::nanmedian(temp, dims[0], /*keepdim=*/true));
                    return tensor.sub(temp);
                });
                continue;
            }
            if (opcfg.operation == "treshold") {
                m_ops.push_back([scalar](const torch::Tensor& tensor) -> torch::Tensor {
                    return tensor > scalar;
                });
                continue;
            }
            if (opcfg.operation == "slice") {
                m_ops.push_back([dims](const torch::Tensor& tensor) -> torch::Tensor {
                    if (dims.size() != 3) {
                        raise<ValueError>("slice transform requires 3 values in dims config, got %d",
                                          dims.size());
                    }
                    return tensor.slice(dims[0], dims[1], dims[2]);
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

