#include "WireCellSpng/Reduce.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/Torch.h"

using WireCell::HanaJsonCPP::from_json;
using WireCell::HanaJsonCPP::to_json;

namespace WireCell::SPNG {

    Reduce::Reduce()
    {
    }

    void Reduce::fanin_combine(const input_vector& inv, output_pointer& out)
    {
        Configuration md;
        std::vector<torch::Tensor> inputs;
        const size_t nin = inv.size();
        for (size_t ind=0; ind<nin; ++ind) {
            md = update(md, inv[ind]->metadata());
            inputs.push_back(inv[ind]->tensor());
        }
        auto reduced = m_op(inputs);
        out = std::make_shared<SimpleTorchTensor>(reduced, md);
        logit(out, "reduce");
    }

    WireCell::Configuration Reduce::default_configuration() const
    {
        auto cfg = FaninBase<ITorchTensor>::default_configuration();
        update(cfg, HanaJsonCPP::to_json(m_config));
        return cfg;
    }

    // Define the signature for the dynamic operation (the class member)
    using ReduceOperation = std::function<torch::Tensor(const std::vector<torch::Tensor>&)>;

    // Define the signature for the element-wise binary operation (e.g., torch::max(A, B))
    using BinaryOp = std::function<torch::Tensor(const torch::Tensor&, const torch::Tensor&)>;

    // Some torch functions like min(a,b) or max(a,b) do not take a vector of
    // tensors.  This helper fills that gap.
    auto iterative_reduction = [](BinaryOp op) {
        return [op](const std::vector<torch::Tensor>& tensors) {
            // Initialize result with the first tensor
            torch::Tensor result = tensors[0];
            // Iterate from the second tensor
            for (size_t i = 1; i < tensors.size(); ++i) {
                // Apply the passed binary operation (op)
                result = op(result, tensors[i]);
            }
            return result;
        };
    };

    void Reduce::configure(const WireCell::Configuration& jconfig) {
        FaninBase<ITorchTensor>::configure(jconfig);
        from_json(m_config, jconfig);

        // Functions requiring 'dim' (Combination/Stacking)

        if (m_config.operation == "cat") {
            int captured_dim = m_config.dim;
            m_op = [captured_dim](const std::vector<torch::Tensor>& tensors) {
                return torch::cat(tensors, captured_dim);
            };
        }

        else if (m_config.operation == "stack") {
            int captured_dim = m_config.dim;
            m_op = [captured_dim](const std::vector<torch::Tensor>& tensors) {
                return torch::stack(tensors, captured_dim);
            };
        }
        
        // Functions ignoring 'dim' (Element-wise Reduction or specific stacking)
        else if (m_config.operation == "max") {
            m_op = iterative_reduction(
                [](const torch::Tensor& a, const torch::Tensor& b) { return torch::max(a, b); }
            );
        }

        else if (m_config.operation == "min") {
            m_op = iterative_reduction(
                [](const torch::Tensor& a, const torch::Tensor& b) { return torch::min(a, b); }
            );
        }

        else if (m_config.operation == "sum") {
            m_op = [](const std::vector<torch::Tensor>& tensors) {
                auto s = torch::stack(tensors, 0);
                return s.sum(0);
            };
        }

        else if (m_config.operation == "mean") {
            m_op = [](const std::vector<torch::Tensor>& tensors) {
                auto s = torch::stack(tensors, 0);
                s = s.sum(0);
                return s / tensors.size();
            };
        }

        else if (m_config.operation == "mul") {
            m_op = iterative_reduction(
                [](const torch::Tensor& a, const torch::Tensor& b) { return torch::mul(a, b); }
            );
        }

        else if (m_config.operation == "hstack") {
            m_op = [](const std::vector<torch::Tensor>& tensors) {
                return torch::hstack(tensors);
            };
        }

        else if (m_config.operation == "vstack") {
            m_op = [](const std::vector<torch::Tensor>& tensors) {
                return torch::vstack(tensors);
            };
        }

        else if (m_config.operation == "dstack") {
            m_op = [](const std::vector<torch::Tensor>& tensors) {
                return torch::dstack(tensors);
            };
        }

        else {
            raise<ValueError>("unsupported reduce operation: %s", m_config.operation);
        }
    }


}
