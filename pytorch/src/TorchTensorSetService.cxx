#include "WireCellPytorch/TorchTensorSetService.h"

#include "WireCellAux/SimpleTensor.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Persist.h"

#include <cstdint>
#include <cstring>

WIRECELL_FACTORY(TorchTensorSetService,
                 WireCell::Pytorch::TorchTensorSetService,
                 WireCell::ITensorForward,
                 WireCell::IConfigurable)

using namespace WireCell;

torch::Dtype Pytorch::torch_dtype_of(const ITensor::pointer& ten)
{
    const auto& et = ten->element_type();
    if (et == typeid(float)) return torch::kFloat32;
    if (et == typeid(double)) return torch::kFloat64;
    if (et == typeid(int32_t)) return torch::kInt32;
    if (et == typeid(int64_t)) return torch::kInt64;
    THROW(ValueError() << errmsg{"TorchTensorSetService: unsupported element type " + ten->dtype()});
}

std::vector<torch::IValue> Pytorch::tensorset_to_ivalues(const ITensorSet::pointer& in, const torch::Device& dev)
{
    std::vector<torch::IValue> ret;
    for (const auto& ten : *in->tensors()) {
        const auto dt = torch_dtype_of(ten);
        std::vector<int64_t> shape;
        for (auto s : ten->shape()) shape.push_back((int64_t) s);
        // from_blob does not take ownership; the ITensor outlives the forward() call.  The module must not write
        // to its inputs (const_cast only because from_blob takes void*).
        auto t = torch::from_blob(const_cast<std::byte*>(ten->data()), shape, torch::TensorOptions().dtype(dt));
        if (dev.is_cuda()) {
            t = t.to(dev);
        }
        ret.push_back(t);
    }
    return ret;
}

namespace {
    template <typename T>
    ITensor::pointer copy_out(const torch::Tensor& t)
    {
        ITensor::shape_t shape;
        for (auto s : t.sizes()) shape.push_back((size_t) s);
        return std::make_shared<Aux::SimpleTensor>(shape, t.data_ptr<T>());
    }

    ITensor::pointer tensor_out(const torch::Tensor& tin)
    {
        auto t = tin.to(torch::kCPU).contiguous();
        switch (t.scalar_type()) {
            case torch::kFloat32: return copy_out<float>(t);
            case torch::kFloat64: return copy_out<double>(t);
            case torch::kInt32: return copy_out<int32_t>(t);
            case torch::kInt64: return copy_out<int64_t>(t);
            default: break;
        }
        THROW(ValueError() << errmsg{"TorchTensorSetService: unsupported output dtype"});
    }
}  // namespace

ITensor::vector Pytorch::ivalue_to_tensors(const torch::IValue& out)
{
    ITensor::vector ret;
    if (out.isTensor()) {
        ret.push_back(tensor_out(out.toTensor()));
    }
    else if (out.isTuple()) {
        for (const auto& e : out.toTupleRef().elements()) ret.push_back(tensor_out(e.toTensor()));
    }
    else if (out.isList()) {
        for (const auto& e : out.toListRef()) ret.push_back(tensor_out(e.toTensor()));
    }
    else if (out.isTensorList()) {
        for (const auto& t : out.toTensorVector()) ret.push_back(tensor_out(t));
    }
    else {
        THROW(ValueError() << errmsg{"TorchTensorSetService: module returned neither a tensor nor a tuple/list of tensors"});
    }
    return ret;
}

ITensorSet::pointer Pytorch::forward_tensorset(torch::jit::script::Module& module, const ITensorSet::pointer& in,
                                               const torch::Device& dev)
{
    torch::NoGradGuard no_grad;
    auto args = tensorset_to_ivalues(in, dev);
    torch::IValue out = module.forward(args);
    auto tens = ivalue_to_tensors(out);
    auto sv = std::make_shared<ITensor::vector>(tens.begin(), tens.end());
    return std::make_shared<Aux::SimpleTensorSet>(in->ident(), in->metadata(), sv);
}

Pytorch::TorchTensorSetService::TorchTensorSetService()
  : Aux::Logger("TorchTensorSetService", "torch")
{
}

Configuration Pytorch::TorchTensorSetService::default_configuration() const
{
    Configuration cfg;
    cfg["model"] = "model.ts";
    cfg["device"] = "cpu";
    return cfg;
}

void Pytorch::TorchTensorSetService::configure(const WireCell::Configuration& cfg)
{
    auto dev = get<std::string>(cfg, "device", "cpu");
    auto sem = get<std::string>(cfg, "semaphore", "");
    m_ctx.connect(dev, sem);

    m_model_path = Persist::resolve(get<std::string>(cfg, "model", ""));
    if (m_model_path.empty()) {
        log->critical("no TorchScript model file found for \"{}\"", get<std::string>(cfg, "model", ""));
        THROW(ValueError() << errmsg{"TorchTensorSetService: no TorchScript model file"});
    }
    torch::NoGradGuard no_grad;
    try {
        m_module = torch::jit::load(m_model_path, m_ctx.device());
    }
    catch (const c10::Error& e) {
        log->critical("error loading model \"{}\" to device \"{}\": {}", m_model_path, dev, e.what());
        throw;
    }
    m_module.eval();
    log->debug("loaded model \"{}\" to device \"{}\"", m_model_path, m_ctx.devname());
}

ITensorSet::pointer Pytorch::TorchTensorSetService::forward(const ITensorSet::pointer& in) const
{
    TorchSemaphore sem(m_ctx);
    try {
        return forward_tensorset(m_module, in, m_ctx.device());
    }
    catch (const std::runtime_error& err) {
        log->error("error running model \"{}\" on device \"{}\": {}", m_model_path, m_ctx.devname(), err.what());
        return nullptr;
    }
}
