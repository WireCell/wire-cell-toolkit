#include "WireCellSpng/TensorToTorch.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/SimpleTorchTensor.h"

WIRECELL_FACTORY(TensorToTorch, WireCell::SPNG::TensorToTorch,
                 WireCell::INamed,
                 WireCell::SPNG::ITensorToTorchSet,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::SPNG;

TensorToTorch::TensorToTorch()
  : Aux::Logger("TensorToTorch", "spng")
{}

TensorToTorch::~TensorToTorch()
{}


bool TensorToTorch::operator()(const input_pointer& tensorset_in, output_pointer& its)
{
    //Setup output + basic check
    its = nullptr;
    if (! tensorset_in) {
        // eos
        return true;
    }

    // log->debug("call={} ", m_count++);

    //Go through each input tensor
    ITorchTensor::vector torch_tensor_vec;
    for (const auto & tensor_ptr : (*tensorset_in->tensors())) {

        //Make a torch::Tensor from the tensor data
        auto torch_tensor = torch::from_blob(
            (void*)tensor_ptr->data(), tensor_ptr->size()
        );

        //Make shared ptr using the tensor's metadata
        torch_tensor_vec.emplace_back(std::make_shared<SimpleTorchTensor>(
            torch_tensor, tensor_ptr->metadata()
        ));
    }

    //Turn into shared TorchTensor vector
    auto sv = std::make_shared<ITorchTensor::vector>(
        torch_tensor_vec.begin(), torch_tensor_vec.end());
    
    //Finally, turn into TorchTensor set w/ ident + metadata from input set
    its = std::make_shared<SimpleTorchTensorSet>(
        tensorset_in->ident(), tensorset_in->metadata(), sv);

    return true;
}

WireCell::Configuration TensorToTorch::default_configuration() const
{
    Configuration cfg;
 
    return cfg;
}
void TensorToTorch::configure(const WireCell::Configuration& cfg)
{
    // m_datapath = get(cfg, "datapath", m_datapath);
    // m_baseline = get(cfg, "baseline", m_baseline);
    // m_offset = get(cfg, "offset", m_offset);
    // m_scale = get(cfg, "scale", m_scale);
    // m_transform = [&](float x) -> float {
    //     return (x + m_baseline) * m_scale + m_offset;
    // };

    // std::map<std::string, TensorToTorchMode> modes = {
    //     {"unified",  TensorToTorchMode::unified},
    //     {"tagged", TensorToTorchMode::tagged },
    //     {"sparse", TensorToTorchMode::sparse }};
    // auto mode = get<std::string>(cfg, "mode", "tagged");
    // m_mode = modes[mode];
    // m_digitize = get(cfg, "digitize", m_digitize);
}
