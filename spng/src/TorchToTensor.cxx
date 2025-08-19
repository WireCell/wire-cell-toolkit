#include "WireCellSpng/TorchToTensor.h"
#include "WireCellUtil/String.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellAux/SimpleTensorSet.h"
#include "WireCellAux/SimpleTensor.h"

WIRECELL_FACTORY(TorchToTensor, WireCell::SPNG::TorchToTensor,
                 WireCell::INamed,
                 WireCell::SPNG::ITorchToTensorSet,
                 WireCell::IConfigurable)

using namespace WireCell;
using namespace WireCell::Aux;
using namespace WireCell::SPNG;

TorchToTensor::TorchToTensor()
  : Aux::Logger("TorchToTensor", "spng")
{}

TorchToTensor::~TorchToTensor()
{}


bool TorchToTensor::operator()(const input_pointer& torchset_in, output_pointer& its)
{
    //Setup output + basic check
    its = nullptr;
    if (! torchset_in) {
        // eos
        return true;
    }

    // log->debug("call={} ", m_count++);

    //Go through each input tensor
    ITensor::vector tensor_vec;
    for (const auto & tensor_ptr : (*torchset_in->tensors())) {
        //Make sure to move to CPU
        auto tensor_clone = tensor_ptr->tensor().clone().to(torch::kCPU).contiguous();
        // const auto shape = ten->shape();
    
        std::vector<double> as_vec(
            tensor_clone.data_ptr<double>(),
            tensor_clone.data_ptr<double>() + tensor_clone.numel());

        //Make a torch::Tensor from the tensor data
        auto sizes = tensor_ptr->tensor().sizes();

        //Make shared ptr using the tensor's metadata
        tensor_vec.emplace_back(std::make_shared<SimpleTensor>(
            WireCell::ITensor::shape_t(sizes.begin(), sizes.end()),
            &(as_vec[0]), tensor_ptr->metadata()
        ));
    }

    //Turn into shared TorchTensor vector
    auto sv = std::make_shared<ITensor::vector>(
        tensor_vec.begin(), tensor_vec.end());
    
    //Finally, turn into TorchTensor set w/ ident + metadata from input set
    its = std::make_shared<SimpleTensorSet>(
        torchset_in->ident(), torchset_in->metadata(), sv);

    return true;
}

WireCell::Configuration TorchToTensor::default_configuration() const
{
    Configuration cfg;
 
    return cfg;
}
void TorchToTensor::configure(const WireCell::Configuration& cfg)
{
    // m_datapath = get(cfg, "datapath", m_datapath);
    // m_baseline = get(cfg, "baseline", m_baseline);
    // m_offset = get(cfg, "offset", m_offset);
    // m_scale = get(cfg, "scale", m_scale);
    // m_transform = [&](float x) -> float {
    //     return (x + m_baseline) * m_scale + m_offset;
    // };

    // std::map<std::string, TorchToTensorMode> modes = {
    //     {"unified",  TorchToTensorMode::unified},
    //     {"tagged", TorchToTensorMode::tagged },
    //     {"sparse", TorchToTensorMode::sparse }};
    // auto mode = get<std::string>(cfg, "mode", "tagged");
    // m_mode = modes[mode];
    // m_digitize = get(cfg, "digitize", m_digitize);
}
