#include "WireCellSpng/ThresholdROIs.h"
#include "WireCellUtil/NamedFactory.h"
#include "WireCellUtil/Exceptions.h"
#include "WireCellSpng/SimpleTorchTensor.h"
#include "WireCellSpng/SimpleTorchTensorSet.h"
#include "WireCellSpng/Util.h"

// #include "WireCellSpng/ITorchFieldResponse.h"

// #include "WireCellSpng/ITorchColdElecResponse.h"

WIRECELL_FACTORY(SPNGThresholdROIs, WireCell::SPNG::ThresholdROIs,
                 WireCell::INamed,
                 WireCell::ITorchTensorSetFilter, WireCell::IConfigurable)

WireCell::SPNG::ThresholdROIs::ThresholdROIs()
  : Aux::Logger("SPNGThresholdROIs", "spng") {

}

WireCell::SPNG::ThresholdROIs::~ThresholdROIs() {};


void WireCell::SPNG::ThresholdROIs::configure(const WireCell::Configuration& config) {

    m_passthrough.append("channel_map");

    m_output_set_tag = get(config, "output_set_tag", m_output_set_tag);
    m_output_tensor_tag = get(config, "output_tensor_tag", m_output_tensor_tag);
    m_unsqueeze_input = get(config, "unsqueeze_input", m_unsqueeze_input);
    m_threshold_rms_factor = get(config, "threshold_rms_factor", m_threshold_rms_factor);
    log->debug("Will tag with Set:{} Tensor:{}", m_output_set_tag.asString(),
               m_output_tensor_tag.asString());
    log->debug("Will find ROIs with threshold_rms_factor {}", m_threshold_rms_factor);
}

bool WireCell::SPNG::ThresholdROIs::operator()(const input_pointer& in, output_pointer& out) {
    out = nullptr;
    if (!in) {
        log->debug("EOS ");
        return true;
    }
    log->debug("Running ThresholdROIs");

    // //Get the cloned tensor from the input
    // size_t ntensors = in->tensors()->size();
    // log->debug("Got {} tensors", ntensors);
    // std::vector<int64_t> prev_shape = {-999, -999};
    // for (size_t i = 0; i < ntensors; ++i) {
    //     auto this_sizes = in->tensors()->at(i)->tensor().sizes();
    //     if (i == 0) {
    //         prev_shape = this_sizes;
    //     }
    //     // if (prev_shape[0] != this_sizes[0] || prev_shape[1] != this_sizes[1]) {
    //     //     log->debug("Size mismatch in tensors");
    //     //     //TODO -- throw error
    //     // } 
    //     log->debug("Tensor {} has shape {} {}", i, this_sizes[0], this_sizes[1]);
    // }

    auto tensor_clone = in->tensors()->at(0)->tensor().clone();

    if (m_unsqueeze_input)
        tensor_clone = torch::unsqueeze(tensor_clone, 0);

    auto sizes = tensor_clone.sizes();
    // std::vector<int64_t> shape;
    for (const auto & s : sizes) {
        log->debug("shape {}", s);
    }


    //Get baseline + subtract for each channel
    auto median = std::get<0>(torch::median(tensor_clone, 2, true));
    auto nantensor = torch::zeros({1});
    nantensor[0] = NAN;
    torch::Device device((
        (torch::cuda::is_available() && !m_debug_force_cpu) ? torch::kCUDA : torch::kCPU
    ));
    nantensor = nantensor.to(device);
    // log->debug("nan device: {} input device: {}", nantensor.device(), tensor_clone.device());
    auto no_outliers = torch::where(torch::abs(tensor_clone - median) < 500., tensor_clone, nantensor);
    auto baseline = std::get<0>(torch::nanmedian(
        no_outliers, 2, true
    ));
    // std::cout << baseline << std::endl;

    tensor_clone = tensor_clone - baseline;
    for (const auto & s : tensor_clone.sizes()) {
        log->debug("shape {}", s);
    }
    auto rms_vals = torch::std(tensor_clone, 2, std::nullopt, true);
    for (const auto & s : rms_vals.sizes()) {
        log->debug("shape {}", s);
    }
    tensor_clone = torch::where(
        (tensor_clone > m_threshold_rms_factor*rms_vals),
        tensor_clone,
        torch::zeros({1}).to(device)
    );

    if (m_unsqueeze_input)
        tensor_clone = torch::squeeze(tensor_clone, 0);


    Configuration set_md, tensor_md;
    metadata_passthrough(in->tensors()->at(0)->metadata(), tensor_md, m_passthrough);
    // log->debug("Passed channel_map\n{}", tensor_md["channel_map"]);
    tensor_md["tag"] = m_output_tensor_tag;
    set_md["tag"] = m_output_set_tag;
    auto output_tensor = std::make_shared<SimpleTorchTensor>(tensor_clone, tensor_md);
    std::vector<ITorchTensor::pointer> itv{
        output_tensor
    };
    out = std::make_shared<SimpleTorchTensorSet>(
        in->ident(), set_md,
        std::make_shared<std::vector<ITorchTensor::pointer>>(itv)
    );
    log->debug("Tagged output with Set:{} Tensor:{}",
               out->metadata()["tag"],
               out->tensors()->at(0)->metadata()["tag"]);

    return true;
}